"""Finite XYZ: explicit precision, independent energy, failures and exports."""
import csv
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, available = sys.argv[1:]
enabled = available.upper() in ("ON", "TRUE", "1")
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")
base = ["4", "--eta", "0.4", "--t", "0.7"]


def run(args, status=0):
    result = subprocess.run([program, *args], text=True, capture_output=True, env=env, timeout=60)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    return result.stdout


def tables(args, status=0):
    document = json.loads(run([*args, "--format", "json"], status))
    assert document["status"] == "complete"
    for table in document["tables"].values():
        assert table["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    return document["tables"]


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]],
                     (v["value"] if isinstance(v, dict) else v for v in row))) for row in table["rows"]]


def check_precision(result, selector, bits):
    for table in result.values():
        assert int(table["metadata"]["Working precision (bits)"]) == bits
        assert table["metadata"]["Precision"] == selector
        for row in table["rows"]:
            for column, cell in zip(table["columns"], row):
                if column["type"] == "real" and cell is not None:
                    assert isinstance(cell, dict) and cell["precision"] == bits, (column, cell)


with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    sentinel = root / "keep.json"
    invalid = ["mpfr", "bits-0", "digits--1"] + ([] if enabled else ["bits-256", "digits-100"])
    for precision in invalid:
        sentinel.write_text("preserve")
        run([*base, "--precision", precision, "--json", str(sentinel), "--force"], 1)
        assert sentinel.read_text() == "preserve"
    if not enabled:
        sys.exit(0)
    with localcontext() as ctx:
        ctx.prec = 160
        # Full spin Hamiltonian, evaluated independently at 110 and 140 decimal digits.
        exact = Decimal("-1.78084826309600006130342170385401606120750417562005774531984592688108743371478564411325843404138151192307371375133700476524497621961684")
        previous = Decimal(1)
        for selector, bits in [("bits-80", 80), ("bits-160", 160), ("bits-256", 256), ("digits-100", 333)]:
            args = [*base, "--precision", selector, "--roots"]
            result = tables(args)
            check_precision(result, selector, bits)
            state = rows(result["states"])[0]
            assert state["converged"] and state["status"] == "converged"
            energy = Decimal(state["energy"])
            error = abs(energy-exact)
            assert error < Decimal(2)**(-bits+13), (selector, error)
            assert error < previous*Decimal("1e-15"), (selector, error, previous)
            if bits >= 256:
                assert error < Decimal("1e-70")
            previous = error
            assert abs(Decimal(state["energy_per_site"])*4-energy) < Decimal(2)**(-bits+13)
            assert Decimal(state["momentum"]) == 0
            roots = rows(result["roots"])
            assert [r["I"] for r in roots] == [-0.5, 0.5]
            assert all(Decimal(r["lambda_real"]) == 0 for r in roots)
            assert Decimal(roots[0]["lambda_imag"]) == -Decimal(roots[1]["lambda_imag"])
            no_retain = tables([*args, "--no-retain"])
            assert all(t["rows"] == no_retain[name]["rows"] for name, t in result.items())

            failed = tables([*args, "--max-iterations", "0"], 2)
            check_precision(failed, selector, bits)
            failure = rows(failed["states"])[0]
            assert not failure["converged"] and failure["status"] == "iteration_limit"
            assert all(failure[key] is None for key in ["energy", "energy_per_site", "momentum"])
            assert all(r["lambda_real"] is None and r["lambda_imag"] is None for r in rows(failed["roots"]))

            files = ["--json", str(root/"state.json")]
            for suffix in ["csv", "tsv"]:
                for name in ["states", "roots"]:
                    files += [f"--{suffix}-table", f"{name}={root/f'{name}.{suffix}'}"]
            # The roots export requests that table without selecting screen output.
            run([*base, "--precision", selector, *files, "--no-retain", "--quiet", "--force"])
            exported = json.loads((root/"state.json").read_text())["tables"]
            assert all(t["rows"] == exported[name]["rows"] for name, t in result.items())
            for suffix, separator in [("csv", ","), ("tsv", "\t")]:
                for name, table in result.items():
                    text = (root/f"{name}.{suffix}").read_text()
                    assert "Working precision (bits)" in text
                    records = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=separator))
                    expected = [{k: "" if v is None else str(v).lower() if isinstance(v, bool) else str(v)
                                 for k, v in row.items()} for row in rows(table)]
                    assert records == expected

        invalid = [["3", "--eta", "0.4", "--t", "0.7"],
                   ["4", "--eta", "nan", "--t", "0.7"],
                   ["4", "--eta", "0.4", "--t", "inf"],
                   [*base, "--tolerance", "nan"], [*base, "--tolerance", "0"]]
        for args in invalid:
            sentinel.write_text("preserve")
            run([*args, "--precision", "bits-256", "--json", str(sentinel), "--force"], 1)
            assert sentinel.read_text() == "preserve"

print("Finite XYZ MPFR contracts passed")
