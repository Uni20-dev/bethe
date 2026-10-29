"""Explicit Potts precision, ordered Async rows, and complete export checks."""
import csv
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, mpfr = sys.argv[1:]
enabled = mpfr.upper() in ("ON", "TRUE", "1")
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")


def run(args, status=0):
    p = subprocess.run([program, *args], capture_output=True, text=True, timeout=120, env=env)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def unwrap(value):
    return value["value"] if isinstance(value, dict) else value


def records(table):
    names = [c["id"] for c in table["columns"]]
    return [dict(zip(names, map(unwrap, row))) for row in table["rows"]]


def result(args, bits, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete"
    assert set(doc["tables"]) == {"levels"}
    table = doc["tables"]["levels"]
    assert int(table["metadata"]["Working precision (bits)"]) == bits
    assert table["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    for row in table["rows"]:
        for value in row:
            if isinstance(value, dict):
                assert value["precision"] == bits
    return table


if enabled:
    with localcontext() as context:
        context.prec = 130
        exact = -2 - 2 * Decimal(3).sqrt()
        charged = -(1 + Decimal(57).sqrt()) / 2
        errors = []
        for bits in (128, 256):
            args = ["2", "--momentum-index", "0", "--precision", f"bits-{bits}"]
            t = result(args, bits)
            rows = records(t)
            assert len(rows) == 3
            assert all(row["status"] == "converged" for row in rows)
            errors.append(abs(Decimal(rows[0]["energy"]) - exact))
            assert errors[-1] < 8192 * Decimal(2) ** (1-bits)
            assert abs(Decimal(rows[1]["energy"]) - charged) < 8192 * Decimal(2) ** (1-bits)
            assert rows[1]["energy"] == rows[2]["energy"]
            assert rows[0]["momentum"] == rows[0]["gap"] == rows[0]["x_scaled"] == "0"
            all_args = ["3", "--precision", f"bits-{bits}"]
            direct = result(all_args, bits)
            parallel = result([*all_args, "--threads", "4", "--no-retain"], bits)
            assert direct["rows"] == parallel["rows"]
            assert parallel["metadata"]["Effective concurrency limit"] == "3"
            assert [(row["charge"], int(row["k"])) for row in records(parallel)] == [
                (0, 0), (1, 0), (-1, 0), (1, 1), (-1, 1), (1, 2), (-1, 2)]
            failed = result([*all_args, "--max-iterations", "0"], bits, 2)
            failed_parallel = result([*all_args, "--max-iterations", "0", "--threads", "4"], bits, 2)
            assert failed["rows"] == failed_parallel["rows"]
            assert all(row["energy"] is None and row["gap"] is None and row["x_scaled"] is None
                       and row["status"] == "iteration_limit" and row["iterations"] == "0"
                       for row in records(failed))
        assert errors[1] < errors[0] * Decimal("1e-30")
        tolerance = "1.0000000000000000000000000000000000000000000000001e-30"
        t = result(["3", "--precision", "digits-100", "--tolerance", tolerance], 333)
        assert abs(Decimal(t["metadata"]["Residual tolerance"]) - Decimal(tolerance)) < Decimal("1e-125")

    with tempfile.TemporaryDirectory() as directory:
        args = ["3", "--precision", "bits-160"]
        expected = [{k: "" if v is None else str(v) for k, v in row.items()}
                    for row in records(result(args, 160))]
        for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
            path = Path(directory) / ("levels." + suffix)
            run([*args, f"--{suffix}-table", f"levels={path}", "--threads", "4", "--quiet"])
            text = path.read_text()
            assert "Hamiltonian" in text and "Momentum convention" in text
            actual = list(csv.DictReader([line for line in text.splitlines() if not line.startswith("#")],
                                        delimiter=delimiter))
            assert actual == expected

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "preserve.json"
    path.write_text("sentinel")
    invalid = [["--precision", value] for value in ("bits-0", "bits-128junk", "digits--1")]
    if enabled:
        invalid += [["--precision", "bits-128", "--tolerance", value] for value in ("0", "-1", "nan", "inf")]
    else:
        invalid += [["--precision", value] for value in ("bits-128", "digits-100")]
    for options in invalid:
        run(["3", *options, "--json", str(path), "--force"], 1)
        assert path.read_text() == "sentinel"

print("Potts explicit precision and export checks passed")
