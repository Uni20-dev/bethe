"""Explicit precision in every biquadratic CLI family and dependent output."""
import csv
from decimal import Decimal, localcontext, InvalidOperation
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, available = sys.argv[1:]
enabled = available.upper() in ("ON", "TRUE", "1")
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")


def run(args, code=0):
    result = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=120)
    assert result.returncode == code, (args, result.returncode, result.stdout, result.stderr)
    return result.stdout


def tables(args, code=0, bits=160):
    document = json.loads(run([*args, "--precision", f"bits-{bits}", "--format", "json"], code))
    assert document["status"] == "complete"
    for table in document["tables"].values():
        assert int(table["metadata"]["Working precision (bits)"]) == bits
        for row in table["rows"]:
            for value in row:
                if isinstance(value, dict):
                    assert value["precision"] == bits, value
        assert table["summary"]["Outcome"] == ("success" if code == 0 else "partial")
    return document["tables"]


def records(table):
    return [dict(zip([column["id"] for column in table["columns"]],
                     [value["value"] if isinstance(value, dict) else value for value in row]))
            for row in table["rows"]]


def equal_tables(a, b):
    assert a.keys() == b.keys()
    for name in a:
        for field in ("columns", "rows"):
            assert a[name][field] == b[name][field], (name, field)


if not enabled:
    assert "UNI20_ENABLE_MPREAL=ON" in subprocess.run(
        [program, "4", "--precision", "bits-160"], capture_output=True, text=True, env=env, timeout=30).stderr
else:
    # Independent exact four-site spin Hamiltonian / one-defect eigenvalues.
    with localcontext() as context:
        context.prec = 120
        af = (-15 - Decimal(17).sqrt()) / 2
        singlet = (-15 + Decimal(17).sqrt()) / 2
        ferro = 6 - Decimal(2).sqrt()
        previous = Decimal(1)
        for bits, bound in ((80, "1e-18"), (160, "1e-42"), (256, "1e-70")):
            errors = []
            for args, expected in ((["4", "--roots"], af),
                                   (["4", "--singlet-excitation", "--roots"], singlet),
                                   (["4", "--q-seed", "1.5,-1.3333333333333333333333333333333333333", "--roots"], singlet),
                                   (["4", "--ferromagnetic", "--one-defect"], ferro)):
                data = tables(args, bits=bits)
                errors.append(abs(Decimal(records(data["states"])[0]["energy"]) - expected))
            error = max(errors)
            assert error < Decimal(bound), (bits, errors)
            assert error < previous * Decimal("1e-12"), (bits, error, previous)
            previous = error
        spectrum = tables(["4", "--q-spectrum", "--roots"])
        energies = sorted(Decimal(row["energy"]) for row in records(spectrum["states"]))
        assert len(energies) == 2
        assert max(abs(a - b) for a, b in zip(energies, [af, singlet])) < Decimal("1e-42")

        # Parsing is tested without Newton updates so a rounded binary64 seed
        # cannot be hidden by eventual convergence to the same solution.
        seed = ("1.500000000000000000000000000000000000000000000000000000000000001",
                "-1.333333333333333333333333333333333333333333333333333333333333333")
        data = tables(["4", "--q-seed", ",".join(seed), "--max-iterations", "0", "--roots"], 2, 333)
        coefficients = records(data["q_coefficients"])
        assert len(coefficients) == 2
        for row, value in zip(coefficients, seed):
            assert abs(Decimal(row["coefficient"]) - Decimal(value)) < Decimal("1e-95")
        assert records(data["states"])[0]["gap"] is None

    scans = [
        ["8", "--sectors"],
        ["8", "--through-lines", "4", "--excitations", "2"],
        ["9", "--ferromagnetic", "--through-lines", "5", "--excitations", "all", "--real-window", "3"],
        ["10", "--ferromagnetic", "--bound-pairs", "2"],
        ["12", "--ferromagnetic", "--bound-triples", "2"],
        ["12", "--ferromagnetic", "--bound-quartets", "2"],
        ["10", "--ferromagnetic", "--pair-defects", "2", "--mixed-window", "2"],
        ["10", "--ferromagnetic", "--triple-defects", "1", "--mixed-window", "1"],
        ["10", "--ferromagnetic", "--two-pair-states", "1", "--pair-window", "2"],
    ]
    for case in scans:
        for budget, code in (([], 0), (["--max-iterations", "0"], 2)):
            args = [*case, *budget, "--roots", "--spin-content"]
            expected = tables(args, code)
            equal_tables(expected, tables([*args, "--threads", "2"], code))
            equal_tables(expected, tables([*args, "--threads", "2", "--no-retain"], code))
            if code:
                assert any(not row["converged"] for table in expected.values() for row in records(table)
                           if "converged" in row)

    # Analytic ground-space endpoints also carry explicit precision.
    data = tables(["7", "--ferromagnetic", "--spin-content"])
    row, = records(data["states"])
    assert Decimal(row["energy"]) == 6 and Decimal(row["gap"]) == 0
    data = tables(["7", "--ferromagnetic", "--through-lines", "5"])
    assert Decimal(records(data["states"])[0]["gap"]) > 1

    args = ["4", "--singlet-excitation", "--roots"]
    expected = tables(args)
    for name, flag, delimiter in (("states", "--csv-table", ","), ("roots", "--tsv-table", "\t")):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "table.dat"
            run([*args, "--precision", "bits-160", flag, f"{name}={path}", "--quiet"])
            content = "\n".join(line for line in path.read_text().splitlines() if not line.startswith("#"))
            exported = list(csv.DictReader(io.StringIO(content), delimiter=delimiter))
            reference = records(expected[name])
            assert len(exported) == len(reference)
            for actual, wanted in zip(exported, reference):
                for key, value in wanted.items():
                    if isinstance(value, bool):
                        assert actual[key].lower() == str(value).lower()
                    elif value is None:
                        assert actual[key] == ""
                    else:
                        try:
                            assert Decimal(actual[key]) == Decimal(str(value)), (key, actual[key], value)
                        except InvalidOperation:
                            assert actual[key] == str(value), (key, actual[key], value)

    # digits-N and tolerance strings bypass native floating-point parsing.
    data = json.loads(run(["4", "--precision", "digits-100", "--tolerance", "1e-80", "--format", "json"]))
    metadata = data["tables"]["states"]["metadata"]
    assert int(metadata["Working precision (bits)"]) == 333
    with localcontext() as context:
        context.prec = 120
        assert abs(Decimal(metadata["Residual tolerance"]) / Decimal("1e-80") - 1) < Decimal("1e-95")

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "preserve.json"
    invalid = [["--precision", p] for p in ("bits-0", "bits-160junk", "digits--1")]
    if enabled:
        invalid += [["--precision", "bits-160", "--tolerance", t] for t in ("0", "-1", "nan", "inf")]
        invalid += [["--precision", "bits-160", "--q-seed", "nan,1"]]
    else:
        invalid += [["--precision", p] for p in ("bits-160", "digits-100")]
    for options in invalid:
        path.write_text("sentinel")
        run(["4", *options, "--json", str(path), "--force"], 1)
        assert path.read_text() == "sentinel"
print("Biquadratic explicit precision, full mode coverage and ordered exports passed")
