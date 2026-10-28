"""TB explicit precision, complex-root exports and provider rejection."""
import csv
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, mpfr, provider = sys.argv[1:]
enabled = all(flag.upper() in ("ON", "TRUE", "1") for flag in (mpfr, provider))
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")


def run(args, status=0):
    p = subprocess.run([program, *args], capture_output=True, text=True, timeout=120, env=env)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def unwrap(value):
    return value["value"] if isinstance(value, dict) else value


def rows(table):
    names = [column["id"] for column in table["columns"]]
    return [dict(zip(names, map(unwrap, row))) for row in table["rows"]]


def result(args, bits, status=0):
    doc = json.loads(run([*args, "--roots", "--format", "json"], status))
    assert doc["status"] == "complete"
    tables = doc["tables"]
    assert set(tables) == {"states", "strings", "roots"}
    for table in tables.values():
        assert int(table["metadata"]["Working precision (bits)"]) == bits
        for row in table["rows"]:
            for value in row:
                if isinstance(value, dict):
                    assert value["precision"] == bits
    return tables


if enabled:
    with localcontext() as context:
        context.prec = 130
        exact = -11 - Decimal(41).sqrt()
        errors = []
        for bits in (128, 256):
            args = ["4", "--precision", f"bits-{bits}"]
            tables = result(args, bits)
            state = rows(tables["states"])[0]
            assert state["converged"] is True
            assert len(tables["strings"]["rows"]) == 2
            assert len(tables["roots"]["rows"]) == 4
            error = abs(Decimal(state["energy"]) - exact)
            assert error < Decimal(2048) * Decimal(2) ** (1-bits)
            errors.append(error)
            failed = result([*args, "--max-iterations", "0"], bits, 2)
            estimate = rows(failed["states"])[0]
            assert estimate["converged"] is False
            assert int(estimate["iterations"]) == 0
            assert Decimal(estimate["energy"]).is_finite()
        assert errors[1] < Decimal("1e-30") * errors[0]
        tolerance = "1.0000000000000000000000000000000000000000000000001e-30"
        table = result(["6", "--precision", "digits-100", "--tolerance", tolerance], 333)["states"]
        assert abs(Decimal(table["metadata"]["Residual tolerance"]) - Decimal(tolerance)) < Decimal("1e-125")

    with tempfile.TemporaryDirectory() as directory:
        for sites in (4, 6):
            args = [str(sites), "--precision", "bits-160"]
            expected = result(args, 160)
            for name, suffix, delimiter in (("roots", "csv", ","), ("strings", "tsv", "\t")):
                path = Path(directory) / f"{sites}-{name}.{suffix}"
                run([*args, f"--{suffix}-table", f"{name}={path}", "--quiet"])
                exported = list(csv.DictReader([line for line in path.read_text().splitlines()
                                               if not line.startswith("#")], delimiter=delimiter))
                # Exact half-integer labels are JSON numbers, but tabular files
                # spell them as decimal text; MPFR coordinates remain strings.
                expected_text = [{key: str(value) for key, value in row.items()}
                                 for row in rows(expected[name])]
                assert exported == expected_text, (exported, expected_text)

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "preserve.json"
    path.write_text("sentinel")
    invalid = [["--precision", spelling] for spelling in ("bits-0", "bits-128junk", "digits--1")]
    if enabled:
        invalid += [["--precision", "bits-128", "--tolerance", tolerance] for tolerance in ("0", "-1", "nan", "inf")]
    else:
        invalid += [["--precision", spelling] for spelling in ("bits-128", "digits-100")]
    for args in invalid:
        run(["4", *args, "--json", str(path), "--force"], 1)
        assert path.read_text() == "sentinel"

print("TB explicit precision and export checks passed")
