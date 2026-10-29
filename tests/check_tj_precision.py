"""t-J explicit precision, branch exports and provider rejection."""
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


def result(args, bits, status=0, free=False):
    doc = json.loads(run([*args, "--roots", "--format", "json"], status))
    assert doc["status"] == "complete"
    tables = doc["tables"]
    assert set(tables) == ({"states", "free_modes"} if free else {"states", "first_roots", "second_roots"})
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
        exact = -3 - Decimal(5).sqrt()
        errors = []
        for bits in (128, 256):
            args = ["5", "--particles", "4", "--sz", "1", "--precision", f"bits-{bits}"]
            tables = result(args, bits)
            state = rows(tables["states"])[0]
            assert state["converged"] is True
            error = abs(Decimal(state["energy"]) - exact)
            assert error < Decimal(2048) * Decimal(2) ** (1-bits)
            errors.append(error)
            # This small sector starts at an analytic solution, so a zero
            # update budget is still a success. Use an interacting nontrivial
            # seed to check the incomplete-state export contract separately.
            seeded = rows(result([*args, "--max-iterations", "0"], bits)["states"])[0]
            assert seeded["converged"] is True
            assert int(seeded["iterations"]) == 0
            failed = result(["8", "--particles", "6", "--precision", f"bits-{bits}",
                             "--max-iterations", "0"], bits, 2)
            estimate = rows(failed["states"])[0]
            assert estimate["converged"] is False
            assert int(estimate["iterations"]) == 0
            assert Decimal(estimate["energy"]).is_finite()
        assert errors[1] < Decimal("1e-30") * errors[0]
        tolerance = "1.0000000000000000000000000000000000000000000000001e-30"
        table = result(["5", "--particles", "4", "--sz", "1", "--precision", "digits-100",
                        "--tolerance", tolerance], 333)["states"]
        assert abs(Decimal(table["metadata"]["Residual tolerance"]) - Decimal(tolerance)) < Decimal("1e-125")
        full = rows(result(["6", "--precision", "bits-256"], 256)["states"])[0]
        assert full["converged"] is True
        assert abs(Decimal(full["energy"]) + 5 + Decimal(13).sqrt()) < Decimal("1e-65")
        vacuum = rows(result(["7", "--particles", "0", "--precision", "bits-128"], 128, free=True)["states"])[0]
        assert vacuum["converged"] is True and Decimal(vacuum["energy"]) == 0

    with tempfile.TemporaryDirectory() as directory:
        for case, (arguments, free, names) in enumerate((
            (["3", "--particles", "2"], False, ("first_roots", "second_roots")),
            (["8", "--particles", "6"], False, ("first_roots", "second_roots")),
            (["7", "--particles", "4", "--sz", "2"], True, ("free_modes",)),
        )):
            args = [*arguments, "--precision", "bits-160"]
            expected = result(args, 160, free=free)
            for name in names:
                for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
                    path = Path(directory) / f"{case}-{name}.{suffix}"
                    run([*args, f"--{suffix}-table", f"{name}={path}", "--quiet"])
                    exported = list(csv.DictReader([line for line in path.read_text().splitlines()
                                                   if not line.startswith("#")], delimiter=delimiter))
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
        run(["5", "--particles", "4", "--sz", "1", *args, "--json", str(path), "--force"], 1)
        assert path.read_text() == "sentinel"

print("t-J explicit precision and export checks passed")
