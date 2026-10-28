"""Ladder precision, descendant energies, branch budgets and root exports."""
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
    assert set(tables) == {"states", "representations", "roots"}
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
        exact = -Decimal("0.5") - Decimal(5).sqrt()
        errors = []
        for bits in (128, 256):
            args = ["6", "--rung", "0", "--singlets", "4", "--precision", f"bits-{bits}"]
            tables = result(args, bits)
            state = rows(tables["states"])[0]
            assert state["converged"] is True and state["descendant"] is True
            error = abs(Decimal(state["energy"]) - exact)
            assert error < Decimal(4096) * Decimal(2) ** (1-bits)
            errors.append(error)
            shifted = rows(result([*args, "--field", "1", "--sz", "1"], bits)["states"])[0]
            assert abs(Decimal(shifted["energy"])+1-Decimal(state["energy"])) < Decimal(8192)*Decimal(2)**(1-bits)
        assert errors[1] < Decimal("1e-30") * errors[0]
        tolerance = "1.0000000000000000000000000000000000000000000000001e-30"
        table = result(["4", "--rung", "1", "--precision", "digits-100", "--tolerance", tolerance], 333)["states"]
        assert abs(Decimal(table["metadata"]["Residual tolerance"]) - Decimal(tolerance)) < Decimal("1e-125")
        for budget in (0, 1):
            tables = result(["6", "--rung", "0", "--sectors", "--max-branches", str(budget),
                             "--precision", "bits-160"], 160, 2)
            states = rows(tables["states"])
            assert all(not s["converged"] and int(s["branches"]) == budget for s in states)
            assert any(s["energy"] is not None for s in states) == bool(budget)
        product = rows(result(["7", "--rung", "6", "--max-branches", "0", "--max-iterations", "0",
                               "--precision", "bits-160"], 160)["states"])[0]
        assert product["converged"] is True and Decimal(product["energy"]) == Decimal(-105)/4

    with tempfile.TemporaryDirectory() as directory:
        for case, arguments in enumerate((["4", "--rung", "1"], ["6", "--rung", "0", "--singlets", "4"])):
            args = [*arguments, "--precision", "bits-160"]
            expected = result(args, 160)
            for name, suffix, delimiter in (("roots", "csv", ","), ("representations", "tsv", "\t")):
                path = Path(directory) / f"{case}-{name}.{suffix}"
                run([*args, f"--{suffix}-table", f"{name}={path}", "--quiet"])
                exported = list(csv.DictReader([line for line in path.read_text().splitlines()
                                               if not line.startswith("#")], delimiter=delimiter))
                expected_text = [{key: str(value) for key, value in row.items()} for row in rows(expected[name])]
                assert exported == expected_text, (exported, expected_text)

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "preserve.json"
    path.write_text("sentinel")
    invalid = [["--precision", spelling] for spelling in ("bits-0", "bits-128junk", "digits--1")]
    if enabled:
        invalid += [["--precision", "bits-128", "--tolerance", t] for t in ("0", "-1", "nan", "inf")]
        invalid += [["--precision", "bits-128", "--field", t] for t in ("nan", "inf")]
    else:
        invalid += [["--precision", spelling] for spelling in ("bits-128", "digits-100")]
    for args in invalid:
        run(["4", "--rung", "1", *args, "--json", str(path), "--force"], 1)
        assert path.read_text() == "sentinel"

print("Ladder explicit precision and export checks passed")
