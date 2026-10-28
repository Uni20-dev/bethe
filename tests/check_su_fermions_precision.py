"""Explicit SU(n) precision, accepted-stage semantics and complete root exports."""
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
pi = "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651"


def run(args, status=0):
    p = subprocess.run([program, *args], capture_output=True, text=True, timeout=120, env=env)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def args(pop="1,1,0", c=pi, precision="bits-160"):
    return ["--populations", pop, "--length", "1", "--c", c, "--precision", precision]


def unwrap(value):
    return value["value"] if isinstance(value, dict) else value


def rows(table):
    names = [column["id"] for column in table["columns"]]
    return [dict(zip(names, map(unwrap, row))) for row in table["rows"]]


def result(arguments, bits, status=0):
    doc = json.loads(run([*arguments, "--roots", "--format", "json"], status))
    assert doc["status"] == "complete"
    tables = doc["tables"]
    assert "states" in tables and "components" in tables
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
        exact = Decimal(pi)**2 / 2
        errors = []
        for bits in (128, 256):
            tables = result(args(precision=f"bits-{bits}"), bits)
            state = rows(tables["states"])[0]
            assert state["converged"] is True
            errors.append(abs(Decimal(state["energy"])-exact))
            assert errors[-1] < Decimal(4096)*Decimal(2)**(1-bits)
        assert errors[1] < Decimal("1e-30")*errors[0]
        tolerance = "1.0000000000000000000000000000000000000000000000001e-30"
        table = result([*args(c="1", precision="digits-100"), "--tolerance", tolerance], 333)["states"]
        assert abs(Decimal(table["metadata"]["Residual tolerance"])-Decimal(tolerance)) < Decimal("1e-125")

    for stages in (0, 1):
        tables = result([*args(pop="3,1,1", c="1"), "--max-stages", str(stages)], 160, 2)
        state = rows(tables["states"])[0]
        assert state["converged"] is False
        assert int(state["stages"]) == stages
        if stages:
            assert Decimal(state["reached_c"]) > 1
            assert Decimal(state["energy"]).is_finite()
            assert tables["roots"]["rows"]
        else:
            assert state["energy"] is None and state["reached_c"] is None
            assert not tables["roots"]["rows"]
    vacuum = rows(result(args(pop="0,0,0", c="0"), 160)["states"])[0]
    assert vacuum["converged"] is True and Decimal(vacuum["energy"]) == 0

    with tempfile.TemporaryDirectory() as directory:
        for pop, coupling, name, suffix, delimiter in (
            ("3,1,1", "1", "roots", "csv", ","),
            ("2,1,0", "0", "free_modes", "tsv", "\t"),
        ):
            arguments = args(pop=pop, c=coupling)
            expected = result(arguments, 160)
            path = Path(directory) / f"{name}.{suffix}"
            run([*arguments, f"--{suffix}-table", f"{name}={path}", "--quiet"])
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
    else:
        invalid += [["--precision", spelling] for spelling in ("bits-128", "digits-100")]
    for extra in invalid:
        run(["--populations", "1,1", "--length", "1", "--c", "1", *extra,
             "--json", str(path), "--force"], 1)
        assert path.read_text() == "sentinel"
print("SU(n) fermion explicit precision and complete exports passed")
