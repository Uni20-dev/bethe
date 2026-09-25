"""TASEP frontend: decay/frequency conventions, native rates and output contracts."""
import csv
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, fp128 = sys.argv[1:]
precisions = ["fp64", "long-double"] + (["fp128"] if fp128.upper() in ("ON", "TRUE", "1") else [])
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")
base = ["5", "--particles", "2"]


def run(args, status=0):
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=30, env=env)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]], r)) for r in table["rows"]]


def tables(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete"
    for table in doc["tables"].values():
        assert table["summary"]["Outcome"] == ("success" if status == 0 else "partial")
        assert Decimal(table["summary"]["Run CPU seconds"]) >= 0
        assert "not quantum energies" in table["metadata"]["Units"]
        assert "columns sum to zero" in table["metadata"]["Generator"]
        assert "energy" not in [c["id"] for c in table["columns"]]
    return doc["tables"]


for precision in precisions:
    flags = ["--precision", precision, "--roots"]
    with localcontext() as ctx:
        ctx.prec = 100
        tol = {"fp64": Decimal("1e-12"), "long-double": Decimal("1e-15"), "fp128": Decimal("1e-29")}[precision]
        expected_gap = Decimal("0.710290069123028672493762532365387684613821535715147285599991234743117333842877636225644555")
        expected_frequency = Decimal("0.326895513054695736451443863914497675039691838878110784103129825870014387978645791366909038")
        for n in (2, 3):
            for rate in (Decimal("1"), Decimal("3.125"), Decimal("1e-40")):
                result = tables(["5", "--particles", str(n), "--rate", str(rate), *flags])
                r = rows(result["relaxation"])[0]
                assert r["has_mode"] and r["converged"] and r["status"] == "converged"
                assert abs(Decimal(r["gap"])/rate-expected_gap) < tol
                assert abs(Decimal(r["frequency"])/rate-expected_frequency) < tol
                assert Decimal(r["lambda_real"]) == -Decimal(r["gap"])
                assert Decimal(r["lambda_imag"]) == Decimal(r["frequency"])
                assert len(rows(result["roots"])) == 2
        for n in (1, 3):
            analytic = tables(["4", "--particles", str(n), "--max-iterations", "0",
                               "--max-seed-iterations", "0", *flags])
            r = rows(analytic["relaxation"])[0]
            assert abs(Decimal(r["gap"])-1) < tol and abs(Decimal(r["frequency"])-1) < tol
            assert int(r["iterations"]) == 0 and int(r["seed_iterations"]) == 0
    half = tables(["8", "--particles", "4", *flags])
    assert Decimal(rows(half["relaxation"])[0]["frequency"]) == 0
    for n in (0, 8):
        empty = tables(["8", "--particles", str(n), *flags])
        r = rows(empty["relaxation"])[0]
        assert not r["has_mode"] and r["status"] == "stationary_only" and r["converged"]
        assert all(r[k] is None for k in ("gap", "lambda_real", "lambda_imag", "frequency"))
        assert rows(empty["roots"]) == []
    for budget, status in (("--max-iterations", "iteration_limit"), ("--max-seed-iterations", "seed_limit")):
        failed = tables([*base, *flags, budget, "0"], 2)
        r = rows(failed["relaxation"])[0]
        assert r["has_mode"] and not r["converged"] and r["status"] == status
        assert all(r[k] is None for k in ("gap", "lambda_real", "lambda_imag", "frequency"))
        assert all(r["Z_real"] is None and r["Z_imag"] is None for r in rows(failed["roots"]))
    huge = "1e308" if precision == "fp64" else "1e4932"
    overflow = tables(["2", "--particles", "1", "--rate", huge, *flags], 2)
    assert rows(overflow["relaxation"])[0]["status"] == "precision_limit"
    normal = tables([*base, *flags])
    assert normal["relaxation"]["rows"] == tables([*base, *flags, "--no-retain"])["relaxation"]["rows"]

invalid = [["4"], ["1", "--particles", "0"], ["4", "--particles", "5"], ["4", "--particles", "-1"],
           [*base, "--rate", "0"], [*base, "--rate", "-1"], [*base, "--rate", "nan"],
           [*base, "--tolerance", "0"], [*base, "--max-iterations", "-1"],
           [*base, "--max-seed-iterations", "-1"], [*base, "--max-sites", "4"],
           [*base, "--left-rate", "1"], [*base, "--excitations", "all"]]
for args in invalid:
    run(args, 1)
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)
    flags = ["--json", str(path/"result.json")]
    names = ("relaxation", "roots")
    for suffix in ("csv", "tsv"):
        for name in names:
            flags += [f"--{suffix}-table", f"{name}={path/f'{name}.{suffix}'}"]
    screen = run([*base, "--roots", *flags, "--format", "plain"])
    assert "CPU time" in screen and "Used for:" not in screen
    doc = json.loads((path/"result.json").read_text())["tables"]
    for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
        for name in names:
            text = (path/f"{name}.{suffix}").read_text()
            assert "not quantum energies" in text
            exported = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
            expected = [{k: "" if v is None else str(v).lower() if isinstance(v, bool) else str(v) for k,v in r.items()} for r in rows(doc[name])]
            assert exported == expected
    old = (path/"result.json").read_text()
    run([*base, "--json", str(path/"result.json")], 1)
    for args in invalid:
        run([*args, "--json", str(path/"result.json"), "--force"], 1)
        assert (path/"result.json").read_text() == old
    run([*base, "--roots", *flags, "--force", "--no-retain", "--quiet"])
    refreshed = json.loads((path/"result.json").read_text())["tables"]
    assert all(refreshed[name]["rows"] == doc[name]["rows"] for name in names)
    run([*base, "--roots", "--stream", "--no-retain", "--format", "plain"])
    run([*base, "--format", "pretty"])
    failure = run([*base, "--max-iterations", "0", "--format", "csv", "--no-preamble"], 2)
    assert next(csv.DictReader(failure.splitlines()))["gap"] == ""

print("TASEP CLI contracts passed")
