"""Bidirectional exclusion CLI: native rates, limits, budgets and export contracts."""
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
base = ["5", "--particles", "2", "--left-rate", "0.5"]


def run(args, status=0):
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=30, env=env)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]


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
        gap = Decimal("1.04039558904641505636142345412607413799326810221158347669691888089548565042267332323361234")
        frequency = Decimal("0.159168821610180253202429215239399086838149271486788653791649970802963326547367324815530809")
        for n in (2, 3):
            for scale in (Decimal("1"), Decimal("3.125"), Decimal("1e-40")):
                for reverse in (False, True):
                    right, left = (scale / 2, scale) if reverse else (scale, scale / 2)
                    result = tables(["5", "--particles", str(n), "--right-rate", str(right), "--left-rate", str(left), *flags])
                    r = rows(result["relaxation"])[0]
                    assert r["has_mode"] and r["converged"] and r["status"] == "converged"
                    assert abs(Decimal(r["gap"])/scale-gap) < tol
                    assert abs(Decimal(r["frequency"])/scale-frequency) < tol
                    assert Decimal(r["lambda_real"]) == -Decimal(r["gap"])
                    assert Decimal(r["lambda_imag"]) == Decimal(r["frequency"])
                    assert len(rows(result["roots"])) == 2
                    meta = result["roots"]["metadata"]
                    assert int(meta["Wave root index"]) in (0, 1)
                    assert Decimal(meta["Last reached rate ratio"]) == Decimal("0.5")
        for n in (1, 3):
            analytic = tables(["4", "--particles", str(n), "--right-rate", "2", "--left-rate", "1",
                               "--max-iterations", "0", "--max-seed-iterations", "0", *flags])
            r = rows(analytic["relaxation"])[0]
            assert abs(Decimal(r["gap"])-3) < tol and abs(Decimal(r["frequency"])-1) < tol
            assert rows(analytic["roots"]) == []
        symmetric = tables(["4", "--particles", "2", "--right-rate", "2", "--left-rate", "2", *flags])
        r = rows(symmetric["relaxation"])[0]
        assert abs(Decimal(r["gap"])-4) < tol and Decimal(r["frequency"]) == 0
        assert rows(symmetric["roots"]) == []
        # Bias retained in native parsing, far below double resolution in fp128.
        delta = {"fp64": Decimal("1e-12"), "long-double": Decimal("1e-16"), "fp128": Decimal("1e-30")}[precision]
        near = tables(["5", "--particles", "2", "--left-rate", str(1-delta), *flags])
        m = near["relaxation"]["metadata"]
        actual_delta = Decimal(m["Right-hop rate"]) - Decimal(m["Left-hop rate"])
        drift = Decimal(rows(near["relaxation"])[0]["frequency"])/actual_delta
        assert abs(drift-Decimal("0.3170188387650511907")) < Decimal("0.001")
    tasep = tables(["5", "--particles", "2", *flags])
    assert abs(Decimal(rows(tasep["relaxation"])[0]["gap"])-Decimal("0.71029006912302867249")) < Decimal("1e-12")
    left_only = tables(["5", "--particles", "2", "--right-rate", "0", "--left-rate", "1", *flags])
    assert left_only["relaxation"]["rows"] == tasep["relaxation"]["rows"]
    for n in (0, 8):
        empty = tables(["8", "--particles", str(n), *flags])
        r = rows(empty["relaxation"])[0]
        assert not r["has_mode"] and r["status"] == "stationary_only" and r["converged"]
        assert all(r[k] is None for k in ("gap", "lambda_real", "lambda_imag", "frequency"))
        assert rows(empty["roots"]) == []
    for budget, status in (("--max-iterations", "iteration_limit"), ("--max-continuation-steps", "continuation_limit"),
                           ("--max-seed-iterations", "seed_limit"), ("--max-seed-newton-iterations", "seed_limit")):
        failed = tables([*base, *flags, budget, "0"], 2)
        r = rows(failed["relaxation"])[0]
        assert r["has_mode"] and not r["converged"] and r["status"] == status
        assert all(r[k] is None for k in ("gap", "lambda_real", "lambda_imag", "frequency"))
        assert all(r["v_real"] is None and r["v_imag"] is None for r in rows(failed["roots"]))
    huge = "1e308" if precision == "fp64" else "1e4932"
    overflow = tables(["2", "--particles", "1", "--right-rate", huge, *flags], 2)
    assert rows(overflow["relaxation"])[0]["status"] == "precision_limit"
    assert tables([*base, *flags])["relaxation"]["rows"] == tables([*base, *flags, "--no-retain"])["relaxation"]["rows"]

invalid = [["4"], ["1", "--particles", "0"], ["4", "--particles", "5"], ["4", "--particles", "-1"],
           [*base, "--right-rate", "-1"], [*base, "--right-rate", "nan"], [*base, "--left-rate", "-1"],
           ["4", "--particles", "2", "--right-rate", "0", "--left-rate", "0"],
           [*base, "--tolerance", "0"], [*base, "--seed-tolerance", "0"],
           ["4", "--particles", "1", "--seed-tolerance", "nan"],
           [*base, "--max-iterations", "-1"], [*base, "--max-continuation-steps", "-1"],
           [*base, "--max-sites", "4"], [*base, "--excitations", "all"]]
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
    failure = run([*base, "--max-continuation-steps", "0", "--format", "csv", "--no-preamble"], 2)
    assert next(csv.DictReader(failure.splitlines()))["gap"] == ""

print("ASEP CLI contracts passed")
