"""Bose-Fermi frontend: native precision, limits, equations, failure and exports."""
import csv
from decimal import Decimal, localcontext
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, fp128 = sys.argv[1:]
precisions = ["fp64", "long-double"] + (["fp128"] if fp128.upper() in ("ON", "TRUE", "1") else [])
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")
base = ["--bosons", "2", "--fermions", "3", "--length", "5", "--c", "1"]


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
        assert "BB/BF" in table["metadata"]["Hamiltonian"]
    return doc["tables"]


for precision in precisions:
    flags = ["--precision", precision, "--roots"]
    with localcontext() as ctx:
        ctx.prec = 100
        pi = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628")
        tol = {"fp64": Decimal("1e-13"), "long-double": Decimal("1e-16"), "fp128": Decimal("1e-30")}[precision]
        for b,f in ((1,1), (2,0)):
            exact = tables(["--bosons", str(b), "--fermions", str(f), "--length", "1", "--c", str(pi), *flags])
            assert abs(Decimal(rows(exact["states"])[0]["energy"])-pi*pi/2) < tol*10
            assert len(rows(exact["auxiliary_roots"])) == f
        weak = tables(["--bosons", "1", "--fermions", "1", "--length", "1", "--c", "1e-40", *flags])
        assert abs(Decimal(rows(weak["states"])[0]["energy"])/Decimal("2e-40")-1) < tol
        for b,f,c in ((0,0,"1"), (3,2,"0"), (0,2,"7"), (2,1,"0")):
            free = tables(["--bosons", str(b), "--fermions", str(f), "--length", "2", "--c", c,
                           "--max-iterations", "0", *flags])
            r = rows(free["states"])[0]
            energy = pi*pi*f*(f*f+(-1 if f%2 else 2))/12
            assert abs(Decimal(r["energy"])-energy) < tol*(1+energy)
            assert int(r["momentum_index"]) == (0 if f%2 else f//2)
            assert len(rows(free["free_modes"])) == b+f
    normal = tables([*base, *flags])
    charge, auxiliary = rows(normal["charge_roots"]), rows(normal["auxiliary_roots"])
    assert len(charge) == 5 and len(auxiliary) == 2
    for r in charge:
        k = float(r["k"])
        phase = 5*k + sum(2*math.atan(2*(k-float(a["lambda"]))) for a in auxiliary)
        assert abs(phase-2*math.pi*r["I"]) < 1e-12
    for r in auxiliary:
        phase = sum(2*math.atan(2*(float(r["lambda"])-float(k["k"]))) for k in charge)
        assert abs(phase-2*math.pi*r["J"]) < 1e-12
    failed = tables([*base, *flags, "--max-iterations", "0"], 2)
    r = rows(failed["states"])[0]
    assert r["status"] == "iteration_limit" and r["energy"] is None and r["momentum"] is None
    assert all(r["k"] is None for r in rows(failed["charge_roots"]))
    assert all(r["lambda"] is None for r in rows(failed["auxiliary_roots"]))
    pure_failed = tables(["--bosons", "3", "--fermions", "0", "--length", "1", "--c", "1",
                          "--max-iterations", "0", *flags], 2)
    assert rows(pure_failed["states"])[0]["energy"] is None
    assert rows(pure_failed["auxiliary_roots"]) == []
    assert all(r["k"] is None for r in rows(pure_failed["charge_roots"]))
    short_length = "1e-300" if precision == "fp64" else "1e-4900"
    overflow = tables(["--bosons", "0", "--fermions", "3", "--length", short_length, "--c", "0", *flags], 2)
    assert rows(overflow["states"])[0]["status"] == "precision_limit"
    assert rows(overflow["states"])[0]["energy"] is None
    assert all(r["k"] is None for r in rows(overflow["free_modes"]))
    assert normal["states"]["rows"] == tables([*base, *flags, "--no-retain"])["states"]["rows"]

invalid = [[], ["--bosons", "1"], ["--bosons", "2", "--fermions", "2", "--length", "5", "--c", "1"],
           ["--bosons", "-1", "--fermions", "3", "--length", "5", "--c", "1"],
           ["--bosons", "2", "--fermions", "3", "--length", "0", "--c", "1"],
           ["--bosons", "2", "--fermions", "3", "--length", "5", "--c", "-1"],
           ["--bosons", "2", "--fermions", "3", "--length", "5", "--c", "nan"],
           [*base, "--tolerance", "0"], [*base, "--max-iterations", "-1"]]
# No options intentionally shows help successfully.
for args in invalid[1:]:
    run(args, 1)
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)
    flags = ["--json", str(path/"state.json")]
    names = ("states", "charge_roots", "auxiliary_roots")
    for suffix in ("csv", "tsv"):
        for name in names:
            flags += [f"--{suffix}-table", f"{name}={path/f'{name}.{suffix}'}"]
    screen = run([*base, "--roots", *flags, "--format", "plain"])
    assert "CPU time" in screen and "Used for:" not in screen
    doc = json.loads((path/"state.json").read_text())["tables"]
    for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
        for name in names:
            text = (path/f"{name}.{suffix}").read_text()
            assert "BB/BF" in text
            exported = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
            expected = [{k: "" if v is None else str(v).lower() if isinstance(v, bool) else str(v) for k,v in r.items()} for r in rows(doc[name])]
            assert exported == expected
    old = (path/"state.json").read_text()
    run([*base, "--json", str(path/"state.json")], 1)
    for args in invalid[1:]:
        run([*args, "--json", str(path/"state.json"), "--force"], 1)
        assert (path/"state.json").read_text() == old
    run([*base, "--roots", *flags, "--force", "--no-retain", "--quiet"])
    refreshed = json.loads((path/"state.json").read_text())["tables"]
    assert all(refreshed[name]["rows"] == doc[name]["rows"] for name in names)
    run([*base, "--stream", "--no-retain", "--format", "plain"])
    run([*base, "--format", "pretty"])
    failure = run([*base, "--max-iterations", "0", "--format", "csv", "--no-preamble"], 2)
    assert next(csv.DictReader(failure.splitlines()))["energy"] == ""

print("Bose-Fermi CLI contracts passed")
