"""q-boson CLI precision, normalization, limits and shared export contracts."""
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
base = ["4", "--particles", "3", "--eta", "1"]


def run(args, status=0):
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=30, env=env)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def rows(t):
    return [dict(zip([c["id"] for c in t["columns"]], r)) for r in t["rows"]]


def tables(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete"
    for t in doc["tables"].values():
        assert t["summary"]["Outcome"] == ("success" if status == 0 else "partial")
        assert Decimal(t["summary"]["Run CPU seconds"]) >= 0
        assert "+2N" in t["metadata"]["Energy reference"]
    return doc["tables"]


for precision in precisions:
    flags = ["--precision", precision, "--roots"]
    with localcontext() as ctx:
        ctx.prec = 100
        tol = {"fp64": Decimal("1e-13"), "long-double": Decimal("1e-16"), "fp128": Decimal("1e-30")}[precision]
        for eta in ("1", "1e-60"):
            t = tables(["2", "--particles", "2", "--eta", eta, *flags])
            r = rows(t["states"])[0]
            loss = 1-(-2*Decimal(eta)).exp()
            exact = 2*loss/(1+(1-loss/2).sqrt())
            assert abs(Decimal(r["energy"])-exact) < tol*exact
            assert r["converged"] and r["momentum"] == "0"
            assert [x["I"] for x in rows(t["roots"])] == [-0.5, 0.5]
        phase = rows(tables(["2", "--particles", "2", "--phase", "--max-iterations", "0", *flags])["states"])[0]
        assert abs(Decimal(phase["energy"])-(4-2*Decimal(2).sqrt())) < tol
        assert phase["iterations"] == "0"
    for n in (0, 1, 4):
        free = tables(["4", "--particles", str(n), "--eta", "0", "--max-iterations", "0", *flags])
        assert rows(free["states"])[0]["energy"] == "0"
        assert len(rows(free["roots"])) == n and all(r["k"] == "0" for r in rows(free["roots"]))
    normal = tables([*base, *flags])
    root_rows = rows(normal["roots"])
    roots = [float(r["k"]) for r in root_rows]
    assert abs(sum(roots)) < 1e-12
    for j, k in enumerate(roots):
        phase = 4*k + sum(2*math.atan2(math.sin((k-q)/2), math.tanh(1)*math.cos((k-q)/2))
                         for i,q in enumerate(roots) if i != j)
        assert abs(phase - 2*math.pi*root_rows[j]["I"]) < 1e-12
    failed = tables([*base, *flags, "--max-iterations", "0"], 2)
    r = rows(failed["states"])[0]
    assert r["status"] == "iteration_limit" and r["energy"] is None and r["momentum"] is None
    assert all(x["k"] is None for x in rows(failed["roots"]))
    assert normal["states"]["rows"] == tables([*base, *flags, "--no-retain"])["states"]["rows"]
    for deformation in (["--eta", "0"], ["--eta", "1"], ["--phase"]):
        args = ["2", "--particles", "2", *deformation, "--excitations", "all", *flags]
        spectrum = tables(args)
        levels = rows(spectrum["states"])
        assert len(levels) == 3 and len(rows(spectrum["reference"])) == 1
        assert spectrum["states"]["metadata"]["Candidate count"] == "3"
        with localcontext() as ctx:
            ctx.prec = 100
            q2 = Decimal(1) if deformation[-1] == "0" else ((-Decimal(2)).exp() if deformation[-1] == "1" else Decimal(0))
            split = 2*(2*(1+q2)).sqrt()
            for level, exact in zip(levels, [4-split, Decimal(4), 4+split]):
                assert abs(Decimal(level["energy"])-exact) < tol*8
                assert abs(Decimal(level["gap"])-(exact-(4-split))) < tol*16
        assert [int(r["momentum_index"]) for r in levels] == [0, 1, 0]
        root_table = rows(spectrum["roots"])
        assert len(root_table) == 8
        for level in levels:
            rr = [r for r in root_table if r["state_id"] == level["state_id"]]
            assert sum(int(r["m"]) for r in rr) % 2 == int(level["momentum_index"])
            assert all(abs(float(r["k"])-math.pi*int(r["m"])-float(r["deviation"])) < 1e-12 for r in rr)
        assert spectrum["states"]["rows"] == tables([*args, "--no-retain"])["states"]["rows"]
    scan_args = [*base, *flags, "--excitations", "all"]
    complete = tables(scan_args)
    assert len(rows(complete["states"])) == math.comb(6, 3)
    limited = tables([*base, *flags, "--excitations", "2"])
    assert rows(limited["states"]) == rows(complete["states"])[:2]
    partial = tables([*scan_args, "--max-iterations", "0"], 2)
    assert "failed" in partial and rows(partial["reference"])[0]["energy"] is None
    assert all(r["gap"] is None and r["converged"] for r in rows(partial["states"]))
    assert rows(partial["failed"])[0]["energy"] is None
    vacuum = tables(["3", "--particles", "0", "--eta", "1", "--excitations", "all", *flags])
    assert len(rows(vacuum["states"])) == 1 and rows(vacuum["roots"]) == []

for args in (["4", "--particles", "3"], [*base, "--phase"], ["4", "--particles", "3", "--phase=false"],
             ["1", "--particles", "2", "--eta", "1"], ["4", "--particles", "-1", "--eta", "1"],
             ["4", "--particles", "3", "--eta", "-1"], ["4", "--particles", "3", "--eta", "inf"],
             ["4", "--particles", "3", "--eta", "nan"], [*base, "--tolerance", "0"],
             [*base, "--max-iterations", "-1"], [*base, "--max-candidates", "10"],
             [*base, "--excitations", "0"], [*base, "--excitations", "all", "--max-candidates", "1"],
             [*base, "--excitations", "all", "--max-candidates", "0"]):
    run(args, 1)
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)
    export_base = [*base, "--excitations", "all"]
    flags = ["--json", str(path/"state.json")]
    for suffix in ("csv", "tsv"):
        for name in ("states", "reference", "roots"):
            flags += [f"--{suffix}-table", f"{name}={path/f'{name}.{suffix}'}"]
    screen = run([*export_base, "--roots", *flags, "--format", "plain"])
    assert "CPU time" in screen and "Used for:" not in screen
    doc = json.loads((path/"state.json").read_text())["tables"]
    for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
        for name in ("states", "reference", "roots"):
            text = (path/f"{name}.{suffix}").read_text()
            assert "+2N" in text and "q-boson" in text
            exported = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
            expected = [{k: "" if v is None else str(v).lower() if isinstance(v, bool) else str(v) for k,v in r.items()} for r in rows(doc[name])]
            assert exported == expected
    old = (path/"state.json").read_text()
    run([*base, "--json", str(path/"state.json")], 1)
    assert (path/"state.json").read_text() == old
    run([*base, "--excitations", "all", "--max-candidates", "1", "--json", str(path/"state.json"), "--force"], 1)
    assert (path/"state.json").read_text() == old
    run(["1", "--particles", "2", "--eta", "1", "--json", str(path/"state.json"), "--force"], 1)
    assert (path/"state.json").read_text() == old
    run([*export_base, "--roots", *flags, "--force", "--no-retain", "--quiet"])
    refreshed = json.loads((path/"state.json").read_text())["tables"]
    assert all(refreshed[name]["rows"] == doc[name]["rows"] for name in doc)
    run([*export_base, "--stream", "--no-retain", "--format", "plain"])
    run([*export_base, "--format", "pretty"])
    failure = run([*base, "--max-iterations", "0", "--format", "csv", "--no-preamble"], 2)
    assert next(csv.DictReader(failure.splitlines()))["energy"] == ""

print("q-boson CLI contracts passed")
