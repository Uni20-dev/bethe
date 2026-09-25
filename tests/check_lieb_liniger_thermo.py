"""Native thermodynamic curves, failure honesty, and shared output contracts."""
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


def run(args, status=0):
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=60, env=env)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]


def table(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete"
    t = doc["tables"]["dispersion"]
    assert t["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    assert Decimal(t["summary"]["Run CPU seconds"]) >= 0
    assert t["metadata"]["Energy reference"] == "fixed-N excitation gap above the bulk ground state"
    return t


for precision in precisions:
    base = ["--c", "4", "--points", "9", "--precision", precision]
    t = table(base)
    r = rows(t)
    assert len(r) == 18 and all(x["status"] == "converged" for x in r)
    for a, b in zip(r[9:], r[:8:-1]):
        # Independently rounded grid momenta need not be bitwise mirror images.
        assert abs(Decimal(a["energy"])-Decimal(b["energy"])) < Decimal("1e-13")
    assert Decimal(r[8]["energy"]) > Decimal(r[1]["energy"]) > 0
    assert r[0]["energy"] == r[9]["energy"] == r[-1]["energy"] == "0"
    assert t["rows"] == table([*base, "--no-retain"])["rows"]
    for branch in ("type-i", "type-ii"):
        single = ["--c", "4", "--momentum", "1", "--branch", branch, "--precision", precision]
        a = rows(table(single))[0]
        b = rows(table(["--c", "8", "--density", "2", "--momentum", "2",
                        "--branch", branch, "--precision", precision]))[0]
        with localcontext() as ctx:
            ctx.prec = 70
            assert abs(Decimal(b["energy"])/Decimal(a["energy"])-4) < Decimal("1e-14")
        tiny = rows(table(["--c", "4", "--momentum", "1e-60", "--branch", branch,
                           "--precision", precision]))[0]
        assert Decimal(tiny["energy"]) > 0 and Decimal(tiny["edge_distance"]) > 0
        failed = rows(table([*single, "--max-iterations", "0"], 2))[0]
        assert failed["energy"] is None and failed["rapidity"] is None
        assert failed["status"] == "momentum_limit"
        strong = rows(table(["--c", "1e40", "--momentum", "1", "--branch", branch,
                             "--precision", precision]))[0]
        with localcontext() as ctx:
            ctx.prec = 70
            pi = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406")
            expected = 2*pi + (1 if branch == "type-i" else -1)
            tol = {"fp64": Decimal("1e-12"), "long-double": Decimal("1e-15"), "fp128": Decimal("1e-28")}[precision]
            assert abs(Decimal(strong["energy"])-expected) < tol
    for option, value, expected in (("--max-background-iterations", "0", "density_limit"),
                                     ("--max-nodes", "16", "mesh_limit")):
        failed = rows(table([*base, option, value], 2))
        assert all(x["energy"] is None and x["status"] == expected for x in failed)

base = ["--c", "4", "--points", "5"]
assert len(rows(table(["--c", "4", "--branch", "type-i", "--p-max", "20", "--points", "3"]))) == 3
for args in (["--c", "0"], ["--c", "-1"], ["--c", "nan"], ["--c", "1", "--density", "0"],
             ["--c", "1", "--points", "1"], ["--c", "1", "--points", "-1"],
             ["--c", "1", "--tolerance", "0"], ["--c", "1", "--max-nodes", "513"],
             ["--c", "1", "--momentum", "7"], ["--c", "1", "--momentum", "-1"],
             ["--c", "1", "--branch", "type-ii", "--p-max", "1"],
             ["--c", "1", "--momentum", "1", "--points", "3"],
             ["--c", "1", "--momentum", "1", "--p-max", "2"]):
    run(args, 1)
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)
    js, cs, ts = [path / ("curves."+suffix) for suffix in ("json", "csv", "tsv")]
    screen = run([*base, "--json", str(js), "--csv", str(cs), "--tsv", str(ts), "--format", "plain"])
    assert "CPU time" in screen and "Used for:" not in screen
    t = json.loads(js.read_text())["tables"]["dispersion"]
    expected = rows(t)
    for file, delimiter in ((cs, ","), (ts, "\t")):
        text = file.read_text()
        assert "#" in text and "fixed-N" in text
        exported = list(csv.DictReader([line for line in text.splitlines() if not line.startswith("#")], delimiter=delimiter))
        assert exported == expected
    old = js.read_text()
    run([*base, "--json", str(js)], 1)
    assert js.read_text() == old
    run(["--c", "-1", "--json", str(js), "--force"], 1)
    assert js.read_text() == old
    run([*base, "--json", str(js), "--force", "--quiet", "--no-retain"])
    assert json.loads(js.read_text())["tables"]["dispersion"]["rows"] == t["rows"]
    run([*base, "--stream", "--no-retain", "--format", "plain"])
    run([*base, "--format", "pretty"])
    failed_csv = run([*base, "--max-background-iterations", "0", "--format", "csv", "--no-preamble"], 2)
    assert all(x["energy"] == "" for x in csv.DictReader(failed_csv.splitlines()))

print("Lieb-Liniger thermodynamic CLI contracts passed")
