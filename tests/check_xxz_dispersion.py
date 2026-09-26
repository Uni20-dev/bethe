"""XXZ thermodynamic physics, native precision, failure and export contracts."""
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
    p = subprocess.run([program, *args], capture_output=True, text=True, timeout=30, env=env)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def records(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]


def table(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete"
    t = doc["tables"]["dispersion"]
    assert t["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    assert Decimal(t["summary"]["Run CPU seconds"]) >= 0
    assert "zero-field thermodynamic" in t["metadata"]["Energy reference"]
    return t


for precision in precisions:
    base = ["--delta", "2", "--points", "5", "--precision", precision]
    t = table(base)
    r = records(t)
    assert len(r) == 10 and all(row["status"] == "converged" for row in r)
    assert all(row["lower"] is None and row["upper"] is None for row in r[:5])
    assert all(row["energy"] is None for row in r[5:])
    with localcontext() as ctx:
        ctx.prec = 80
        tol = {"fp64": Decimal("5e-14"), "long-double": Decimal("5e-17"), "fp128": Decimal("5e-31")}[precision]
        gap = Decimal("0.1949011357100147072770116351591698798800292182717762239260031033")
        maximum = Decimal("2.070496154807134434659160021557335658575303075001046084961949819")
        assert abs(Decimal(r[0]["energy"]) - gap) < tol
        assert r[0]["energy"] == r[4]["energy"]
        assert abs(Decimal(r[2]["energy"]) - maximum) < tol
        assert abs(Decimal(r[5]["lower"]) - 2*gap) < tol
        assert r[5]["lower"] == r[5]["upper"]
        assert abs(Decimal(r[7]["upper"]) - 2*maximum) < tol
        assert abs(Decimal(r[6]["lower"]) - (maximum+gap)) < tol
        folded = records(table([*base, "--folded"]))
        assert abs(Decimal(folded[5]["upper"]) - 2*maximum) < tol
        assert folded[5]["upper"] == folded[7]["upper"]
        assert folded[5]["lower"] == folded[7]["lower"]
        scaled = records(table([*base, "--exchange", "3"]))
        for x, y in zip(r, scaled):
            for key in ("energy", "lower", "upper"):
                if x[key] is not None:
                    assert abs(Decimal(y[key])-3*Decimal(x[key])) < 10*tol
        near = records(table(["--delta", "1.0078125", "--branch", "spinon", "--momentum", "0", "--precision", precision]))[0]
        expected = Decimal("4.3944033654951077815566705015087735721807786658387465951809322194e-17")
        assert abs(Decimal(near["energy"])/expected-1) < 100*tol
    assert r[0]["cell_momentum"] == r[4]["cell_momentum"] == r[7]["cell_momentum"] == "0"
    assert t["rows"] == table([*base, "--no-retain"])["rows"]
    for branch in ("spinon", "two-spinon"):
        single = records(table(["--delta", "0.5", "--branch", branch, "--momentum", "0", "--precision", precision]))
        assert len(single) == 1
        assert single[0]["energy" if branch == "spinon" else "lower"] == "0"
    failed_bulk = table([*base, "--max-evaluations", "0"], 2)
    assert "Bulk energy/site" not in failed_bulk["metadata"]
    assert failed_bulk["metadata"]["Bulk status"] == "evaluation_limit"
    assert failed_bulk["rows"] == t["rows"]  # analytic dispersions are independently valid
    underflow = table(["--delta", "1.00000001", "--points", "2", "--precision", precision], 2)
    assert "Single-spinon gap" not in underflow["metadata"]
    assert "Bulk energy/site" in underflow["metadata"]
    assert all(x["energy"] is None and x["lower"] is None and x["upper"] is None and x["status"] == "precision_limit"
               for x in records(underflow))

for args in (["--delta", "-1"], ["--delta", "nan"], ["--delta", "inf"], ["--delta", "2", "--exchange", "0"],
             ["--delta", "2", "--points", "1"], ["--delta", "2", "--points", "-1"],
             ["--delta", "2", "--momentum", "-1"], ["--delta", "2", "--momentum", "4"],
             ["--delta", "2", "--momentum", "7", "--branch", "two-spinon"],
             ["--delta", "2", "--momentum", "1", "--points", "4"],
             ["--delta", "2", "--folded", "--branch", "spinon"],
             ["--delta", "2", "--max-levels", "31"], ["--delta", "2", "--tolerance", "0"]):
    run(args, 1)
assert len(records(table(["--delta", "2", "--branch", "two-spinon", "--momentum", "4"]))) == 1

base = ["--delta", "2", "--points", "3"]
with tempfile.TemporaryDirectory() as directory:
    folder = Path(directory)
    js, cs, ts = [folder / ("curves."+suffix) for suffix in ("json", "csv", "tsv")]
    screen = run([*base, "--json", str(js), "--csv", str(cs), "--tsv", str(ts), "--format", "plain"])
    assert "CPU time" in screen and "Used for:" not in screen
    t = json.loads(js.read_text())["tables"]["dispersion"]
    expected = [{key: value if value is not None else "" for key, value in row.items()} for row in records(t)]
    for file, delimiter in ((cs, ","), (ts, "\t")):
        text = file.read_text()
        assert "#" in text and "Hamiltonian" in text and "Spinon sector" in text
        rows = list(csv.DictReader([line for line in text.splitlines() if not line.startswith("#")], delimiter=delimiter))
        assert rows == expected
    old = js.read_text()
    run([*base, "--json", str(js)], 1)
    assert js.read_text() == old
    run(["--delta", "-1", "--json", str(js), "--force"], 1)
    assert js.read_text() == old
    run([*base, "--json", str(js), "--force", "--quiet", "--no-retain"])
    assert json.loads(js.read_text())["tables"]["dispersion"]["rows"] == t["rows"]
    run([*base, "--stream", "--no-retain", "--format", "plain"])
    run([*base, "--format", "pretty"])
    failed = run(["--delta", "1.00000001", "--points", "2", "--format", "csv", "--no-preamble"], 2)
    assert all(x["energy"] == x["lower"] == x["upper"] == "" for x in csv.DictReader(failed.splitlines()))

print("XXZ thermodynamic CLI contracts passed")
