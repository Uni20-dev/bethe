"""XYZ thermodynamic native precision, momentum, failure and export contracts."""
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
    assert "S=sigma/2" in t["metadata"]["Hamiltonian"]
    return t


for precision in precisions:
    base = ["--eta", "0.75", "--t", "1", "--points", "3", "--precision", precision]
    t = table(base)
    rows = records(t)
    assert len(rows) == 18 and all(r["status"] == "converged" for r in rows)
    assert [r["branch"] for r in rows] == ["spinon"]*3 + ["two-spinon"]*3 + ["bound"]*12
    assert all(r["s"] is None and r["delta_rx"] is None and r["delta_rz"] is None for r in rows[:6])
    with localcontext() as ctx:
        ctx.prec = 90
        tol = {"fp64": Decimal("5e-14"), "long-double": Decimal("5e-17"), "fp128": Decimal("5e-31")}[precision]
        gap = Decimal("0.233183237167837415865617249661776631695819568738375906775664238")
        maximum = Decimal("0.502312502487677101434967797501115423981948185227108626838861918")
        b1 = Decimal("0.2437867651629989290195459944981206910241294413764639662718940275")
        b2 = Decimal("0.4098243055858140121240503172543673264480114334945959030163822254")
        assert abs(Decimal(rows[0]["energy"])-gap) < tol
        assert abs(Decimal(rows[1]["energy"])-maximum) < tol
        for r in rows[3:6]:
            assert r["energy"] is None
            assert abs(Decimal(r["lower"])-2*gap) < tol
            assert abs(Decimal(r["upper"])-2*maximum) < tol
        for offset, s, expected in ((6, 1, b1), (12, 2, b2)):
            copy = rows[offset:offset+6]
            assert all(r["s"] == s and r["delta_rx"] == (-1)**s for r in copy)
            assert [r["delta_rz"] for r in copy] == [1]*3+[-1]*3
            assert abs(Decimal(copy[0]["energy"])-expected) < tol
            assert copy[0]["energy"] == copy[2]["energy"] == copy[4]["energy"]
            assert copy[1]["energy"] == copy[3]["energy"] == copy[5]["energy"]
            assert all(r["cell_momentum"] == "0" for r in copy)
            assert copy[0]["reduced_q"] == copy[4]["reduced_q"] == "0"
        scaled = records(table([*base, "--exchange", "3"]))
        for r, s in zip(rows, scaled):
            for key in ("energy", "lower", "upper"):
                if r[key] is not None:
                    assert abs(Decimal(s[key])-3*Decimal(r[key])) < 10*tol
    assert rows == records(table([*base, "--no-retain"]))
    selected = records(table([*base, "--bound-state", "2"]))
    assert selected == rows[:6]+rows[12:]
    point = records(table(["--eta", "0.75", "--t", "1", "--branch", "bound", "--bound-state", "1",
                           "--momentum", "0", "--precision", precision]))
    assert len(point) == 2 and point[0] == rows[6] and point[1] == rows[9]
    repulsive = records(table(["--eta", "0.25", "--t", "1", "--points", "2", "--precision", precision]))
    assert len(repulsive) == 4 and all(r["branch"] != "bound" for r in repulsive)
    failed = table(["--eta", "0.25", "--t", "10000", "--points", "2", "--precision", precision], 2)
    assert "Single-spinon gap" not in failed["metadata"]
    assert all(r["energy"] is None and r["lower"] is None and r["upper"] is None and r["status"] == "precision_limit"
               for r in records(failed))

invalid = [["--eta", "0", "--t", "1"], ["--eta", "1", "--t", "1"], ["--eta", "nan", "--t", "1"],
           ["--eta", "0.75", "--t", "0"], ["--eta", "0.75", "--t", "inf"],
           ["--eta", "0.25", "--t", "1", "--branch", "bound"],
           ["--eta", "0.9999", "--t", "1"]]
base = ["--eta", "0.75", "--t", "1"]
invalid += [[*base, *flags] for flags in (["--bound-state", "0"], ["--bound-state", "3"],
              ["--bound-state", "-1"], ["--bound-state", "4294967296"], ["--exchange", "0"],
              ["--points", "1"], ["--points", "-1"], ["--points", "1000000"], ["--momentum", "-1"],
              ["--momentum", "4"], ["--branch", "bound", "--momentum", "7"],
              ["--momentum", "0", "--points", "2"], ["--branch", "spinon", "--bound-state", "1"])]
for args in invalid:
    run(args, 1)
assert len(records(table([*base, "--branch", "bound", "--momentum", "4"]))) == 4
assert len(records(table(["--eta", "0.9999", "--t", "1", "--bound-state", "1", "--points", "2"]))) == 8

with tempfile.TemporaryDirectory() as directory:
    paths = [Path(directory)/("curves."+suffix) for suffix in ("json", "csv", "tsv")]
    js, cs, ts = paths
    args = [*base, "--points", "3"]
    screen = run([*args, "--json", str(js), "--csv", str(cs), "--tsv", str(ts), "--format", "plain"])
    assert "CPU time" in screen and "Used for:" not in screen
    t = json.loads(js.read_text())["tables"]["dispersion"]
    expected = [{k: "" if v is None else str(v) for k, v in row.items()} for row in records(t)]
    for path, delimiter in ((cs, ","), (ts, "\t")):
        text = path.read_text()
        assert "Hamiltonian" in text and "Bound momentum" in text
        actual = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
        assert actual == expected
    original = js.read_text()
    run([*args, "--json", str(js)], 1)
    for bad in invalid:
        run([*bad, "--json", str(js), "--force"], 1)
        assert js.read_text() == original
    run([*args, "--json", str(js), "--force", "--quiet", "--no-retain"])
    assert json.loads(js.read_text())["tables"]["dispersion"]["rows"] == t["rows"]
    run([*args, "--stream", "--no-retain", "--format", "plain"])
    run([*args, "--format", "pretty"])
    failure = run(["--eta", "0.25", "--t", "10000", "--points", "2", "--format", "csv", "--no-preamble"], 2)
    assert all(r["energy"] == r["lower"] == r["upper"] == "" for r in csv.DictReader(failure.splitlines()))

print("XYZ thermodynamic CLI contracts passed")
