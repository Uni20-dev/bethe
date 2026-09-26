"""Exact two-hole NLIE frontend: level/gap distinction and output contracts."""
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
base = ["--length", "1", "--p", "1"]


def run(args, status=0):
    result = subprocess.run([program, *args], text=True, capture_output=True, timeout=300, env=env)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    return result.stdout


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]], r)) for r in table["rows"]]


def tables(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete" and list(doc["tables"]) == ["levels", "source", "gap"]
    for table in doc["tables"].values():
        assert table["summary"]["Outcome"] == ("success" if status == 0 else "partial")
        assert Decimal(table["summary"]["Run CPU seconds"]) >= 0
        assert "bulk-subtracted" in table["metadata"]["Energy convention"]
        assert "exact continuum" in table["metadata"]["Calculation"]
        assert table["metadata"]["Dimensionless observables"] == "Y=L*E_C; scaled_gap=L*gap/(2*pi)"
        assert int(table["metadata"]["Max total nonlinear updates per state"]) == 32768
    return doc["tables"]


for precision in precisions:
    with localcontext() as ctx:
        ctx.prec = 100
        pi = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406286")
        tol = {"fp64": Decimal("1e-8"), "long-double": Decimal("1e-11"), "fp128": Decimal("1e-25")}[precision]
        for number in ("0.5", "1.5"):
            t = tables([*base, "--precision", precision, "--number", number, "--charge", "-2"])
            source = rows(t["source"])[0]
            gap = rows(t["gap"])[0]
            levels = rows(t["levels"])
            reference = 2 * (1 + (2*pi*Decimal(number))**2).sqrt()
            assert abs(Decimal(gap["gap"]) - reference) < tol
            assert abs(Decimal(gap["scaled_gap"]) - reference/(2*pi)) < tol
            assert source["charge"] == -2 and Decimal(source["number"]) == Decimal(number)
            assert levels[0]["state"] == "vacuum" and levels[1]["state"] == "two_soliton"
            assert abs(Decimal(levels[1]["casimir_energy"])-Decimal(levels[0]["casimir_energy"])-reference) < tol
            assert int(source["source_evaluations"]) == 0
            assert Decimal(t["source"]["metadata"]["UV weight Delta+=Delta-"]) == Decimal(number)

# Distinguish the exact interacting level from both BY and the vacuum-relative gap.
t = tables(["--p", "2", "--length", "1", "--tolerance", "1e-7"])
levels, gap = rows(t["levels"]), rows(t["gap"])[0]
assert abs(float(levels[1]["casimir_energy"])-5.11055662981005) < 1e-7
assert abs(float(gap["gap"])-5.44401859555224) < 1e-7

for flags, expected in ((["--max-root-iterations", "0"], "iteration_limit"),
                         (["--max-intervals", "64"], "mesh_limit"),
                         (["--max-cutoffs", "1"], "cutoff_limit")):
    t = tables([*base, *flags], 2)
    assert rows(t["source"])[0]["status"] == expected
    assert rows(t["source"])[0]["rapidity"] is None
    assert rows(t["levels"])[1]["casimir_energy"] is None
    assert rows(t["gap"])[0]["gap"] is None
    if flags[0] == "--max-root-iterations":
        assert rows(t["levels"])[0]["converged"]  # successful vacuum cannot rescue a failed excited level

invalid = [["--length", "1"], ["--p", "1"], [*base, "--mass", "0"],
           ["--length", "1", "--p", "0.5"], [*base, "--number", "1"], [*base, "--number", "2.5"],
           [*base, "--charge", "0"], [*base, "--tolerance", "0"], [*base, "--max-root-iterations", "-1"],
           [*base, "--initial-intervals", "9"], [*base, "--contour-shift", "2"]]
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)
    flags = ["--json", str(path/"result.json"), "--csv-table", "gap="+str(path/"gap.csv"),
             "--tsv-table", "source="+str(path/"source.tsv")]
    screen = run([*base, *flags, "--format", "plain"])
    assert "CPU time" in screen and "Used for:" not in screen
    expected = json.loads((path/"result.json").read_text())["tables"]
    for table, suffix, delimiter in (("gap", "csv", ","), ("source", "tsv", "\t")):
        text = (path/f"{table}.{suffix}").read_text()
        assert "bulk-subtracted" in text
        actual = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
        normalized = [{k: "" if v is None else str(v).lower() if isinstance(v, bool) else str(v)
                       for k, v in row.items()} for row in rows(expected[table])]
        assert actual == normalized
    old = (path/"result.json").read_text()
    run([*base, "--json", str(path/"result.json")], 1)
    for args in invalid:
        run([*args, "--json", str(path/"result.json"), "--force"], 1)
        assert (path/"result.json").read_text() == old
    run([*base, *flags, "--force", "--no-retain", "--quiet"])
    actual = json.loads((path/"result.json").read_text())["tables"]
    assert all(actual[k]["rows"] == expected[k]["rows"] for k in expected)
    run([*base, "--stream", "--no-retain", "--format", "plain"])
    run([*base, "--format", "pretty"])

print("Sine-Gordon exact excitation CLI contracts passed")
