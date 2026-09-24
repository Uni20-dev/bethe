"""Two-pair CLI: exhaustive labels, bounded scans, tables and failure contracts."""
from decimal import Decimal
import csv
import itertools
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
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=30, env=env)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def tables(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete"
    for table in doc["tables"].values():
        assert table["summary"]["Outcome"] == ("success" if status == 0 else "partial")
        assert "NOT the full" in table["metadata"]["Coverage"]
    return doc["tables"]


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]


def labels(t):
    result = {}
    for s in rows(t["string"]):
        result.setdefault(s["state_id"], []).append(int(s["string_label"]))
    return {sid: tuple(js) for sid, js in result.items()}


for precision in precisions:
    for n in (8, 9, 10, 16):
        base = [str(n), "--ferromagnetic", "--precision", precision]
        full = tables([*base, "--two-pair-states", "all", "--roots"])
        s = rows(full["states"])
        keys = labels(full)
        assert set(keys.values()) == set(itertools.combinations(range(1, n-5), 2))
        assert len(s) == (n-6)*(n-7)//2
        gaps = [Decimal(x["gap"]) for x in s]
        assert gaps == sorted(gaps) and all(x["converged"] for x in s)
        assert len(rows(full["roots"])) == 4*len(s)
        assert all(int(x["through_lines"]) == n-8 for x in s)
        for x in rows(full["string"]):
            expected = -1 if (n-int(x["string_label"])-int(x["pair_index"])-1) % 2 else 1
            assert x["deviation_sign"] == expected
        for x in s:
            js = keys[x["state_id"]]
            selected = tables([*base, "--two-pairs", f"{js[0]},{js[1]}"])
            assert rows(selected["states"])[0]["gap"] == x["gap"]
            assert abs(Decimal(x["energy"])-Decimal(x["gap"])-(n-1)) < Decimal("1e-12")
        width = min(3, n-6)
        window_args = [*base, "--two-pair-states", "all", "--pair-window", str(width), "--roots"]
        window = tables(window_args)
        assert set(labels(window).values()) == set(itertools.combinations(range(n-5-width, n-5), 2))
        first = tables([*base, "--two-pair-states", "1", "--pair-window", str(width)])
        assert rows(first["states"])[0] == rows(window["states"])[0]
        assert first["states"]["metadata"]["Scanned candidates"] == str(width*(width-1)//2)
        streamed = tables([*window_args, "--no-retain"])
        assert all(streamed[name]["rows"] == t["rows"] for name, t in window.items())
        r = rows(window["roots"])
        for i in range(0, len(r), 2):
            assert r[i]["u_real"] == r[i+1]["u_real"]
            assert Decimal(r[i]["u_imag"]) == Decimal(r[i+1]["u_imag"]).copy_negate()
    large = tables(["100000", "--ferromagnetic", "--two-pair-states", "all", "--pair-window", "2",
                    "--precision", precision])
    assert all(x["deviation"] is None for x in rows(large["string"]))
    assert rows(large["states"])[0]["multiplicity"] is None
    assert abs(Decimal(rows(large["states"])[0]["gap"])-Decimal(10)/3) < Decimal("1e-7")
    for budget in (0, 1):
        failed = tables(["16", "--ferromagnetic", "--two-pair-states", "1", "--pair-window", "3",
                         "--max-iterations", str(budget), "--precision", precision, "--roots"], 2)
        s = rows(failed["states"])
        assert len(s) == 1 and s[0]["gap"] is None and not s[0]["converged"]
        assert int(s[0]["iterations"]) == budget
        assert failed["states"]["metadata"]["Failed candidates"] == "3"
        assert len(rows(failed["string"])) == 2 and len(rows(failed["roots"])) == 4

base = ["16", "--ferromagnetic"]
partial = tables([*base, "--two-pair-states", "1", "--max-iterations", "4", "--roots"], 2)
metadata = partial["states"]["metadata"]
assert 0 < int(metadata["Converged candidates"]) < 45
assert int(metadata["Converged candidates"]) + int(metadata["Failed candidates"]) == 45
s = rows(partial["states"])
assert len(s) == 2 and s[0]["converged"] and not s[1]["converged"]
assert s[0]["gap"] is not None and s[1]["gap"] is None
assert len(rows(partial["string"])) == 4 and len(rows(partial["roots"])) == 8
assert rows(partial["reference"])[0]["state_id"] == "2"
invalid = [
    ["7", "--ferromagnetic", "--two-pair-states", "all"],
    ["8", "--two-pair-states", "all"],
    [*base, "--two-pair-states", "0"],
    [*base, "--pair-window", "2"],
    [*base, "--two-pair-states", "all", "--pair-window", "1"],
    [*base, "--two-pair-states", "all", "--pair-window", "11"],
    [*base, "--two-pair-states", "1", "--max-candidates", "44"],
    ["1000000000", "--ferromagnetic", "--two-pair-states", "1"],
]
for text in ("1", "1,2,3", "2,1", "1,1", "0,2", "1,11", "x,2", "1,", "-1,2"):
    invalid.append([*base, "--two-pairs", text])
for option in (["--two-pair-states", "1"], ["--bound-pairs", "1"], ["--bound-triples", "1"],
               ["--pair-defects", "all"], ["--pair-defect", "1,1"], ["--through-lines", "8"],
               ["--max-candidates", "1"], ["--real-defects", "1"], ["--mixed-window", "2"],
               ["--excitations", "1"], ["--one-defect"], ["--q-spectrum"], ["--sectors"]):
    invalid.append([*base, "--two-pairs", "9,10", *option])
for args in invalid:
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=10, env=env)
    assert p.returncode != 0, (args, p.stdout)

with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp)
    args = [*base, "--two-pair-states", "all", "--pair-window", "3", "--roots"]
    expected = tables(args)
    run([*args, "--quiet", "--no-retain", "--json", str(path/"out.json"),
         "--csv", str(path/"states.csv"), "--tsv-table", f"string={path/'strings.tsv'}"])
    doc = json.loads((path/"out.json").read_text())
    assert all(doc["tables"][name]["rows"] == t["rows"] for name, t in expected.items())
    for name, delimiter, count in (("states.csv", ",", 3), ("strings.tsv", "\t", 6)):
        lines = [line for line in (path/name).read_text().splitlines() if not line.startswith("#")]
        assert len(list(csv.DictReader(lines, delimiter=delimiter))) == count
