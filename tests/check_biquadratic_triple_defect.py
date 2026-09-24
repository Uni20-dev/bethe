"""Triple-plus-defect selections and scans through the shared cluster CLI."""
from decimal import Decimal, localcontext
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
        assert "NOT the full four-defect" in table["metadata"]["Coverage"]
        assert "not a global excitation rank" in table["metadata"]["Ordering"]
    return doc["tables"]


def rows(t):
    return [dict(zip([c["id"] for c in t["columns"]], row)) for row in t["rows"]]


def labels(t):
    return [(int(x["real_label"]), int(x["string_label"])) for x in rows(t["labels"])]


for precision in precisions:
    for n in (8, 9, 10, 16):
        base = [str(n), "--ferromagnetic", "--precision", precision]
        full = tables([*base, "--triple-defects", "all", "--roots"])
        s = rows(full["states"])
        assert len(s) == (n-3)*(n-7)
        assert set(labels(full)) == set(itertools.product(range(1, n-2), range(1, n-6)))
        assert all(x["converged"] and int(x["through_lines"]) == n-8 for x in s)
        gaps = [Decimal(x["gap"]) for x in s]
        assert gaps == sorted(gaps)
        assert len(rows(full["roots"])) == 4*len(s)
        assert len(rows(full["string"])) == len(s)
        # Selected solves must reproduce the scan, including the finite deviation.
        for idx in (0, len(s)//2, len(s)-1):
            i, j = labels(full)[idx]
            selected = tables([*base, "--triple-defect", f"{i},{j}"])
            for key in ("gap", "energy", "reference_energy", "residual"):
                assert rows(selected["states"])[0][key] == s[idx][key]
            for key in ("center", "log_deviation", "deviation_phase"):
                assert rows(selected["string"])[0][key] == rows(full["string"])[idx][key]
        with localcontext() as ctx:
            ctx.prec = 70
            tol = Decimal("1e-29" if precision == "fp128" else "1e-12")
            for x in s:
                assert abs(Decimal(x["energy"])-Decimal(x["gap"])-(n-1)) < tol
                assert Decimal(x["gap"]) == Decimal(x["tl_energy"])
        width = min(3, n-7)
        args = [*base, "--triple-defects", "all", "--mixed-window", str(width), "--roots"]
        window = tables(args)
        assert len(rows(window["states"])) == width*width
        assert set(labels(window)) == set(itertools.product(range(n-2-width, n-2), range(n-6-width, n-6)))
        lookup = dict(zip(labels(full), [x["gap"] for x in s]))
        for key, x in zip(labels(window), rows(window["states"])):
            assert x["gap"] == lookup[key]
        first = tables([*base, "--triple-defects", "1", "--mixed-window", str(width)])
        assert rows(first["states"])[0] == rows(window["states"])[0]
        assert first["states"]["metadata"]["Scanned candidates"] == str(width*width)
        streamed = tables([*args, "--no-retain"])
        assert all(streamed[name]["rows"] == t["rows"] for name, t in window.items())
        roots = rows(window["roots"])
        for idx, label in enumerate(rows(window["labels"])):
            central, plus, minus, real = roots[4*idx:4*idx+4]
            assert central["state_id"] == plus["state_id"] == minus["state_id"] == real["state_id"] == str(idx)
            assert Decimal(central["u_real"]) == Decimal(real["u_real"]) == 0
            assert plus["u_real"] == minus["u_real"]
            assert Decimal(plus["u_imag"]) == Decimal(minus["u_imag"]).copy_negate()
            with localcontext() as ctx:
                ctx.prec = 70
                assert abs(2*Decimal(real["u_imag"])-Decimal(label["real_rapidity"])) < Decimal("1e-12")
    large = tables(["100000", "--ferromagnetic", "--triple-defects", "all", "--mixed-window", "2",
                    "--precision", precision, "--roots"])
    assert all(x["deviation_real"] is None and x["deviation_imag"] is None for x in rows(large["string"]))
    assert all(x["multiplicity"] is None for x in rows(large["states"]))
    assert abs(Decimal(rows(large["states"])[0]["gap"])-3) < Decimal("1e-7")
    for budget in (0, 1):
        failed = tables(["16", "--ferromagnetic", "--triple-defects", "2", "--mixed-window", "2",
                         "--precision", precision, "--max-iterations", str(budget), "--roots"], 2)
        s = rows(failed["states"])
        assert len(s) == 1 and s[0]["gap"] is None and not s[0]["converged"]
        assert int(s[0]["iterations"]) == budget
        assert failed["states"]["metadata"]["Failed candidates"] == "4"
        assert len(rows(failed["roots"])) == 4
        assert len(rows(failed["labels"])) == len(rows(failed["string"])) == 1

base = ["16", "--ferromagnetic"]
partial = tables([*base, "--triple-defects", "1", "--max-iterations", "4", "--roots"], 2)
metadata = partial["states"]["metadata"]
assert 0 < int(metadata["Converged candidates"]) < 117
assert int(metadata["Converged candidates"]) + int(metadata["Failed candidates"]) == 117
s = rows(partial["states"])
assert len(s) == 2 and s[0]["converged"] and not s[1]["converged"]
assert s[0]["gap"] is not None and s[1]["gap"] is None
assert len(rows(partial["roots"])) == 8 and len(rows(partial["labels"])) == 2
assert rows(partial["reference"])[0]["state_id"] == "2"
invalid = [["7", "--ferromagnetic", "--triple-defects", "all"], ["8", "--triple-defects", "all"],
           [*base, "--triple-defects", "0"], [*base, "--triple-defects", "all", "--mixed-window", "0"],
           [*base, "--triple-defects", "all", "--mixed-window", "10"],
           [*base, "--triple-defects", "1", "--mixed-window", "3", "--max-candidates", "8"],
           ["1000000000", "--ferromagnetic", "--triple-defects", "1"]]
for text in ("1", "1,2,3", "0,1", "1,0", "14,1", "1,10", "1,x", "1,", "-1,1"):
    invalid.append([*base, "--triple-defect", text])
for option in (["--triple-defects", "all"], ["--pair-defects", "all"], ["--pair-defect", "1,1"],
               ["--bound-pairs", "1"], ["--bound-triples", "1"], ["--two-pairs", "9,10"],
               ["--two-pair-states", "1"], ["--real-defects", "1"], ["--mixed-window", "2"],
               ["--pair-window", "2"], ["--excitations", "1"], ["--one-defect"], ["--q-spectrum"],
               ["--through-lines", "8"], ["--max-candidates", "1"], ["--sectors"]):
    invalid.append([*base, "--triple-defect", "13,9", *option])
for args in invalid:
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=10, env=env)
    assert p.returncode != 0, (args, p.stdout)

with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp)
    args = [*base, "--triple-defects", "3", "--mixed-window", "2", "--roots"]
    expected = tables(args)
    run([*args, "--quiet", "--no-retain", "--json", str(path/"out.json"),
         "--csv", str(path/"states.csv"), "--tsv-table", f"string={path/'strings.tsv'}"])
    doc = json.loads((path/"out.json").read_text())
    assert all(doc["tables"][name]["rows"] == t["rows"] for name, t in expected.items())
    for name, delimiter in (("states.csv", ","), ("strings.tsv", "\t")):
        lines = [line for line in (path/name).read_text().splitlines() if not line.startswith("#")]
        assert len(list(csv.DictReader(lines, delimiter=delimiter))) == 3
