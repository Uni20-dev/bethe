"""Mixed-family selection/scans through the shared cluster CLI and tables."""
from decimal import Decimal, localcontext
import csv
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
    for t in doc["tables"].values():
        assert t["summary"]["Outcome"] == ("success" if status == 0 else "partial")
        assert "NOT the full three-defect" in t["metadata"]["Coverage"]
        assert "not a global excitation rank" in t["metadata"]["Ordering"]
    return doc["tables"]


def rows(t):
    return [dict(zip([c["id"] for c in t["columns"]], row)) for row in t["rows"]]


for precision in precisions:
    base = ["6", "--ferromagnetic", "--precision", precision]
    t = tables([*base, "--pair-defects", "all", "--roots"])
    s = rows(t["states"])
    assert len(s) == 3 and all(x["converged"] for x in s)
    assert [int(x["real_label"]) for x in rows(t["labels"])] == [3, 2, 1]
    assert all(x["string_label"] == "1" for x in rows(t["labels"]))
    assert len(rows(t["roots"])) == 9
    with localcontext() as ctx:
        ctx.prec = 70
        tol = Decimal("1e-29" if precision == "fp128" else "1e-12")
        # Independent six-site module characteristic polynomial.
        g = Decimal(s[0]["gap"])
        assert abs(((g-17)*g+80)*g-106) < 100 * tol
        assert abs(Decimal(s[1]["gap"])-6) < tol
        assert abs(Decimal(s[2]["gap"])-7) < tol
        for x in s:
            assert abs(Decimal(x["energy"]) - Decimal(x["gap"]) - 5) < tol
            assert Decimal(x["gap"]) == Decimal(x["tl_energy"])
    selected = tables([*base, "--pair-defect", "3,1", "--roots"])
    assert rows(selected["states"])[0] == s[0]
    assert rows(selected["labels"])[0] == rows(t["labels"])[0]
    for n in (7, 8, 16):
        scan = [str(n), "--ferromagnetic", "--precision", precision, "--pair-defects", "all"]
        full = tables(scan)
        assert len(rows(full["states"])) == (n-3)*(n-5)
        energies = [Decimal(x["gap"]) for x in rows(full["states"])]
        assert energies == sorted(energies)
        window = tables([*scan, "--mixed-window", "2", "--roots"])
        labels = rows(window["labels"])
        assert {(int(x["real_label"]), int(x["string_label"])) for x in labels} == {
            (i, j) for i in (n-4, n-3) for j in (n-6, n-5)}
        lookup = {(x["real_label"], x["string_label"]): state["gap"]
                  for x, state in zip(rows(full["labels"]), rows(full["states"]))}
        for label, state in zip(labels, rows(window["states"])):
            assert state["gap"] == lookup[label["real_label"], label["string_label"]]
        first = tables([str(n), "--ferromagnetic", "--precision", precision,
                        "--pair-defects", "1", "--mixed-window", "2"])
        assert rows(first["states"])[0] == rows(window["states"])[0]
        assert first["states"]["metadata"]["Scanned candidates"] == "4"
        streamed = tables([*scan, "--mixed-window", "2", "--roots", "--no-retain"])
        assert all(streamed[name]["rows"] == value["rows"] for name, value in window.items())
        roots = rows(window["roots"])
        for idx in range(4):
            real, plus, minus = roots[3*idx:3*idx+3]
            assert real["state_id"] == plus["state_id"] == minus["state_id"] == str(idx)
            assert Decimal(real["u_real"]) == 0
            assert plus["u_real"] == minus["u_real"]
            assert Decimal(plus["u_imag"]) == Decimal(minus["u_imag"]).copy_negate()
    large = tables(["100000", "--ferromagnetic", "--pair-defects", "all", "--mixed-window", "2",
                    "--precision", precision, "--roots"])
    assert len(rows(large["states"])) == 4
    assert all(s["multiplicity"] is None for s in rows(large["states"]))
    assert all(s["deviation"] is None for s in rows(large["string"]))
    assert abs(Decimal(rows(large["states"])[0]["gap"]) - Decimal(8)/3) < Decimal("1e-7")
    for budget in (0, 1):
        failed = tables(["16", "--ferromagnetic", "--pair-defects", "2", "--mixed-window", "2",
                         "--precision", precision, "--max-iterations", str(budget), "--roots"], 2)
        assert failed["states"]["metadata"]["Failed candidates"] == "4"
        assert failed["states"]["metadata"]["Retained converged levels"] == "0"
        s = rows(failed["states"])
        assert len(s) == 1 and not s[0]["converged"] and s[0]["gap"] is None
        assert int(s[0]["iterations"]) == budget
        assert len(rows(failed["labels"])) == 1 and len(rows(failed["roots"])) == 3
        assert rows(failed["reference"])[0]["state_id"] == "1"

partial = tables(["16", "--ferromagnetic", "--pair-defects", "2", "--mixed-window", "4",
                  "--max-iterations", "5", "--roots"], 2)
meta = partial["states"]["metadata"]
assert 0 < int(meta["Converged candidates"]) < 16
assert int(meta["Failed candidates"]) + int(meta["Converged candidates"]) == 16
partial_rows = rows(partial["states"])
assert len(partial_rows) == 3 and all(s["converged"] for s in partial_rows[:2])
assert not partial_rows[-1]["converged"] and partial_rows[-1]["gap"] is None
assert rows(partial["reference"])[0]["state_id"] == "3"

rounded = rows(tables(["100000000", "--ferromagnetic", "--pair-defects", "all", "--mixed-window", "4"])["states"])
assert len({r["energy"] for r in rounded}) == 1
gaps = [Decimal(r["gap"]) for r in rounded]
assert len(set(gaps)) > 1 and gaps == sorted(gaps)

for extra in (["--pair-defects", "0"], ["--pair-defects", "all", "--mixed-window", "0"],
              ["--pair-defects", "all", "--mixed-window", "4"], ["--mixed-window", "2"],
              ["--pair-defect", "1"], ["--pair-defect", "1,2,3"], ["--pair-defect", "-1,2"],
              ["--pair-defect", "1/2,1"], ["--pair-defect", "0,1"], ["--pair-defect", "6,1"],
              ["--pair-defect", "1,4"], ["--pair-defect", "1,1", "--max-candidates", "2"],
              ["--pair-defects", "1", "--max-candidates", "14"],
              ["--pair-defects", "all", "--mixed-window", "2", "--max-candidates", "3"]):
    run(["8", "--ferromagnetic", *extra], 1)
for option in ("--pair-defect", "--pair-defects"):
    selector = [option, "1,1" if option == "--pair-defect" else "all"]
    run(["8", *selector], 1)
    run(["5", "--ferromagnetic", *selector], 1)
    for conflict in (["--bound-pairs", "1"], ["--bound-triples", "1"], ["--excitations", "1"],
                     ["--through-lines", "2"], ["--one-defect"], ["--sectors"], ["--q-spectrum"],
                     ["--quantum-numbers", "1,2,3"], ["--singlet-excitation"]):
        run(["8", "--ferromagnetic", *selector, *conflict], 1)
run(["8", "--ferromagnetic", "--pair-defect", "1,1", "--pair-defects", "all"], 1)
run(["100000", "--ferromagnetic", "--pair-defects", "1"], 1)  # Budget covers scan, not retained count.

with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp)
    args = ["8", "--ferromagnetic", "--pair-defects", "all", "--mixed-window", "2", "--roots"]
    run([*args, "--json", str(path/"mixed.json"), "--csv", str(path/"states.csv"),
         "--tsv-table", f"labels={path/'labels.tsv'}", "--csv-table", f"string={path/'string.csv'}"])
    doc = json.loads((path/"mixed.json").read_text())
    assert len(doc["tables"]["labels"]["rows"]) == 4
    for name, sep in (("states.csv", ","), ("labels.tsv", "\t"), ("string.csv", ",")):
        content = (path/name).read_text().splitlines()
        assert any(line.startswith("#") for line in content)
        data = list(csv.DictReader([line for line in content if not line.startswith("#")], delimiter=sep))
        assert len(data) == 4
    # Invalid labels fail before creating output files.
    run(["8", "--ferromagnetic", "--pair-defect", "6,1", "--json", str(path/"invalid.json")], 1)
    assert not (path/"invalid.json").exists()

assert "--pair-defect" in run(["--help"])
