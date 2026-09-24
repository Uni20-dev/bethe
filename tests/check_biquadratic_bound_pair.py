"""Targeted signed two-string modes, native precision and output contracts."""
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
    assert doc["status"] == "complete"  # Transport completion, not scientific outcome.
    for table in doc["tables"].values():
        assert table["summary"]["Outcome"] == ("success" if status == 0 else "partial")
        assert "NOT the full two-defect" in table["metadata"]["Coverage"]
        assert "not a global excitation rank" in table["metadata"]["Ordering"]
        assert "exact degenerate ferro" in table["metadata"]["Gap reference"]
        assert Decimal(table["summary"]["Run CPU seconds"]) >= 0
    return doc["tables"]


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]


for precision in precisions:
    base = ["4", "--ferromagnetic", "--bound-pairs", "all", "--precision", precision]
    t = tables([*base, "--roots"])
    s = rows(t["states"])[0]
    with localcontext() as ctx:
        ctx.prec = 70
        tol = Decimal("1e-29" if precision == "fp128" else "1e-12")
        exact = (15 - Decimal(17).sqrt()) / 2
        assert abs(Decimal(s["energy"]) - exact) < tol
        assert abs(Decimal(s["gap"]) - exact + 3) < tol
        assert abs(Decimal(s["tl_energy"]) - Decimal(s["gap"])) < tol
        assert abs(Decimal(s["energy"]) + 2 * Decimal(s["reference_energy"]) - Decimal(21) / 4) < tol
    assert s["through_lines"] == "0" and s["multiplicity"] == "1" and s["mode"] == "1"
    assert s["converged"] and len(rows(t["roots"])) == 2
    assert rows(t["reference"])[0]["multiplicity"] == "55"
    assert int(rows(t["string"])[0]["deviation_sign"]) == 1
    for n in (5, 8, 16):
        args = [str(n), "--ferromagnetic", "--bound-pairs", "all", "--precision", precision, "--roots"]
        full = tables(args)
        states = rows(full["states"])
        assert len(states) == n - 3
        assert [int(s["mode"]) for s in states] == list(range(1, n - 2))
        assert all(s["converged"] and int(s["through_lines"]) == n - 4 for s in states)
        energies = [Decimal(s["energy"]) for s in states]
        assert energies == sorted(energies) and len(set(energies)) == len(energies)
        for i, string in enumerate(rows(full["string"])):
            assert int(string["state_id"]) == i and int(string["string_label"]) == n - 3 - i
            assert int(string["deviation_sign"]) == (1 if i % 2 == 0 else -1)
            assert Decimal(string["deviation"]) * int(string["deviation_sign"]) > 0
        roots = rows(full["roots"])
        assert len(roots) == 2 * (n - 3)
        for i in range(n - 3):
            assert roots[2*i]["state_id"] == roots[2*i+1]["state_id"] == str(i)
            assert Decimal(roots[2*i]["u_imag"]) == Decimal(roots[2*i+1]["u_imag"]).copy_negate()
        first = tables([str(n), "--ferromagnetic", "--bound-pairs", "1", "--precision", precision])
        assert rows(first["states"])[0] == states[0]
        streamed = tables([*args, "--no-retain"])
        assert all(streamed[name]["rows"] == table["rows"] for name, table in full.items())
    large = tables(["100000", "--ferromagnetic", "--bound-pairs", "2", "--precision", precision, "--roots"])
    assert all(s["multiplicity"] is None and s["converged"] for s in rows(large["states"]))
    assert all(s["deviation"] is None and float(s["log_deviation"]) > 80000 for s in rows(large["string"]))
    assert abs(float(rows(large["states"])[0]["gap"]) - 5/3) < 4e-10
    failed = tables(["8", "--ferromagnetic", "--bound-pairs", "2", "--max-iterations", "0",
                     "--precision", precision, "--roots"], 2)
    assert all(not s["converged"] and s["gap"] is None and s["iterations"] == "0" for s in rows(failed["states"]))
    assert rows(failed["reference"])[0]["energy"] == "7"

for args in (["8", "--bound-pairs", "1"], ["3", "--ferromagnetic", "--bound-pairs", "1"],
             ["8", "--ferromagnetic", "--bound-pairs", "0"],
             ["8", "--ferromagnetic", "--bound-pairs", "all", "--max-candidates", "4"],
             ["8", "--ferromagnetic", "--bound-pairs", "1", "--max-candidates", "0"],
             ["8", "--ferromagnetic", "--bound-pairs", "1", "--tolerance", "nan"]):
    run(args, 1)
for selection in (["--one-defect"], ["--sectors"], ["--singlet-excitation"], ["--through-lines", "4"],
                  ["--excitations", "all"], ["--quantum-numbers", "1,2"], ["--q-spectrum"],
                  ["--q-seed", "1,2"], ["--max-attempts", "5"]):
    run(["8", "--ferromagnetic", "--bound-pairs", "1", *selection], 1)
assert len(rows(tables(["8", "--ferromagnetic", "--bound-pairs", "100"])["states"])) == 5
assert len(rows(tables(["100000", "--ferromagnetic", "--bound-pairs", "1", "--max-candidates", "1"])["states"])) == 1

with tempfile.TemporaryDirectory() as tmp:
    base = Path(tmp)
    args = ["8", "--ferromagnetic", "--bound-pairs", "2", "--roots"]
    exports = ["--json", str(base / "pairs.json")]
    for extension in ("csv", "tsv"):
        for name in ("states", "reference", "string", "roots"):
            exports += [f"--{extension}-table", f"{name}={base / f'pairs.{name}.{extension}'}"]
    run([*args, *exports])
    doc = json.loads((base / "pairs.json").read_text())
    assert len(rows(doc["tables"]["states"])) == 2
    for extension, delimiter in (("csv", ","), ("tsv", "\t")):
        for name, count in (("states", 2), ("reference", 1), ("string", 2), ("roots", 4)):
            path = base / f"pairs.{name}.{extension}"
            data = list(csv.DictReader((line for line in path.read_text().splitlines() if not line.startswith("#")), delimiter=delimiter))
            assert len(data) == count, (path, data)
assert "--bound-pairs" in run(["--help"])
assert "nachtergaele-spitzer-starr-2007" in run(["--references"])
