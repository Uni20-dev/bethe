"""Shared two-/three-string CLI, native precision and output contracts."""
from decimal import Decimal, localcontext
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, fp128, defects = sys.argv[1:]
defects = int(defects)
assert defects in (2, 3)
option = "--bound-pairs" if defects == 2 else "--bound-triples"
other_option = "--bound-triples" if defects == 2 else "--bound-pairs"
family = "two" if defects == 2 else "three"
minimum_n = 2 * defects
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
        assert f"NOT the full {family}-defect" in table["metadata"]["Coverage"]
        assert "not a global excitation rank" in table["metadata"]["Ordering"]
        assert "exact degenerate ferro" in table["metadata"]["Gap reference"]
        assert Decimal(table["summary"]["Run CPU seconds"]) >= 0
    return doc["tables"]


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]


for precision in precisions:
    base = [str(minimum_n), "--ferromagnetic", option, "all", "--precision", precision]
    t = tables([*base, "--roots"])
    s = rows(t["states"])[0]
    with localcontext() as ctx:
        ctx.prec = 70
        tol = Decimal("1e-29" if precision == "fp128" else "1e-12")
        if defects == 2:
            exact = (15 - Decimal(17).sqrt()) / 2
        else:
            gap = Decimal("2.3")
            for _ in range(15):
                gap -= (((gap-17)*gap+80)*gap-106) / ((3*gap-34)*gap+80)
            exact = 5 + gap
        assert abs(Decimal(s["energy"]) - exact) < tol
        assert abs(Decimal(s["gap"]) - exact + minimum_n - 1) < tol
        assert abs(Decimal(s["tl_energy"]) - Decimal(s["gap"])) < tol
        assert abs(Decimal(s["energy"]) + 2 * Decimal(s["reference_energy"]) - Decimal(7*(minimum_n-1)) / 4) < tol
    assert s["through_lines"] == "0" and s["multiplicity"] == "1" and s["mode"] == "1"
    assert s["converged"] and len(rows(t["roots"])) == defects
    assert rows(t["reference"])[0]["multiplicity"] == ("55" if defects == 2 else "377")
    for n in (minimum_n + 1, 8, 16, 65):
        count = n - (minimum_n - 1)
        args = [str(n), "--ferromagnetic", option, "all", "--precision", precision, "--roots"]
        full = tables(args)
        states = rows(full["states"])
        assert len(states) == count
        assert [int(s["mode"]) for s in states] == list(range(1, count + 1))
        assert all(s["converged"] and int(s["through_lines"]) == n - minimum_n for s in states)
        energies = [Decimal(s["energy"]) for s in states]
        assert energies == sorted(energies) and len(set(energies)) == len(energies)
        for i, string in enumerate(rows(full["string"])):
            assert int(string["state_id"]) == i and int(string["string_label"]) == count - i
            if defects == 2:
                assert int(string["deviation_sign"]) == (1 if i % 2 == 0 else -1)
                assert Decimal(string["deviation"]) * int(string["deviation_sign"]) > 0
            else:
                assert -3.142 < float(string["deviation_phase"]) < 3.142
                assert float(string["log_deviation"]) > 0
                assert string["deviation_real"] is not None and string["deviation_imag"] is not None
        roots = rows(full["roots"])
        assert len(roots) == defects * count
        for i in range(count):
            cluster = roots[defects*i:defects*(i+1)]
            assert all(r["state_id"] == str(i) for r in cluster)
            assert [int(r["index"]) for r in cluster] == list(range(defects))
            assert Decimal(cluster[-2]["u_imag"]) == Decimal(cluster[-1]["u_imag"]).copy_negate()
            assert cluster[-2]["u_real"] == cluster[-1]["u_real"]
            if defects == 3:
                assert Decimal(cluster[0]["u_real"]) == 0
        first = tables([str(n), "--ferromagnetic", option, "1", "--precision", precision])
        assert rows(first["states"])[0] == states[0]
        streamed = tables([*args, "--no-retain"])
        assert all(streamed[name]["rows"] == table["rows"] for name, table in full.items())
    large = tables(["100000", "--ferromagnetic", option, "2", "--precision", precision, "--roots"])
    assert all(s["multiplicity"] is None and s["converged"] for s in rows(large["states"]))
    components = ("deviation",) if defects == 2 else ("deviation_real", "deviation_imag")
    assert all(all(s[key] is None for key in components) and float(s["log_deviation"]) > 80000
               for s in rows(large["string"]))
    assert abs(float(rows(large["states"])[0]["gap"]) - (5/3 if defects == 2 else 2)) < 4e-10
    failed = tables(["8", "--ferromagnetic", option, "2", "--max-iterations", "0",
                     "--precision", precision, "--roots"], 2)
    assert all(not s["converged"] and s["gap"] is None and s["iterations"] == "0" for s in rows(failed["states"]))
    assert rows(failed["reference"])[0]["energy"] == "7"

for args in (["8", option, "1"], [str(minimum_n - 1), "--ferromagnetic", option, "1"],
             ["8", "--ferromagnetic", option, "0"],
             ["8", "--ferromagnetic", option, "all", "--max-candidates", "2"],
             ["8", "--ferromagnetic", option, "1", "--max-candidates", "0"],
             ["8", "--ferromagnetic", option, "1", "--tolerance", "nan"]):
    run(args, 1)
for selection in (["--one-defect"], ["--sectors"], ["--singlet-excitation"], ["--through-lines", "4"],
                  ["--excitations", "all"], ["--quantum-numbers", "1,2"], ["--q-spectrum"],
                  ["--q-seed", "1,2"], ["--max-attempts", "5"], [other_option, "1"]):
    run(["8", "--ferromagnetic", option, "1", *selection], 1)
assert len(rows(tables(["8", "--ferromagnetic", option, "100"])["states"])) == 9 - minimum_n
assert len(rows(tables(["100000", "--ferromagnetic", option, "1", "--max-candidates", "1"])["states"])) == 1

with tempfile.TemporaryDirectory() as tmp:
    base = Path(tmp)
    args = ["8", "--ferromagnetic", option, "2", "--roots"]
    exports = ["--json", str(base / "clusters.json")]
    for extension in ("csv", "tsv"):
        for name in ("states", "reference", "string", "roots"):
            exports += [f"--{extension}-table", f"{name}={base / f'clusters.{name}.{extension}'}"]
    run([*args, *exports])
    doc = json.loads((base / "clusters.json").read_text())
    assert len(rows(doc["tables"]["states"])) == 2
    for extension, delimiter in (("csv", ","), ("tsv", "\t")):
        for name, count in (("states", 2), ("reference", 1), ("string", 2), ("roots", 2 * defects)):
            path = base / f"clusters.{name}.{extension}"
            data = list(csv.DictReader((line for line in path.read_text().splitlines() if not line.startswith("#")), delimiter=delimiter))
            assert len(data) == count, (path, data)
assert option in run(["--help"])
assert "nachtergaele-spitzer-starr-2007" in run(["--references"])
