"""Ferromagnetic excitation ordering, coverage and export contracts."""
from decimal import Decimal, localcontext
import csv
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


def run(args, status=0):
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=60, env=env)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def tables(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete"
    for table in doc["tables"].values():
        assert table["summary"]["Outcome"] == ("success" if status == 0 else "partial")
        assert table["metadata"]["Hamiltonian"] == "H=+sum_i (S_i.S_(i+1))^2"
        assert "exact degenerate ferro" in table["metadata"]["Gap reference"]
        assert Decimal(table["summary"]["Run CPU seconds"]) >= 0
    return doc["tables"]


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]


for precision in precisions:
    base = ["4", "--ferromagnetic", "--precision", precision]
    band = tables([*base, "--one-defect"])
    states = rows(band["states"])
    assert [r["mode"] for r in states] == ["3", "2", "1"]
    assert all(r["through_lines"] == "2" and r["defects"] == "1" and r["multiplicity"] == "8" for r in states)
    assert rows(band["reference"])[0]["multiplicity"] == "55"
    assert "NOT the full excited spectrum" in band["states"]["metadata"]["Coverage"]
    assert "not lattice momentum" in band["states"]["metadata"]["Wave number"]
    real = rows(tables([*base, "--excitations", "all", "--roots"])["states"])
    truncated = rows(tables([*base, "--excitations", "1"])["states"])
    selected = rows(tables([*base, "--quantum-numbers", "3"])["states"])
    minimum = rows(tables([*base, "--through-lines", "2"])["states"])
    with localcontext() as context:
        context.prec = 70
        tol = Decimal("1e-29" if precision == "fp128" else "1e-12")
        energies = [6 - Decimal(2).sqrt(), Decimal(6), 6 + Decimal(2).sqrt()]
        for a, b, energy in zip(states, real, energies):
            assert abs(Decimal(a["energy"]) - energy) < tol
            assert abs(Decimal(a["gap"]) - energy + 3) < tol
            assert abs(Decimal(b["energy"]) - energy) < tol
            assert abs(Decimal(b["gap"]) - energy + 3) < tol
            assert abs(Decimal(b["tl_energy"]) - energy + 3) < tol
            assert abs(Decimal(21)/4 - 2*Decimal(b["reference_energy"]) - energy) < tol
        for variant in (truncated, selected, minimum):
            assert len(variant) == 1 and abs(Decimal(variant[0]["energy"]) - energies[0]) < tol
        # The real family misses the lowest ferro singlet; it is a complex-root level.
        module = tables([*base, "--q-spectrum", "--through-lines", "0", "--roots"])
        full = rows(module["states"])
        assert len(full) == 2 and all(r["converged"] for r in full)
        singlets = [(15 - Decimal(17).sqrt())/2, (15 + Decimal(17).sqrt())/2]
        for r, energy in zip(full, singlets):
            assert abs(Decimal(r["energy"]) - energy) < tol
            assert abs(Decimal(r["gap"]) - energy + 3) < tol
        complex_roots = [r for r in rows(module["roots"]) if r["state_id"] == "0"]
        assert len(complex_roots) == 2 and all(abs(Decimal(r["x_imag"])) > 1 for r in complex_roots)
        coefficients = [r for r in rows(module["q_coefficients"]) if r["state_id"] == "0"]
        seed = ",".join(r["coefficient"] for r in coefficients)
        replay = rows(tables([*base, "--q-seed", seed])["states"])[0]
        assert abs(Decimal(replay["energy"]) - singlets[0]) < tol
        restricted = rows(tables([*base, "--excitations", "all", "--through-lines", "0"])["states"])
        assert len(restricted) == 1 and abs(Decimal(restricted[0]["energy"]) - singlets[1]) < tol
    for args in ([*base], [*base, "--quantum-numbers", "none"], [*base, "--q-seed", "none"],
                 [*base, "--through-lines", "4", "--excitations", "all", "--max-iterations", "0"]):
        s = rows(tables(args)["states"])[0]
        assert Decimal(s["energy"]) == 3 and Decimal(s["gap"]) == 0 and s["multiplicity"] == "55"
    assert len(rows(tables([*base, "--q-spectrum"])["states"])) == 3  # Default ferro ell=N-2.
    empty = tables([*base, "--q-spectrum", "--max-attempts", "0"], 2)
    assert not rows(empty["states"])
    assert "NOT guaranteed lowest" in empty["states"]["metadata"]["Ordering"]
    failed = tables([*base, "--q-seed", "1.5,-1.333", "--max-iterations", "0"], 2)
    assert not rows(failed["states"])[0]["converged"] and rows(failed["states"])[0]["gap"] is None
    assert rows(failed["reference"])[0]["converged"]  # Ground space does not require a numerical solve.
    pole = tables([*base, "--q-seed=-1.5", "--max-iterations", "0"], 2)
    assert rows(pole["states"])[0]["energy"] is None
    no_retain = tables([*base, "--one-defect", "--no-retain"])
    assert all(no_retain[name]["rows"] == table["rows"] for name, table in band.items())

for n in (2, 3, 5, 64, 129):
    band = rows(tables([str(n), "--ferromagnetic", "--one-defect"])["states"])
    assert len(band) == n-1
    assert abs(float(band[0]["gap"]) - (3-2*math.cos(math.pi/n))) < 2e-13
    assert abs(float(band[-1]["gap"]) - (3+2*math.cos(math.pi/n))) < 2e-13
overflow = tables(["48", "--ferromagnetic", "--one-defect"])
assert all(r["multiplicity"] is None for r in rows(overflow["states"]))
partial = tables(["8", "--ferromagnetic", "--excitations", "all", "--through-lines", "2", "--max-iterations", "0"], 2)
assert "failed" in partial and rows(partial["reference"])[0]["converged"]
failed = rows(tables(["8", "--ferromagnetic", "--quantum-numbers", "1,2,3", "--max-iterations", "0"], 2)["states"])[0]
assert failed["gap"] is None and not failed["converged"]

for args in (["4", "--one-defect"], ["0", "--ferromagnetic"], ["1", "--ferromagnetic"],
             ["4", "--ferromagnetic", "--sectors"], ["4", "--ferromagnetic", "--singlet-excitation"],
             ["4", "--ferromagnetic", "--through-lines", "0"],
             ["4", "--ferromagnetic", "--one-defect", "--roots"],
             ["4", "--ferromagnetic", "--one-defect", "--excitations", "all"],
             ["4", "--ferromagnetic", "--one-defect", "--q-spectrum"],
             ["4", "--ferromagnetic", "--one-defect", "--quantum-numbers", "1"],
             ["4", "--ferromagnetic", "--one-defect", "--through-lines", "2"],
             ["4", "--ferromagnetic", "--one-defect", "--max-candidates", "2"],
             ["4", "--ferromagnetic", "--one-defect", "--max-candidates", "0"],
             ["4", "--ferromagnetic", "--tolerance", "nan"],
             ["5", "--ferromagnetic", "--quantum-numbers", "1"],
             ["5", "--ferromagnetic", "--q-spectrum"],
             ["5", "--ferromagnetic", "--excitations", "all"]):
    run(args, 1)

with tempfile.TemporaryDirectory(prefix="bethe-ferro-") as folder:
    root = Path(folder)
    js, cv, tv = root/"ferro.json", root/"band.csv", root/"reference.tsv"
    run(["5", "--ferromagnetic", "--one-defect", "--quiet", "--json", str(js),
         "--csv-table", f"states={cv}", "--tsv-table", f"reference={tv}"])
    assert len(json.loads(js.read_text())["tables"]["states"]["rows"]) == 4
    for path, delimiter, count in ((cv, ",", 4), (tv, "\t", 1)):
        content = path.read_text().splitlines()
        assert any(line.startswith("#") for line in content)
        records = list(csv.reader((line for line in content if line and not line.startswith("#")), delimiter=delimiter))
        assert len(records) == count+1

assert "--ferromagnetic" in run(["--help"]) and "--one-defect" in run(["--help"])
assert "[koma-nachtergaele-1997]" in run(["--references"])
assert "[zhou-2025-biquadratic]" in run(["--references"])
