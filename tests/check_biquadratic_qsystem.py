"""Q-system scientific and typed-output contracts; standard-library only."""
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


def run(args, status=0):
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=60,
                       env=dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never"))
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def tables(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete"  # Document transport, not scientific convergence.
    return doc["tables"]


def records(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]


for precision in precisions:
    args = ["4", "--q-spectrum", "--roots", "--precision", precision]
    doc = tables(args)
    assert set(doc) == {"states", "reference", "q_coefficients", "roots"}
    rows = records(doc["states"])
    assert len(rows) == 2 and all(r["converged"] for r in rows)
    assert all(r["multiplicity"] == "1" and r["through_lines"] == "0" for r in rows)
    assert doc["states"]["metadata"]["Expected module dimension"] == "2"
    with localcontext() as context:
        context.prec = 70
        tol = Decimal("1e-29" if precision == "fp128" else "1e-12")
        expected = [(-15 - Decimal(17).sqrt()) / 2, (-15 + Decimal(17).sqrt()) / 2]
        ground = Decimal(records(doc["reference"])[0]["energy"])
        for row, energy in zip(rows, expected):
            assert abs(Decimal(row["energy"]) - energy) < tol, (precision, row)
            assert abs(Decimal(row["gap"]) - energy + ground) < tol
            assert abs(Decimal(row["tl_energy"]) - energy - 3) < tol
            assert abs(2 * Decimal(row["reference_energy"]) - Decimal(21) / 4 - energy) < tol
        roots = [r for r in records(doc["roots"]) if r["state_id"] == "1"]
        assert len(roots) == 2 and all(abs(Decimal(r["x_imag"])) > 1 for r in roots)
        # Feed the exported native-precision polynomial back to the selected solver.
        seed = ",".join(r["coefficient"] for r in records(doc["q_coefficients"]) if r["state_id"] == "1")
        selected = records(tables(["4", "--q-seed", seed, "--precision", precision])["states"])[0]
        assert selected["converged"] and abs(Decimal(selected["energy"]) - expected[1]) < tol
    replay = tables([*args, "--no-retain"])
    assert all(doc[key]["rows"] == replay[key]["rows"] for key in doc)
    failed = tables(["4", "--q-seed", "1.5,-1.333", "--max-iterations", "0", "--precision", precision], 2)
    row, = records(failed["states"])
    assert not row["converged"] and row["gap"] is None and row["bethe_residual"] is None
    pole = tables(["4", "--q-seed=-1.5", "--max-iterations", "0", "--precision", precision], 2)
    assert records(pole["states"])[0]["energy"] is None
    vacuum = tables(["4", "--q-seed", "none", "--precision", precision])
    assert records(vacuum["states"])[0]["multiplicity"] == "55"
    assert Decimal(records(vacuum["states"])[0]["energy"]) == -3
    empty = tables(["4", "--q-spectrum", "--max-attempts", "0", "--precision", precision], 2)
    assert not records(empty["states"])
    assert "NOT guaranteed lowest" in empty["states"]["metadata"]["Ordering"]

module = tables(["6", "--q-spectrum", "--through-lines", "2"])
assert len(records(module["states"])) == 9
assert all(r["multiplicity"] == "8" for r in records(module["states"]))

for args in (["4", "--q-spectrum", "--excitations", "all"], ["4", "--q-spectrum", "--sectors"],
             ["4", "--q-spectrum", "--quantum-numbers", "1,2"], ["4", "--q-seed", "1,2", "--through-lines", "0"],
             ["4", "--q-seed", "1,2", "--q-spectrum"], ["4", "--max-attempts", "2"],
             ["10", "--q-spectrum"], ["34", "--q-seed", "none"], ["4", "--q-seed", "1,,2"],
             ["4", "--q-seed", "nan"], ["4", "--q-spectrum", "--through-lines", "1"],
             ["4", "--q-spectrum", "--max-attempts", "-1"]):
    run(args, 1)

with tempfile.TemporaryDirectory(prefix="bethe-qsystem-") as folder:
    root = Path(folder)
    target = root / "levels.json"
    csv_path, tsv_path = root / "levels.csv", root / "coefficients.tsv"
    run(["4", "--q-spectrum", "--roots", "--json", str(target),
         "--csv-table", f"states={csv_path}", "--tsv-table", f"q_coefficients={tsv_path}"])
    assert len(json.loads(target.read_text())["tables"]["states"]["rows"]) == 2
    for path, delimiter, count in ((csv_path, ",", 2), (tsv_path, "\t", 4)):
        content = path.read_text().splitlines()
        assert any(line.startswith("#") for line in content)
        parsed = list(csv.reader((line for line in content if not line.startswith("#") and line), delimiter=delimiter))
        assert len(parsed) == count + 1, parsed

assert "[bajnok-2020]" in run(["--references"])
assert "--q-spectrum" in run(["--help"])
