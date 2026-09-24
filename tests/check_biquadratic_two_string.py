"""Targeted singlet: precision, finite deviations, failed estimates and exports."""
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
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=120,
                       env=dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never"))
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def tables(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete"
    assert all(t["summary"]["Outcome"] == ("success" if status == 0 else "partial")
               for t in doc["tables"].values())
    return doc["tables"]


def records(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]


large = []
for precision in precisions:
    args = ["4", "--singlet-excitation", "--roots", "--precision", precision]
    doc = tables(args)
    assert set(doc) == {"states", "reference", "string", "roots"}
    row, = records(doc["states"])
    assert row["converged"] and row["multiplicity"] == "1" and row["through_lines"] == "0"
    assert "not an exhaustive" in doc["states"]["metadata"]["Ordering"]
    with localcontext() as context:
        context.prec = 70
        tol = Decimal("1e-29" if precision == "fp128" else "1e-12")
        expected = (-15 + Decimal(17).sqrt()) / 2
        assert abs(Decimal(row["energy"]) - expected) < tol
        assert abs(Decimal(row["gap"]) - Decimal(17).sqrt()) < 2 * tol
        assert abs(Decimal(row["tl_energy"]) - expected - 3) < tol
        assert abs(2 * Decimal(row["reference_energy"]) - Decimal(21) / 4 - expected) < tol
        string, = records(doc["string"])
        d = Decimal(string["deviation"])
        assert 0 < d < Decimal("0.01")
        assert abs(d.ln() + Decimal(string["log_deviation"])) < tol
    roots = records(doc["roots"])
    assert len(roots) == 2 and all(r["kind"] == "two-string" for r in roots)
    assert Decimal(roots[0]["u_imag"]) == Decimal(roots[1]["u_imag"]).copy_negate()
    replay = tables([*args, "--no-retain"])
    assert all(doc[key]["rows"] == replay[key]["rows"] for key in doc)
    failed = tables(["16", "--singlet-excitation", "--max-iterations", "0", "--precision", precision], 2)
    row, = records(failed["states"])
    assert not row["converged"] and row["gap"] is None and row["iterations"] == "0"
    assert row["energy"] is not None and "unconverged" in row["status"]
    assert not records(failed["reference"])[0]["converged"]
    big = tables(["128", "--singlet-excitation", "--precision", precision])
    large.append((records(big["states"])[0], records(big["string"])[0]))

for state, string in large:
    assert state["converged"] and Decimal(state["gap"]) > 0
    assert abs(Decimal(state["energy"]) - Decimal(large[-1][0]["energy"])) < Decimal("1e-10")
    assert abs(Decimal(string["log_deviation"]) - Decimal(large[-1][1]["log_deviation"])) < Decimal("1e-10")

big = tables(["512", "--singlet-excitation"])
assert "unresolved in rounded u" in big["states"]["metadata"]["String deviation"]
assert Decimal(records(big["string"])[0]["log_deviation"]) > 50
assert Decimal(records(big["string"])[0]["deviation"]) > 0

for extra in (["--q-spectrum"], ["--q-seed", "1,2"], ["--sectors"], ["--through-lines", "0"],
              ["--excitations", "all"], ["--quantum-numbers", "1,2"], ["--max-candidates", "1"],
              ["--max-attempts", "1"]):
    run(["4", "--singlet-excitation", *extra], 1)
for n in ("2", "3", "5"):
    run([n, "--singlet-excitation"], 1)

with tempfile.TemporaryDirectory(prefix="bethe-two-string-") as folder:
    root = Path(folder)
    target, csv_path, tsv_path = root / "singlet.json", root / "levels.csv", root / "string.tsv"
    run(["4", "--singlet-excitation", "--roots", "--json", str(target),
         "--csv-table", f"states={csv_path}", "--tsv-table", f"string={tsv_path}"])
    assert len(json.loads(target.read_text())["tables"]["string"]["rows"]) == 1
    for path, delimiter in ((csv_path, ","), (tsv_path, "\t")):
        lines = path.read_text().splitlines()
        assert any(line.startswith("#") for line in lines)
        rows = list(csv.reader((line for line in lines if line and not line.startswith("#")), delimiter=delimiter))
        assert len(rows) == 2

assert "--singlet-excitation" in run(["--help"])
assert "[bajnok-2020]" in run(["--references"])
