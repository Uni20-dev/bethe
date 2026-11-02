# SPDX-License-Identifier: GPL-3.0-or-later
"""Boundary normalization, signed response, failure records and precision exports."""
import csv
from decimal import Decimal, localcontext
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, mpfr = sys.argv[1:]
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")


def run(args, code=0):
    result = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=60)
    assert result.returncode == code, (args, result.stdout, result.stderr)
    return result


def rows(table):
    return [dict(zip((c["id"] for c in table["columns"]),
                     (v["value"] if isinstance(v, dict) else v for v in row))) for row in table["rows"]]


args = ["2", "--time-max", "2", "--time-points", "3", "--format", "json", "--spectrum", "--density"]
tables = json.loads(run(args).stdout)["tables"]
for t, row in enumerate(rows(tables["time"])):
    assert abs(float(row["local_real"]) - math.cos(t) / 4) < 2e-14
    assert abs(float(row["local_imag"]) + math.sin(t) / 4) < 2e-14
    assert float(row["end_to_end_real"]) == -float(row["local_real"])
    assert float(row["end_to_end_imag"]) == -float(row["local_imag"])
line = rows(tables["spectrum"])[0]
assert line["candidate_id"] == "N2-obc-I"
assert abs(float(line["gap"]) - 1) < 2e-14
assert abs(float(line["zz_weight"]) - .25) < 2e-14
assert float(line["raising_weight"]) == 2 * float(line["zz_weight"])
assert float(line["end_to_end_coefficient"]) == -float(line["zz_weight"])
moment = rows(tables["moments"])[0]
assert float(moment["full_local_weight"]) == .25
assert abs(float(moment["full_local_first_moment"]) - .25) < 2e-14
rho = {(int(r["row"]), int(r["column"])): float(r["value"]) for r in rows(tables["density"])}
assert abs(sum(rho[(i, i)] for i in range(4)) - 1) < 2e-14
assert abs(rho[(1, 2)] + .5) < 2e-14
assert tables["time"]["metadata"]["Channel"] == "zz"

small = ["8", "--time-points", "2", "--format", "json", "--diagnostics", "--spectrum"]
serial = json.loads(run(small).stdout)["tables"]
parallel = json.loads(run([*small, "--threads", "2"]).stdout)["tables"]
for name in ("time", "spectrum", "diagnostics", "moments"):
    assert rows(serial[name]) == rows(parallel[name]), name
assert len(rows(serial["diagnostics"])) == 10
assert float(rows(serial["moments"])[0]["omitted_weight"]) > 0

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    diagnostics = root / "failed.csv"
    failed = json.loads(run(["8", "--max-iterations", "0", "--format", "json", "--csv-table",
                             f"diagnostics={diagnostics}"], 2).stdout)["tables"]
    assert rows(failed["time"]) == [] and rows(failed["moments"]) == []
    assert "diagnostics" not in failed
    text = diagnostics.read_text()
    assert "# Hamiltonian:" in text and "# Precision:" in text
    failed_rows = list(csv.DictReader(line for line in text.splitlines() if not line.startswith("#")))
    assert len(failed_rows) == 10
    assert all(r["root_status"] == r["factor_status"] == "not_attempted" for r in failed_rows)
    assert all(r["residual"] == r["pivot_margin"] == "" for r in failed_rows)
    exported = root / "roots.tsv"
    doc = json.loads(run(["2", "--time-points", "2", "--format", "json", "--tsv-table",
                          f"roots={exported}" ]).stdout)
    assert "roots" not in doc["tables"] and "quantum_number" in exported.read_text()
    for options in (["0"], ["3"], ["8", "--max-candidates", "1"], ["2", "--time-points", "1"],
                    ["2", "--time-points", "1000001"], ["2", "--time-max", "-1"],
                    ["2", "--time-max", "inf"], ["2", "--tolerance", "0"]):
        target = root / "must-not-exist.csv"
        result = subprocess.run([program, *options, "--csv", str(target)], capture_output=True, text=True, env=env)
        assert result.returncode != 0 and not target.exists(), (options, result)
    target = root / "unknown.csv"
    bad = subprocess.run([program, "2", "--csv-table", f"failed={target}"], capture_output=True, text=True, env=env)
    assert bad.returncode != 0 and not target.exists()

help_text = run(["--help"]).stdout
for name in ("time", "moments", "spectrum", "diagnostics", "roots", "density"):
    assert name in help_text
references = run(["--references"]).stdout
assert "Grijalva" in references and "Wang" in references
if mpfr.upper() in ("ON", "TRUE", "1"):
    table = json.loads(run(["2", "--precision", "bits-128", "--threads", "2", "--time-points", "2",
                            "--format", "json", "--spectrum"]).stdout)["tables"]["spectrum"]
    assert table["metadata"]["Working precision (bits)"] == "128"
    with localcontext() as context:
        context.prec = 50
        assert abs(Decimal(rows(table)[0]["zz_weight"]) - Decimal(".25")) < Decimal("1e-32")
print("boundary CLI: dimer time/signs, density, parallel identity, failures, exports and precision passed")
