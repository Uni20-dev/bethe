# SPDX-License-Identifier: GPL-3.0-or-later
"""Finite-four frontend normalization, failure accounting and export contract."""
from decimal import Decimal, localcontext
import csv
import json
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


selected = ["12", "--holes=-3,-1,0,2", "--branch", "1"]
zz = json.loads(run([*selected, "--format", "json", "--diagnostics"]).stdout)["tables"]
raising = json.loads(run([*selected, "--channel", "raising", "--format", "json"]).stdout)["tables"]
line = rows(zz["spectrum"])[0]
assert line["candidate_id"] == "N12-s4-h-3,-1,0,2-b1"
assert int(line["momentum_index"]) == 3 and line["multiplicity"] == 1
assert float(rows(raising["spectrum"])[0]["weight"]) == 2 * float(line["weight"])
assert abs(float(line["weight"])/0.00000600720973266099492641108831817 - 1) < 1e-8
assert rows(zz["diagnostics"])[0]["accepted"] is True
assert "Ground energy" in zz["spectrum"]["metadata"]
assert zz["spectrum"]["metadata"]["Channel"] == "zz"
assert rows(zz["moments"])[3]["retained_weight"] == line["weight"]
permuted = json.loads(run(["12", "--holes=2,-3,0,-1", "--branch", "1", "--format", "json",
                           "--diagnostics"]).stdout)["tables"]
assert rows(permuted["spectrum"]) == rows(zz["spectrum"])
assert rows(permuted["diagnostics"]) == rows(zz["diagnostics"])

zero = json.loads(run(["12", "--holes=-2,-1,0,2", "--branch", "-1", "--momentum-index", "0",
                       "--format", "json", "--diagnostics"]).stdout)["tables"]
assert len(rows(zero["spectrum"])) == 1
assert Decimal(str(rows(zero["spectrum"])[0]["weight"])) == 0
assert rows(zero["diagnostics"])[0]["accepted"] is True
assert rows(zero["diagnostics"])[0]["factor_status"] == "converged"
assert Decimal(str(rows(zero["moments"])[0]["retained_weight"])) == 0

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    rejected = root / "rejected.csv"
    result = run(["12", "--holes=-2,-1,1,2", "--branch", "0", "--format", "json",
                  "--csv-table", f"diagnostics={rejected}"], 2)
    doc = json.loads(result.stdout)
    assert rows(doc["tables"]["spectrum"]) == []
    assert "diagnostics" not in doc["tables"]  # A file does not request screen output.
    text = rejected.read_text()
    assert "# Hamiltonian:" in text and "# Precision:" in text
    records = list(csv.DictReader(line for line in text.splitlines() if not line.startswith("#")))
    assert len(records) == 1
    assert records[0]["root_status"] == "unsupported_geometry"
    assert records[0]["factor_status"] == "not_attempted"
    assert records[0]["deviation"] == ""
    assert records[0]["diagnostic"]
    for args in (["12", "--max-candidates", "1"],
                 ["12", "--max-candidates", "104"],
                 ["36"],
                 ["12", "--momentum-index", "3", "--max-candidates", "1"]):
        target = root / "over-budget.csv"
        result = run([*args, "--csv", str(target)], 1)
        assert "four-spinon catalogue exceeds max_candidates" in result.stderr, result
        assert not target.exists(), result
    for args in ([*selected, "--momentum-index", "9"],
                 ["12", "--branch", "0"],
                 ["12", "--holes=-2,-2,1,2", "--branch", "0"]):
        target = root / "must-not-exist.csv"
        result = subprocess.run([program, *args, "--csv", str(target)], capture_output=True, text=True, env=env)
        assert result.returncode != 0 and not target.exists(), result
    target = root / "unknown.csv"
    result = subprocess.run([program, *selected, "--csv-table", f"unknown={target}"],
                            capture_output=True, text=True, env=env)
    assert result.returncode != 0 and not target.exists()

help_text = run(["--help"]).stdout
for name in ("spectrum", "diagnostics", "moments", "roots"):
    assert name in help_text
assert run(["--references"]).stdout
if mpfr.upper() in ("ON", "TRUE", "1"):
    doc = json.loads(run([*selected, "--precision", "bits-128", "--threads", "2", "--format", "json"]).stdout)
    table = doc["tables"]["spectrum"]
    assert table["metadata"]["Working precision (bits)"] == "128"
    with localcontext() as context:
        context.prec = 50
        expected = Decimal("0.000012014419465321989852822176636341757165614650731119761969") / 2
        assert abs(Decimal(rows(table)[0]["weight"])/expected - 1) < Decimal("1e-31")
print("finite-four CLI normalization, rejection records, strict exports and precision passed")
