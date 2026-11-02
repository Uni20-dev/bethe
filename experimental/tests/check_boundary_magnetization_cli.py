# SPDX-License-Identifier: GPL-3.0-or-later
"""Signed magnetization exports, reconstruction, failure and precision contracts."""
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


base = ["2", "--time-max", "100", "--time-points", "3", "--format", "json", "--pairs", "--states"]
tables = json.loads(run(base).stdout)["tables"]
for row in rows(tables["time"]):
    assert abs(float(row["left_magnetization"]) - .5) < 2e-14
    assert row["left_magnetization"] == row["right_magnetization"]
assert len(rows(tables["pairs"])) == 1
pair = rows(tables["pairs"])[0]
assert float(pair["matrix_element"]) == .5
assert float(pair["frequency"]) == 0
assert abs(float(pair["coefficient"]) - .5) < 2e-14
assert pair["left_candidate"] == pair["right_candidate"] == "N2-obc-I"
summary = rows(tables["summary"])[0]
assert abs(float(summary["projected_norm"]) - 1) < 2e-14
assert summary["matrix_ready"] is True

small = ["8", "--time-points", "3", "--format", "json", "--states", "--pairs"]
serial = json.loads(run(small).stdout)["tables"]
parallel = json.loads(run([*small, "--threads", "2"]).stdout)["tables"]
for name in ("time", "states", "pairs", "summary"):
    assert rows(serial[name]) == rows(parallel[name]), name
assert len(rows(serial["states"])) == 10
assert len(rows(serial["pairs"])) == 55
assert float(rows(serial["summary"])[0]["omitted_norm"]) > 0
for row in rows(serial["time"]):
    t = float(row["time"])
    for endpoint, parity in (("left_magnetization", False), ("right_magnetization", True)):
        value = math.fsum(float(p["coefficient"]) * math.cos(float(p["frequency"]) * t)
                          * (int(p["opposite_end_parity"]) if parity else 1) for p in rows(serial["pairs"]))
        assert abs(value - float(row[endpoint])) < 1e-13

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    target = root / "pairs.csv"
    failed = json.loads(run(["8", "--max-iterations", "0", "--format", "json", "--csv-table",
                             f"pairs={target}"], 2).stdout)["tables"]
    assert rows(failed["time"]) == [] and "pairs" not in failed
    assert rows(failed["summary"])[0]["projected_norm"] is None
    text = target.read_text()
    assert "# Hamiltonian:" in text and "# Precision:" in text
    pairs = list(csv.DictReader(line for line in text.splitlines() if not line.startswith("#")))
    assert len(pairs) == 55
    assert all(p["status"] == "not_attempted" and p["matrix_element"] == p["coefficient"] == "" for p in pairs)
    exported = root / "roots.tsv"
    doc = json.loads(run(["2", "--time-points", "2", "--format", "json", "--tsv-table",
                          f"roots={exported}" ]).stdout)
    assert "roots" not in doc["tables"] and "quantum_number" in exported.read_text()
    for options in (["0"], ["3"], ["64", "--max-pairs", "139655"], ["2", "--time-points", "1"],
                    ["2", "--time-points", "1000001"], ["2", "--time-max", "-1"],
                    ["2", "--time-max", "inf"], ["2", "--tolerance", "0"]):
        target = root / "must-not-exist.csv"
        bad = subprocess.run([program, *options, "--csv", str(target)], capture_output=True, text=True, env=env)
        assert bad.returncode != 0 and not target.exists(), (options, bad)
    bad = subprocess.run([program, "2", "--csv-table", f"failed={root/'unknown.csv'}"], capture_output=True, text=True, env=env)
    assert bad.returncode != 0 and not (root/'unknown.csv').exists()

help_text = run(["--help"]).stdout
for name in ("time", "summary", "states", "pairs", "roots"):
    assert name in help_text
assert "Grijalva" in run(["--references"]).stdout
if mpfr.upper() in ("ON", "TRUE", "1"):
    table = json.loads(run(["2", "--precision", "bits-128", "--threads", "2", "--time-points", "2",
                            "--format", "json"]).stdout)["tables"]["time"]
    assert table["metadata"]["Working precision (bits)"] == "128"
    with localcontext() as context:
        context.prec = 50
        assert abs(Decimal(rows(table)[0]["left_magnetization"]) - Decimal(".5")) < Decimal("1e-32")
print("boundary magnetization CLI: signed terms, reconstruction, concurrency, failures and precision passed")
