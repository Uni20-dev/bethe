"""SU(3) dispersion conventions, native precision and export contracts."""
import csv
from decimal import Decimal, localcontext
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
    result = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=30)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    return result.stdout


def table(args, status=0):
    document = json.loads(run([*args, "--format", "json"], status))
    assert document["status"] == "complete"
    result = document["tables"]["dispersion"]
    assert result["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    assert Decimal(result["summary"]["Run CPU seconds"]) >= 0
    assert "H=J sum P" in result["metadata"]["Hamiltonian"]
    return result


def rows(t):
    return [dict(zip([c["id"] for c in t["columns"]], row)) for row in t["rows"]]


for precision in precisions:
    base = ["--points", "3", "--precision", precision]
    t = table(base)
    records = rows(t)
    assert len(records) == 12
    assert [r["branch"] for r in records] == [b for b in ("3", "bar3", "two-soliton", "four-soliton") for _ in range(3)]
    assert all(r["status"] == "converged" for r in records)
    assert all(r["lower"] is None and r["upper"] is None for r in records[:6])
    assert all(r["energy"] is None for r in records[6:])
    for i in (0, 2, 3, 5):
        assert records[i]["energy"] == "0"
    assert rows(table([*base, "--no-retain"])) == records
    with localcontext() as context:
        context.prec = 85
        tolerance = {"fp64": Decimal("4e-14"), "long-double": Decimal("4e-17"), "fp128": Decimal("4e-32")}[precision]
        for branch, reference in (
            ("3", "0.672943171908697176654766974738581598626465808498253705740688258859362234"),
            ("bar3", "0.564928975070772007354510973939731453793808737285571166075921275577481670"),
        ):
            r = rows(table(["--branch", branch, "--momentum", "0.3", "--precision", precision]))[0]
            assert abs(Decimal(r["energy"]) - Decimal(reference)) < tolerance
        two, four = records[7], records[10]
        assert abs(Decimal(two["lower"]) - 2 * Decimal(four["lower"])) < tolerance
        folded = rows(table([*base, "--folded"]))
        assert folded[:6] == records[:6]
        assert abs(Decimal(folded[7]["lower"]) - Decimal(four["lower"])) < tolerance
        assert Decimal(folded[6]["upper"]) > 0  # Q=0 has nonzero image envelope
    scaled = rows(table([*base, "--exchange", "2"]))
    assert abs(float(scaled[7]["upper"]) - 2 * float(records[7]["upper"])) < 1e-12

for args in (["--exchange", "0"], ["--exchange", "nan"], ["--exchange", "inf"],
             ["--points", "1"], ["--points", "-1"], ["--points", "1000001"],
             ["--momentum", "3"], ["--branch", "bar3", "--momentum", "3"],
             ["--momentum", "nan"], ["--momentum", "-1"], ["--branch", "unknown"],
             ["--branch", "3", "--folded"], ["--branch", "bar3", "--folded"],
             ["--points", "3", "--momentum", "0"]):
    run(args, 1)

overflow = rows(table(["--branch", "four-soliton", "--points", "2", "--exchange", "2e307"], 2))
assert all(r["lower"] is None and r["upper"] is None and r["status"] == "precision_limit" for r in overflow)

with tempfile.TemporaryDirectory() as directory:
    js, cs, ts = [Path(directory) / ("uls." + suffix) for suffix in ("json", "csv", "tsv")]
    args = ["--points", "3", "--json", str(js), "--csv", str(cs), "--tsv", str(ts), "--format", "plain"]
    screen = run(args)
    assert "CPU time" in screen and "Used for:" not in screen
    expected = [{k: "" if v is None else str(v) for k, v in r.items()}
                for r in rows(json.loads(js.read_text())["tables"]["dispersion"])]
    for path, delimiter in ((cs, ","), (ts, "\t")):
        text = path.read_text()
        assert "Hamiltonian" in text and "Continuum momentum" in text
        actual = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
        assert actual == expected
    original = js.read_bytes()
    run(["--points", "3", "--json", str(js)], 1)
    assert js.read_bytes() == original
    run(["--exchange", "-1", "--force", "--json", str(js)], 1)
    assert js.read_bytes() == original
