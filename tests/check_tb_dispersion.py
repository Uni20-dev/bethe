"""TB dispersion normalization, half-integer labels and export contracts."""
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
    assert "S.S-(S.S)^2" in result["metadata"]["Hamiltonian"]
    assert result["metadata"]["Spinon spin"] == "0.5"
    return result


def rows(t):
    return [dict(zip([c["id"] for c in t["columns"]], row)) for row in t["rows"]]


for precision in precisions:
    base = ["--points", "3", "--precision", precision]
    t = table(base)
    records = rows(t)
    assert len(records) == 9
    assert [r["branch"] for r in records] == [b for b in ("spinon", "two-spinon", "four-spinon") for _ in range(3)]
    assert all(r["status"] == "converged" for r in records)
    assert all(r["spin"] == 0.5 and r["lower"] is None and r["upper"] is None for r in records[:3])
    assert all(r["spin"] is None and r["energy"] is None for r in records[3:])
    assert records[3]["sectors"] == "S=0,1" and records[6]["sectors"] == "S=0,1,2"
    assert records[0]["energy"] == records[2]["energy"] == "0"
    assert records[3]["upper"] == "0" and Decimal(records[6]["upper"]) > 0
    assert rows(table([*base, "--no-retain"])) == records
    with localcontext() as context:
        context.prec = 85
        tolerance = {"fp64": Decimal("8e-14"), "long-double": Decimal("8e-17"), "fp128": Decimal("8e-32")}[precision]
        reference = Decimal("1.8568082204692037760139169230174695786304118186757373077249143016552658568")
        r = rows(table(["--branch", "spinon", "--momentum", "0.3", "--precision", precision]))[0]
        assert abs(Decimal(r["energy"]) - reference) < tolerance
        velocity = Decimal(t["metadata"]["Spinon velocity"])
        assert abs(Decimal(records[1]["energy"]) - velocity) < tolerance
        assert abs(Decimal(records[4]["upper"]) - 2 * velocity) < tolerance
        assert abs(Decimal(records[6]["upper"]) - 4 * velocity) < tolerance
        folded = rows(table([*base, "--folded"]))
        assert folded[:3] == records[:3]
        assert abs(Decimal(folded[3]["upper"]) - 2 * velocity) < tolerance
        assert abs(Decimal(folded[7]["upper"]) - 4 * velocity) < tolerance

for args in (["--exchange", "0"], ["--exchange", "-1"], ["--exchange", "nan"], ["--exchange", "inf"],
             ["--points", "1"], ["--points", "-1"], ["--points", "1000001"],
             ["--momentum", "4"], ["--branch", "two-spinon", "--momentum", "7"],
             ["--momentum", "nan"], ["--momentum", "-1"], ["--branch", "unknown"],
             ["--branch", "spinon", "--folded"], ["--points", "3", "--momentum", "0"]):
    run(args, 1)

overflow = rows(table(["--branch", "four-spinon", "--points", "2", "--exchange", "1e307"], 2))
assert all(r["lower"] is None and r["upper"] is None and r["status"] == "precision_limit" for r in overflow)
underflow = rows(table(["--branch", "spinon", "--momentum", "1e-100", "--exchange", "1e-300"], 2))
assert underflow[0]["energy"] is None and underflow[0]["status"] == "precision_limit"

with tempfile.TemporaryDirectory() as directory:
    js, cs, ts = [Path(directory) / ("tb." + suffix) for suffix in ("json", "csv", "tsv")]
    screen = run(["--points", "3", "--json", str(js), "--csv", str(cs), "--tsv", str(ts), "--format", "plain"])
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
