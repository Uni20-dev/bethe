"""Ordered independent Bethe-Yang levels, complete diagnostics and bounded output."""
import csv
from decimal import Decimal
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, wide = sys.argv[1:]
precisions = ["fp64", "long-double"] + (["fp128"] if wide.upper() in ("ON", "TRUE", "1") else [])
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")


def run(args, status=0):
    result = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=60)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    return result.stdout


def table(args, status=0):
    document = json.loads(run([*args, "--format", "json"], status))
    assert document["status"] == "complete"
    result = document["tables"]["levels"]
    assert result["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    assert Decimal(result["summary"]["Run CPU seconds"]) >= 0
    assert result["metadata"]["CPU timing scope"] == "level grid including ordered row output"
    parallel = result["metadata"]["Execution mode"] == "bounded Async level grid"
    assert result["metadata"]["Maximum levels in flight"] == ("128" if parallel else "1")
    return result


def rows(result):
    return [dict(zip([c["id"] for c in result["columns"]], row)) for row in result["rows"]]


def compare(args, status=0, expected=None):
    expected = expected or table(args, status)
    for threads in (1, 4):
        actual = table([*args, "--threads", str(threads), "--no-retain"], status)
        assert actual["columns"] == expected["columns"]
        assert actual["rows"] == expected["rows"]  # every status, error and work count
        assert actual["metadata"]["Scheduler concurrency limit"] == str(threads)
        parallel = (threads > 1 and len(expected["rows"]) > 1
                    and Decimal(actual["metadata"]["Coupling p"]) != 1)
        assert actual["metadata"]["Execution mode"] == ("bounded Async level grid" if parallel else "direct serial")
    return expected


for precision in precisions:
    # Arithmetic retains its selected precision; use a loose numerical budget
    # for scheduling checks. The existing physics test uses default tolerances.
    base = ["--p", "2", "--length", "10", "--levels", "7", "--number", "1.5",
            "--tolerance", "1e-7", "--precision", precision]
    expected = compare(base)
    assert [r["number"] for r in rows(expected)] == [i + 1.5 for i in range(7)]
    assert all(r["status"] == "converged" for r in rows(expected))
    updates = [int(r["iterations"]) for r in rows(expected)]
    assert min(updates) < max(updates)
    partial = compare([*base, "--max-iterations", str(min(updates))], 2)
    assert {r["status"] for r in rows(partial)} == {"converged", "iteration_limit"}
    for flag in ("--max-iterations", "--max-phase-evaluations", "--max-fourier-cutoffs"):
        failed = compare([*base, flag, "0"], 2)
        assert all(r["energy"] is None and r["rapidity"] is None for r in rows(failed))
    free = compare(["--p", "1", "--length", "10", "--levels", "129", "--charge", "-2",
                    "--precision", precision])
    assert len(free["rows"]) == 129
    assert all(r["charge"] == -2 and int(r["evaluations"]) == 0 for r in rows(free))
    compare(["--p", "2", "--length", "10", "--tolerance", "1e-7", "--precision", precision])

# Interacting scan crosses the default admission window, with deliberately
# modest quadrature accuracy to keep a scheduler regression inexpensive in CI.
large = compare(["--p", "2", "--length", "500", "--levels", "129", "--tolerance", "1e-6"])
assert len(large["rows"]) == 129

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    js, cs, ts = [root / ("levels." + ext) for ext in ("json", "csv", "tsv")]
    base = ["--p", "2", "--length", "10", "--levels", "7", "--tolerance", "1e-7"]
    args = [*base, "--threads", "4", "--stream", "--no-retain", "--quiet",
            "--json", str(js), "--csv", str(cs), "--tsv", str(ts)]
    run(args)
    exported = json.loads(js.read_text())["tables"]["levels"]
    assert exported["rows"] == table(base)["rows"]
    expected = [{k: "" if v is None else str(v) for k, v in row.items()} for row in rows(exported)]
    for path, delimiter in ((cs, ","), (ts, "\t")):
        actual = list(csv.DictReader((line for line in path.read_text().splitlines()
                                     if not line.startswith("#")), delimiter=delimiter))
        assert actual == expected
    original = js.read_bytes()
    for extra in (["--threads", "0"], ["--threads", "-2"], ["--threads", "2", "--number", "1"],
                  ["--threads", "2", "--levels", "10001"], ["--threads", "2", "--tolerance", "0"]):
        # Avoid duplicate numeric options so this reaches the intended check.
        run(["--p", "2", "--length", "10", *extra, "--json", str(js), "--force"], 1)
        assert js.read_bytes() == original

assert "--threads" in run(["--help"])
print("Sine-Gordon Bethe-Yang parallel level contracts passed")
