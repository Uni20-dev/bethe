"""Complete ordered outputs, partial scans and export preflight with concurrency."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, fp128 = sys.argv[1:]
precisions = ["fp64", "long-double"] + (["fp128"] if fp128.upper() in ("ON", "TRUE", "1") else [])
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")
variable = {"Command", "Date", "CPU time", "Run CPU seconds", "Elapsed seconds",
            "Scheduler concurrency limit", "Solver concurrency"}


def stable(value):
    if isinstance(value, dict):
        return {k: stable(v) for k, v in value.items() if k not in variable}
    if isinstance(value, list):
        return [stable(v) for v in value]
    return value


def run(args, allowed=(0, 2)):
    p = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=30)
    assert p.returncode in allowed, (args, p.returncode, p.stdout, p.stderr)
    return p


cases = [
    ["8", "--sectors"],
    ["10", "--ferromagnetic", "--bound-pairs", "all"],
    ["10", "--ferromagnetic", "--bound-triples", "all"],
    ["10", "--ferromagnetic", "--bound-quartets", "all"],
    ["10", "--ferromagnetic", "--pair-defects", "all", "--mixed-window", "3"],
    ["10", "--ferromagnetic", "--pair-defects", "all", "--mixed-window", "3", "--real-defects", "2"],
    ["9", "--ferromagnetic", "--pair-defects", "1", "--mixed-window", "3"],
    ["12", "--ferromagnetic", "--triple-defects", "all", "--mixed-window", "3"],
    ["10", "--ferromagnetic", "--two-pair-states", "all"],
    ["10", "--ferromagnetic", "--two-pair-states", "1", "--pair-window", "3"],
    ["8", "--excitations", "all"],  # Existing continuous real-root path.
]
for precision in precisions:
    for case in cases:
        for budget in ([], ["--max-iterations", "0"], ["--max-iterations", "1"]):
            args = [*case, *budget, "--precision", precision, "--roots", "--format", "json"]
            serial = run(args)
            expected = stable(json.loads(serial.stdout))
            for threads in (2, 4):
                concurrent = run([*args, "--threads", str(threads)])
                assert concurrent.returncode == serial.returncode, (args, threads)
                actual = json.loads(concurrent.stdout)
                assert stable(actual) == expected, (args, threads)
                assert actual["tables"]["states"]["metadata"]["Scheduler concurrency limit"] == str(threads)

# A scan of one candidate is direct, even when more workers are requested.
single = json.loads(run(["8", "--ferromagnetic", "--bound-quartets", "1",
                        "--threads", "8", "--format", "json"]).stdout)
assert single["tables"]["states"]["metadata"]["Solver concurrency"] == "1"

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    output, csv, tsv = root / "levels.json", root / "states.csv", root / "roots.tsv"
    args = ["10", "--ferromagnetic", "--pair-defects", "all", "--mixed-window", "3", "--threads", "4"]
    expected = json.loads(run([*args, "--roots", "--format", "json"]).stdout)
    run([*args, "--json", str(output), "--csv-table", f"states={csv}", "--tsv-table", f"roots={tsv}",
         "--roots", "--no-retain", "--quiet"])
    assert stable(json.loads(output.read_text())) == stable(expected)
    assert "state" in csv.read_text() and "state" in tsv.read_text()
    invalid = [
        ["8", "--threads", "1"],
        ["8", "--ferromagnetic", "--one-defect", "--threads", "4"],
        ["8", "--ferromagnetic", "--pair-defect", "4,3", "--threads", "4"],
        ["8", "--ferromagnetic", "--two-pairs", "1,2", "--threads", "4"],
        ["8", "--ferromagnetic", "--triple-defect", "1,1", "--threads", "4"],
        ["8", "--q-spectrum", "--threads", "4"],
        ["8", "--singlet-excitation", "--threads", "4"],
        ["8", "--sectors", "--threads", "0"],
        ["8", "--sectors", "--threads", "-2"],
        ["8", "--sectors", "--threads", "1.5"],
        ["7", "--sectors", "--threads", "4"],
        ["1000000000000000000", "--sectors", "--threads", "4"],
        ["10", "--ferromagnetic", "--pair-defects", "all", "--max-candidates", "1", "--threads", "4"],
    ]
    for args in invalid:
        output.write_text("preserve json")
        csv.write_text("preserve csv")
        run([*args, "--json", str(output), "--csv-table", f"states={csv}"], (1,))
        assert output.read_text() == "preserve json" and csv.read_text() == "preserve csv", args
