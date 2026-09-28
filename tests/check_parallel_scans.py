"""XXX/XXZ frontends: scheduling changes execution, not table contents."""
import json
import os
from pathlib import Path
import subprocess
import sys

env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")

def run(program, args, code=0):
    p = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=60)
    assert p.returncode == code, (program, args, p.stdout, p.stderr)
    return p.stdout

for program in sys.argv[1:]:
    base = ["10"]
    if "xxz" in Path(program).name:
        base += ["--delta", "0.7"]
    assert "--threads" in run(program, ["--help"])
    run(program, [*base, "--threads", "2"], 1)  # no ignored option outside a scan
    for extra, code in [([], 0), (["--max-iterations", "0"], 2)]:
        args = [*base, "--excitations", "all", "--roots", "--format", "json", *extra]
        serial = json.loads(run(program, [*args, "--threads", "1"], code))["tables"]
        parallel = json.loads(run(program, [*args, "--threads", "4"], code))["tables"]
        assert serial.keys() == parallel.keys()
        for name in serial:
            assert serial[name]["rows"] == parallel[name]["rows"], (program, name)
            assert parallel[name]["metadata"]["Scheduler concurrency limit"] == "4"
print("Parallel spin-chain CLI contracts passed")
