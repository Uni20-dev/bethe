"""Excitation frontends: scheduling changes execution, not table contents."""
from decimal import Decimal
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")

def run(program, args, code=0):
    p = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=60)
    assert p.returncode == code, (program, args, p.stdout, p.stderr)
    return p.stdout

cases = {
    "bethe-xxx-pbc": (["10"], []),
    "bethe-xxx-obc": (["10"], []),
    "bethe-xxz-pbc": (["10", "--delta", "0.7"], []),
    "bethe-xxz-obc": (["10", "--delta", "0.7"], []),
    "bethe-lieb-liniger-pbc": (["3", "--length", "4", "--c", "2"], ["--padding", "1"]),
    "bethe-lieb-liniger-obc": (["3", "--length", "4", "--c", "2"], ["--padding", "3"]),
    "bethe-q-boson-pbc": (["4", "--particles", "3", "--eta", "1"], []),
    "bethe-xxz-qg-obc": (["8", "--delta", "0.6", "--through-lines", "4"], []),
    "bethe-biquadratic-obc": (["8", "--through-lines", "2"], []),
}


def tables(program, args, code=0):
    return json.loads(run(program, [*args, "--roots", "--format", "json"], code),
                      parse_float=Decimal)["tables"]


def same_tables(program, args, code=0):
    serial = tables(program, [*args, "--threads", "1"], code)
    for threads in (2, 4):
        parallel = tables(program, [*args, "--threads", str(threads)], code)
        assert serial.keys() == parallel.keys()
        for name in serial:
            assert serial[name]["columns"] == parallel[name]["columns"], (program, name)
            assert serial[name]["rows"] == parallel[name]["rows"], (program, name, args)
            assert serial[name]["metadata"]["Scheduler concurrency limit"] == "1"
            assert parallel[name]["metadata"]["Scheduler concurrency limit"] == str(threads)
        # Retention must not change any streamed row or partial-failure record.
        streamed = tables(program, [*args, "--threads", str(threads), "--no-retain"], code)
        assert all(streamed[name]["rows"] == parallel[name]["rows"] for name in parallel)
    return serial


for program in sys.argv[1:]:
    name = Path(program).stem
    base, scan_options = cases[name]
    assert "--threads" in run(program, ["--help"])
    args = [*base, *scan_options, "--excitations", "all"]
    full = same_tables(program, args)
    implicit = tables(program, args)  # default is a recorded serial limit
    assert all(implicit[key]["rows"] == full[key]["rows"] for key in full)
    assert all(t["metadata"]["Scheduler concurrency limit"] == "1" for t in implicit.values())
    for count in ("1", "3"):
        limited = same_tables(program, [*base, *scan_options, "--excitations", count])
        primary = "levels" if name == "bethe-xxz-qg-obc" else "states"
        assert limited[primary]["rows"] == full[primary]["rows"][:int(count)]
    same_tables(program, [*args, "--max-iterations", "0"], 2)
    # Argument rejection must happen before overwriting an export, even with --force.
    with tempfile.TemporaryDirectory() as tmp:
        output = Path(tmp) / "keep.json"
        output.write_text("do not overwrite\n")
        invalid = [[*base, "--threads", str(t)] for t in (1, 4)]
        invalid += [[*args, "--threads", t] for t in ("0", "-1", "1.5", "invalid")]
        invalid += [[*args, "--threads", "4", "--max-candidates", "1"]]
        for bad in invalid:
            run(program, [*bad, "--json", str(output), "--force"], 1)
            assert output.read_text() == "do not overwrite\n"

    if name.startswith("bethe-lieb-liniger-"):
        same_tables(program, ["0", *base[1:], "--padding", "2", "--excitations", "all"])
    elif name == "bethe-q-boson-pbc":
        for deformation in (["--eta", "0"], ["--phase"]):
            same_tables(program, ["4", "--particles", "3", *deformation, "--excitations", "all"])
        same_tables(program, ["4", "--particles", "0", "--eta", "1", "--excitations", "all"])
    elif name == "bethe-xxz-qg-obc":
        same_tables(program, ["8", "--delta", "0.6", "--through-lines", "8", "--excitations", "all"])
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "keep.json"
            output.write_text("do not overwrite\n")
            run(program, ["8", "--delta", "0", "--excitations", "all", "--threads", "4",
                          "--json", str(output), "--force"], 1)
            assert output.read_text() == "do not overwrite\n"
    elif name == "bethe-biquadratic-obc":
        for window in ([], ["--real-window", "3"]):
            for extra, code in (([], 0), (["--max-iterations", "0"], 2)):
                same_tables(program, ["9", "--ferromagnetic", "--through-lines", "5",
                                      "--excitations", "3", *window, *extra], code)
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "keep.json"
            output.write_text("do not overwrite\n")
            for mode in (["--q-spectrum"], ["--singlet-excitation"],
                         ["--ferromagnetic", "--one-defect"],
                         ["--ferromagnetic", "--bound-pairs", "all"]):
                for extra in ([], ["--excitations", "all"]):
                    run(program, ["8", *mode, *extra, "--threads", "4",
                                  "--json", str(output), "--force"], 1)
                    assert output.read_text() == "do not overwrite\n"

print("Parallel excitation CLI contracts passed")
