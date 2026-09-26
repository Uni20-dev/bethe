"""Every executable obeys the same informational and exact-option lifecycle."""
import os
from pathlib import Path
import subprocess
import sys

env = dict(os.environ, UNI20_COLOR="never", COLUMNS="4096")
base_arguments = {
    "bethe-xxz-dispersion": ["--delta", "2", "--points", "3"],
    "bethe-xyz-dispersion": ["--eta", "0.75", "--t", "1", "--points", "3"],
    "bethe-potts-pbc": ["4", "--momentum-index", "0"],
    "bethe-hubbard-continuum": ["--u", "4", "--momentum", "0"],
    "bethe-xxz-qg-obc": ["4", "--delta", "0.25"],
    "bethe-kondo-response": ["--field", "0", "--scale", "1"],
    "bethe-sine-gordon-vacuum": ["--length", "1", "--p", "1"],
    "bethe-lee-yang-vacuum": ["--length", "1"],
    "bethe-lee-yang-excited": ["--length", "5"],
    "bethe-asep-pbc": ["5", "--particles", "2", "--left-rate", "0.5"],
    "bethe-tasep-pbc": ["5", "--particles", "2"],
    "bethe-xyz-pbc": ["4", "--eta", "0.4", "--t", "0.7"],
    "bethe-xxz-pbc": ["4", "--delta", "0.5"],
    "bethe-xxz-obc": ["4", "--delta", "0.5"],
    "bethe-hubbard-pbc": ["4", "--u", "4"],
    "bethe-hubbard-obc": ["4", "--u", "4"],
    "bethe-hubbard-dispersion": ["--u", "4", "--points", "3"],
    "bethe-lieb-liniger-pbc": ["4", "--length", "4", "--c", "1"],
    "bethe-lieb-liniger-obc": ["4", "--length", "4", "--c", "1"],
    "bethe-lieb-liniger-dispersion": ["--c", "4", "--points", "3"],
    "bethe-lieb-liniger-thermal": ["--c", "4", "--temperature", "1", "--mu", "-1"],
    "bethe-q-boson-pbc": ["4", "--particles", "3", "--eta", "1"],
    "bethe-bose-fermi-pbc": ["--bosons", "2", "--fermions", "3", "--length", "5", "--c", "1"],
    "bethe-sutherland-pbc": ["4", "--length", "4", "--lambda", "2"],
    "bethe-su3-pbc": ["6"],
    "bethe-gaudin-yang-pbc": ["6", "--length", "6", "--c", "1"],
    "bethe-richardson": ["--levels", "0,1,2,3", "--pairs", "2", "--g", "1"],
    "bethe-central-spin": ["--couplings", "1,0.7,0.3", "--field", "1", "--sz", "0"],
    "bethe-sun-fermions-pbc": ["--populations", "3,1,1", "--length", "5", "--c", "1"],
    "bethe-ladder-pbc": ["4", "--rung", "1"],
}


def run(program, args, status=0):
    p = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=30)
    assert p.returncode == status, (program, args, p.returncode, p.stdout, p.stderr)
    if status == 0:
        assert not p.stderr, (program, args, p.stderr)
    return p


for program in sys.argv[1:]:
    name = Path(program).name
    base = base_arguments.get(name, ["4"])
    help_text = run(program, ["--help"]).stdout
    assert name in help_text and "Usage" in help_text
    assert "--references" in help_text and "Used for:" not in help_text
    assert run(program, ["-h"]).stdout == help_text
    empty = run(program, [], 1)
    assert not empty.stdout and empty.stderr == help_text
    assert run(program, ["--references", "--help"]).stdout == help_text
    for flag in ["--help", "--references", "--version", "--build-info"]:
        result = run(program, ["--precision", "invalid", "--max-iterations", "bad", flag]).stdout
        assert name in result
        if flag == "--references":
            assert "Used for:" in result
        else:
            assert "Used for:" not in result
    unknown = run(program, [*base, "--not-a-bethe-option"], 1)
    assert not unknown.stdout and "--not-a-bethe-option" in unknown.stderr
    duplicate = run(program, [*base, "--precision", "fp64", "--precision", "long-double"], 1)
    assert not duplicate.stdout and "--precision" in duplicate.stderr
    assert "received 2" in duplicate.stderr
print(f"Shared argument contracts passed for {len(sys.argv) - 1} frontends")
