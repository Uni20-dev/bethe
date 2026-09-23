"""Uni20 CLI integration: help lifecycle, exact tokens, exports, and diagnostics."""
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

program = str(Path(sys.argv[1]).resolve())
fp128 = sys.argv[2].upper() in ("ON", "TRUE", "1")


def run(*args, status=0, cwd=None, **environment):
    env = dict(os.environ, UNI20_COLOR="never", COLUMNS="160")
    env.update(environment)
    result = subprocess.run([program, *map(str, args)], capture_output=True, text=True,
                            timeout=60, cwd=cwd, env=env)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    if status == 0:
        assert not result.stderr, result.stderr
    return result


help_text = run("--help").stdout
assert run("-h").stdout == help_text
empty = run(status=1)
assert not empty.stdout and empty.stderr == help_text
for text in ("bethe-hubbard-dispersion", "Copyright", "GPL-3.0-or-later", "--references",
             "Examples", "default: 33", "--u REAL", "required", "--momentum", "excludes: --points",
             "--csv FILE", "--no-preamble", "--no-retain", "fp128", "MPLAPACK", "CITATIONS.md"):
    assert text in help_text, (text, help_text)
references = run("--references").stdout
assert "\nReferences\n" in references and "Used for:" in references
assert "\nReferences\n" not in help_text and "Used for:" not in help_text
assert "\nExamples\n" not in references and "--max-evaluations" not in references
# Ordinary help wins when both information flags are supplied, independent of order.
assert run("--references", "--help").stdout == help_text
assert run("--help", "--references").stdout == help_text
version = run("--version").stdout
assert re.fullmatch(r"bethe-hubbard-dispersion \d+\.\d+\.\d+ \([^\n]+\)", version.strip()), version
build = run("--build-info").stdout
for text in ("Compiler", "UNI20_BUILD_CLI", "CLI11"):
    assert text in build, (text, build)
assert "References" not in version + build

for width in (40, 16):
    narrow = run("--help", COLUMNS=str(width)).stdout
    for token in ("bethe-hubbard-dispersion", "--precision=long-double"):
        assert token in narrow, narrow
    assert "--u=4.000000000000000001" in narrow
assert "\x1b[" in run("--help", UNI20_COLOR="always").stdout
assert "\x1b[" in run("--references", UNI20_COLOR="always").stdout
assert "https://arxiv.org/abs/cond-mat/9808018" in run("--references", COLUMNS="16").stdout
assert "\x1b" not in help_text

if Path("/dev/full").exists():
    with open("/dev/full", "w") as full:
        failed = subprocess.run([program, "--references"], stdout=full, stderr=subprocess.PIPE,
                                text=True, timeout=10)
    assert failed.returncode == 1 and failed.stderr, (failed.returncode, failed.stderr)

with tempfile.TemporaryDirectory(prefix="bethe-arguments-") as directory:
    root = Path(directory)
    existing = root / "existing.csv"
    absent = root / "absent.json"
    existing.write_text("sentinel")
    # Information takes precedence over required parameters, invalid option values,
    # solver setup, and destructive --force output requests.
    for info in ("--help", "-h", "--version", "--build-info", "--references"):
        result = run("--csv", existing, "--json", absent, "--force", "--precision=invalid",
                     "--max-nodes=-1", "--density=nan", info)
        assert result.stdout and existing.read_text() == "sentinel" and not absent.exists()

    # Errors have a concise stderr report, never a bibliography or partial output.
    for args in (("--unknown",), ("--u",), ("--points=2",), ("--u=--help",), ("--", "--help"),
                 ("--u=--references",), ("--", "--references"),
                 ("--u=0",), ("--u=nan",), ("--u=4", "--density=0"),
                 ("--u=4", "--precision=invalid"), ("--u=4", "--points=-1"),
                 ("--u=4", "--points=18446744073709551616"), ("--u=4", "--points=2x"),
                 ("--u=4", "--points=1"), ("--u=4", "--points=2", "--momentum=0"),
                 ("--u=4", "--u=8"), ("--u=4", "--points=2", "--points=3"),
                 ("--u=4", "--precision=fp64", "--precision=long-double"),
                 ("--u=4", "--tolerance=1e-6", "--tolerance=1e-8")):
        result = run(*args, "--csv", existing, "--json", absent, "--force", status=1)
        assert not result.stdout and result.stderr and "References" not in result.stderr
        assert "--help" in result.stderr
        assert existing.read_text() == "sentinel" and not absent.exists()

    # Both spellings of an option value work, negative momenta stay values, and
    # each repeated export consumes exactly one path, including spaces and '='.
    first, second = root / "one = first.csv", root / "two.csv"
    result = run("--u=4", "--branch=holon", "--momentum=-1", "--format=json",
                 "--csv", first, "--csv=" + str(second), UNI20_COLOR="always")
    assert first.read_text() == second.read_text()
    doc = json.loads(result.stdout)
    assert len(doc["rows"]) == 1 and doc["metadata"]["Requested momentum"] == "-1"
    assert "\x1b" not in result.stdout and "References" not in result.stdout

    # A literal help token in an option value is a filename, not an information request.
    run("--u=4", "--branch=spinon", "--points=2", "--quiet", "--no-retain", "--csv=--help", cwd=root)
    assert (root / "--help").read_text().startswith("# ")
    run("--u=4", "--branch=spinon", "--points=2", "--quiet", "--csv=--references", cwd=root)
    assert (root / "--references").read_text().startswith("# ")
    run("--u=4", "--branch=spinon", "--points=2", "--quiet", "--csv=strict.csv",
        "--no-preamble", cwd=root)
    assert (root / "strict.csv").read_text().startswith("branch,p,")

if fp128:
    # Well below long-double resolution: either narrowing or order-dependent
    # conversion would erase the nonzero difference from four.
    token = "4.000000000000000000000000000001"
    with localcontext() as context:
        context.prec = 60
        for arguments in (("--precision=fp128", "--u=" + token),
                          ("--u=" + token, "--precision=fp128")):
            result = run(*arguments, "--branch=spinon", "--momentum=0", "--format=json")
            value = Decimal(json.loads(result.stdout)["metadata"]["U (t=1)"])
            assert value > Decimal(4) and abs(value - Decimal(token)) < Decimal("1e-32"), value

print("Hubbard Uni20 argument/help contracts passed")
