"""Sutherland CLI contracts, including native fp128 parsing and formatting."""
import csv
from decimal import Decimal, localcontext
import io
import subprocess
import sys

program, fp128 = sys.argv[1:]


def run(*args, status=0):
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=30)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    if status == 0:
        assert not p.stderr, p.stderr
    return p.stdout, p.stderr


def rows(text, delimiter=","):
    assert "# CPU time:" in text and "References:" not in text
    lines = "\n".join(x for x in text.splitlines() if not x.startswith("#"))
    return list(csv.DictReader(io.StringIO(lines), delimiter=delimiter))


precisions = ["fp64", "long-double"]
if fp128.upper() in ("ON", "TRUE", "1"):
    precisions.append("fp128")
else:
    run("2", "--length", "4", "--lambda", "2", "--precision", "fp128", status=1)
with localcontext() as ctx:
    ctx.prec = 65
    pi = Decimal("3.1415926535897932384626433832795028841971693993751058209749445923")
    for precision in precisions:
        base = ["2", "--length", "4", "--lambda", "2", "--precision", precision, "--format", "csv"]
        ground = rows(run(*base)[0])
        assert len(ground) == 1 and ground[0]["labels"] == "0 0" and ground[0]["gap"] == "0"
        tol = Decimal("1e-31" if precision == "fp128" else "1e-13")
        assert abs(Decimal(ground[0]["energy"]) - pi * pi / 2) < tol, ground
        excited = rows(run(*base, "--labels", "0,1", "--pseudomomenta")[0])[0]
        assert excited["momentum_index"] == "1"
        assert abs(Decimal(excited["gap"]) - 3 * pi * pi / 4) < tol, excited
        assert abs(Decimal(excited["p"]) - pi / 2) < tol
        ks = [Decimal(k) for k in excited["pseudomomenta"].split()]
        assert len(ks) == 2 and abs(ks[0] + pi / 2) < tol and abs(ks[1] - pi) < tol
        if precision == "fp128":
            length = "4.000000000000000000000000002"
            coupling = "2.000000000000000000000000001"
            text = run("2", "--length", length, "--lambda", coupling, "--labels", "0,1",
                       "--precision", precision, "--format", "csv")[0]
            data = rows(text)[0]
            q2 = (2 * pi / Decimal(length)) ** 2
            expected = q2 * (Decimal(coupling) ** 2 / 2 + 1 + Decimal(coupling))
            assert abs(Decimal(data["energy"]) - expected) < tol, data
            metadata = dict(line[2:].split(": ", 1) for line in text.splitlines() if line.startswith("# "))
            assert abs(Decimal(metadata["Lambda"]) - Decimal(coupling)) < tol
            assert abs(Decimal(metadata["Length"]) - Decimal(length)) < tol

base = ["3", "--length", "4", "--lambda", "2", "--format", "csv"]
data = rows(run(*base, "--levels", "all", "--window", "2", "--max-states", "35")[0])
assert len(data) == 35 and data[0]["labels"] == "0 0 0"
assert len({r["labels"] for r in data}) == 35
assert len(rows(run(*base, "--levels", "4", "--window", "2")[0])) == 4
assert len(rows(run(*base, "--levels", "all", "--window", "0")[0])) == 1
assert not rows(run(*base, "--levels", "all", "--window", "2", "--max-states", "34", status=2)[0])
run("4", "--length", "4", "--lambda", "0.25", "--format", "pretty")
data = rows(run("2", "--length", "4", "--lambda", "0", "--labels", "-2,3", "--format", "tsv")[0], "\t")
assert data[0]["labels"] == "-2 3" and data[0]["momentum_index"] == "1"
data = rows(run("4", "--length", "4", "--lambda", "2", "--labels",
                "-9223372036854775808,-9223372036854775808,9223372036854775807,9223372036854775807",
                "--format", "csv")[0])
assert data[0]["momentum_index"] == "-2"
for options in [("--levels", "all"), ("--window", "2"), ("--labels", ""),
                ("--labels", "0,1"), ("--labels", "0,2,1"), ("--labels", "0,1,2,"),
                ("--labels", "0,1,1/2"), ("--labels", "0,0,9223372036854775808"),
                ("--labels", "0,1,2", "--levels", "all", "--window", "2"),
                ("--levels", "0", "--window", "2"), ("--levels", "all", "--window", "1000001"),
                ("--lambda", "-1"), ("--lambda", "nan"), ("--length", "0"), ("--length", "inf"),
                ("--format", "unknown"), ("--precision", "unknown"), ("--window",)]:
    run(*base, *options, status=1)
run("0", "--length", "4", "--lambda", "2", status=1)
run("1000001", "--length", "4", "--lambda", "2", status=1)
run("3", "--length", "4", status=1)
run("3", "--lambda", "2", status=1)
print("Sutherland CLI contracts passed")
