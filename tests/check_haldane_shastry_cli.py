"""Native scalar output and exact spectral-rule CLI contracts."""
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

def rows(text):
    assert "# CPU time:" in text and "References:" not in text
    return list(csv.DictReader(io.StringIO("\n".join(x for x in text.splitlines() if not x.startswith("#")))) )

precisions = ["fp64", "long-double"]
if fp128.upper() in ("ON", "TRUE", "1"):
    precisions.append("fp128")
else:
    run("4", "--precision", "fp128", status=1)
with localcontext() as ctx:
    ctx.prec = 60
    for precision in precisions:
        data = rows(run("4", "--precision", precision, "--format", "csv")[0])
        assert len(data) == 1 and data[0]["motif"] == "1 3" and data[0]["gap"] == "0"
        ref = Decimal("-2.15897596273829719787004490622290806085")
        tol = Decimal("1e-32" if precision == "fp128" else "1e-14")
        assert abs(Decimal(data[0]["energy"]) - ref) < tol, data
        data = rows(run("5", "--precision", precision, "--format", "csv")[0])
        assert len(data) == 2 and {r["momentum_index"] for r in data} == {"1", "4"}
        assert all(r["s_max"] == "1/2" and r["degeneracy"] == "2" for r in data)
data = rows(run("6", "--levels", "all", "--format", "csv")[0])
assert len(data) == 13 and sum(int(r["degeneracy"]) for r in data) == 64
assert len(rows(run("6", "--levels", "3", "--format", "csv")[0])) == 3
assert rows(run("4", "--motif", "", "--format", "csv")[0])[0]["degeneracy"] == "5"
assert len(rows(run("6", "--levels", "all", "--max-motifs", "12", "--format", "csv", status=2)[0])) == 0
run("8", "--sz", "-1", "--format", "pretty")
for args in [("1",), ("1000001",), ("6", "--sz", "1/2"), ("6", "--motif", "1,2"),
             ("6", "--motif", "6"), ("6", "--levels", "0"), ("6", "--levels", "all", "--sz", "0"),
             ("6", "--format", "unknown"), ("6", "--motif", "1,"), ("6", "--levels")]:
    run(*args, status=1)
print("Haldane-Shastry CLI contracts passed")
