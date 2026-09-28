"""XXX two-spinon CLI: normalization, completeness diagnostics and exports."""
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

def run(args, code=0):
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=180, env=env)
    assert p.returncode == code, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout

def document(args, code=0):
    data = json.loads(run([*args, "--format", "json"], code))
    assert data["status"] == "complete"
    assert data["tables"]["spectrum"]["summary"]["Outcome"] == ("success" if code == 0 else "partial")
    return data["tables"]

def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]

for precision in precisions:
    args = ["8", "--precision", precision]
    z = document(args)
    parallel = document([*args, "--threads", "4"])
    for name in z:
        assert parallel[name]["rows"] == z[name]["rows"]
    r = document([*args, "--channel", "raising"])
    assert set(z) == {"spectrum", "moments"}
    assert len(rows(z["spectrum"])) == 10
    with localcontext() as ctx:
        ctx.prec = 100
        tolerance = {"fp64": Decimal("1e-12"), "long-double": Decimal("1e-15"), "fp128": Decimal("1e-29")}[precision]
        for zz, raising in zip(rows(z["spectrum"]), rows(r["spectrum"])):
            assert zz["state_id"] == raising["state_id"]
            assert abs(2*Decimal(zz["weight"])-Decimal(raising["weight"])) < tolerance
        weight = sum(Decimal(x["weight"]) for x in rows(z["spectrum"]))/8
        assert Decimal(".249") < weight < Decimal(".25")
    for table in z.values():
        assert Decimal(table["summary"]["Run CPU seconds"]) >= 0
        assert "partial DSF" in table["metadata"]["Family"]
    moment0 = rows(z["moments"])[0]
    assert Decimal(moment0["weight"]) == Decimal(moment0["full_first_moment"]) == 0
    assert moment0["first_moment_fraction"] is None
    failed = document([*args, "--max-iterations", "0", "--diagnostics"], 2)
    assert not rows(failed["spectrum"])
    assert any(not row["roots_converged"] for row in rows(failed["diagnostics"]))
    assert all(row["full_first_moment"] is None for row in rows(failed["moments"]))
    # A root solve accepted under an intentionally loose tolerance must still
    # fail the independent on-shell gate in the form-factor layer.
    loose = document([*args, "--tolerance", "0.1", "--diagnostics"], 2)
    assert not rows(loose["spectrum"])
    assert any(row["form_factor_status"] == "precision_limit" for row in rows(loose["diagnostics"]))

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)
    screen = run(["6", "--csv-table", f"roots={path/'roots.csv'}", "--format", "plain"])
    assert "Bethe roots (ground=0" not in screen
    text = (path/"roots.csv").read_text()
    roots = list(csv.DictReader(line for line in text.splitlines() if not line.startswith("#")))
    assert len(roots) == 15
    assert "1/2" not in "\n".join(line for line in text.splitlines() if not line.startswith("#"))
    out = path/"spectrum.json"
    run(["6", "--json", str(out), "--quiet", "--no-retain"])
    previous = out.read_text()
    invalid = [["8", "--threads", "0"], ["8", "--threads", "-1"], ["8", "--threads", "1.5"],
               ["8", "--threads", "2147483648"],
               ["5"], ["0"], ["8", "--max-candidates", "9"], ["8", "--channel", "+-"],
               ["8", "--csv-table", f"failed={path/'bad.csv'}"], ["8", "--table", "failed"],
               ["8", "--tolerance", "nan"], ["8", "--tolerance", "0"]]
    for args in invalid:
        run([*args, "--json", str(out), "--force"], 1)
        assert out.read_text() == previous
    assert not (path/"bad.csv").exists()
    for flag in ("--csv", "--tsv"):
        run(["6", flag, str(path/flag[2:]), "--stream", "--no-retain", "--quiet"])
    run(["6", "--format", "pretty"])
    help_text = run(["--help"])
    for table in ("spectrum", "moments", "diagnostics", "roots"):
        assert table in help_text
    assert "Used for:" not in help_text
    assert "Caux" in run(["--references"])

print("XXX structure factor CLI contracts passed")
