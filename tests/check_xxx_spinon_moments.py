"""Fresh direct moments, physical normalization, failure and export contracts."""
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, mpfr = sys.argv[1:]
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")

def invoke(args, expected=0):
    p = subprocess.run([program, *map(str, args)], capture_output=True, text=True, env=env, timeout=240)
    if p.returncode != expected:
        detail = p.stdout[-4000:]
        try:
            doc = json.loads(p.stdout)
            detail = {name: rows(doc, name) for name in ("moments", "diagnostics") if name in doc["tables"]}
        except (ValueError, KeyError, TypeError):
            pass
        raise AssertionError((args, p.returncode, detail, p.stderr))
    return p.stdout

def document(args, expected=0):
    return json.loads(invoke([*args, "--format", "json"], expected))

def rows(doc, table="moments"):
    t = doc["tables"][table]
    names = [x["id"] for x in t["columns"]]
    return [dict(zip(names, [x["value"] if isinstance(x, dict) else x for x in row])) for row in t["rows"]]

def number(x): return Decimal(str(x))
def close(a, b, tolerance="3e-13"):
    with localcontext() as c:
        c.prec = 130
        a, b = number(a), number(b)
        assert abs(a-b) <= Decimal(tolerance)*max(abs(a), abs(b), Decimal("1e-100")), (a,b)

def indexed(doc): return {(r["sector"], int(r["moment"])):r for r in rows(doc)}
def sums(doc, sectors):
    r = indexed(doc)
    label = "sum:" + "+".join(map(str, sectors))
    for moment in (0,1):
        assert all(r[(str(n),moment)]["value"] is not None for n in sectors)
        with localcontext() as c:
            c.prec = 130
            close(r[(label,moment)]["value"], sum(number(r[(str(n),moment)]["value"]) for n in sectors))

# Refinement must reject bad coordinates and deferred controls before --force
# can open either whole-document or named-table exports. Select only sector 4:
# another sector's eager validation must not conceal a deferred-query bug.
with tempfile.TemporaryDirectory() as tmp:
    protected = Path(tmp)/"protected.json"
    table = Path(tmp)/"protected.csv"
    absent = Path(tmp)/"not-created.tsv"
    base = ["--sectors", "4", "--precision", "long-double", "--four-order", "2", "--four-refine-order", "3"]
    fixed_controls = ["--four-cutoff", "2", "--four-extended-cutoff", "3"]
    global_controls = ["--global", "--four-gap-cutoffs", "1,1,1", "--four-extended-gap-cutoffs", "2,2,2"]
    bad_requests = [["--momentum", q, *fixed_controls] for q in ("-1", "7", "nan", "inf")]
    bad_requests += [["--momentum", "1", *fixed_controls, "--four-max-nodes", "0"],
                     [*global_controls, "--four-max-nodes", "0"],
                     [*global_controls, "--four-max-nodes", "20"],
                     ["--global", "--four-gap-cutoffs", "1,1,1", "--four-extended-gap-cutoffs", "8193,2,2"],
                     ["--global", "--four-gap-cutoffs", "1,1,1", "--four-extended-gap-cutoffs", "3000,3000,3000"]]
    for bad in bad_requests:
        protected.write_bytes(b"preserve prior JSON\n")
        table.write_bytes(b"preserve prior CSV\n")
        invoke([*base, *bad, "--json", protected, "--csv-table", f"moments={table}",
                "--tsv-table", f"diagnostics={absent}", "--force"], 1)
        assert protected.read_bytes() == b"preserve prior JSON\n", bad
        assert table.read_bytes() == b"preserve prior CSV\n", bad
        assert not absent.exists(), bad

g = document(["--global", "--diagnostics"])
r = indexed(g)
assert set(r) == {("2",0), ("2",1)}
close(r[("2",0)]["sum_rule"], ".25")
for moment in (0,1):
    close(r[("2",moment)]["fraction"], number(r[("2",moment)]["value"])/number(r[("2",moment)]["sum_rule"]))
    assert r[("2",moment)]["status"] == "estimated"
assert 0 < number(r[("2",0)]["fraction"]) < 1
assert rows(g,"diagnostics")[0]["sampling_error"] is None

fixed = ["--momentum","1.5","--sectors","2,4","--precision","long-double",
         "--two-order","4","--two-gap-cutoff","4","--four-order","2","--four-cutoff","2","--diagnostics",
         "--density-relative-tolerance","1e-4","--density-absolute-tolerance","1e-24"]
f = document(fixed); sums(f,(2,4))
assert all(row["sum_rule"] is None for row in rows(f) if int(row["moment"])==0)
scaled = document([*fixed,"--exchange","3","--channel","raising","--threads","2"])
fr,sr=indexed(f),indexed(scaled)
for key,row in fr.items():
    close(sr[key]["value"], (2 if key[1]==0 else 6)*number(row["value"]))
    if row["fraction"] is not None: close(sr[key]["fraction"],row["fraction"])

global_four=["--global","--sectors","2,4","--precision","long-double","--two-order","4",
             "--two-gap-cutoff","4","--four-order","2","--four-gap-cutoffs","1,1,1","--diagnostics",
             "--density-relative-tolerance","1e-4","--density-absolute-tolerance","1e-24"]
sums(document(global_four),(2,4))
refined=document([*global_four,"--four-refine-order","3","--four-extended-gap-cutoffs","2,2,2"])
diag=[r for r in rows(refined,"diagnostics") if r["sector"]=="4"]
assert all(r["quadrature_difference"] is not None and r["cutoff_change"] is not None for r in diag)
failed=document([*fixed,"--four-max-nodes","1"],2)
assert all(row["value"] is None for row in rows(failed) if row["sector"] in ("4","sum:2+4"))
zero=document(["--momentum","0","--sectors","2,4,6","--precision","long-double","--replicate-output","--covariance"])
assert all(number(r["value"])==0 and r["fraction"] is None for r in rows(zero))
assert rows(zero,"replicates")==[]
assert all(number(r["covariance"])==0 for r in rows(zero,"covariance"))

with tempfile.TemporaryDirectory() as tmp:
    path=Path(tmp)/"protected.json";path.write_text("keep this file\n");original=path.read_bytes()
    for bad in [["--sectors","2,2"],["--csv-table","unknown=x.csv"],["--csv-table","replicates=x.csv"],
                ["--two-gap-cutoff","-1"],["--four-refine-order","4"],["--covariance"]]:
        invoke(["--global",*bad,"--json",path,"--force"],1)
        assert path.read_bytes()==original
    invoke(["--global","--momentum","1.5"],1)
    invoke(["--sectors","2"],1)
    csv=Path(tmp)/"checks.csv"
    invoke([*fixed,"--quiet","--no-retain","--csv-table",f"diagnostics={csv}"])
    assert "shift_tail_bound" in csv.read_text()

if mpfr.upper() in ("ON","TRUE","1"):
    # This fixture checks 3e-13 normalization/covariance identities, not precision
    # convergence. 192-bit geometry retains these observations; an explicit
    # kernel tolerance avoids requesting 4096*epsilon accuracy for sector 2.
    mpfr_controls=["--precision","bits-192","--tolerance","1e-20"]
    args=["--momentum","1.57079632679489661923132169163975","--sectors","2,4,6",*mpfr_controls,
          "--two-order","4","--two-gap-cutoff","4","--four-order","2","--four-cutoff","2",
          "--replicates","2","--samples","4","--seed","42","--density-relative-tolerance","1e-5",
          "--density-absolute-tolerance","1e-30",
          "--threads","2","--diagnostics","--replicate-output","--covariance"]
    doc=document(args);sums(doc,(2,4,6))
    assert len(rows(doc,"replicates"))==2
    assert all(r["failed_sample"] is None for r in rows(doc,"replicates"))
    assert number(indexed(doc)[("6",0)]["value"])>0
    covariance=rows(doc,"covariance")
    assert len(covariance)==4 and number(covariance[1]["covariance"])>0
    raised=document([*args,"--channel","raising"])
    for a,b in zip(covariance,rows(raised,"covariance")): close(b["covariance"],4*number(a["covariance"]))
    for a,b in zip(rows(doc,"replicates"),rows(raised,"replicates")):
        close(b["m0"],2*number(a["m0"]));close(b["m1"],2*number(a["m1"]))
    glob=document(["--global","--sectors","2,4,6",*mpfr_controls,"--two-order","4","--two-gap-cutoff","4",
                   "--four-order","2","--four-gap-cutoffs","1,1,1","--replicates","2","--samples","4","--seed","42",
                   "--density-relative-tolerance","1e-5","--density-absolute-tolerance","1e-30",
                   "--threads","2","--covariance"])
    sums(glob,(2,4,6))
print("Direct spinon moments, normalization, covariance and export contracts passed")
