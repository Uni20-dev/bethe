"""Lee-Yang one-particle levels, vacuum-relative gaps and partial output."""
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
base = ["--length", "5"]


def run(args, status=0):
    result = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=60)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    if status == 0:
        assert not result.stderr
    if status == 2:
        assert "gap unavailable" in result.stderr
    return result.stdout


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]


def result(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete"
    tables = doc["tables"]
    assert list(tables) == ["levels", "source", "gap"]
    for table in tables.values():
        meta, summary = table["metadata"], table["summary"]
        assert "bulk-subtracted" in meta["Energy convention"]
        assert "E1_C-E0_C" in meta["Gap convention"]
        assert "no UV continuation" in meta["Scope"]
        assert Decimal(summary["Run CPU seconds"]) >= 0
        assert summary["Outcome"] == ("success" if status == 0 else "partial")
        for row in rows(table):
            for key, value in row.items():
                if key not in ("state", "status", "converged") and value is not None:
                    assert Decimal(value).is_finite(), (key, value)
    assert len(tables["levels"]["rows"]) == 2
    assert len(tables["source"]["rows"]) == len(tables["gap"]["rows"]) == 1
    return tables


for precision in precisions:
    flags = ["--precision", precision]
    tables = result([*base, *flags])
    vacuum, excited = rows(tables["levels"])
    gap = rows(tables["gap"])[0]
    source = rows(tables["source"])[0]
    assert vacuum["state"] == "vacuum" and excited["state"] == "one_particle"
    assert vacuum["converged"] and excited["converged"] and gap["converged"] and source["converged"]
    assert "effective_central_charge" not in excited
    with localcontext() as ctx:
        ctx.prec = 80
        tol = Decimal(tables["levels"]["metadata"]["Absolute Y tolerance per state"])
        y0, y1 = Decimal(vacuum["scaling_function"]), Decimal(excited["scaling_function"])
        assert abs(y1-Decimal("5.146781165267865")) < Decimal("2e-12")
        assert abs(Decimal(gap["gap"])-Decimal("1.0306379212700556")) < Decimal("2e-12")
        assert abs(Decimal(gap["scaled_gap"])-(y1-y0)) < tol
        assert abs(Decimal(gap["gap"])-(y1-y0)/5) < tol
        assert Decimal(gap["gap_error"]) >= (Decimal(vacuum["cutoff_error"])+Decimal(excited["cutoff_error"]))/5
        assert Decimal(gap["gap_error"]) < tol
        assert Decimal(source["quantization_residual"]) <= tol
    # Loose but explicit native tolerance keeps repeated fp128 CLI checks cheap;
    # default native precision is exercised above and in the C++ mesh tests.
    flags += ["--tolerance", "1e-10"]
    scaled = result(["--mass", "2", "--length", "2.5", *flags])
    assert abs(Decimal(rows(scaled["gap"])[0]["gap"])-2*Decimal("1.0306379212700556")) < Decimal("1e-10")
    partial = result([*base, *flags, "--max-root-iterations", "0"], 2)
    v, e = rows(partial["levels"])
    assert v["converged"] and v["casimir_energy"] is not None
    assert not e["converged"] and e["casimir_energy"] is None and e["scaling_function"] is None
    assert e["status"] == "iteration_limit"
    s = rows(partial["source"])[0]
    assert s["beta"] is None and s["pole_displacement"] is None
    g = rows(partial["gap"])[0]
    assert all(g[k] is None for k in ("gap", "scaled_gap", "gap_error"))
    assert g["status"] == "excited_iteration_limit"
    for controls, status in [(["--max-iterations", "0"], "iteration_limit"),
                             (["--max-intervals", "32"], "mesh_limit"),
                             (["--max-cutoffs", "0"], "cutoff_limit"),
                             (["--max-kernel-products", "0"], "work_limit")]:
        failed = result([*base, *flags, *controls], 2)
        for level in rows(failed["levels"]):
            assert not level["converged"] and level["casimir_energy"] is None
            assert level["status"] == status
        assert rows(failed["gap"])[0]["gap"] is None
    failed = result([*base, "--precision", precision, "--tolerance", "1e-40"], 2)
    assert rows(failed["levels"])[1]["status"] == "precision_limit"

# Both levels are representable, but the slightly larger vacuum-relative
# gap is not. This must not leak an infinity or claim scientific success.
overflow = result(["--mass", "1.7745e308", "--length", repr(6/1.7745e308)], 2)
assert all(row["converged"] for row in rows(overflow["levels"]))
assert rows(overflow["source"])[0]["converged"]
g = rows(overflow["gap"])[0]
assert g["status"] == "gap_precision_limit"
assert all(g[key] is None for key in ("gap", "scaled_gap", "gap_error"))

invalid = [[], ["--length", "4"], ["--length", "31"], ["--length", "nan"],
           [*base, "--mass", "0"], [*base, "--mass", "inf"], [*base, "--mass", "bad"],
           [*base, "--initial-cutoff", "0"], [*base, "--initial-intervals", "1"],
           [*base, "--max-intervals", "16384"], [*base, "--max-root-iterations", "-1"],
           [*base, "--tolerance", "0"], [*base, "--tolerance", "nan"]]
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    flags = ["--json", str(root/"result.json")]
    for suffix in ("csv", "tsv"):
        for name in ("levels", "source", "gap"):
            flags += [f"--{suffix}-table", f"{name}={root/f'{name}.{suffix}'}"]
    screen = run([*base, *flags, "--format", "plain"])
    assert "CPU time" in screen and "Used for:" not in screen
    expected = json.loads((root/"result.json").read_text())["tables"]
    for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
        for name, table in expected.items():
            text = (root/f"{name}.{suffix}").read_text()
            assert "bulk-subtracted" in text
            records = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
            want = [{k: str(v).lower() if isinstance(v, bool) else str(v) for k,v in row.items()} for row in rows(table)]
            assert records == want
    old = (root/"result.json").read_text()
    for args in invalid:
        run([*args, "--json", str(root/"result.json"), "--force"], 1)
        assert (root/"result.json").read_text() == old
    run([*base, *flags, "--force", "--quiet", "--no-retain"])
    actual = json.loads((root/"result.json").read_text())["tables"]
    assert all(actual[name]["rows"] == expected[name]["rows"] for name in actual)
    run([*base, "--stream", "--no-retain", "--format", "plain"])
    run([*base, "--format", "pretty"])
    run([*base, "--max-root-iterations", "0", "--quiet", "--csv-table", f"gap={root/'failed.csv'}"], 2)
    text = (root/"failed.csv").read_text()
    row = next(csv.DictReader(line for line in text.splitlines() if not line.startswith("#")))
    assert row["gap"] == row["scaled_gap"] == row["gap_error"] == ""

print("Lee-Yang excited CLI contracts passed")
