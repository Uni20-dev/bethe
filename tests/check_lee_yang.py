"""Periodic Lee-Yang frontend: conventions, precision, failures and exports."""
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
base = ["--length", "1"]
observables = ("casimir_energy", "scaling_function", "effective_central_charge")
errors = ("nonlinear_residual", "nonlinear_error", "mesh_error", "cutoff_error", "direct_tail_bound")


def run(args, status=0):
    result = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=30)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    if status == 0:
        assert not result.stderr
    return result.stdout


def record(table):
    assert len(table["rows"]) == 1
    return dict(zip([c["id"] for c in table["columns"]], table["rows"][0]))


def result(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete" and list(doc["tables"]) == ["vacuum"]
    table = doc["tables"]["vacuum"]
    meta, summary = table["metadata"], table["summary"]
    assert "bulk-subtracted" in meta["Energy convention"]
    assert "not the CFT central charge" in meta["Effective charge convention"]
    assert meta["CFT central charge (theory)"] == "-22/5"
    assert meta["Lowest conformal weight (theory)"] == "-1/5"
    assert meta["UV effective charge (theory)"] == "2/5"
    assert Decimal(summary["Run CPU seconds"]) >= 0
    assert summary["Outcome"] == ("success" if status == 0 else "partial")
    row = record(table)
    assert row["converged"] == (status == 0)
    for key in (*observables, *errors):
        assert row[key] is None or Decimal(row[key]).is_finite()
    if status:
        assert all(row[key] is None for key in observables)
    return table, row


failures = [(["--max-iterations", "0"], "iteration_limit"),
            (["--max-intervals", "32"], "mesh_limit"),
            (["--max-cutoffs", "0"], "cutoff_limit"),
            (["--max-cutoffs", "1"], "cutoff_limit"),
            (["--max-kernel-products", "0"], "work_limit"),
            (["--tolerance", "1e-40"], "precision_limit")]
for precision in precisions:
    flags = ["--precision", precision]
    table, row = result([*base, *flags])
    assert table["metadata"]["Precision"] == precision and row["status"] == "converged"
    with localcontext() as ctx:
        ctx.prec = 80
        y = Decimal(row["scaling_function"])
        assert abs(y-Decimal("-0.15320688011013006")) < Decimal("1e-11")
        assert Decimal(row["casimir_energy"]) == y
        pi = Decimal("3.1415926535897932384626433832795028841971693993751058209749445923078164")
        tol = Decimal(table["metadata"]["Absolute Y tolerance"])
        assert abs(Decimal(row["effective_central_charge"])+6*y/pi) < tol
        scaled = result(["--mass", "2", "--length", "0.5", *flags])[1]
        assert scaled["scaling_function"] == row["scaling_function"]
        assert scaled["effective_central_charge"] == row["effective_central_charge"]
        assert abs(Decimal(scaled["casimir_energy"])-2*y) < tol
    for controls, status in failures:
        failed = result([*base, *flags, *controls], 2)[1]
        assert failed["status"] == status
        if status == "work_limit":
            assert int(failed["kernel_products"]) == 0
            assert all(failed[key] is None for key in errors)
    uv = result(["--length", "1e-5", "--tolerance", "1e-10", *flags])[1]
    assert abs(Decimal(uv["effective_central_charge"])-Decimal("0.4")) < Decimal("1e-8")

invalid = [[], ["--length", "0"], ["--length", "-1"], ["--length", "nan"], ["--length", "inf"],
           [*base, "--mass", "0"], [*base, "--mass", "-1"], [*base, "--mass", "bad"],
           [*base, "--tolerance", "0"], [*base, "--tolerance", "1"], [*base, "--tolerance", "nan"],
           [*base, "--initial-cutoff", "0"], [*base, "--initial-cutoff", "inf"],
           [*base, "--initial-intervals", "1"], [*base, "--max-intervals", "16384"],
           [*base, "--max-intervals", "16"], [*base, "--max-iterations", "-1"],
           [*base, "--max-kernel-products", "-1"], [*base, "--p", "1"]]
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)
    flags = ["--json", str(path/"vacuum.json"), "--csv", str(path/"vacuum.csv"), "--tsv", str(path/"vacuum.tsv")]
    plain = run([*base, *flags, "--format", "plain"])
    assert "CPU time" in plain and "Used for:" not in plain
    table = json.loads((path/"vacuum.json").read_text())["tables"]["vacuum"]
    expected = {k: str(v).lower() if isinstance(v, bool) else str(v) for k,v in record(table).items()}
    for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
        content = (path/f"vacuum.{suffix}").read_text()
        assert "bulk-subtracted" in content
        rows = list(csv.DictReader((line for line in content.splitlines() if not line.startswith("#")), delimiter=delimiter))
        assert rows == [expected]
    old = (path/"vacuum.json").read_text()
    for args in invalid:
        run([*args, "--json", str(path/"vacuum.json"), "--force"], 1)
        assert (path/"vacuum.json").read_text() == old
    run([*base, *flags, "--quiet", "--force", "--no-retain"])
    assert json.loads((path/"vacuum.json").read_text())["tables"]["vacuum"]["rows"] == table["rows"]
    run([*base, "--stream", "--no-retain", "--format", "plain"])
    run([*base, "--format", "pretty"])
    text = run([*base, "--max-kernel-products", "0", "--format", "csv", "--no-preamble"], 2)
    row = next(csv.DictReader(text.splitlines()))
    assert all(row[key] == "" for key in (*observables, *errors))

print("Lee-Yang CLI contracts passed")
