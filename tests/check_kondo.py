"""Kondo CLI: native response, budget forwarding and shared export contracts."""
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
base = ["--field", "2", "--scale", "1"]


def run(args, status=0):
    result = subprocess.run([program, *args], text=True, capture_output=True, timeout=120, env=env)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    return result.stdout


def row(table):
    return dict(zip([c["id"] for c in table["columns"]], table["rows"][0]))


def table(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete" and list(doc["tables"]) == ["response"]
    result = doc["tables"]["response"]
    assert result["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    assert Decimal(result["summary"]["Run CPU seconds"]) >= 0
    assert "E_imp(b)-E_imp(0)" in result["metadata"]["Energy convention"]
    assert "T_B=2*T1" in result["metadata"]["Scale convention"]
    return result


for precision in precisions:
    flags = ["--precision", precision]
    with localcontext() as ctx:
        ctx.prec = 100
        m = Decimal("0.27516363791067304687136724858153892615232627968484052219868642339")
        e = Decimal("-0.34629170027457867780906807893015178590470848543255021208916271878")
        tol = {"fp64": Decimal("3e-10"), "long-double": Decimal("2e-13"), "fp128": Decimal("3e-28")}[precision]
        positive = row(table([*base, *flags]))
        negative = row(table(["--field", "-8", "--scale", "4", *flags]))
        assert positive["converged"] and positive["status"] == "converged"
        assert abs(Decimal(positive["magnetization"])-m) < tol
        assert abs(Decimal(positive["energy_change"])-e) < 2*tol
        assert abs(Decimal(negative["magnetization"])+m) < tol
        assert abs(Decimal(negative["energy_change"])/4-e) < 2*tol
        assert abs(Decimal(negative["zero_field_susceptibility"])*4
                   - Decimal(positive["zero_field_susceptibility"])) < tol
    custom = table([*base, *flags, "--tolerance", "1e-6"])
    assert abs(Decimal(custom["metadata"]["Absolute M and Delta E/|b| tolerance"])
               - Decimal("1e-6")) < Decimal("1e-21")
    zero = row(table(["--field", "0", "--scale", "1", "--max-series-terms", "0", *flags]))
    assert Decimal(zero["energy_change"]) == Decimal(zero["magnetization"]) == 0
    assert int(zero["series_terms"]) == int(zero["evaluations"]) == 0
    for option, expected in (("--max-series-terms", "series_limit"), ("--max-lobes", "tail_limit"),
                             ("--max-evaluations", "quadrature_limit"),
                             ("--max-quadrature-levels", "quadrature_limit")):
        r = row(table([*base, *flags, option, "0"], 2))
        assert not r["converged"] and r["status"] == expected
        assert all(r[k] is None for k in ("energy_change", "magnetization", "zero_field_susceptibility"))

underflow = row(table(["--field", "1e-300", "--scale", "1e300"], 2))
assert underflow["status"] == "precision_limit" and underflow["energy_change"] is None

invalid = [["--field", "1"], ["--scale", "1"], ["--field", "1", "--scale", "0"],
           ["--field", "nan", "--scale", "1"], [*base, "--tolerance", "0"],
           [*base, "--max-lobes", "4097"], [*base, "--max-series-terms", "4097"],
           [*base, "--max-evaluations", "-1"], [*base, "--roots"]]
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)
    flags = ["--json", str(path/"result.json"), "--csv", str(path/"response.csv"), "--tsv", str(path/"response.tsv")]
    screen = run([*base, *flags, "--format", "plain"])
    assert "CPU time" in screen and "Used for:" not in screen
    exported = json.loads((path/"result.json").read_text())["tables"]["response"]
    expected = {k: "" if v is None else str(v).lower() if isinstance(v, bool) else str(v)
                for k, v in row(exported).items()}
    for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
        text = (path/f"response.{suffix}").read_text()
        assert "T_B=2*T1" in text
        actual = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
        assert actual == [expected]
    old = (path/"result.json").read_text()
    run([*base, "--json", str(path/"result.json")], 1)
    for args in invalid:
        run([*args, "--json", str(path/"result.json"), "--force"], 1)
        assert (path/"result.json").read_text() == old
    run([*base, *flags, "--force", "--no-retain", "--quiet"])
    assert json.loads((path/"result.json").read_text())["tables"]["response"]["rows"] == exported["rows"]
    run([*base, "--stream", "--no-retain", "--format", "plain"])
    run([*base, "--format", "pretty"])
    failure = run([*base, "--max-series-terms", "0", "--format", "csv", "--no-preamble"], 2)
    assert next(csv.DictReader(failure.splitlines()))["energy_change"] == ""

print("Kondo CLI contracts passed")
