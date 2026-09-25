"""Sine-Gordon vacuum CLI: native precision, budgets and shared export contracts."""
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
base = ["--length", "1", "--p", "1"]


def run(args, status=0):
    result = subprocess.run([program, *args], text=True, capture_output=True, timeout=300, env=env)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    return result.stdout


def row(table):
    return dict(zip([c["id"] for c in table["columns"]], table["rows"][0]))


def table(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete" and list(doc["tables"]) == ["vacuum"]
    result = doc["tables"]["vacuum"]
    assert result["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    assert Decimal(result["summary"]["Run CPU seconds"]) >= 0
    assert "bulk-subtracted" in result["metadata"]["Energy convention"]
    assert "soliton mass" in result["metadata"]["Units"]
    return result


for precision in precisions:
    flags = ["--precision", precision]
    with localcontext() as ctx:
        ctx.prec = 100
        reference = Decimal("-0.3456042161410258359590522072055135625147457700998595063575857703631793172742704659519")
        pi = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406286")
        tol = {"fp64": Decimal("3e-10"), "long-double": Decimal("2e-13"), "fp128": Decimal("3e-28")}[precision]
        for mass, length in (("1", "1"), ("4", "0.25")):
            t = table(["--mass", mass, "--length", length, "--p", "1", *flags])
            r = row(t)
            assert r["converged"] and r["status"] == "converged"
            assert abs(Decimal(r["scaling_function"]) - reference) < tol
            assert abs(Decimal(r["casimir_energy"]) * Decimal(length) - reference) < tol
            assert abs(Decimal(r["effective_central_charge"]) + 6*reference/pi) < 2*tol
            assert int(r["iterations"]) == 0 and int(r["kernel_evaluations"]) == 0
            assert Decimal(t["metadata"]["Verification contour shift"]) < Decimal(t["metadata"]["Contour shift"])
        # Library tests independently validate interacting fp128 at its default tolerance.
        # A looser explicit tolerance keeps repeated CLI/export checks affordable.
        r = row(table(["--length", "1", "--p", "2", "--tolerance", "1e-8", *flags]))
        assert abs(Decimal(r["scaling_function"]) - Decimal("-0.3334619657421919836183854")) < Decimal("1e-8")
    for args, expected in (([*base, "--max-intervals", "64"], "mesh_limit"),
                           ([*base, "--max-cutoffs", "1"], "cutoff_limit"),
                           (["--length", "1", "--p", "2", "--max-iterations", "0"], "iteration_limit"),
                           (["--length", "1", "--p", "2", "--max-kernel-evaluations", "0"], "kernel_limit"),
                           (["--length", "1", "--p", "2", "--max-kernel-levels", "0"], "kernel_limit"),
                           (["--length", "1", "--p", "2", "--max-fourier-cutoffs", "0"], "kernel_limit")):
        # No expensive high-precision kernel construction is needed to test budget forwarding.
        r = row(table([*args, *flags, "--tolerance", "1e-8"], 2))
        assert not r["converged"] and r["status"] == expected
        assert all(r[k] is None for k in ("casimir_energy", "scaling_function", "effective_central_charge"))

attractive = row(table(["--length", "1", "--p", "0.5", "--contour-shift", "0.3", "--tolerance", "1e-8"]))
assert abs(Decimal(attractive["scaling_function"])-Decimal("-0.3860018985655546")) < Decimal("1e-8")
short_cutoff = row(table([*base, "--initial-cutoff", "0.5", "--max-cutoffs", "2", "--tolerance", "1e-4"], 2))
assert short_cutoff["status"] == "cutoff_limit" and short_cutoff["scaling_function"] is None
overflow = row(table(["--mass", "1e308", "--length", "2", "--p", "1"], 2))
assert overflow["status"] == "precision_limit" and overflow["casimir_energy"] is None

invalid = [["--length", "1"], ["--p", "1"], [*base, "--mass", "0"], [*base, "--mass", "nan"],
           ["--length", "-1", "--p", "1"], ["--length", "1", "--p", "0"],
           [*base, "--tolerance", "0"], [*base, "--contour-shift", "2"], [*base, "--initial-cutoff", "-1"],
           [*base, "--initial-intervals", "9"], [*base, "--max-intervals", "16385"],
           [*base, "--max-iterations", "-1"], [*base, "--roots"]]
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)
    flags = ["--json", str(path/"result.json"), "--csv", str(path/"vacuum.csv"), "--tsv", str(path/"vacuum.tsv")]
    screen = run([*base, *flags, "--format", "plain"])
    assert "CPU time" in screen and "Used for:" not in screen
    exported_table = json.loads((path/"result.json").read_text())["tables"]["vacuum"]
    expected = {k: "" if v is None else str(v).lower() if isinstance(v, bool) else str(v)
                for k, v in row(exported_table).items()}
    for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
        text = (path/f"vacuum.{suffix}").read_text()
        assert "bulk-subtracted" in text
        actual = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
        assert actual == [expected]
    old = (path/"result.json").read_text()
    run([*base, "--json", str(path/"result.json")], 1)
    for args in invalid:
        run([*args, "--json", str(path/"result.json"), "--force"], 1)
        assert (path/"result.json").read_text() == old
    run([*base, *flags, "--force", "--no-retain", "--quiet"])
    assert json.loads((path/"result.json").read_text())["tables"]["vacuum"]["rows"] == exported_table["rows"]
    run([*base, "--stream", "--no-retain", "--format", "plain"])
    run([*base, "--format", "pretty"])
    failure = run([*base, "--max-intervals", "64", "--format", "csv", "--no-preamble"], 2)
    assert next(csv.DictReader(failure.splitlines()))["casimir_energy"] == ""

print("Sine-Gordon CLI contracts passed")
