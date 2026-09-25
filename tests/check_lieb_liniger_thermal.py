"""Thermal ensembles, native precision, scan and shared export contracts."""
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
base = ["--c", "4", "--temperature", "1", "--mu", "-1"]


def run(args, status=0):
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=120, env=env)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def table(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete"
    t = doc["tables"]["thermodynamics"]
    assert t["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    assert Decimal(t["summary"]["Run CPU seconds"]) >= 0
    assert "k_B=1" in t["metadata"]["Hamiltonian"]
    return t


def rows(t):
    return [dict(zip([c["id"] for c in t["columns"]], row)) for row in t["rows"]]


with localcontext() as ctx:
    ctx.prec = 70
    pi = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406")
    z = Decimal(-1).exp()
    p = n = Decimal(0)
    power = z
    for l in range(1, 200):
        term = power / Decimal(l).sqrt()
        n += term
        p += term / l
        power *= -z
    n /= 2 * pi.sqrt()
    p /= 2 * pi.sqrt()
    for precision in precisions:
        tol = {"fp64": Decimal("1e-12"), "long-double": Decimal("1e-15"), "fp128": Decimal("1e-28")}[precision]
        tonks = rows(table(["--c", "1e40", "--temperature", "1", "--mu", "-1", "--precision", precision]))[0]
        assert abs(Decimal(tonks["pressure"]) - p) < tol * p
        assert abs(Decimal(tonks["density"]) - n) < tol * n
        assert abs(Decimal(tonks["energy_per_length"]) - p/2) < tol * p
        assert abs(Decimal(tonks["entropy_per_length"]) - (3*p/2+n)) < tol
        normal = rows(table([*base, "--precision", precision]))[0]
        canonical = rows(table(["--c", "4", "--temperature", "1", "--density", normal["density"],
                                "--precision", precision]))[0]
        assert abs(Decimal(canonical["chemical_potential"]) + 1) < 10 * tol
        assert canonical["status"] == "converged"
        assert Decimal(canonical["density_error"]) < 10 * tol

scan = [*base, "--temperature-end", "2", "--points", "3"]
t = table(scan)
assert [r["temperature"] for r in rows(t)] == ["1", "1.5", "2"]
assert t["rows"] == table([*scan, "--no-retain"])["rows"]
reverse = rows(table(["--c", "4", "--temperature", "2", "--temperature-end", "1", "--points", "3", "--density", "0.1"]))
assert [r["temperature"] for r in reverse] == ["2", "1.5", "1"]
assert all(abs(Decimal(r["density"]) - Decimal("0.1")) < Decimal("1e-12") for r in reverse)
for option, value, expected in (("--max-iterations", "0", "iteration_limit"),
                                ("--max-nodes", "16", "mesh_limit"),
                                ("--max-cutoffs", "0", "cutoff_limit")):
    r = rows(table([*base, option, value], 2))[0]
    assert r["status"] == expected and all(r[k] is None for k in ("density", "pressure", "energy_per_length", "entropy_per_length"))
failed = rows(table(["--c", "4", "--temperature", "1", "--density", "1", "--max-evaluations", "0"], 2))[0]
assert failed["status"] == "density_limit" and failed["chemical_potential"] is None
mixed = rows(table(["--c", "4", "--temperature", "0.0001", "--temperature-end", "1",
                    "--points", "2", "--mu", "-1"], 2))
assert mixed[0]["status"] == "precision_limit" and mixed[0]["pressure"] is None
assert mixed[1]["status"] == "converged" and Decimal(mixed[1]["pressure"]) > 0
# The documented lower-temperature scan needs a larger mesh than the default.
run(["--c", "4", "--temperature", "0.5", "--temperature-end", "2", "--points", "9",
     "--density", "1", "--max-nodes", "512", "--quiet"])

for args in (["--c", "4", "--temperature", "1"], [*base, "--density", "1"],
             [*base, "--points", "3"], [*base, "--temperature-end", "0"],
             [*base, "--temperature-end", "2", "--points", "1"],
             [*base, "--max-nodes", "513"], [*base, "--tolerance", "nan"],
             [*base, "--density-tolerance", "1e-8"], [*base, "--max-evaluations", "1"],
             ["--c", "-1", "--temperature", "1", "--mu", "0"],
             ["--c", "4", "--temperature", "0", "--mu", "0"],
             ["--c", "4", "--temperature", "1", "--density", "0"]):
    run(args, 1)

with tempfile.TemporaryDirectory() as directory:
    js, cs, ts = [Path(directory) / ("thermal." + suffix) for suffix in ("json", "csv", "tsv")]
    screen = run([*scan, "--json", str(js), "--csv", str(cs), "--tsv", str(ts), "--format", "plain"])
    assert "CPU time" in screen and "Used for:" not in screen
    expected = rows(json.loads(js.read_text())["tables"]["thermodynamics"])
    for file, delimiter in ((cs, ","), (ts, "\t")):
        text = file.read_text()
        assert "#" in text and "Hamiltonian" in text
        exported = list(csv.DictReader([line for line in text.splitlines() if not line.startswith("#")], delimiter=delimiter))
        assert exported == [{k: "" if v is None else v for k, v in r.items()} for r in expected]
    old = js.read_text()
    run([*base, "--json", str(js)], 1)
    assert js.read_text() == old
    run(["--c", "-1", "--temperature", "1", "--mu", "0", "--json", str(js), "--force"], 1)
    assert js.read_text() == old
    run([*scan, "--json", str(js), "--force", "--quiet", "--no-retain"])
    assert rows(json.loads(js.read_text())["tables"]["thermodynamics"]) == expected
    run([*base, "--stream", "--no-retain", "--format", "plain"])
    run([*base, "--format", "pretty"])
    failed = run([*base, "--max-iterations", "0", "--format", "csv", "--no-preamble"], 2)
    assert next(csv.DictReader(failed.splitlines()))["pressure"] == ""

print("Lieb-Liniger thermal CLI contracts passed")
