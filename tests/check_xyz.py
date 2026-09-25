"""XYZ frontend: normalization, native precision, failure, metadata and exports."""
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
base = ["8", "--eta", "0.4", "--t", "0.7"]


def run(args, status=0):
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=30, env=env)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]], r)) for r in table["rows"]]


def tables(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete"
    for table in doc["tables"].values():
        assert table["summary"]["Outcome"] == ("success" if status == 0 else "partial")
        assert Decimal(table["summary"]["Run CPU seconds"]) >= 0
        assert "S=sigma/2" in table["metadata"]["Hamiltonian"]
    return doc["tables"]


for precision in precisions:
    flags = ["--precision", precision, "--roots"]
    with localcontext() as ctx:
        ctx.prec = 100
        tol = {"fp64": Decimal("1e-13"), "long-double": Decimal("1e-16"), "fp128": Decimal("1e-30")}[precision]
        # At eta=1/3, t=100: Jx=Jy=1 and Jz=1/2 up to exponentially small terms.
        eta = str(Decimal(1)/3)
        exact = tables(["2", "--eta", eta, "--t", "100", "--max-iterations", "0", *flags])
        r = rows(exact["states"])[0]
        assert abs(Decimal(r["energy"])+Decimal("1.25")) < tol
        assert abs(Decimal(r["energy_per_site"])+Decimal("0.625")) < tol
        assert Decimal(rows(exact["roots"])[0]["lambda_imag"]) == 0
        for eta in ("0.2", "0.5", "0.8"):
            two = tables(["2", "--eta", eta, "--t", "0.7", *flags])
            metadata = two["states"]["metadata"]
            energy = -sum(Decimal(metadata[key]) for key in ("Jx", "Jy", "Jz"))/2
            assert abs(Decimal(rows(two["states"])[0]["energy"])-energy) < tol*10
    normal = tables([*base, *flags])
    roots = rows(normal["roots"])
    assert [r["I"] for r in roots] == [-1.5, -0.5, 0.5, 1.5]
    assert all(Decimal(r["lambda_real"]) == 0 for r in roots)
    assert [Decimal(r["lambda_imag"]) for r in roots] == [Decimal(r["lambda_imag"]).copy_negate() for r in reversed(roots)]
    assert rows(normal["states"])[0]["status"] == "converged"
    failed = tables([*base, *flags, "--max-iterations", "0"], 2)
    r = rows(failed["states"])[0]
    assert r["status"] == "iteration_limit" and r["energy"] is None and r["momentum"] is None
    assert all(r["lambda_real"] is None and r["lambda_imag"] is None for r in rows(failed["roots"]))
    overflow = tables(["4", "--eta", "0.4", "--t", "0.000001", *flags], 2)
    assert rows(overflow["states"])[0]["status"] == "precision_limit"
    assert rows(overflow["states"])[0]["energy"] is None
    assert normal["states"]["rows"] == tables([*base, *flags, "--no-retain"])["states"]["rows"]

invalid = [["4"], ["3", "--eta", "0.4", "--t", "1"], ["-4", "--eta", "0.4", "--t", "1"],
           ["4", "--eta", "0", "--t", "1"], ["4", "--eta", "1", "--t", "1"],
           ["4", "--eta", "nan", "--t", "1"], ["4", "--eta", "0.4", "--t", "0"],
           ["4", "--eta", "0.4", "--t", "inf"], [*base, "--tolerance", "0"],
           [*base, "--max-iterations", "-1"], [*base, "--sz", "0"], [*base, "--excitations", "all"]]
for args in invalid:
    run(args, 1)
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)
    flags = ["--json", str(path/"state.json")]
    names = ("states", "roots")
    for suffix in ("csv", "tsv"):
        for name in names:
            flags += [f"--{suffix}-table", f"{name}={path/f'{name}.{suffix}'}"]
    screen = run([*base, "--roots", *flags, "--format", "plain"])
    assert "CPU time" in screen and "Used for:" not in screen
    doc = json.loads((path/"state.json").read_text())["tables"]
    for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
        for name in names:
            text = (path/f"{name}.{suffix}").read_text()
            assert "S=sigma/2" in text
            exported = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
            expected = [{k: "" if v is None else str(v).lower() if isinstance(v, bool) else str(v) for k,v in r.items()} for r in rows(doc[name])]
            assert exported == expected
    old = (path/"state.json").read_text()
    run([*base, "--json", str(path/"state.json")], 1)
    for args in invalid:
        run([*args, "--json", str(path/"state.json"), "--force"], 1)
        assert (path/"state.json").read_text() == old
    run([*base, "--roots", *flags, "--force", "--no-retain", "--quiet"])
    refreshed = json.loads((path/"state.json").read_text())["tables"]
    assert all(refreshed[name]["rows"] == doc[name]["rows"] for name in names)
    run([*base, "--roots", "--stream", "--no-retain", "--format", "plain"])
    run([*base, "--format", "pretty"])
    failure = run([*base, "--max-iterations", "0", "--format", "csv", "--no-preamble"], 2)
    assert next(csv.DictReader(failure.splitlines()))["energy"] == ""

print("XYZ CLI contracts passed")
