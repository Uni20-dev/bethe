"""Particle/Bethe-Yang physics, asymptotic labels, failure and export contracts."""
import csv
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

dispersion, bethe_yang, fp128 = sys.argv[1:]
precisions = ["fp64", "long-double"] + (["fp128"] if fp128.upper() in ("ON", "TRUE", "1") else [])
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")


def run(program, args, status=0):
    p = subprocess.run([program, *args], capture_output=True, text=True, timeout=60, env=env)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def table(program, args, status=0):
    doc = json.loads(run(program, [*args, "--format", "json"], status))
    assert doc["status"] == "complete"
    result = doc["tables"]["dispersion" if program == dispersion else "levels"]
    assert result["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    assert Decimal(result["summary"]["Run CPU seconds"]) >= 0
    return result


def rows(t):
    return [dict(zip([c["id"] for c in t["columns"]], row)) for row in t["rows"]]


for precision in precisions:
    base = ["--p", "0.4", "--momentum", "0.7", "--precision", precision]
    t = table(dispersion, base)
    data = {r["branch"]: r for r in rows(t)}
    assert len(data) == 14  # four species and ten unordered pairs
    assert data["soliton"]["charge"] == 1 and data["antisoliton"]["charge"] == -1
    assert data["B1"]["charge"] == 0 and data["soliton+antisoliton"]["kind"] == "threshold"
    assert t["rows"] == table(dispersion, [*base, "--no-retain"])["rows"]
    assert t["rows"] == table(dispersion, [*base, "--stream"])["rows"]
    with localcontext() as ctx:
        ctx.prec = 85
        tol = {"fp64": Decimal("8e-12"), "long-double": Decimal("8e-15"), "fp128": Decimal("8e-30")}[precision]
        reference = Decimal("1.368198089185226552659329247836824968488111471138021993646805247075555208")
        assert abs(Decimal(data["B1"]["energy"]) - reference) < tol / 100
        assert abs(Decimal(data["soliton+antisoliton"]["energy"]) - Decimal("4.49").sqrt()) < tol / 100
        pair = table(bethe_yang, ["--p", "2", "--length", "10", "--precision", precision])
        assert "wrapping corrections omitted" in pair["metadata"]["Approximation"]
        r = rows(pair)[0]
        assert r["number"] == 0.5 and r["charge"] == 2
        reference = Decimal("2.085735618785140488030977932237507957090038265403091568663327979382012603")
        assert abs(Decimal(r["energy"]) - reference) < tol
        assert Decimal(r["residual"]) + Decimal(r["quadrature_error"]) + Decimal(r["tail_bound"]) <= Decimal(pair["metadata"]["Counting tolerance"])
    free = rows(table(bethe_yang, ["--p", "1", "--length", "10", "--levels", "2", "--charge", "-2", "--precision", precision]))
    assert [r["number"] for r in free] == [0.5, 1.5]
    assert all(r["charge"] == -2 and int(r["evaluations"]) == 0 for r in free)
    assert Decimal(free[1]["energy"]) > Decimal(free[0]["energy"])

assert len(rows(table(dispersion, ["--p", "0.5", "--branch", "particles", "--momentum", "0"]))) == 3
assert len(rows(table(dispersion, ["--p", "1", "--branch", "particles", "--momentum", "0"]))) == 2
positive = rows(table(dispersion, ["--p", "2", "--momentum", "1"]))
negative = rows(table(dispersion, ["--p", "2", "--momentum", "-1"]))
assert [r["energy"] for r in positive] == [r["energy"] for r in negative]
for flag, expected, phase in [("--max-iterations", "iteration_limit", "converged"),
                              ("--max-phase-evaluations", "phase_limit", "quadrature_limit"),
                              ("--max-fourier-cutoffs", "phase_limit", "cutoff_limit")]:
    r = rows(table(bethe_yang, ["--p", "2", "--length", "10", flag, "0"], 2))[0]
    assert r["energy"] is None and r["rapidity"] is None
    assert r["status"] == expected and r["phase_status"] == phase
for args in (["--p", "0"], ["--p", "nan"], ["--p", "0.0001"], ["--p", "1", "--points", "1"],
             ["--p", "1", "--mass", "0"], ["--p", "1", "--momentum", "nan"],
             ["--p", "1", "--momentum", "0", "--points", "3"]):
    run(dispersion, args, 1)
for args in (["--p", "0.5", "--length", "10"], ["--p", "2", "--length", "0"],
             ["--p", "2", "--length", "10", "--number", "1"], ["--p", "2", "--length", "10", "--number", "-0.5"],
             ["--p", "2", "--length", "10", "--levels", "0"], ["--p", "2", "--length", "10", "--tolerance", "0"]):
    run(bethe_yang, args, 1)

for program, base in [(dispersion, ["--p", "0.4", "--points", "3"]),
                      (bethe_yang, ["--p", "1", "--length", "10", "--levels", "2"])]:
    with tempfile.TemporaryDirectory() as directory:
        js, cs, ts = [Path(directory) / ("data." + ext) for ext in ("json", "csv", "tsv")]
        screen = run(program, [*base, "--json", str(js), "--csv", str(cs), "--tsv", str(ts), "--format", "plain"])
        assert "CPU time" in screen and "Used for:" not in screen
        t = next(iter(json.loads(js.read_text())["tables"].values()))
        expected = [{k: "" if v is None else str(v) for k, v in row.items()} for row in rows(t)]
        for file, sep in ((cs, ","), (ts, "\t")):
            content = file.read_text()
            assert "Energy reference" in content
            actual = list(csv.DictReader((l for l in content.splitlines() if not l.startswith("#")), delimiter=sep))
            assert actual == expected
        old = js.read_bytes()
        run(program, [*base, "--json", str(js)], 1)
        assert js.read_bytes() == old
        run(program, [*base, "--mass", "0", "--force", "--json", str(js)], 1)
        assert js.read_bytes() == old
