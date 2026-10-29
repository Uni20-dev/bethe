"""Runtime precision, non-dyadic scans, scalar gating and thermal exports."""
import csv
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, mpfr = sys.argv[1:3]
long_running = sys.argv[3:] == ["--long-running"]
# Routine runs exercise precision, scans and export contracts at a modest
# spatial tolerance. The long-running variant retains the original tighter
# numerical controls, including both canonical round trips.
equilibrium_tolerance = "1e-20" if long_running else "1e-12"
canonical_tolerance = "1e-18" if long_running else "1e-10"
density_tolerance = "1e-16" if long_running else "1e-8"
density_error_bound = Decimal("1e-17" if long_running else "1e-9")
scan_tolerance = "1e-14" if long_running else "1e-12"
enabled = mpfr.upper() in ("ON", "TRUE", "1")
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")
base = ["--c", "4", "--temperature", "1", "--mu", "-1"]


def run(args, status=0):
    p = subprocess.run([program, *args], text=True, capture_output=True, env=env, timeout=1200 if long_running else 180)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def table(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete"
    t = doc["tables"]["thermodynamics"]
    assert t["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    return t


def rows(t):
    names = [c["id"] for c in t["columns"]]
    return [dict(zip(names, [v["value"] if isinstance(v, dict) else v for v in row])) for row in t["rows"]]


def precision(t, bits):
    assert t["metadata"]["Arithmetic format"] == "arbitrary binary"
    assert int(t["metadata"]["Working precision (bits)"]) == bits
    for row in t["rows"]:
        for value in row:
            if isinstance(value, dict):
                assert value["precision"] == bits, value


if enabled:
    for spelling, bits in (("bits-128", 128), ("bits-256", 256)):
        args = [*base, "--precision", spelling, "--tolerance", equilibrium_tolerance]
        state = table(args)
        precision(state, bits)
        assert rows(state)[0]["status"] == "converged"
        for budget, status in (("--max-cutoffs", "cutoff_limit"), ("--max-iterations", "iteration_limit")):
            failed = table([*args, budget, "0"], 2)
            precision(failed, bits)
            r = rows(failed)[0]
            assert r["status"] == status
            assert all(r[key] is None for key in ("pressure", "density", "energy_per_length", "entropy_per_length"))
        canonical = table(["--c", "4", "--temperature", "1", "--density", "0.1", "--precision", spelling,
                           "--tolerance", canonical_tolerance, "--density-tolerance", density_tolerance])
        precision(canonical, bits)
        assert abs(Decimal(rows(canonical)[0]["density"]) - Decimal("0.1")) < density_error_bound
        failed = table(["--c", "4", "--temperature", "1", "--density", "0.1", "--precision", spelling,
                        "--max-evaluations", "0"], 2)
        precision(failed, bits)
        assert rows(failed)[0]["status"] == "density_limit"

    with localcontext() as ctx:
        ctx.prec = 100
        for count in (4, 6):
            scan = table([*base, "--temperature-end", "2", "--points", str(count),
                          "--precision", "bits-256", "--tolerance", scan_tolerance, "--no-retain"])
            precision(scan, 256)
            for i, row in enumerate(rows(scan)):
                assert row["status"] == "converged"
                assert abs(Decimal(row["temperature"]) - (1 + Decimal(i) / Decimal(count-1))) < Decimal("1e-75")

        # Parse controls after dispatch, preserving a change below binary128.
        c = "4.0000000000000000000000000000000000000000000000000000001"
        tolerance = "1.0000000000000000000000000000000000000000000000001e-20"
        parsing = ["--c", c, "--temperature", "1", "--mu", "-1", "--tolerance", tolerance,
                   "--precision", "digits-100"]
        # Zero budget isolates literal parsing from expensive refinement.
        t = table(parsing if long_running else [*parsing, "--max-cutoffs", "0"], 0 if long_running else 2)
        precision(t, 333)
        assert abs(Decimal(t["metadata"]["c"])-Decimal(c)) < Decimal("1e-95")
        assert abs(Decimal(t["metadata"]["Requested inner tolerance"])/Decimal(tolerance)-1) < Decimal("1e-95")

    # The merged frontend's scheduler must preserve every numerical column,
    # work count and failure row at the explicitly selected precision.
    for ensemble in (("--mu", "-1"), ("--density", "0.1")):
        for fail in (False, True):
            args = ["--c", "4", "--temperature", "1", "--temperature-end", "2", "--points", "3",
                    *ensemble, "--precision", "bits-160", "--tolerance", "1e-12", "--no-retain"]
            if ensemble[0] == "--density":
                args += ["--density-tolerance", "1e-10"]
            if fail:
                args += ["--max-iterations" if ensemble[0] == "--mu" else "--max-evaluations", "0"]
            serial = table([*args, "--threads", "1"], 2 if fail else 0)
            parallel = table([*args, "--threads", "2"], 2 if fail else 0)
            precision(serial, 160)
            precision(parallel, 160)
            assert serial["columns"] == parallel["columns"]
            assert serial["rows"] == parallel["rows"]

    with tempfile.TemporaryDirectory() as folder:
        args = [*base, "--precision", "bits-160", "--tolerance", equilibrium_tolerance]
        expected = rows(table(args))
        for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
            path = Path(folder) / f"thermal.{suffix}"
            run([*args, f"--{suffix}", str(path), "--quiet", "--no-retain"])
            text = path.read_text()
            exported = list(csv.DictReader([line for line in text.splitlines() if not line.startswith("#")],
                                           delimiter=delimiter))
            assert exported == [{key: "" if value is None else value for key, value in row.items()}
                                for row in expected], (suffix, exported, expected)

with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "preserve.json"
    path.write_text("sentinel")
    cases = [["--precision", p] for p in ("bits-0", "digits--1", "bits-128junk")]
    if enabled:
        cases += [["--precision", "bits-128", "--tolerance", t] for t in ("0", "-1", "inf", "nan")]
    else:
        cases += [["--precision", p] for p in ("bits-128", "digits-100")]
    for extra in cases:
        run([*base, *extra, "--json", str(path), "--force"], 1)
        assert path.read_text() == "sentinel"

print("Thermal runtime precision and export checks passed")
