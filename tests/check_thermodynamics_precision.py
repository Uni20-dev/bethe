"""Explicit precision, physical corrections and exports for Lieb-Liniger thermodynamics."""
from decimal import Decimal, localcontext
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, mpfr, mplapack = sys.argv[1:]
enabled = all(flag.upper() in ("ON", "TRUE", "1") for flag in (mpfr, mplapack))
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")


def run(args, status=0, error=None):
    result = subprocess.run([program, *args], text=True, capture_output=True, env=env, timeout=60)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    if error is not None:
        assert not result.stdout and error in result.stderr, (args, result.stdout, result.stderr)
    return result.stdout


def table(args, status=0):
    document = json.loads(run([*args, "--format", "json"], status))
    assert document["status"] == "complete"
    result = document["tables"]["dispersion"]
    assert result["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    return result


def rows(t):
    columns = [column["id"] for column in t["columns"]]
    return [dict(zip(columns, [v["value"] if isinstance(v, dict) else v for v in row])) for row in t["rows"]]


def check_precision(t, precision, bits):
    assert t["metadata"]["Precision"] == precision
    assert t["metadata"]["Arithmetic format"] == "arbitrary binary"
    assert int(t["metadata"]["Working precision (bits)"]) == bits
    for row in t["rows"]:
        for value in row:
            if isinstance(value, dict):
                assert value["precision"] == bits, value


if enabled:
    for precision, bits in (("bits-128", 128), ("bits-256", 256), ("digits-100", 333)):
        base = ["--c", "4", "--points", "3", "--precision", precision]
        serial = table([*base, "--threads", "1"])
        parallel = table([*base, "--threads", "3", "--no-retain"])
        check_precision(serial, precision, bits)
        check_precision(parallel, precision, bits)
        assert serial["rows"] == parallel["rows"]
        r = rows(serial)
        assert len(r) == 6 and all(row["status"] == "converged" for row in r)
        assert Decimal(r[0]["energy"]) == Decimal(r[3]["energy"]) == Decimal(r[5]["energy"]) == 0
        for args, expected in ((["--max-background-iterations", "0"], "density_limit"),
                               (["--max-nodes", "16"], "mesh_limit")):
            failed = table([*base, *args], 2)
            check_precision(failed, precision, bits)
            assert all(row["energy"] is None and row["status"] == expected for row in rows(failed))
        failed = table(["--c", "4", "--momentum", "1", "--precision", precision, "--max-iterations", "0"], 2)
        assert all(row["energy"] is None and row["status"] == "momentum_limit" for row in rows(failed))

    with localcontext() as ctx:
        ctx.prec = 120
        # Independently known strong-coupling asymptotics expose a real energy
        # correction that is erased by native-double/binary128 arithmetic.
        pi = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647")
        strong = table(["--c", "1e50", "--momentum", "0", "--precision", "digits-100"])
        energy = Decimal(strong["metadata"]["Ground energy/length"])
        assert abs((energy / (pi*pi/3)-1)*Decimal("1e50")+4) < Decimal("1e-43")

        base = ["--c", "4", "--momentum", "0", "--precision", "digits-100"]
        reference = table(base)
        density = Decimal(1) + Decimal("1e-50")
        scaled = table(["--c", str(4*density), "--density", str(density), "--momentum", "0",
                        "--precision", "digits-100"])
        expected = Decimal(reference["metadata"]["Ground energy/length"]) * density**3
        assert abs(Decimal(scaled["metadata"]["Ground energy/length"])-expected) < Decimal("1e-95")

        tolerance = "1.00000000000000000000000000000000000000001e-60"
        first = table(["--c", "4", "--momentum", "1", "--tolerance", tolerance, "--precision", "digits-100"])
        last = table(["--c", "4", "--momentum", "1", "--precision", "bits-333", "--tolerance", tolerance])
        assert first["rows"] == last["rows"]
        assert abs(Decimal(first["metadata"]["Relative tolerance"])/Decimal(tolerance)-1) < Decimal("1e-95")

    with tempfile.TemporaryDirectory() as directory:
        base = ["--c", "4", "--momentum", "1", "--precision", "bits-256"]
        expected = rows(table(base))
        for extension, delimiter in (("csv", ","), ("tsv", "\t")):
            path = Path(directory) / f"curves.{extension}"
            run([*base, f"--{extension}", str(path), "--quiet"])
            text = path.read_text()
            assert "256" in text
            exported = list(csv.DictReader([line for line in text.splitlines() if not line.startswith("#")],
                                           delimiter=delimiter))
            for a, b in zip(expected, exported):
                for name in ("p", "energy", "rapidity", "edge_distance", "energy_error", "momentum_error"):
                    assert a[name] == b[name], (name, a[name], b[name])

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "preserve.json"
    path.write_text("sentinel")
    invalid = [["--precision", p] for p in ("bits-0", "digits--1", "bits-128junk")]
    if enabled:
        invalid += [["--precision", "bits-128", "--tolerance", t] for t in ("0", "-1", "inf", "nan")]
    else:
        invalid += [["--precision", p] for p in ("bits-128", "digits-100")]
    for args in invalid:
        missing_provider = not enabled and args in (["--precision", "bits-128"], ["--precision", "digits-100"])
        run(["--c", "4", "--momentum", "0", *args, "--json", str(path), "--force"], 1,
            error="arbitrary-precision linear solves require" if missing_provider else None)
        assert path.read_text() == "sentinel"

print("Thermodynamic precision, physical corrections and exports passed")
