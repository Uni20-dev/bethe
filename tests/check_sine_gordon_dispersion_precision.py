"""MPFR particle masses, strict poles, arbitrary range and complete exports."""
import csv
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, mpfr = sys.argv[1:]
enabled = mpfr.upper() in ("ON", "TRUE", "1")
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")


def run(args, status=0):
    result = subprocess.run([program, *args], text=True, capture_output=True, timeout=30, env=env)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    return result.stdout


def table(args):
    document = json.loads(run([*args, "--format", "json"]))
    assert document["status"] == "complete"
    result = document["tables"]["dispersion"]
    assert result["summary"]["Outcome"] == "success"
    return result


def rows(t):
    columns = [c["id"] for c in t["columns"]]
    return [dict(zip(columns, [v["value"] if isinstance(v, dict) else v for v in row])) for row in t["rows"]]


def check_precision(t, bits):
    assert t["metadata"]["Arithmetic format"] == "arbitrary binary"
    assert int(t["metadata"]["Working precision (bits)"]) == bits
    # Each real cell, including zero momentum and rest energy, must carry
    # the requested precision rather than remaining an exact rational.
    for row in t["rows"]:
        for i in (3, 4, 5):
            assert isinstance(row[i], dict) and row[i]["precision"] == bits, row


if enabled:
    with localcontext() as ctx:
        ctx.prec = 150
        previous = {"B1": Decimal(1), "B2": Decimal(1)}
        mass = "1.00000000000000000000000000000000000000000000000001"
        momentum = "0.70000000000000000000000000000000000000000000000001"
        for precision, bits in (("bits-128", 128), ("bits-256", 256), ("digits-100", 333)):
            tol = Decimal(2) ** (10-bits)
            t = table(["--p", "0.4", "--mass", mass, "--momentum", momentum, "--precision", precision])
            check_precision(t, bits)
            data = {r["branch"]: r for r in rows(t)}
            assert len(data) == 14 and all(r["status"] == "converged" for r in data.values())
            # Independent algebraic mass formula at p=2/5, no sine/hypot.
            for branch, sign in (("B1", -1), ("B2", 1)):
                square = (5+sign*Decimal(5).sqrt())/2
                expected = (Decimal(mass)**2*square + Decimal(momentum)**2).sqrt()
                error = abs(Decimal(data[branch]["energy"])-expected)
                assert error < tol and error < previous[branch]/Decimal("1e15"), (bits, branch, error)
                previous[branch] = error
            assert abs(Decimal(t["metadata"]["Soliton mass M"])-Decimal(mass)) < tol
            assert abs(Decimal(data["B1"]["momentum"])-Decimal(momentum)) < tol
            for count in (4, 6):
                base = ["--p", "0.4", "--points", str(count), "--precision", precision]
                t = table(base)
                check_precision(t, bits)
                records = rows(t)
                assert len(records) == 14*count
                assert records == rows(table([*base, "--no-retain"]))
                assert records == rows(table([*base, "--stream"]))
                for start in range(0, len(records), count):
                    for j, r in enumerate(records[start:start+count]):
                        assert abs(Decimal(r["momentum"])-5*Decimal(j)/(count-1)) < tol
                        e = (Decimal(r["rest_energy"])**2+Decimal(r["momentum"])**2).sqrt()
                        assert abs(Decimal(r["energy"])-e) < tol
            for scale in ("1e-7000", "1e7000"):
                records = rows(table(["--p", "0.4", "--mass", scale, "--momentum", "-"+scale, "--precision", precision]))
                first = next(r for r in records if r["branch"] == "soliton")
                assert abs(Decimal(first["energy"])/Decimal(scale)-Decimal(2).sqrt()) < tol

        half = Decimal("0.5")
        for coupling, count in ((half-Decimal("1e-60"), 4), (half, 3), (half+Decimal("1e-60"), 3), (Decimal(1), 2)):
            t = table(["--p", str(coupling), "--branch", "particles", "--momentum", "0", "--precision", "digits-100"])
            check_precision(t, 333)
            assert len(rows(t)) == count
            assert all(Decimal(r["momentum"]) == 0 for r in rows(t))

    with tempfile.TemporaryDirectory() as directory:
        base = ["--p", "0.4", "--points", "6", "--precision", "digits-100"]
        expected = [{k: str(v) for k, v in r.items()} for r in rows(table(base))]
        for extension, delimiter in (("csv", ","), ("tsv", "\t")):
            path = Path(directory)/f"sg.{extension}"
            run([*base, f"--{extension}-table", f"dispersion={path}", "--quiet"])
            text = path.read_text()
            assert "333" in text and "Energy reference" in text
            actual = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
            assert actual == expected

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)/"preserve.json"
    path.write_text("sentinel")
    invalid = [["--precision", p] for p in ("bits-0", "digits--1", "bits-128junk")]
    if enabled:
        invalid += [["--precision", "bits-128", "--mass", x] for x in ("0", "-1", "nan", "inf")]
        invalid += [["--precision", "bits-128", "--momentum", x] for x in ("nan", "inf")]
        invalid += [["--precision", "bits-128", "--max-momentum", x] for x in ("0", "-1", "nan", "inf")]
    else:
        invalid += [["--precision", p] for p in ("bits-128", "digits-100")]
    for args in invalid:
        assert not run(["--p", "0.4", *args, "--json", str(path), "--force"], 1)
        assert path.read_text() == "sentinel"
    if enabled:
        for coupling in ("0", "-1", "nan", "inf", "1e-7000"):
            assert not run(["--p", coupling, "--precision", "digits-100", "--json", str(path), "--force"], 1)
            assert path.read_text() == "sentinel"
