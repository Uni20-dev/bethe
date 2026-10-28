"""Explicit MPFR SU(3) bands, non-dyadic grids, precision growth and exports."""
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
    result = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=30)
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
    for row in t["rows"]:
        for value in row:
            if isinstance(value, dict):
                assert value["precision"] == bits, value


def cosine(x):
    # Independent Decimal Taylor series, including the cancellation-prone
    # cosine difference deliberately avoided by the production band formula.
    value = term = Decimal(1)
    for n in range(1, 250):
        term *= -x*x / Decimal((2*n-1)*(2*n))
        last, value = value, value + term
        if value == last:
            return value
    raise AssertionError("reference series did not converge")


if enabled:
    with localcontext() as ctx:
        ctx.prec = 150
        pi = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128")
        exchange = "1.00000000000000000000000000000000000000000000000001"
        momentum = "0.30000000000000000000000000000000000000000000000001"
        amplitude = 4*pi / (3*Decimal(3).sqrt())
        previous = {"3": Decimal(1), "bar3": Decimal(1)}
        for precision, bits in (("bits-128", 128), ("bits-256", 256), ("digits-100", 333)):
            tolerance = Decimal(2) ** (10-bits)
            for count in (4, 6):
                base = ["--precision", precision, "--points", str(count)]
                t = table(base)
                check_precision(t, bits)
                records = rows(t)
                assert len(records) == 4*count and all(r["status"] == "converged" for r in records)
                assert rows(table([*base, "--no-retain"])) == records
                for offset in (0, count):
                    assert Decimal(records[offset]["energy"]) == Decimal(records[offset+count-1]["energy"]) == 0
                for offset in (2*count, 3*count):
                    assert Decimal(records[offset]["lower"]) == Decimal(records[offset+count-1]["lower"]) == 0
                    assert Decimal(records[offset+count-1]["cell_momentum"]) == 0
                for index, multiplier in enumerate((Decimal(4)/3, Decimal(2)/3, Decimal(2), Decimal(2))):
                    for i, r in enumerate(records[index*count:(index+1)*count]):
                        assert abs(Decimal(r["p"])-multiplier*pi*i/(count-1)) < tolerance
                folded = table([*base, "--folded"])
                check_precision(folded, bits)
                assert rows(folded)[:2*count] == records[:2*count]
            for branch, centre in (("3", 2*pi/3), ("bar3", pi/3)):
                expected = Decimal(exchange)*amplitude*(cosine(Decimal(momentum)-centre)-cosine(centre))
                t = table(["--branch", branch, "--exchange", exchange, "--momentum", momentum, "--precision", precision])
                check_precision(t, bits)
                r = rows(t)[0]
                error = abs(Decimal(r["energy"])-expected)
                assert error < tolerance and error < previous[branch]/Decimal("1e15"), (bits, branch, error)
                previous[branch] = error
                assert abs(Decimal(t["metadata"]["Exchange J"])-Decimal(exchange)) < tolerance
                tiny = rows(table(["--branch", branch, "--momentum", "1e-7000", "--precision", precision]))[0]
                assert abs(Decimal(tiny["energy"])/Decimal("1e-7000")-2*pi/3) < tolerance

    with tempfile.TemporaryDirectory() as directory:
        base = ["--points", "6", "--folded", "--precision", "digits-100"]
        expected = [{k: "" if v is None else str(v) for k, v in r.items()} for r in rows(table(base))]
        for extension, delimiter in (("csv", ","), ("tsv", "\t")):
            path = Path(directory) / f"uls.{extension}"
            run([*base, f"--{extension}", str(path), "--quiet"])
            text = path.read_text()
            assert "333" in text and "Hamiltonian" in text
            actual = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
            assert actual == expected

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "preserve.json"
    path.write_text("sentinel")
    invalid = [["--precision", p] for p in ("bits-0", "digits--1", "bits-128junk")]
    if enabled:
        invalid += [["--precision", "bits-128", "--exchange", j] for j in ("0", "-1", "nan", "inf")]
        invalid += [["--precision", "digits-100", "--momentum", q] for q in ("-1", "nan", "inf", "3")]
    else:
        invalid += [["--precision", p] for p in ("bits-128", "digits-100")]
    for args in invalid:
        assert not run([*args, "--force", "--json", str(path)], 1)
        assert path.read_text() == "sentinel"
