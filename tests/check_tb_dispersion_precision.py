"""MPFR TB bands: non-dyadic grids, independent series references and exports."""
import csv
from decimal import Decimal as D, localcontext
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
    doc = json.loads(run([*args, "--format", "json"]))
    assert doc["status"] == "complete"
    t = doc["tables"]["dispersion"]
    assert t["summary"]["Outcome"] == "success"
    return t


def rows(t):
    names = [c["id"] for c in t["columns"]]
    return [dict(zip(names, [x["value"] if isinstance(x, dict) else x for x in row])) for row in t["rows"]]


def precision(t, bits):
    assert t["metadata"]["Arithmetic format"] == "arbitrary binary"
    assert int(t["metadata"]["Working precision (bits)"]) == bits
    for row in t["rows"]:
        for x in row:
            if isinstance(x, dict): assert x["precision"] == bits


def sine(x):
    term = total = x
    for n in range(1, 250):
        term *= -x*x / D((2*n)*(2*n+1))
        previous, total = total, total + term
        if total == previous: return total
    raise AssertionError("reference series did not converge")


if enabled:
    with localcontext() as ctx:
        ctx.prec = 150
        pi = D("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128")
        exchange = "1.00000000000000000000000000000000000000000000000001"
        momentum = "0.30000000000000000000000000000000000000000000000001"
        previous = D(1)
        for option, bits in (("bits-128", 128), ("bits-256", 256), ("digits-100", 333)):
            eps = D(2) ** (12-bits)
            for count in (4, 6):
                base = ["--precision", option, "--points", str(count)]
                t = table(base)
                precision(t, bits)
                records = rows(t)
                assert len(records) == 3*count and all(r["status"] == "converged" for r in records)
                assert rows(table([*base, "--no-retain"])) == records
                folded = table([*base, "--folded"])
                precision(folded, bits)
                folded_rows = rows(folded)
                assert folded_rows[:count] == records[:count]
                for index, branch in enumerate(("spinon", "two-spinon", "four-spinon")):
                    for i, r in enumerate(records[index*count:(index+1)*count]):
                        assert r["branch"] == branch
                        q = (pi if index == 0 else 2*pi)*i/(count-1)
                        assert abs(D(r["p"])-q) < eps
                        reduced = min(q, 2*pi-q)
                        v = 2*pi
                        if index == 0:
                            assert r["spin"] == 0.5 and r["lower"] is None and r["upper"] is None
                            assert abs(D(r["energy"])-v*sine(q)) < eps
                        else:
                            assert r["spin"] is None and r["energy"] is None
                            assert abs(D(r["lower"])-v*abs(sine(q))) < eps
                            expected = 2*v*sine(q/2) if index == 1 else 4*v*sine(pi/2-reduced/4)
                            assert abs(D(r["upper"])-expected) < eps
                            f = folded_rows[index*count+i]
                            folded_q = min(reduced, pi-reduced)
                            expected = 2*v*sine(pi/2-folded_q/2) if index == 1 else 4*v*sine(pi/2-folded_q/4)
                            assert abs(D(f["upper"])-expected) < eps
                    last = records[(index+1)*count-1]
                    assert D(last["energy"] if index == 0 else last["lower"]) == 0
                    assert D(last["cell_momentum"]) == 0
            t = table(["--branch", "spinon", "--exchange", exchange, "--momentum", momentum, "--precision", option])
            precision(t, bits)
            expected = 2*pi*D(exchange)*sine(D(momentum))
            error = abs(D(rows(t)[0]["energy"])-expected)
            assert error < eps and error < previous/D("1e15"), (bits, error)
            previous = error
            assert abs(D(t["metadata"]["Exchange J"])-D(exchange)) < eps
            tiny = table(["--momentum", "1e-7000", "--precision", option])
            precision(tiny, bits)
            r = rows(tiny)
            assert abs(D(r[0]["energy"])/D("1e-7000")-2*pi) < eps
            assert abs(D(r[1]["upper"])/D("1e-7000")-2*pi) < eps
            assert r[1]["lower"] == r[2]["lower"]
            huge = table(["--branch", "spinon", "--exchange", "1e7000", "--momentum", "1e-7000", "--precision", option])
            assert abs(D(rows(huge)[0]["energy"])-2*pi) < eps

    with tempfile.TemporaryDirectory() as directory:
        base = ["--points", "6", "--folded", "--precision", "digits-100"]
        expected = [{k: "" if v is None else str(v) for k, v in r.items()} for r in rows(table(base))]
        for extension, delimiter in (("csv", ","), ("tsv", "\t")):
            path = Path(directory)/f"tb.{extension}"
            run([*base, f"--{extension}", str(path), "--quiet"])
            text = path.read_text()
            assert "333" in text and "Hamiltonian" in text
            actual = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
            assert actual == expected

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)/"preserve.json"
    path.write_text("sentinel")
    invalid = [["--precision", p] for p in ("bits-0", "digits--1", "bits-128junk")]
    if enabled:
        invalid += [["--precision", "bits-128", "--exchange", j] for j in ("0", "-1", "nan", "inf")]
        invalid += [["--precision", "digits-100", "--momentum", q] for q in ("-1", "nan", "inf", "4")]
    else:
        invalid += [["--precision", p] for p in ("bits-128", "digits-100")]
    for args in invalid:
        assert not run([*args, "--force", "--json", str(path)], 1)
        assert path.read_text() == "sentinel"
assert "bits-N" in run(["--help"]) and "digits-N" in run(["--help"])
