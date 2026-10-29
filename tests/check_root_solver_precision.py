"""Runtime precision in open XXX and finite-size Lieb-Liniger frontends."""
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

xxx, pbc, obc, mpfr = sys.argv[1:]
enabled = mpfr.upper() in ("ON", "TRUE", "1")
supported = {xxx: enabled, pbc: enabled, obc: enabled}
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")
pi = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647")
base = {xxx: ["4"], pbc: ["2", "--length", "1", "--c", str(pi)],
        obc: ["2", "--length", "1", "--c", "1"]}


def run(program, args, code=0):
    result = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=60)
    assert result.returncode == code, (program, args, result.stdout, result.stderr)
    return result


def document(program, args, code=0):
    return json.loads(run(program, [*args, "--format", "json"], code).stdout)["tables"]


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]],
                     [v["value"] if isinstance(v, dict) else v for v in row]))
            for row in table["rows"]]


def energy(tables):
    return Decimal(rows(tables["states"])[0]["energy"])


def check_precision(tables, precision, bits):
    for table in tables.values():
        assert table["metadata"]["Precision"] == precision
        assert table["metadata"]["Arithmetic format"] == "arbitrary binary"
        assert int(table["metadata"]["Working precision (bits)"]) == bits
        for row in table["rows"]:
            for value in row:
                if isinstance(value, dict):
                    assert value["precision"] == bits, value


def same_rows(a, b):
    assert set(a) == set(b)
    for name in a:
        assert a[name]["rows"] == b[name]["rows"], name


if enabled:
    for precision, bits, error in (("bits-128", 128, "1e-33"), ("bits-256", 256, "1e-71"),
                                   ("digits-100", 333, "1e-95")):
        opts = ["--precision", precision, "--roots"]
        chain = document(xxx, [*base[xxx], *opts])
        explicit = document(xxx, ["4", "--quantum-numbers", "1,2", *opts])
        same_rows(chain, explicit)
        sectors = document(xxx, ["5", "--sectors", *opts])
        negative = document(xxx, ["5", "--sz", "-1.5", *opts])
        vacuum = document(xxx, ["4", "--quantum-numbers", "", *opts])
        assert energy(vacuum) == Decimal("0.75")
        assert rows(negative["states"])[0]["spin_reversed"]
        sector_rows = rows(sectors["states"])
        assert len(sector_rows) == 6
        assert sector_rows[1]["energy"] == rows(negative["states"])[0]["energy"]
        for a, b in zip(sector_rows, reversed(sector_rows)):
            assert a["energy"] == b["energy"]
        with localcontext() as ctx:
            ctx.prec = 120
            tol = Decimal(error)
            assert abs(energy(chain)+Decimal("0.75")+Decimal(3).sqrt()/2) < tol
            if enabled:
                ring = document(pbc, [*base[pbc], *opts])
                wall = document(obc, ["1", "--length", "3", "--c", "1", *opts])
                assert abs(energy(ring)-pi*pi/2) < tol
                assert abs(energy(wall)-pi*pi/9) < tol
                check_precision(ring, precision, bits)
                check_precision(wall, precision, bits)
        for tables in (chain, sectors, negative, vacuum):
            check_precision(tables, precision, bits)
        for program in (xxx, pbc, obc):
            if not supported[program]:
                continue
            scan_args = [*base[program], "--excitations", "all", *opts]
            if program != xxx:
                scan_args += ["--padding", "1"]
                empty = document(program, ["0", "--length", "1", "--c", "1", *opts])
                assert energy(empty) == 0 and not empty["roots"]["rows"]
                check_precision(empty, precision, bits)
            serial = document(program, [*scan_args, "--threads", "1"])
            parallel = document(program, [*scan_args, "--threads", "2"])
            same_rows(serial, parallel)
            check_precision(parallel, precision, bits)
            assert all(row["gap"] is not None for row in rows(serial["states"]))
            failed_scan = document(program, [*scan_args, "--max-iterations", "0", "--threads", "2"], 2)
            check_precision(failed_scan, precision, bits)
            assert not rows(failed_scan["reference"])[0]["converged"]
            assert not rows(failed_scan["failed"])[0]["converged"]
            assert all(row["gap"] is None for row in rows(failed_scan["states"]))
            failed = document(program, [*base[program], "--max-iterations", "0", *opts], 2)
            check_precision(failed, precision, bits)
            assert not rows(failed["states"])[0]["converged"]

if enabled:
    with localcontext() as ctx:
        ctx.prec = 120
        # Distinct box roots coalesce in fp64, and their O(c) energy shift
        # cannot be retained even in binary128. This checks a physical result.
        weak_args = ["2", "--length", "1", "--c", "1e-50"]
        weak = document(obc, [*weak_args, "--precision", "digits-100", "--roots"])
        assert abs((energy(weak)-2*pi*pi)/Decimal("1e-50")-3) < Decimal("1e-45")
        run(obc, [*weak_args, "--precision", "fp64"], 1)
        # Parsing length through double would erase this small energy change.
        length = Decimal(1)+Decimal("1e-50")
        one = document(obc, ["1", "--length", str(length), "--c", "1", "--precision", "digits-100"])
        assert abs(energy(one)-pi*pi/length**2) < Decimal("1e-95")

if enabled:
    tolerance = "1.00000000000000000000000000000000000000001e-80"
    for program in (xxx, pbc, obc):
        if not supported[program]:
            continue
        first = document(program, [*base[program], "--tolerance", tolerance, "--precision", "digits-100"])
        last = document(program, [*base[program], "--precision", "bits-333", "--tolerance", tolerance])
        same_rows(first, last)
        with localcontext() as ctx:
            ctx.prec = 120
            assert abs(Decimal(first["states"]["metadata"]["Residual tolerance"])/Decimal(tolerance)-1) < Decimal("1e-95")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"roots.tsv"
            run(program, [*base[program], "--precision", "digits-100", "--tsv-table", f"roots={path}"])
            text = path.read_text()
            assert "333" in text
            # Full-precision exported roots agree with the JSON table, not
            # merely the shorter screen formatting.
            root_table = document(program, [*base[program], "--precision", "digits-100", "--roots"])["roots"]
            key = "rapidity" if program == xxx else "k"
            assert rows(root_table)[0][key] in text

with tempfile.TemporaryDirectory() as directory:
    for program in (xxx, pbc, obc):
        path = Path(directory)/"preserve.json"
        path.write_text("sentinel")
        invalid = [["--precision", value] for value in ("bits-0", "digits--1", "mpfr", "bits-128junk")]
        if supported[program]:
            invalid += [["--precision", "bits-128", "--tolerance", value] for value in ("0", "-1", "nan", "inf")]
        else:
            invalid += [["--precision", value] for value in ("bits-128", "digits-100")]
        for args in invalid:
            result = run(program, [*base[program], *args, "--json", str(path), "--force"], 1)
            assert path.read_text() == "sentinel"

print("Open XXX and finite Lieb-Liniger precision contracts passed")
