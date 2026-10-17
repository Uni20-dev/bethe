"""Explicit runtime precision across every periodic XXX frontend mode."""
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


def run(args, code=0):
    result = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=60)
    assert result.returncode == code, (args, result.stdout, result.stderr)
    return result


def document(args, code=0):
    return json.loads(run([*args, "--format", "json"], code).stdout)["tables"]


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]],
                     [v["value"] if isinstance(v, dict) else v for v in row]))
            for row in table["rows"]]


def same_rows(left, right):
    assert set(left) == set(right)
    for name in left:
        assert left[name]["rows"] == right[name]["rows"], name


if enabled:
    for precision, bits, error in (("bits-128", 128, "1e-34"), ("bits-256", 256, "1e-72"),
                                   ("digits-100", 333, "1e-95")):
        args = ["--precision", precision, "--roots"]
        ground = document(["6", *args])
        sectors = document(["5", "--sectors", *args])
        negative = document(["5", "--sz", "-1.5", *args])
        spinons = document(["5", "--spinons", *args])
        explicit = document(["4", "--quantum-numbers=-1/2,1/2", *args])
        vacuum = document(["4", "--quantum-numbers", "", *args])
        scan = document(["8", "--excitations", "all", "--threads", "1", *args])
        parallel = document(["8", "--excitations", "all", "--threads", "2", *args])
        same_rows(scan, parallel)
        assert len(rows(scan["states"])) == 10
        assert all(row["gap"] is not None for row in rows(scan["states"]))
        assert not vacuum["roots"]["rows"]
        assert Decimal(rows(vacuum["states"])[0]["energy"]) == 1
        assert rows(negative["states"])[0]["spin_reversed"]
        with localcontext() as ctx:
            ctx.prec = 120
            tol = Decimal(error)
            assert abs(Decimal(rows(ground["states"])[0]["energy"]) + (2+Decimal(13).sqrt())/2) < tol
            assert abs(Decimal(rows(explicit["states"])[0]["energy"]) + 2) < tol
            sector_rows = rows(sectors["states"])
            assert len(sector_rows) == 6
            assert abs(Decimal(sector_rows[1]["energy"]) + Decimal(5).sqrt()/4) < tol
            assert sector_rows[1]["energy"] == rows(negative["states"])[0]["energy"]
            for a, b in zip(sector_rows, reversed(sector_rows)):
                assert a["energy"] == b["energy"]
            center = rows(spinons["spinons"])[1]
            assert abs(Decimal(center["bulk_subtracted_energy"]) - (5*Decimal(2).ln()-2)) < tol
            assert abs(Decimal(center["k"])-Decimal(center["epsilon_inf"])) < tol
            # Tutorial regression: a long-wavelength magnon needs both a
            # well-resolved root and enough digits to retain its tiny shift
            # in the extensive total energy. No half-filled large-N solve.
            n, shift_tolerance, root_tolerance = {
                128: (10**8, "1e-14", "1e-25"),
                256: (10**16, "1e-28", "1e-58"),
                333: (10**16, "1e-50", "1e-80"),
            }[bits]
            magnon = document([str(n), "--quantum-numbers", str(n//2-1), *args])
            pi = Decimal("3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865")
            x = pi/n
            # The omitted sine term contributes relative O(N^-6), safely
            # below the requested checks at these sizes.
            exact_shift = 2*(x-x**3/6+x**5/120)**2
            shift = Decimal(n)/4-Decimal(rows(magnon["states"])[0]["energy"])
            z = Decimal(rows(magnon["roots"])[0]["rapidity"])
            assert abs(shift/exact_shift-1) < Decimal(shift_tolerance)
            assert abs(2/(1+z*z)/exact_shift-1) < Decimal(root_tolerance)
        for tables in (ground, sectors, negative, spinons, explicit, vacuum, scan):
            for table in tables.values():
                assert table["metadata"]["Precision"] == precision
                assert table["metadata"]["Arithmetic format"] == "arbitrary binary"
                assert int(table["metadata"]["Working precision (bits)"]) == bits
                for row in table["rows"]:
                    for value in row:
                        if isinstance(value, dict):
                            assert value["precision"] == bits

    # Option order cannot determine the precision used to parse a real control.
    tolerance = "1.00000000000000000000000000000000000000001e-90"
    first = document(["6", "--tolerance", tolerance, "--precision", "digits-100"])
    last = document(["6", "--precision", "bits-333", "--tolerance", tolerance])
    same_rows(first, last)
    with localcontext() as ctx:
        ctx.prec = 120
        assert abs(Decimal(first["states"]["metadata"]["Residual tolerance"])/Decimal(tolerance)-1) < Decimal("1e-95")
    same_rows(document(["5", "--spinons", "--precision", "digits-100"]),
              document(["5", "--spinons", "--precision", "bits-333"]))

    failed = document(["8", "--excitations", "all", "--max-iterations", "0", "--precision", "bits-128"], 2)
    assert not rows(failed["reference"])[0]["converged"]
    assert not rows(failed["failed"])[0]["converged"]
    assert all(row["gap"] is None for row in rows(failed["states"]))
    failed_spinons = document(["7", "--spinons", "--max-iterations", "0", "--precision", "bits-128"], 2)
    assert any(not row["converged"] for row in rows(failed_spinons["states"]))

    with tempfile.TemporaryDirectory() as directory:
        # A roots file requests root production, not an extra screen table.
        path = Path(directory)/"roots.tsv"
        screen = run(["4", "--precision", "digits-100", "--tsv-table", f"roots={path}"]).stdout
        assert "Rapidities" not in screen
        assert "0.577350269189625764509148780501957" in path.read_text()

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)/"preserve.json"
    path.write_text("sentinel")
    invalid = [["4", "--precision", value] for value in
               ("bits-0", "digits--1", "bits-128junk", "mpfr", "bits-99999999999999999999999999")]
    if enabled:
        invalid += [["6", "--precision", "bits-128", "--tolerance", value] for value in ("nan", "0", "-1")]
        invalid += [["6", "--spinons", "--precision", "bits-128"],
                    ["8", "--excitations", "all", "--max-candidates", "9", "--precision", "bits-128"]]
    else:
        invalid += [["4", "--precision", value] for value in ("bits-128", "digits-100")]
    for args in invalid:
        run([*args, "--json", str(path), "--force"], 1)
        assert path.read_text() == "sentinel"

print("Periodic XXX precision contracts passed")
