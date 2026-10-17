"""End-to-end arbitrary precision for Richardson pairing and central spin."""
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

richardson, central, mpfr = sys.argv[1:]
enabled = mpfr.upper() in ("ON", "TRUE", "1")
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")
base = {
    richardson: ["--levels", "0,1,2,3", "--pairs", "2", "--g", "1"],
    central: ["--couplings", "1,0.7,0.3", "--sz", "0", "--field", "1"],
}


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


if enabled:
    for precision, bits, error in (("bits-128", 128, "1e-33"), ("bits-256", 256, "1e-71"),
                                   ("digits-100", 333, "1e-95")):
        opts = ["--precision", precision, "--variables"]
        pair = document(richardson, ["--levels", "0,1", "--pairs", "1", "--g", "1", *opts])
        spin = document(central, ["--couplings", "1", "--field", "1", "--sz", "0", *opts])
        isolated = document(central, ["--couplings", "", "--field", "3", "--sz", "-0.5", *opts])
        blocked = document(richardson, ["--levels", "0,1,2,3", "--pairs", "1", "--g", "0",
                                        "--blocked", "0,2", *opts])
        vacuum = document(richardson, ["--levels", "0,1", "--pairs", "0", "--g", "1", *opts])
        with localcontext() as ctx:
            ctx.prec = 120
            sqrt2, tolerance = Decimal(2).sqrt(), Decimal(error)
            assert abs(energy(pair) + sqrt2) < tolerance
            assert abs(energy(spin) + Decimal("0.25") + sqrt2/2) < tolerance
            assert energy(isolated) == Decimal("-1.5")
            assert energy(blocked) == 4
            assert energy(vacuum) == 0
        for tables in (pair, spin, isolated, blocked, vacuum):
            check_precision(tables, precision, bits)
        assert rows(blocked["variables"])[0]["y"] is None
        assert rows(isolated["variables"])[0]["v"] is None

        for program in (richardson, central):
            full = document(program, [*base[program], *opts])
            check_precision(full, precision, bits)
            for budget in ("--max-iterations", "--max-stages"):
                failed = document(program, [*base[program], budget, "0", *opts], 2)
                check_precision(failed, precision, bits)
                state = rows(failed["states"])[0]
                assert not state["converged"] and int(state["iterations"]) == 0
                if program == richardson:
                    assert Decimal(state["reached_g"]) == 0 and energy(failed) == 2
                else:
                    assert state["reached_field"] is None and state["energy"] is None

    with localcontext() as ctx:
        ctx.prec = 120
        # Distinct physical inputs collapse to identical values if any parser
        # or solver stage narrows to a native type.
        delta = Decimal("1e-50")
        level = str(1+delta)
        close = document(richardson, ["--levels", f"1,{level}", "--pairs", "1", "--g", str(delta),
                                      "--precision", "digits-100"])
        assert abs(energy(close)-(2-Decimal(2).sqrt()*delta)) < Decimal("1e-95")
        run(richardson, ["--levels", f"1,{level}", "--pairs", "1", "--g", str(delta),
                         "--precision", "fp64"], 1)
        spin = document(central, ["--couplings", "1", "--field", level, "--sz", "0",
                                   "--precision", "digits-100"])
        assert abs(energy(spin)+Decimal("0.25")+(1+Decimal(level)**2).sqrt()/2) < Decimal("1e-95")

    tolerance = "1.00000000000000000000000000000000000000001e-80"
    for program in (richardson, central):
        first = document(program, [*base[program], "--tolerance", tolerance, "--precision", "digits-100"])
        last = document(program, [*base[program], "--precision", "bits-333", "--tolerance", tolerance])
        assert first["states"]["rows"] == last["states"]["rows"]
        with localcontext() as ctx:
            ctx.prec = 120
            got = Decimal(first["states"]["metadata"]["Residual tolerance"])
            assert abs(got/Decimal(tolerance)-1) < Decimal("1e-95")

    with tempfile.TemporaryDirectory() as directory:
        for i, program in enumerate((richardson, central)):
            path = Path(directory)/f"variables-{i}.tsv"
            screen = run(program, [*base[program], "--precision", "digits-100", "--format", "plain",
                                    "--tsv-table", f"variables={path}"]).stdout
            assert "not occupations" not in screen
            assert "333" in path.read_text()  # metadata, plus high-precision variables

with tempfile.TemporaryDirectory() as directory:
    for program in (richardson, central):
        path = Path(directory)/"preserve.json"
        path.write_text("sentinel")
        invalid = [["--precision", value] for value in ("bits-0", "digits--1", "mpfr", "bits-128junk")]
        if enabled:
            invalid += [["--precision", "bits-128", "--tolerance", value] for value in ("0", "-1", "nan", "inf")]
        else:
            invalid += [["--precision", value] for value in ("bits-128", "digits-100")]
        for args in invalid:
            run(program, [*base[program], *args, "--json", str(path), "--force"], 1)
            assert path.read_text() == "sentinel"

print("Richardson and central-spin precision contracts passed")
