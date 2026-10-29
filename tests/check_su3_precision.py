"""Explicit precision and independent analytic references for finite SU(3)."""
import csv
from decimal import Decimal, localcontext
import io
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
    result = subprocess.run([program, *args], env=env, capture_output=True, text=True, timeout=60)
    assert result.returncode == code, (args, result.stdout, result.stderr)
    return result


def document(args, code=0):
    return json.loads(run([*args, "--format", "json"], code).stdout)["tables"]


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]],
                     [v["value"] if isinstance(v, dict) else v for v in row]))
            for row in table["rows"]]


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
    for precision, bits in (("bits-128", 128), ("bits-256", 256), ("digits-100", 333)):
        opts = ["--precision", precision, "--roots"]
        three, six = (document([str(n), *opts]) for n in (3, 6))
        check_precision(three, precision, bits)
        check_precision(six, precision, bits)
        with localcontext() as ctx:
            ctx.prec = 140
            bound = Decimal(2)**(12-bits)
            assert abs(Decimal(rows(three["states"])[0]["energy"])+3) < bound
            assert abs(Decimal(rows(six["states"])[0]["energy"])+1+Decimal(13).sqrt()) < bound
            assert abs(Decimal(rows(three["first_roots"])[1]["rapidity"])-1/Decimal(12).sqrt()) < bound
        assert Decimal(rows(three["second_roots"])[0]["rapidity"]) == 0
        assert Decimal(rows(six["states"])[0]["p"]) == 0
        failed = document(["9", *opts, "--max-iterations", "0"], 2)
        check_precision(failed, precision, bits)
        state = rows(failed["states"])[0]
        assert not state["converged"] and int(state["iterations"]) == 0

    tolerance = "1.00000000000000000000000000000000000000001e-80"
    first = document(["6", "--tolerance", tolerance, "--precision", "digits-100"])
    last = document(["6", "--precision", "bits-333", "--tolerance", tolerance])
    assert first["states"]["rows"] == last["states"]["rows"]
    with localcontext() as ctx:
        ctx.prec = 140
        got = Decimal(first["states"]["metadata"]["Residual tolerance"])
        assert abs(got/Decimal(tolerance)-1) < Decimal("1e-95")

    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory)
        args = ["6", "--precision", "digits-100", "--no-preamble", "--format", "plain",
                "--csv", str(path/"states.csv"), "--json", str(path/"states.json"),
                "--tsv-table", f"first_roots={path/'roots.tsv'}"]
        run(args)
        exported = json.loads((path/"states.json").read_text())
        if "tables" in exported:
            exported = exported["tables"]["states"]
        ref = rows(exported)[0]
        csv_rows = list(csv.DictReader(io.StringIO((path/"states.csv").read_text())))
        assert Decimal(csv_rows[0]["energy"]) == Decimal(ref["energy"])
        root_rows = list(csv.DictReader(io.StringIO((path/"roots.tsv").read_text()), delimiter="\t"))
        expected = document(["6", "--precision", "digits-100", "--roots"])
        assert [Decimal(r["rapidity"]) for r in root_rows] == [
            Decimal(r["rapidity"]) for r in rows(expected["first_roots"])]

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)/"preserve.csv"
    invalid = [["--precision", p] for p in ("bits-0", "digits--2")]
    if not enabled:
        invalid += [["--precision", p] for p in ("bits-128", "digits-100")]
    else:
        invalid += [["--precision", "digits-100", "--tolerance", t] for t in ("0", "-1", "nan")]
    for args in invalid:
        path.write_text("preserve\n")
        result = run(["6", *args, "--csv", str(path), "--force"], 1)
        assert not result.stdout
        assert path.read_text() == "preserve\n"
print("SU(3) precision, analytic references and exports passed")
