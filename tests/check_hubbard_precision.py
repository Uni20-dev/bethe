"""Finite Hubbard precision, symmetry mappings, provider and export contracts."""
import csv
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

pbc, obc, mpfr = sys.argv[1:]
enabled = mpfr.upper() in ("ON", "TRUE", "1")
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")


def run(program, args, code=0):
    r = subprocess.run([program, *args], env=env, capture_output=True, text=True, timeout=60)
    assert r.returncode == code, (program, args, r.stdout, r.stderr)
    return r


def tables(program, args, code=0):
    return json.loads(run(program, [*args, "--format", "json"], code).stdout)["tables"]


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]],
                     [v["value"] if isinstance(v, dict) else v for v in row])) for row in table["rows"]]


def check_precision(document, precision, bits):
    for table in document.values():
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
        for program, hopping in ((pbc, 64), (obc, 16)):
            for interaction in ("-4", "0", "0.25", "4"):
                document = tables(program, ["2", "--u", interaction, *opts])
                check_precision(document, precision, bits)
                row = rows(document["states"])[0]
                assert row["converged"]
                with localcontext() as ctx:
                    ctx.prec = 140
                    u = Decimal(interaction)
                    expected = (u-(u*u+hopping).sqrt())/2
                    assert abs(Decimal(row["energy"])-expected) < Decimal(2)**(14-bits)
            for n in (0, 6, 8):
                document = tables(program, ["4", "--particles", str(n), "--u", "-4", *opts])
                check_precision(document, precision, bits)
                row = rows(document["states"])[0]
                assert row["converged"]
                if n in (0, 8): assert Decimal(row["energy"]) == (0 if n == 0 else -16)
                else: assert row["auxiliary_roots"]
            failed = tables(program, ["4", "--u", "-0.25", *opts, "--max-iterations", "0"], 2)
            check_precision(failed, precision, bits)
            assert not rows(failed["states"])[0]["converged"]

    interaction = "4.000000000000000000000000000000000000000000000000000000000001"
    tolerance = "1.000000000000000000000000000000000000000000001e-90"
    for program in (pbc, obc):
        a = tables(program, ["2", "--u", interaction, "--tolerance", tolerance, "--precision", "digits-100"])
        b = tables(program, ["2", "--precision", "bits-333", "--u", interaction, "--tolerance", tolerance])
        assert a["states"]["rows"] == b["states"]["rows"]
        with localcontext() as ctx:
            ctx.prec = 140
            metadata = a["states"]["metadata"]
            assert abs(Decimal(metadata["U"])-Decimal(interaction)) < Decimal("1e-98")
            assert abs(Decimal(metadata["Residual tolerance"])/Decimal(tolerance)-1) < Decimal("1e-98")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            args = ["4", "--u", "4", "--precision", "digits-100", "--no-preamble", "--format", "plain",
                    "--csv", str(path/"states.csv"), "--json", str(path/"states.json"),
                    "--tsv-table", f"charge_roots={path/'roots.tsv'}"]
            run(program, args)
            exported = json.loads((path/"states.json").read_text())
            if "tables" in exported: exported = exported["tables"]["states"]
            expected = [{k: "" if v is None else str(v).lower() if isinstance(v, bool) else str(v)
                         for k, v in row.items()} for row in rows(exported)]
            assert list(csv.DictReader((path/"states.csv").read_text().splitlines())) == expected
            roots = list(csv.DictReader((path/"roots.tsv").read_text().splitlines(), delimiter="\t"))
            document = tables(program, ["4", "--u", "4", "--precision", "digits-100", "--roots"])
            assert [Decimal(r["k"]) for r in roots] == [Decimal(r["k"]) for r in rows(document["charge_roots"])]

for program in (pbc, obc):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory)/"preserve.csv"
        path.write_text("preserve\n")
        invalid = [["--precision", x] for x in ("bits-0", "digits--2")]
        if not enabled: invalid += [["--precision", x] for x in ("bits-128", "digits-100")]
        else: invalid += [["--precision", "digits-100", "--tolerance", x] for x in ("0", "nan")]
        invalid += [["--particles", "9"]]
        if program == pbc: invalid += [["--particles", "3", "--sz", "1/2"]]
        for options in invalid:
            result = run(program, ["4", "--u", "4", *options, "--csv", str(path), "--force"], 1)
            assert not result.stdout and path.read_text() == "preserve\n"
print("Finite Hubbard precision contracts passed")
