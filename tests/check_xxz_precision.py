"""Explicit precision, physical references and exports for finite XXZ."""
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


def document(program, args, code=0):
    return json.loads(run(program, [*args, "--format", "json"], code).stdout)["tables"]


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]],
                     [v["value"] if isinstance(v, dict) else v for v in row])) for row in table["rows"]]


def check_precision(tables, precision, bits):
    for table in tables.values():
        assert table["metadata"]["Precision"] == precision
        assert table["metadata"]["Arithmetic format"] == "arbitrary binary"
        assert int(table["metadata"]["Working precision (bits)"]) == bits
        for row in table["rows"]:
            for v in row:
                if isinstance(v, dict):
                    assert v["precision"] == bits, v


if enabled:
    # A rounded-away continuation increment must terminate without waiting
    # for a Newton iteration to consume the budget. Keep the established
    # unconverged-energy output policy, with an explicit precision diagnostic.
    stalled = document(obc, ["3", "--delta", "2", "--precision", "bits-4"], 2)
    check_precision(stalled, "bits-4", 4)
    state = rows(stalled["states"])[0]
    assert not state["converged"] and "continuation stalled; increase precision" in state["status"], state
    assert state["energy"] is not None and int(state["iterations"]) == 0, state
    for precision, bits in (("bits-128", 128), ("bits-256", 256), ("digits-100", 333)):
        opts = ["--precision", precision, "--roots"]
        for delta in ("-0.5", "0", "0.5", "1", "2"):
            with localcontext() as ctx:
                ctx.prec = 140
                d = Decimal(delta)
                bound = Decimal(2)**(12-bits)
                for program, n, energy in ((pbc, "4", -d/2-(d*d/4+2).sqrt()),
                                            (obc, "3", -d/4-(d*d/16+Decimal("0.5")).sqrt())):
                    tables = document(program, [n, "--delta", delta, *opts])
                    check_precision(tables, precision, bits)
                    assert abs(Decimal(rows(tables["states"])[0]["energy"])-energy) < bound
            # Exercise the even-open boundary root, including a failed iterate.
            tables = document(obc, ["4", "--delta", delta, *opts])
            check_precision(tables, precision, bits)
            failed = document(obc, ["4", "--delta", delta, *opts, "--max-iterations", "0"], 2)
            check_precision(failed, precision, bits)
            assert not rows(failed["states"])[0]["converged"]
        for program in (pbc, obc):
            sectors = document(program, ["4", "--delta", "-0.5", "--sectors", *opts])
            check_precision(sectors, precision, bits)
            assert len(sectors["states"]["rows"]) == 5
            args = ["6", "--delta", "0.5", "--excitations", "all", "--sz", "-1", *opts]
            serial, parallel = [document(program, [*args, "--threads", str(t)]) for t in (1, 3)]
            assert set(serial) == set(parallel)
            check_precision(parallel, precision, bits)
            for name in serial:
                assert serial[name]["rows"] == parallel[name]["rows"], name
            empty = document(program, ["4", "--delta", "0.5", "--quantum-numbers", "", *opts])
            check_precision(empty, precision, bits)
            failed = document(program, [*args, "--max-iterations", "0"], 2)
            check_precision(failed, precision, bits)
            assert failed["failed"]["rows"]

    delta = "1.000000000000000000000000000000000000000000000000000000000001"
    tolerance = "1.000000000000000000000000000000000000000000000001e-90"
    for program in (pbc, obc):
        a = document(program, ["4", "--delta", delta, "--tolerance", tolerance, "--precision", "digits-100"])
        b = document(program, ["4", "--precision", "bits-333", "--delta", delta, "--tolerance", tolerance])
        assert a["states"]["rows"] == b["states"]["rows"]
        with localcontext() as ctx:
            ctx.prec = 140
            metadata = a["states"]["metadata"]
            assert abs(Decimal(metadata["Delta"])-Decimal(delta)) < Decimal("1e-98")
            assert abs(Decimal(metadata["Residual tolerance"])/Decimal(tolerance)-1) < Decimal("1e-98")

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            args = ["4", "--delta", "0.5", "--precision", "digits-100", "--no-preamble", "--format", "plain",
                    "--csv", str(path/"states.csv"), "--json", str(path/"states.json"),
                    "--tsv-table", f"roots={path/'roots.tsv'}"]
            run(program, args)
            exported = json.loads((path/"states.json").read_text())
            if "tables" in exported: exported = exported["tables"]["states"]
            assert list(csv.DictReader((path/"states.csv").read_text().splitlines())) == [
                {k: "" if v is None else str(v).lower() if isinstance(v, bool) else str(v) for k, v in row.items()}
                for row in rows(exported)]
            expected = document(program, ["4", "--delta", "0.5", "--precision", "digits-100", "--roots"])
            root_rows = list(csv.DictReader((path/"roots.tsv").read_text().splitlines(), delimiter="\t"))
            assert len(root_rows) == len(expected["roots"]["rows"])
            for name in ("z", "rapidity"):
                if name in root_rows[0]:
                    assert [Decimal(r[name]) for r in root_rows] == [Decimal(r[name]) for r in rows(expected["roots"])]
                    break
            else:
                raise AssertionError(root_rows)

for program in (pbc, obc):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory)/"preserve.csv"
        path.write_text("preserve\n")
        invalid = [["--precision", x] for x in ("bits-0", "digits--2")]
        if not enabled: invalid += [["--precision", x] for x in ("bits-128", "digits-100")]
        else: invalid += [["--precision", "digits-100", "--tolerance", x] for x in ("0", "nan")]
        for options in invalid:
            result = run(program, ["4", "--delta", "0.5", *options, "--csv", str(path), "--force"], 1)
            assert not result.stdout and path.read_text() == "preserve\n"
        if enabled:
            result = run(program, ["4", "--delta", "-0.5", "--excitations", "all", "--precision", "bits-128",
                                   "--csv", str(path), "--force"], 1)
            assert not result.stdout and path.read_text() == "preserve\n"
            if program == obc:
                result = run(program, ["4", "--delta", "2", "--precision", "bits-2",
                                       "--csv", str(path), "--force"], 1)
                assert "too low to seed continuation" in result.stderr
                assert not result.stdout and path.read_text() == "preserve\n"
print("XXZ precision contracts passed")
