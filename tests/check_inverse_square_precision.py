"""Scalar-only MPFR inverse-square models: accuracy, exports and refusal semantics."""
import csv
from decimal import Decimal as D, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

hs, su, mpfr = sys.argv[1:]
enabled = mpfr.upper() in ("ON", "TRUE", "1")
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")


def run(program, args, status=0):
    result = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=30)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    return result.stdout


def tables(program, args, status=0):
    return json.loads(run(program, [*args, "--format", "json"], status))["tables"]


def rows(table):
    names = [c["id"] for c in table["columns"]]
    return [dict(zip(names, [x["value"] if isinstance(x, dict) else x for x in row])) for row in table["rows"]]


def precision(table, bits):
    assert table["metadata"]["Arithmetic format"] == "arbitrary binary"
    assert int(table["metadata"]["Working precision (bits)"]) == bits
    for row in table["rows"]:
        for x in row:
            if isinstance(x, dict):
                assert x["precision"] == bits


if enabled:
    with localcontext() as ctx:
        ctx.prec = 150
        pi = D("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128")
        previous = [D(1), D(1)]
        length = "7.25000000000000000000000000000000000000000000000001"
        coupling = "2.25000000000000000000000000000000000000000000000001"
        for option, bits in (("bits-128", 128), ("bits-256", 256), ("digits-100", 333)):
            eps = D(2) ** (10-bits)
            t = tables(hs, ["2", "--levels", "all", "--precision", option])["levels"]
            precision(t, bits)
            r = rows(t)
            assert len(r) == 2
            error = abs(D(r[0]["energy"]) + 3*pi*pi/16)
            assert error < eps and error < previous[0]/D("1e15")
            previous[0] = error
            assert D(r[0]["gap"]) == 0 and D(r[1]["p"]) == 0
            assert abs(D(r[1]["gap"])-pi*pi/4) < eps
            odd = tables(hs, ["7", "--sz", "3/2", "--precision", option])["levels"]
            precision(odd, bits)
            assert len(rows(odd)) == 2 and rows(odd)[0]["energy"] == rows(odd)[1]["energy"]
            t = tables(su, ["2", "--length", length, "--lambda", coupling, "--labels=-3,2", "--precision", option, "--pseudomomenta"])
            for table in t.values(): precision(table, bits)
            actual = rows(t["states"])[0]
            expected = 2*(pi/D(length))**2*(1+(5+D(coupling))**2)
            error = abs(D(actual["energy"])-expected)
            assert error < eps and error < previous[1]/D("1e15")
            previous[1] = error
            assert abs(D(t["states"]["metadata"]["Length"])-D(length)) < eps
            assert abs(D(t["states"]["metadata"]["Lambda"])-D(coupling)) < eps
            for n, row in zip((-3, 2), rows(t["pseudomomenta"])):
                expected_k = 2*pi/D(length)*(n + (-1 if n == -3 else 1)*D(coupling)/2)
                assert abs(D(row["k"])-expected_k) < eps
            for program, args, name, count in ((hs, ["8", "--levels", "all"], "levels", 34),
                    (su, ["3", "--length", "4", "--lambda", "2", "--levels", "all", "--window", "2"], "states", 35)):
                base = [*args, "--precision", option]
                t = tables(program, base)[name]
                precision(t, bits)
                assert len(rows(t)) == count
                assert rows(tables(program, [*base, "--no-retain"])[name]) == rows(t)
            # Uniform boosts survive enormous E0, even when energy-E0 is zero.
            t = tables(su, ["4", "--length", "5", "--lambda", "1e1000", "--labels=1,1,1,1", "--precision", option])["states"]
            r = rows(t)[0]
            assert D(r["energy"]) == D(t["metadata"]["Ground energy"])
            assert abs(D(r["gap"])-16*pi*pi/25) < eps
            tiny = tables(su, ["1", "--length", "1e7000", "--lambda", "0", "--labels=1", "--precision", option])["states"]
            precision(tiny, bits)
            assert abs(D(rows(tiny)[0]["energy"])/D("1e-14000")-4*pi*pi) < eps
        # A physically nonzero correction below binary128 resolution.
        t = tables(su, ["2", "--length", "5", "--lambda", "1e-60", "--labels=0,1", "--precision", "digits-100"])["states"]
        slope = (D(rows(t)[0]["gap"])-4*pi*pi/25)/D("1e-60")
        assert abs(slope-4*pi*pi/25) < D("1e-38")

    for program, args, name, cap in ((hs, ["8", "--levels", "all"], "levels", "--max-motifs"),
            (su, ["3", "--length", "4", "--lambda", "2", "--levels", "all", "--window", "2"], "states", "--max-states")):
        base = [*args, "--precision", "digits-100"]
        refused = tables(program, [*base, cap, "0"], 2)[name]
        precision(refused, 333)
        assert not refused["rows"] and "budget exceeded" in refused["summary"]["Status"]
        with tempfile.TemporaryDirectory() as directory:
            expected = [{k: "" if v is None else str(v) for k, v in row.items()} for row in rows(tables(program, base)[name])]
            for extension, delimiter in (("csv", ","), ("tsv", "\t")):
                path = Path(directory)/f"levels.{extension}"
                run(program, [*base, f"--{extension}", str(path), "--quiet"])
                text = path.read_text()
                assert "333" in text and "Hamiltonian" in text
                actual = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
                assert actual == expected
            # Additional table export does not imply screen output.
            extra = "spin_content" if program == hs else "pseudomomenta"
            path = Path(directory)/"extra.csv"
            screen = tables(program, [*base, "--csv-table", f"{extra}={path}"])
            assert extra not in screen and path.exists()
            assert extra in tables(program, [*base, "--spin-content" if program == hs else "--pseudomomenta"])

for program, args in ((hs, ["8"]), (su, ["3", "--length", "4", "--lambda", "2"])):
    help_text = run(program, ["--help"])
    assert "bits-N" in help_text and "digits-N" in help_text
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory)/"preserve.json"
        path.write_text("sentinel")
        invalid = [["--precision", p] for p in ("bits-0", "digits--1", "bits-128junk")]
        if enabled:
            invalid += [["--precision", "bits-128", "--motif", "1,2"]] if program == hs else [
                ["--precision", "bits-128", "--length", "0"], ["--precision", "bits-128", "--lambda", "nan"]]
        else:
            invalid += [["--precision", p] for p in ("bits-128", "digits-100")]
        for case in invalid:
            assert not run(program, [*args, *case, "--force", "--json", str(path)], 1)
            assert path.read_text() == "sentinel"
