"""Explicit q-boson precision, independent energies, scans and export contracts."""
import csv
from decimal import Decimal, localcontext
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, enabled = sys.argv[1:]
enabled = enabled.upper() in ("ON", "TRUE", "1")
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")
base = ["2", "--particles", "2", "--eta", "1"]


def run(args, code=0):
    result = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=60)
    assert result.returncode == code, (args, result.returncode, result.stdout, result.stderr)
    return result


def tables(args, code=0):
    doc = json.loads(run([*args, "--format", "json"], code).stdout)
    assert doc["status"] == "complete"
    for t in doc["tables"].values():
        assert t["summary"]["Outcome"] == ("success" if code == 0 else "partial")
    return doc["tables"]


def rows(table):
    result = []
    bits = int(table["metadata"]["Working precision (bits)"])
    for row in table["rows"]:
        decoded = []
        for value in row:
            if isinstance(value, dict):
                assert value["precision"] == bits
                value = value["value"]
            decoded.append(value)
        result.append(dict(zip([c["id"] for c in table["columns"]], decoded)))
    return result


with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)
    protected = path / "protected.json"
    protected.write_text("sentinel")
    if not enabled:
        for mode in ("bits-128", "digits-100"):
            result = run([*base, "--precision", mode, "--json", str(protected), "--force"], 1)
            assert "UNI20_ENABLE_MPREAL=ON" in result.stderr
            assert "UNI20_BACKEND_MPLAPACK=ON" in result.stderr
            assert protected.read_text() == "sentinel"
    else:
        with localcontext() as ctx:
            ctx.prec = 200
            # Independent eigenvalue of the two-site occupation Hamiltonian.
            eta_text = "1.000000000000000000000000000000000000000000000000000000000000123456789"
            eta = Decimal(eta_text)
            loss = 1 - (-2 * eta).exp()
            exact = 2 * loss / (1 + (1 - loss / 2).sqrt())
            previous = Decimal(1)
            for mode, bits in (("bits-128", 128), ("bits-256", 256), ("digits-100", 333)):
                flags = ["--precision", mode, "--roots"]
                args = ["2", "--particles", "2", "--eta", eta_text, *flags]
                t = tables(args)
                meta = t["states"]["metadata"]
                assert int(meta["Working precision (bits)"]) == bits
                assert abs(Decimal(meta["eta=log(q)"]) - eta) < Decimal(2) ** (2 - bits)
                e = Decimal(rows(t["states"])[0]["energy"])
                error = abs(e - exact)
                assert error < 1024 * Decimal(2) ** (1 - bits)
                assert error < previous * Decimal("1e-15"), (bits, error, previous)
                previous = error
                for threads in (1, 4):
                    for budget, code in ((10000, 0), (0, 2)):
                        scan = ["4", "--particles", "3", "--eta", "1", *flags,
                                "--excitations", "all", "--max-iterations", str(budget)]
                        serial = tables(scan, code)
                        parallel = tables([*scan, "--threads", str(threads), "--no-retain"], code)
                        assert serial.keys() == parallel.keys()
                        for name in serial:
                            assert serial[name]["rows"] == parallel[name]["rows"]
                        assert int(serial["states"]["metadata"]["Candidate count"]) == 20
                        if code:
                            assert "failed" in serial
                            assert rows(serial["reference"])[0]["energy"] is None
                        else:
                            assert len(rows(serial["states"])) == 20
                for n in (0, 1, 4):
                    free = tables(["4", "--particles", str(n), "--eta", "0", *flags])
                    assert rows(free["states"])[0]["energy"] == "0"
                    assert rows(free["states"])[0]["gap"] == "0"
                phase = rows(tables(["2", "--particles", "2", "--phase", "--max-iterations", "0", *flags])["states"])[0]
                assert abs(Decimal(phase["energy"]) - (4 - 2 * Decimal(2).sqrt())) < 1024 * Decimal(2) ** (1 - bits)
                assert phase["iterations"] == "0"
            for eta_text in ("1e-60", "1e-7000"):
                t = tables(["2", "--particles", "2", "--eta", eta_text, "--precision", "digits-100", "--roots"])
                e = Decimal(rows(t["states"])[0]["energy"])
                eta = Decimal(eta_text)
                assert e > 0 and abs(e/(2*eta)-1) < Decimal("1e-55")
                if eta_text == "1e-60":
                    assert abs((e/(2*eta)-1)/eta + Decimal(3)/4) < Decimal("1e-35")
                assert rows(t["roots"])[0]["deviation"] != rows(t["roots"])[1]["deviation"]
        args = [*base, "--excitations", "all", "--precision", "digits-100", "--threads", "4"]
        export = ["--json", str(path/"result.json"), "--quiet"]
        for name in ("states", "reference", "roots"):
            for suffix in ("csv", "tsv"):
                export += [f"--{suffix}-table", f"{name}={path / (name + '.' + suffix)}"]
        run([*args, *export])
        doc = json.loads((path/"result.json").read_text())["tables"]
        for name in doc:
            expected = [{k: "" if v is None else str(v).lower() if isinstance(v, bool) else str(v)
                         for k, v in row.items()} for row in rows(doc[name])]
            for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
                text = (path/(name + "." + suffix)).read_text()
                actual = list(csv.DictReader(io.StringIO("\n".join(line for line in text.splitlines()
                                                    if not line.startswith("#"))), delimiter=delimiter))
                assert actual == expected
        run([*args, *export, "--force", "--no-retain", "--stream"])
        refreshed = json.loads((path/"result.json").read_text())["tables"]
        for name in doc:
            assert doc[name]["rows"] == refreshed[name]["rows"]
        for args in ([*base, "--precision", "digits-100", "--tolerance", "0"],
                     [*base, "--precision", "digits-100", "--excitations", "all", "--max-candidates", "1"],
                     ["2", "--particles", "2", "--eta", "nan", "--precision", "bits-128"]):
            run([*args, "--json", str(protected), "--force"], 1)
            assert protected.read_text() == "sentinel"
print("q-boson explicit-precision CLI contracts passed")
