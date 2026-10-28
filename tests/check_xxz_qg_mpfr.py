"""MPFR regular and free quantum-group XXZ frontend contracts."""
import csv
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, enabled = sys.argv[1:]
enabled = enabled.upper() in ("ON", "TRUE", "1")
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")
base = ["4", "--delta", "0.25"]


def run(args, code=0):
    result = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=60)
    assert result.returncode == code, (args, result.returncode, result.stdout, result.stderr)
    return result


def tables(args, code=0):
    doc = json.loads(run([*args, "--format", "json"], code).stdout)
    assert doc["status"] == "complete"
    for table in doc["tables"].values():
        assert table["summary"]["Outcome"] == ("success" if code == 0 else "partial")
        bits = int(table["metadata"]["Working precision (bits)"])
        for row in table["rows"]:
            for value in row:
                if isinstance(value, dict):
                    assert value["precision"] == bits
    return doc["tables"]


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]],
                     [v["value"] if isinstance(v, dict) else v for v in row])) for row in table["rows"]]


with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)
    protected = path / "protected.json"
    protected.write_text("sentinel")
    if not enabled:
        for mode in ("bits-128", "digits-100"):
            result = run([*base, "--precision", mode, "--json", str(protected), "--force"], 1)
            assert "UNI20_ENABLE_MPREAL=ON" in result.stderr
            assert protected.read_text() == "sentinel"
    else:
        with localcontext() as ctx:
            ctx.prec = 200
            dtext = "0.250000000000000000000000000000000000000000000000000000000000123456789"
            d = Decimal(dtext)
            exact = -3*d/4-(d*d+2).sqrt()/2
            previous = Decimal(1)
            for mode, bits in (("bits-128", 128), ("bits-256", 256), ("digits-100", 333)):
                flags = ["--precision", mode]
                result = tables(["4", "--delta", dtext, "--roots", *flags])
                meta = result["state"]["metadata"]
                assert int(meta["Working precision (bits)"]) == bits
                assert abs(Decimal(meta["Delta"])-d) < Decimal(2)**(2-bits)
                error = abs(Decimal(rows(result["state"])[0]["energy"])-exact)
                assert error < 2048*Decimal(2)**(1-bits)
                assert error < previous*Decimal("1e-15")
                previous = error
                free = rows(tables(["4", "--delta", "0", *flags])["blocks"])
                assert [int(r["block_size"]) for r in free] == [1, 2, 2, 1]
                assert free[0]["energy"] == free[-1]["energy"] == "0"
                assert abs(Decimal(free[2]["energy"])+Decimal(2).sqrt()/2) < 64*Decimal(2)**(1-bits)
                odd = tables(["3", "--delta", "0", "--sz", "-1/2", *flags])["blocks"]
                assert odd["metadata"]["Sector dimension"] == "3"
                vacuum = rows(tables([*base, "--numbers", "none", *flags])["state"])[0]
                assert vacuum["energy_shift"] == vacuum["residual"] == "0"
                scan = ["8", "--delta", "0.6", "--through-lines", "4", "--excitations", "all", "--roots", *flags]
                full = tables(scan)
                assert len(rows(full["levels"])) == 10
                budget = min(int(r["iterations"]) for r in rows(full["levels"]))
                for maximum, code in ((10000, 0), (0, 2), (budget, 2)):
                    args = [*scan, "--max-iterations", str(maximum)]
                    expected = tables(args, code)
                    if maximum == budget:
                        assert 0 < len(rows(expected["levels"])) < 10
                    for workers in (1, 4):
                        actual = tables([*args, "--threads", str(workers), "--no-retain"], code)
                        assert actual.keys() == expected.keys()
                        assert all(actual[name]["rows"] == expected[name]["rows"] for name in expected)
            for d in (Decimal("1e-60"), 1-Decimal("1e-60")):
                result = tables(["4", "--delta", str(d), "--precision", "digits-100", "--roots"])
                assert len(rows(result["roots"])) == 2
                expected = -3*d/4-(d*d+2).sqrt()/2
                assert abs(Decimal(rows(result["state"])[0]["energy"])-expected) < Decimal("1e-95")
            run(["4", "--delta", "1e-60", "--precision", "bits-128"], 1)
        for args, names in (([*base, "--roots"], ("state", "roots")),
                            (["4", "--delta", "0"], ("blocks",))):
            args += ["--precision", "digits-100"]
            export = ["--json", str(path/"result.json"), "--quiet", "--force"]
            for name in names:
                for suffix in ("csv", "tsv"):
                    export += [f"--{suffix}-table", f"{name}={path / (name + '.' + suffix)}"]
            run([*args, *export])
            doc = json.loads((path/"result.json").read_text())["tables"]
            for name in names:
                expected = [{k: "" if v is None else str(v).lower() if isinstance(v, bool) else str(v)
                             for k, v in r.items()} for r in rows(doc[name])]
                for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
                    text = (path/(name+"."+suffix)).read_text()
                    actual = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
                    assert actual == expected
            run([*args, *export, "--stream", "--no-retain"])
            again = json.loads((path/"result.json").read_text())["tables"]
            assert all(again[name]["rows"] == doc[name]["rows"] for name in names)
        for args in ([*base, "--tolerance", "0"], ["4", "--delta", "1"],
                     ["4", "--delta", "0", "--max-blocks", "3"],
                     ["8", "--delta", "0.6", "--through-lines", "4", "--excitations", "all", "--max-candidates", "9"]):
            run([*args, "--precision", "digits-100", "--json", str(protected), "--force"], 1)
            assert protected.read_text() == "sentinel"
print("quantum-group XXZ MPFR frontend contracts passed")
