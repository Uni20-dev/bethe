"""Explicit XYZ precision, complete row/export fields, and strict branch selection."""
import csv
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, available = sys.argv[1:]
enabled = available.upper() in ("ON", "TRUE", "1")
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")
base = ["--eta", "0.75", "--t", "1"]

def run(args, status=0):
    result = subprocess.run([program, *args], text=True, capture_output=True, env=env, timeout=60)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    return result.stdout

def table(args):
    document = json.loads(run([*args, "--format", "json"]))
    assert document["status"] == "complete"
    result = document["tables"]["dispersion"]
    assert result["summary"]["Outcome"] == "success"
    return result

def rows(t):
    return [dict(zip([c["id"] for c in t["columns"]],
                     (v["value"] if isinstance(v, dict) else v for v in r))) for r in t["rows"]]

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    sentinel = root / "keep.json"
    invalid = ["mpfr", "bits-0", "digits--1"] + ([] if enabled else ["bits-256", "digits-100"])
    for precision in invalid:
        sentinel.write_text("preserve")
        run([*base, "--precision", precision, "--json", str(sentinel), "--force"], 1)
        assert sentinel.read_text() == "preserve"
    if not enabled:
        sys.exit(0)
    with localcontext() as ctx:
        ctx.prec = 160
        maximum = Decimal("0.50231250248767710143496779750111542398194818522710862683886191804879928080121848509228955262596002358304688095332672859655675497491330174716")
        gap = Decimal("0.23318323716783741586561724966177663169581956873837590677566423799910344391046032261791952456024744843434643026203040665150922915381773705932")
        # Independent JKM Jacobi reference: 110/140 digit evaluations agree beyond 105 digits.
        previous = Decimal(1)
        for selector, bits in [("bits-80", 80), ("bits-160", 160), ("bits-256", 256), ("digits-100", 333)]:
            args = [*base, "--precision", selector, "--points", "8"]
            t = table(args)
            r = rows(t)
            assert len(r) == 48 and all(x["status"] == "converged" for x in r)
            assert int(t["metadata"]["Working precision (bits)"]) == bits
            assert t["metadata"]["Precision"] == selector
            error = max(abs(Decimal(t["metadata"]["Spinon maximum energy"]) - maximum),
                        abs(Decimal(t["metadata"]["Single-spinon gap"]) - gap))
            assert error < Decimal(2)**(-bits+12), (selector, error)
            assert error < previous * Decimal("1e-15"), (selector, error, previous)
            if bits >= 256:
                assert error < Decimal("1e-70")
            previous = error
            pi = Decimal(r[7]["p"])
            for family in range(6):
                for i, row in enumerate(r[family*8:(family+1)*8]):
                    end = pi if family == 0 else 2*pi
                    expected = end*Decimal(i)/7
                    assert abs(Decimal(row["p"])-expected) < Decimal(2)**(-bits+12)
                    assert abs(Decimal(row["p_over_pi"])-expected/pi) < Decimal(2)**(-bits+12)
                    cell = 2*(expected if expected < pi else expected-pi if expected < 2*pi else 0)
                    assert abs(Decimal(row["cell_momentum"])-cell) < Decimal(2)**(-bits+13)
            assert t["rows"] == table([*args, "--no-retain"])["rows"]
            selected = rows(table([*args, "--bound-state", "2"]))
            assert selected == r[:16]+r[32:]
            out = root/selector
            run([*args, "--json", str(out.with_suffix('.json')), "--csv", str(out.with_suffix('.csv')),
                 "--tsv", str(out.with_suffix('.tsv')), "--no-retain", "--quiet"])
            exported = json.loads(out.with_suffix('.json').read_text())["tables"]["dispersion"]
            assert exported["rows"] == t["rows"]
            for suffix, separator in [('.csv',','),('.tsv','\t')]:
                lines = [line for line in out.with_suffix(suffix).read_text().splitlines() if not line.startswith('#')]
                records = list(csv.DictReader(lines, delimiter=separator))
                assert len(records) == len(r)
                for actual, expected in zip(records,r):
                    for key, value in expected.items():
                        assert actual[key] == ("" if value is None else str(value)), (key,actual[key],value)
        tiny = rows(table([*base,"--precision","bits-256","--branch","spinon","--momentum","1e-80"]))[0]
        assert abs(Decimal(tiny["p"])/Decimal("1e-80")-1) < Decimal("1e-70")
        extreme = table(["--eta","0.25","--t","10000","--precision","bits-256","--points","2"])
        assert Decimal(0) < Decimal(extreme["metadata"]["Single-spinon gap"]) < Decimal("1e-10000")
        assert all(x["status"]=="converged" for x in rows(extreme))
        near = "0." + "5" + "0"*58 + "1"
        strict = rows(table(["--eta",near,"--t","1","--branch","bound","--precision","digits-100","--points","2"]))
        assert len(strict) == 4 and all(x["s"] == 1 for x in strict)
        sentinel.write_text("preserve")
        run(["--eta","0.5","--t","1","--branch","bound","--precision","bits-256","--json",str(sentinel),"--force"],1)
        assert sentinel.read_text() == "preserve"
