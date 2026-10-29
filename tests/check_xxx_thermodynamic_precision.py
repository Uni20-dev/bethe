"""Independent spectral accuracy and precision-preserving CLI grid contracts."""
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
base = ["--momentum", "1.5", "--omega", "1.8"]


def run(args, code=0):
    p = subprocess.run([program, *map(str, args)], text=True, capture_output=True, env=env, timeout=180)
    assert p.returncode == code, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def tables(args, bits, code=0):
    d = json.loads(run([*args, "--format", "json"], code))
    assert d["status"] == "complete"
    assert set(d["tables"]) == {"spectrum", "continuum"}
    for t in d["tables"].values():
        assert t["summary"]["Outcome"] == ("success" if code == 0 else "partial")
        assert int(t["metadata"]["Working precision (bits)"]) == bits
        for row in t["rows"]:
            for cell in row:
                if isinstance(cell, dict):
                    assert cell["precision"] == bits, cell
    return d["tables"]


def rows(t):
    return [dict(zip([c["id"] for c in t["columns"]],
                     [c["value"] if isinstance(c, dict) else c for c in row])) for row in t["rows"]]


if enabled:
    with localcontext() as ctx:
        ctx.prec = 140
        # Independent tanh-sinh oracle, exp(-x) subtraction and unit panels;
        # 105/130 decimal digit evaluations agree beyond 94 decimal places.
        reference = Decimal("0.8731304551163434158090722608166793299789868199172442958199344007169137742788878727691620633517433651166155890938972983891747")
        previous = Decimal(1)
        for bits in (80, 160, 256, 333):
            t = tables([*base, "--precision", f"bits-{bits}"], bits)
            p = rows(t["spectrum"])[0]
            error = abs(Decimal(p["density"]) - reference)
            assert p["status"] == "converged"
            assert error <= Decimal(p["estimated_error"])
            assert error < previous * Decimal("1e-14"), (bits, error, previous)
            previous = error
        assert previous < Decimal("1e-94")
        args = [*base, "--precision", "bits-160"]
        zz = rows(tables(args, 160)["spectrum"])[0]
        raising = rows(tables([*args, "--channel", "raising"], 160)["spectrum"])[0]
        scaled = rows(tables(["--momentum", "1.5", "--omega", "3.6", "--exchange", "2",
                              "--precision", "bits-160"], 160)["spectrum"])[0]
        assert abs(Decimal(raising["density"]) / Decimal(zz["density"]) - 2) < Decimal("1e-46")
        assert abs(Decimal(scaled["density"]) / Decimal(zz["density"]) - Decimal("0.5")) < Decimal("1e-46")
        failure = rows(tables([*args, "--max-evaluations", "0"], 160, 2)["spectrum"])[0]
        assert failure["density"] is None and failure["estimated_error"] is None
        assert failure["status"] == "numerical_failure"
        for omega in ("-1", "5"):
            p = rows(tables(["--momentum", "1.5", "--omega", omega, "--precision", "bits-160"], 160)["spectrum"])[0]
            assert Decimal(p["density"]) == 0 and p["status"] == "outside_continuum"

        # A non-dyadic grid larger than the admission window, including zero,
        # threshold, exterior and interior rows. Grid controls are separate
        # from arithmetic precision so this is a bounded CLI regression.
        grid = ["--momentum-points", "7", "--points", "21", "--omega-max", "3.2",
                "--precision", "bits-128", "--tolerance", "1e-12"]
        serial = tables(grid, 128)
        parallel = tables([*grid, "--threads", "2", "--no-retain"], 128)
        for name in serial:
            assert serial[name]["rows"] == parallel[name]["rows"]
        r = rows(serial["spectrum"])
        assert {"outside_continuum", "lower_threshold", "converged"} <= {p["status"] for p in r}
        assert all(abs(Decimal(p["omega"]) - Decimal("3.2") * (i % 21) / 20) < Decimal("1e-36")
                   for i, p in enumerate(r))
        failed = [*grid, "--max-evaluations", "0"]
        a, b = tables(failed, 128, 2), tables([*failed, "--threads", "2", "--no-retain"], 128, 2)
        for name in a:
            assert a[name]["rows"] == b[name]["rows"]

        q = "1.5000000000000000000000000000000000000000000000000000001"
        tolerance = "1.0000000000000000000000000000000000000000000000000000001e-20"
        t = tables(["--momentum", q, "--omega", "1.8", "--tolerance", tolerance,
                    "--precision", "digits-100"], 333)
        assert abs(Decimal(rows(t["spectrum"])[0]["q"]) - Decimal(q)) < Decimal("1e-95")
        assert abs(Decimal(t["spectrum"]["metadata"]["Kernel tolerance"]) / Decimal(tolerance) - 1) < Decimal("1e-95")

    with tempfile.TemporaryDirectory() as folder:
        args = [*base, "--precision", "bits-160"]
        expected = rows(tables(args, 160)["spectrum"])
        for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
            path = Path(folder) / f"spectral.{suffix}"
            run([*args, f"--{suffix}-table", f"spectrum={path}", "--quiet", "--no-retain"])
            exported = list(csv.DictReader([line for line in path.read_text().splitlines()
                                           if not line.startswith("#")], delimiter=delimiter))
            assert exported == expected, (exported, expected)

with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "preserve.json"
    path.write_text("sentinel")
    invalid = [["--precision", s] for s in ("bits-0", "digits--1", "bits-128junk")]
    if enabled:
        invalid += [["--precision", "bits-128", "--max-quadrature-order", n] for n in (0, 16, 33, 4096)]
        invalid += [["--precision", "bits-128", "--tolerance", t] for t in ("0", "-1", "inf", "nan")]
    else:
        invalid += [["--precision", s] for s in ("bits-128", "digits-100")]
    for args in invalid:
        run([*base, *args, "--json", path, "--force"], 1)
        assert path.read_text() == "sentinel"
print("Thermodynamic XXX runtime precision checks passed")
