"""Observable grids: native values, row order, failures and exports survive scheduling."""
import csv
from decimal import Decimal
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

fp128, hubbard, lieb_liniger = sys.argv[1:]
precisions = ["fp64", "long-double"] + (["fp128"] if fp128.upper() in ("ON", "TRUE", "1") else [])
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")


def run(program, args, code=0):
    p = subprocess.run([program, *args], text=True, capture_output=True, env=env, timeout=60)
    assert p.returncode == code, (program, args, p.returncode, p.stdout, p.stderr)
    if code == 0:
        assert not p.stderr, p.stderr
    return p.stdout


def decode(text):
    doc = json.loads(text, parse_float=Decimal)
    # The Hubbard frontend has a single-table document; Lieb-Liniger has named tables.
    if "tables" in doc:
        assert doc["status"] == "complete"
        return doc["tables"]["dispersion"]
    return doc


def table(program, args, code=0):
    t = decode(run(program, [*args, "--format", "json"], code))
    assert t["summary"]["Outcome"] == ("success" if code == 0 else "partial")
    assert Decimal(t["summary"]["Run CPU seconds"]) >= 0
    assert Decimal(t["summary"]["CPU time"].split()[0]) >= 0
    assert t["metadata"]["CPU timing scope"] == "background + grids (including row output)"
    return t


def same_tables(program, args, code=0):
    serial = table(program, [*args, "--threads", "1"], code)
    for execution, limit in (([], 1), (["--threads", "2"], 2), (["--threads", "4"], 4),
                             (["--threads", "4", "--no-retain"], 4)):
        t = table(program, [*args, *execution], code)
        # Compare emitted native-precision decimal strings, never narrowed binary64 values.
        assert t["columns"] == serial["columns"], (program, args)
        assert t["rows"] == serial["rows"], (program, args, execution)
        assert t["metadata"]["Scheduler concurrency limit"] == str(limit)
    return serial


cases = ((hubbard, ["--u", "4"], "--max-iterations"),
         (hubbard, ["--u", "4", "--density", "0.5"], "--max-iterations"),
         (lieb_liniger, ["--c", "4"], "--max-iterations"))
for program, base, budget in cases:
    assert "--threads" in run(program, ["--help"])
    for precision in precisions:
        args = [*base, "--precision", precision, "--points", "5"]
        t = same_tables(program, args)
        status_column = next(i for i, c in enumerate(t["columns"]) if c["id"] == "status")
        assert all(row[status_column] == "converged" for row in t["rows"])
        failed = same_tables(program, [*args, budget, "0"], 2)
        assert len(failed["rows"]) == len(t["rows"])
        assert any(row[status_column] == "momentum_limit" for row in failed["rows"])
        branch = "type-i" if program == lieb_liniger else "spinon"
        same_tables(program, [*base, "--precision", precision, "--branch", branch, "--momentum", "0.3"])
    if program == lieb_liniger or "--density" in base:
        failed = same_tables(program, [*base, "--points", "3", "--max-background-iterations", "0"], 2)
        column = next(i for i, c in enumerate(failed["columns"]) if c["id"] == "status")
        assert all(row[column] == "density_limit" for row in failed["rows"])

# More points than the admission window, and a custom momentum range.
same_tables(lieb_liniger, ["--c", "4", "--branch", "type-i", "--p-max", "20", "--points", "129"])
same_tables(hubbard, ["--u", "4", "--branch", "holon", "--momentum", "1",
                      "--max-evaluations", "0"], 2)

# One requested momentum can now overlap different Hubbard branches. Several
# branch boundaries also lie inside the same admission window on a short grid.
for base in (["--u", "4"], ["--u", "4", "--density", "0.5"]):
    for precision in precisions:
        same_tables(hubbard, [*base, "--momentum", "1", "--precision", precision])
    args = [*base, "--points", "65", "--max-iterations", "0"]
    failed = same_tables(hubbard, args, 2)
    assert len(failed["rows"]) == 195
    column = next(i for i, c in enumerate(failed["columns"]) if c["id"] == "branch")
    branches = ["spinon", "holon", "charge-particle" if "--density" in base else "antiholon"]
    assert [row[column] for row in failed["rows"]] == [b for b in branches for _ in range(65)]

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)/"preserve.csv"
    path.write_text("preserve\n")
    overflow = subprocess.run([hubbard, "--u", "4", "--points", str((2*sys.maxsize+1)//3+1),
                               "--csv", str(path), "--force"],
                              text=True, capture_output=True, env=env, timeout=60)
    assert overflow.returncode == 1 and "selected branch grid is too large" in overflow.stderr
    assert not overflow.stdout
    assert path.read_text() == "preserve\n"

for program, base, _ in cases:
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory)
        args = [*base, "--points", "5", "--threads", "4"]
        expected = table(program, args)
        # Stream all supported file formats without retaining the table in memory.
        outputs = {fmt: path / ("grid." + fmt) for fmt in ("json", "csv", "tsv")}
        exports = [part for fmt, file in outputs.items() for part in ("--" + fmt, str(file))]
        run(program, [*args, *exports, "--quiet", "--no-retain"])
        assert decode(outputs["json"].read_text())["rows"] == expected["rows"]
        for fmt, delimiter in (("csv", ","), ("tsv", "\t")):
            contents = outputs[fmt].read_text()
            assert "# Scheduler concurrency limit: 4" in contents
            rows = list(csv.reader((line for line in contents.splitlines() if not line.startswith("#")),
                                   delimiter=delimiter))
            serial = run(program, [*base, "--points", "5", "--format", fmt, "--no-preamble"])
            assert rows == list(csv.reader(serial.splitlines(), delimiter=delimiter))

        # Reject invalid requests before touching existing exports, even with --force.
        keep = outputs["json"].read_text()
        invalid = [[*base, "--threads", value] for value in ("0", "-1", "1.5", "invalid")]
        invalid += [[*base, "--threads", "4", *bad] for bad in
                    (["--points", "1"], ["--momentum", "1", "--points", "3"],
                     ["--density", "0"], ["--tolerance", "0"])]
        for bad in invalid:
            run(program, [*bad, "--json", str(outputs["json"]), "--force"], 1)
            assert outputs["json"].read_text() == keep

print("Parallel observable-grid CLI contracts passed")
