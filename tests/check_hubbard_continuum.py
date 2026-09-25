"""Hubbard two-particle CLI: grids, native precision, failures and streaming exports."""
import csv
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, fp128 = sys.argv[1:]
precisions = ["fp64", "long-double"] + (["fp128"] if fp128.upper() in ("ON", "TRUE", "1") else [])
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")
base = ["--u", "4", "--points", "5"]


def run(args, status=0):
    result = subprocess.run([program, *args], text=True, capture_output=True, timeout=180, env=env)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    return result.stdout


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]


def table(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete" and list(doc["tables"]) == ["two_spinon"]
    result = doc["tables"]["two_spinon"]
    assert result["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    assert Decimal(result["summary"]["Run CPU seconds"]) >= 0
    assert "DeltaN=0" in result["metadata"]["Energy convention"]
    return result


for precision in precisions:
    flags = ["--precision", precision]
    with localcontext() as ctx:
        ctx.prec = 100
        tol = {"fp64": Decimal("1e-12"), "long-double": Decimal("1e-15"), "fp128": Decimal("1e-30")}[precision]
        maximum = 2*Decimal("1.24228148941956239984591553194945497655659639")
        grid = rows(table([*base, *flags]))
        assert len(grid) == 5
        for i, row in enumerate(grid):
            assert row["converged"] and row["status"] == "converged"
            assert abs(Decimal(row["momentum_over_pi"])-Decimal(i-2)/2) < tol
            assert Decimal(row["upper_energy"]) >= Decimal(row["lower_energy"]) >= 0
            assert abs(Decimal(row["lower_energy"])-Decimal(grid[4-i]["lower_energy"])) < tol
            assert abs(Decimal(row["upper_energy"])-Decimal(grid[4-i]["upper_energy"])) < tol
        assert Decimal(grid[2]["lower_energy"]) == Decimal(grid[2]["upper_energy"]) == 0
        assert abs(Decimal(grid[0]["upper_energy"])-maximum) < tol
        p = "0.147435530616030589362559903539606986354219786"
        reference = Decimal("0.180195419367165283255675569344877169275747983")
        single = rows(table(["--u", "4", "--momentum", "-"+p, *flags]))[0]
        assert abs(Decimal(single["lower_energy"])-reference) < tol
    partial = rows(table([*base, *flags, "--max-evaluations", "0"], 2))
    assert partial[2]["converged"]
    for i in (0, 1, 3, 4):
        assert partial[i]["lower_energy"] is None and partial[i]["upper_energy"] is None
        assert partial[i]["status"] == "quadrature_limit"
    for option, status in (("--max-levels", "quadrature_limit"), ("--max-iterations", "momentum_limit")):
        failed = rows(table(["--u", "4", "--momentum", "1", *flags, option, "0"], 2))[0]
        assert failed["status"] == status and failed["lower_energy"] is None

invalid = [[], ["--u", "0"], ["--u", "nan"], ["--u", "4", "--momentum", "4"],
           ["--u", "4", "--momentum", "nan"], [*base, "--momentum", "0"],
           ["--u", "4", "--points", "1"], ["--u", "4", "--points", "1000001"],
           [*base, "--max-levels", "25"], [*base, "--tolerance", "0"],
           [*base, "--max-evaluations", "-1"], [*base, "--density", "0.5"]]
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)
    flags = ["--json", str(path/"result.json"), "--csv", str(path/"edges.csv"), "--tsv", str(path/"edges.tsv")]
    screen = run([*base, *flags, "--format", "plain"])
    assert "CPU time" in screen and "Used for:" not in screen
    exported = json.loads((path/"result.json").read_text())["tables"]["two_spinon"]
    expected = [{k: "" if v is None else str(v).lower() if isinstance(v, bool) else str(v)
                 for k, v in row.items()} for row in rows(exported)]
    for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
        text = (path/f"edges.{suffix}").read_text()
        assert "DeltaN=0" in text
        assert list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter)) == expected
    old = (path/"result.json").read_text()
    run([*base, "--json", str(path/"result.json")], 1)
    for args in invalid:
        run([*args, "--json", str(path/"result.json"), "--force"], 1)
        assert (path/"result.json").read_text() == old
    run([*base, *flags, "--force", "--no-retain", "--quiet"])
    assert json.loads((path/"result.json").read_text())["tables"]["two_spinon"]["rows"] == exported["rows"]
    run([*base, "--stream", "--no-retain", "--format", "plain"])
    run([*base, "--format", "pretty"])
    failure = run([*base, "--max-evaluations", "0", "--format", "csv", "--no-preamble"], 2)
    assert next(csv.DictReader(failure.splitlines()))["lower_energy"] == ""

charge = ["--u", "4", "--momentum", "0", "--channel", "spinon-holon",
          "--search-tolerance", "1e-7", "--position-tolerance", "1e-5", "--tolerance", "1e-10"]


def charge_table(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete" and list(doc["tables"]) == ["charge_continuum"]
    result = doc["tables"]["charge_continuum"]
    assert result["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    assert Decimal(result["summary"]["Run CPU seconds"]) >= 0
    return result


for precision in precisions:
    result = charge_table([*charge, "--precision", precision])
    assert result["metadata"]["DeltaN"] == "-1"
    row = rows(result)[0]
    assert row["converged"] and row["search_status"] == row["constituent_status"] == "converged"
    assert abs(float(row["lower_energy"])-1.8856450004260146) < 1e-6
    assert abs(float(row["upper_energy"])-3.499612327587533) < 1e-6
    assert abs(float(row["lower_p1"])+float(row["lower_p2"])) < 1e-10
    assert int(row["quadrature_evaluations"]) > 0 and int(row["objective_evaluations"]) > 0
    for option, search_status, constituent_status in (
            ("--max-total-evaluations", "objective_failure", "quadrature_limit"),
            ("--max-evaluations", "objective_failure", "quadrature_limit"),
            ("--max-search-evaluations", "evaluation_limit", "converged"),
            ("--max-search-iterations", "iteration_limit", "converged"),
            ("--max-levels", "objective_failure", "quadrature_limit"),
            ("--max-iterations", "objective_failure", "momentum_limit")):
        failed = rows(charge_table([*charge, "--precision", precision, option, "0"], 2))[0]
        assert failed["search_status"] == search_status, failed
        assert failed["constituent_status"] == constituent_status, failed
        for key in ("lower_energy", "upper_energy", "lower_symmetric_energy", "upper_symmetric_energy",
                    "lower_p1", "lower_p2", "upper_p1", "upper_p2", "lower_error", "upper_error"):
            assert failed[key] is None, failed

# Convention changes shift energies, not the search or momentum witnesses.
shifted = rows(charge_table([*charge, "--convention", "unshifted"]))[0]
fermi = rows(charge_table([*charge, "--convention", "unshifted", "--reference", "fermi"]))[0]
assert abs(float(shifted["lower_energy"])-float(fermi["lower_energy"])+2) < 1e-12
assert shifted["lower_p1"] == fermi["lower_p1"]
anti = ["--u", "4", "--channel", "spinon-antiholon", "--momentum", "-3.141592653589793",
        "--search-tolerance", "1e-7", "--position-tolerance", "1e-5", "--tolerance", "1e-10"]
anti_row = rows(charge_table([*anti, "--convention", "unshifted"]))[0]
assert abs(float(anti_row["lower_energy"])-float(shifted["lower_energy"])-4) < 1e-6
neutral = ["--u", "4", "--channel", "holon-antiholon", "--momentum", "0"]
neutral_row = rows(charge_table(neutral))[0]  # Default fp64 search controls.
assert abs(float(neutral_row["lower_energy"])-1.2867270220129044) < 1e-10
assert abs(float(neutral_row["upper_energy"])-9.2867270220129044) < 1e-10
mesh_failed = rows(charge_table([*charge, "--max-intervals", "16"], 2))[0]
assert mesh_failed["search_status"] == "mesh_limit" and mesh_failed["lower_energy"] is None
grid_args = charge.copy()
index = grid_args.index("--momentum")
del grid_args[index:index+2]
grid = rows(charge_table([*grid_args, "--points", "3"]))
assert len(grid) == 3 and all(row["converged"] for row in grid)
assert abs(float(grid[0]["lower_energy"])-3.499612327587533) < 1e-6
assert abs(float(grid[0]["upper_energy"])-5.885645000426015) < 1e-6
assert abs(float(grid[0]["lower_energy"])-float(grid[-1]["lower_energy"])) < 1e-10

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)
    flags = ["--json", str(path/"charge.json"), "--csv", str(path/"charge.csv"), "--tsv", str(path/"charge.tsv")]
    run([*charge, *flags, "--quiet", "--no-retain"])
    exported = json.loads((path/"charge.json").read_text())["tables"]["charge_continuum"]
    expected = [{k: "" if v is None else str(v).lower() if isinstance(v, bool) else str(v)
                 for k, v in row.items()} for row in rows(exported)]
    for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
        text = (path/f"charge.{suffix}").read_text()
        assert "spinon-holon" in text
        assert list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter)) == expected
    original = (path/"charge.json").read_text()
    for option, value in (("--search-tolerance", "0"), ("--search-tolerance", "nan"),
                          ("--position-tolerance", "0"), ("--position-tolerance", "inf"),
                          ("--initial-intervals", "3"), ("--max-intervals", "8"),
                          ("--max-intervals", "4097"), ("--max-total-evaluations", "-1"),
                          ("--convention", "bad"), ("--reference", "bad"), ("--channel", "bad")):
        # Do not duplicate an already present option: CLI rejects duplicates independently.
        args = charge.copy()
        if option in args:
            index = args.index(option)
            del args[index:index+2]
        run([*args, option, value, "--json", str(path/"charge.json"), "--force"], 1)
        assert (path/"charge.json").read_text() == original
    run([*base, "--search-tolerance", "1e-7"], 1)
    grid = rows(charge_table(["--u", "4", "--channel", "holon-antiholon", "--points", "3",
                             "--max-total-evaluations", "0", "--no-retain"], 2))
    assert len(grid) == 3 and all(row["lower_p1"] is None for row in grid)
    failed_csv = run([*charge, "--max-total-evaluations", "0", "--format", "csv", "--no-preamble"], 2)
    failed_row = next(csv.DictReader(failed_csv.splitlines()))
    assert failed_row["lower_energy"] == failed_row["lower_p1"] == failed_row["lower_error"] == ""
    run([*charge, "--stream", "--no-retain", "--format", "plain"])
    run([*charge, "--format", "pretty"])

print("Hubbard continuum CLI contracts passed")
