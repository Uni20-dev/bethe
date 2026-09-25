"""Non-Hermitian quantum-group XXZ frontend contracts and native energies."""
import csv
from decimal import Decimal, localcontext
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, fp128 = sys.argv[1:]
precisions = ["fp64", "long-double"] + (["fp128"] if fp128.upper() in ("ON", "TRUE", "1") else [])
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")
base = ["4", "--delta", "0.25"]


def run(args, status=0):
    result = subprocess.run([program, *args], text=True, capture_output=True, env=env, timeout=30)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    if status == 0:
        assert not result.stderr
    return result.stdout


def records(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]


def tables(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete"
    result = doc["tables"]
    assert list(result) == (["state", "roots"] if "--roots" in args else ["state"])
    for table in result.values():
        assert table["summary"]["Outcome"] == ("success" if status == 0 else "partial")
        assert Decimal(table["summary"]["Run CPU seconds"]) >= 0
        assert "Sz_1-Sz_N" in table["metadata"]["Hamiltonian"]
        assert "no completeness" in table["metadata"]["Scope"]
    return result


for precision in precisions:
    flags = ["--precision", precision]
    result = tables([*base, "--roots", *flags])
    state = records(result["state"])[0]
    with localcontext() as ctx:
        ctx.prec = 80
        d = Decimal("0.25")
        exact = -3*d/4-(d*d+2).sqrt()/2
        tol = {"fp64": Decimal("1e-13"), "long-double": Decimal("1e-16"), "fp128": Decimal("1e-31")}[precision]
        assert abs(Decimal(state["energy"])-exact) < tol
        assert abs(Decimal(state["energy"])-Decimal(state["energy_shift"])-3*d/4) < tol
    assert state["converged"] and state["status"] == "converged"
    roots = records(result["roots"])
    assert [r["I"] for r in roots] == [1, 2]
    assert all(r["converged"] and Decimal(r["lambda"]) > 0 for r in roots)
    odd = tables(["3", "--delta", "0.25", *flags])
    assert odd["state"]["metadata"]["ell=N-2M"] == "1"
    assert abs(Decimal(records(odd["state"])[0]["energy"])+Decimal("0.625")) < Decimal("1e-14")
    vacuum = tables([*base, "--numbers", "none", "--roots", "--max-iterations", "0", *flags])
    assert vacuum["roots"]["rows"] == []
    assert Decimal(records(vacuum["state"])[0]["energy"]) == Decimal("0.1875")
    selected = tables(["8", "--delta", "0.6", "--numbers", "1,3", "--roots", *flags])
    assert [r["I"] for r in records(selected["roots"])] == [1, 3]
    sea = tables(["8", "--delta", "0.6", "--through-lines", "4", "--roots", *flags])
    assert [r["I"] for r in records(sea["roots"])] == [1, 2]
    failed = tables([*base, "--roots", "--max-iterations", "0", *flags], 2)
    row = records(failed["state"])[0]
    assert row["energy"] is None and row["energy_shift"] is None and row["status"] == "iteration_limit"
    assert all(not r["converged"] for r in records(failed["roots"]))

def free_tables(args):
    doc = json.loads(run([*args, "--format", "json"]))
    assert doc["status"] == "complete" and list(doc["tables"]) == ["blocks"]
    table = doc["tables"]["blocks"]
    assert table["summary"]["Outcome"] == "success"
    assert Decimal(table["summary"]["Run CPU seconds"]) >= 0
    assert "complete fixed-Sz" in table["metadata"]["Scope"]
    assert "Sz_1-Sz_N" in table["metadata"]["Hamiltonian"]
    return table


endpoint = ["4", "--delta", "0"]
for precision in precisions:
    table = free_tables([*endpoint, "--precision", precision])
    rows = records(table)
    assert [int(r["block_size"]) for r in rows] == [1, 2, 2, 1]
    assert [r["modes"] for r in rows] == ["1,3", "1", "3", ""]
    assert [int(r["zero_occupation"]) for r in rows] == [0, 1, 1, 2]
    with localcontext() as ctx:
        ctx.prec = 80
        exact = -Decimal(2).sqrt()/2
        tol = {"fp64": Decimal("1e-14"), "long-double": Decimal("1e-18"), "fp128": Decimal("1e-32")}[precision]
        assert abs(Decimal(rows[2]["energy"])-exact) < tol
    odd = free_tables(["3", "--delta", "0", "--sz", "-1/2", "--precision", precision])
    assert odd["metadata"]["Down spins"] == "2"
    assert odd["metadata"]["Sector dimension"] == "3"
    assert all(int(r["block_size"]) == 1 for r in records(odd))
    polarized = records(free_tables([*endpoint, "--sz", "-2", "--precision", precision]))
    assert len(polarized) == 1 and Decimal(polarized[0]["energy"]) == 0
    assert int(polarized[0]["zero_occupation"]) == 2 and int(polarized[0]["block_size"]) == 1
for n in range(2, 9):
    for down in range(n+1):
        table = free_tables([str(n), "--delta", "0", "--sz", str((n-2*down)/2)])
        rows = records(table)
        assert sum(int(r["block_size"]) for r in rows) == math.comb(n, down)
        assert int(table["metadata"]["Sector dimension"]) == math.comb(n, down)
        assert all(len(r["modes"].split(","))*(r["modes"] != "")+int(r["zero_occupation"]) == down for r in rows)

invalid = [[], ["1", "--delta", "0.25"], ["4", "--delta", "0", "--roots"], ["4", "--delta", "1"],
           ["4", "--delta", "-0.25"], ["4", "--delta", "nan"], ["4", "--delta", "inf"],
           [*base, "--through-lines", "1"], [*base, "--through-lines", "6"],
           [*base, "--numbers", "1", "--through-lines", "2"], [*base, "--numbers", "1/2"],
           [*base, "--numbers", "2,1"], [*base, "--numbers", "1,1"], [*base, "--numbers", "0"],
           [*base, "--numbers", "3"], [*base, "--numbers", "1,2,3"], [*base, "--tolerance", "0"],
           [*base, "--max-iterations", "-1"], [*base, "--sz", "0"]]
invalid += [[*endpoint, *flags] for flags in (
    ["--numbers", "none"], ["--through-lines", "0"], ["--tolerance", "1e-8"],
    ["--max-iterations", "10000"], ["--sz", "1/2"], ["--sz", "3"], ["--sz", "-3"],
    ["--max-blocks", "3"], ["--max-mode-entries", "3"], ["--max-blocks", "0"],
    ["--max-mode-entries", "-1"])]
invalid += [[*base, "--max-blocks", "10"], [*base, "--max-mode-entries", "10"], ["1", "--delta", "0"]]
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)
    flags = ["--json", str(path/"state.json")]
    for suffix in ("csv", "tsv"):
        for name in ("state", "roots"):
            flags += [f"--{suffix}-table", f"{name}={path/f'{name}.{suffix}'}"]
    screen = run([*base, "--roots", *flags, "--format", "plain"])
    assert "CPU time" in screen and "Used for:" not in screen
    exported = json.loads((path/"state.json").read_text())["tables"]
    for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
        for name in ("state", "roots"):
            text = (path/f"{name}.{suffix}").read_text()
            expected = [{k: "" if v is None else str(v).lower() if isinstance(v, bool) else str(v)
                         for k,v in r.items()} for r in records(exported[name])]
            assert list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter)) == expected
    old = (path/"state.json").read_text()
    for args in invalid:
        run([*args, "--json", str(path/"state.json"), "--force"], 1)
        assert (path/"state.json").read_text() == old
    run([*base, "--roots", *flags, "--force", "--quiet", "--no-retain"])
    refreshed = json.loads((path/"state.json").read_text())["tables"]
    assert all(refreshed[name]["rows"] == exported[name]["rows"] for name in exported)
    run([*base, "--roots", "--stream", "--no-retain", "--format", "plain"])
    run([*base, "--roots", "--format", "pretty"])
    failed = run([*base, "--max-iterations", "0", "--format", "csv", "--no-preamble"], 2)
    assert next(csv.DictReader(failed.splitlines()))["energy"] == ""
    flags = ["--json", str(path/"endpoint.json"), "--csv", str(path/"endpoint.csv"),
             "--tsv", str(path/"endpoint.tsv")]
    screen = run([*endpoint, *flags, "--format", "plain"])
    assert "CPU time" in screen and "Used for:" not in screen
    expected = records(json.loads((path/"endpoint.json").read_text())["tables"]["blocks"])
    for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
        content = (path/f"endpoint.{suffix}").read_text()
        rows = list(csv.DictReader((line for line in content.splitlines() if not line.startswith("#")), delimiter=delimiter))
        assert rows == [{k: str(v) for k,v in r.items()} for r in expected]
    run([*endpoint, *flags, "--quiet", "--force", "--no-retain"])
    assert records(json.loads((path/"endpoint.json").read_text())["tables"]["blocks"]) == expected
    run([*endpoint, "--stream", "--no-retain", "--format", "plain"])
    run([*endpoint, "--format", "pretty"])

print("Critical quantum-group XXZ CLI contracts passed")
