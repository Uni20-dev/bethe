"""Hard-wall CLI: native precision, one-edge windows and shared table exports."""
import csv
from decimal import Decimal, localcontext
import itertools
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


def run(args, status=0):
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=30, env=env)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]


def tables(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete"
    for name, t in doc["tables"].items():
        assert t["summary"]["Outcome"] == ("success" if status == 0 else "partial")
        assert Decimal(t["summary"]["Run CPU seconds"]) >= 0
        assert t["metadata"]["Boundaries"] == "Dirichlet walls at x=0 and x=length"
        assert "no conserved total momentum" in t["metadata"]["Momentum convention"]
        assert not {c["id"] for c in t["columns"]} & {"momentum_index", "p", "momentum"}
    if "roots" in doc["tables"]:
        roots = rows(doc["tables"]["roots"])
        for name in ("states", "reference", "failed"):
            if name not in doc["tables"]:
                continue
            for state in rows(doc["tables"][name]):
                own = [r for r in roots if r["state_id"] == state["state_id"]]
                k = [float(r["k"]) for r in own]
                assert all(x > 0 for x in k) and k == sorted(set(k))
                assert abs(sum(x*x for x in k)-float(state["energy"])) < 1e-11*(1+float(state["energy"]))
                if state["converged"]:
                    meta = doc["tables"][name]["metadata"]
                    length, c = float(meta["Length"]), float(meta["c"])
                    for j, x in enumerate(k):
                        phase = length*x + sum(math.atan((x-y)/c)+math.atan((x+y)/c)
                                               for l, y in enumerate(k) if l != j)
                        assert abs(phase-math.pi*float(own[j]["quantum_number"])) < 1e-10
    return doc["tables"]


for precision in precisions:
    base = ["--length", "1", "--c", "1", "--precision", precision, "--roots"]
    one = rows(tables(["1", *base, "--quantum-numbers", "3", "--max-iterations", "0"])["states"])[0]
    with localcontext() as ctx:
        ctx.prec = 70
        pi = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406")
        tol = {"fp64": Decimal("1e-13"), "long-double": Decimal("1e-16"), "fp128": Decimal("1e-30")}[precision]
        assert abs(Decimal(one["energy"])-9*pi*pi) < tol
    assert one["converged"] and one["iterations"] == "0"
    assert rows(tables(["0", *base, "--quantum-numbers", "none"])["states"])[0]["energy"] == "0"
    for n in (2, 4, 8):
        single = tables([str(n), *base])
        explicit = tables([str(n), *base, "--quantum-numbers", ",".join(map(str, range(1,n+1)))])
        assert single["states"]["rows"] == explicit["states"]["rows"]
    for n, padding in ((0,2), (1,2), (2,2), (4,1)):
        args = [str(n), *base, "--excitations", "all", "--padding", str(padding)]
        full = tables(args)
        states, roots = rows(full["states"]), rows(full["roots"])
        assert len(states) == math.comb(n+padding, n)
        seen = {tuple(int(Decimal(str(r["quantum_number"]))) for r in roots if r["state_id"] == s["state_id"])
                for s in states}
        assert seen == set(itertools.combinations(range(1,n+padding+1), n))
        assert [Decimal(s["energy"]) for s in states] == sorted(Decimal(s["energy"]) for s in states)
        assert all(s["converged"] and Decimal(s["gap"]) >= 0 for s in states)
        assert full["states"]["metadata"]["Padding at upper edge"] == str(padding)
        streamed = tables([*args, "--no-retain"])
        assert all(t["rows"] == streamed[name]["rows"] for name, t in full.items())
    failed = tables(["4", *base, "--max-iterations", "0", "--excitations", "all", "--padding", "1"], 2)
    assert not rows(failed["states"])
    assert not rows(failed["reference"])[0]["converged"]
    assert rows(failed["failed"])[0]["gap"] is None
    assert rows(failed["failed"])[0]["iterations"] == "0"
    tables(["4", *base, "--max-iterations", "1"], 2)

base = ["4", "--length", "4", "--c", "1"]
for extra in (["--excitations", "all"], ["--padding", "1"], ["--max-candidates", "2"],
              ["--excitations", "all", "--padding", "1", "--max-candidates", "4"],
              ["--excitations", "0", "--padding", "1"], ["--quantum-numbers", "1,2"],
              ["--quantum-numbers", "0,1,2,3"], ["--quantum-numbers", "1,1,2,3"],
              ["--quantum-numbers", "1,2,4,3"], ["--quantum-numbers", "1/2,3/2,5/2,7/2"],
              ["--excitations", "all", "--padding", "1", "--quantum-numbers", "1,2,3,4"],
              ["--tolerance", "nan"]):
    run([*base, *extra], 1)
for c in ("0", "-1", "nan", "inf"):
    run(["4", "--length", "4", "--c", c], 1)
for length in ("0", "-1", "nan", "inf"):
    run(["4", "--length", length, "--c", "1"], 1)
help_text = run(["--help"])
assert "bethe-lieb-liniger-obc" in help_text and "N+P" in help_text
assert "[gaudin-1971]" not in help_text
refs = run(["--references"])
assert "[gaudin-1971]" in refs and "[reichert-2019]" in refs
assert "Total energy" in run(base)
assert "Positive Bethe wave numbers" in run([*base, "--roots", "--format", "pretty"])
if "fp128" not in precisions:
    assert "fp128 is unavailable" in subprocess.run([program, *base, "--precision", "fp128"], text=True, capture_output=True).stderr

with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp)
    args = [*base, "--excitations", "2", "--padding", "1", "--roots"]
    flags = ["--json", str(path / "box.json")]
    for suffix in ("csv", "tsv"):
        for name in ("states", "reference", "roots"):
            flags += [f"--{suffix}-table", f"{name}={path / f'{name}.{suffix}'}"]
    run([*args, *flags, "--quiet", "--no-retain"])
    doc = json.loads((path / "box.json").read_text())
    assert len(rows(doc["tables"]["states"])) == 2
    assert doc["tables"]["states"]["metadata"]["Candidate count"] == "5"
    for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
        for name in ("states", "reference", "roots"):
            text = (path / f"{name}.{suffix}").read_text()
            assert "Dirichlet" in text
            data = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
            assert len(data) == len(rows(doc["tables"][name]))
