"""Native scalar output and exact spectral-rule CLI contracts."""
import csv
from decimal import Decimal, localcontext
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile

program, fp128 = sys.argv[1:]

def run(*args, status=0):
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=30)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    if status == 0:
        assert not p.stderr, p.stderr
    return p.stdout, p.stderr

def rows(text):
    assert "# CPU time:" in text and "References:" not in text
    return list(csv.DictReader(io.StringIO("\n".join(x for x in text.splitlines() if not x.startswith("#")))) )

precisions = ["fp64", "long-double"]
if fp128.upper() in ("ON", "TRUE", "1"):
    precisions.append("fp128")
else:
    run("4", "--precision", "fp128", status=1)
with localcontext() as ctx:
    ctx.prec = 60
    for precision in precisions:
        data = rows(run("4", "--precision", precision, "--format", "csv")[0])
        assert len(data) == 1 and data[0]["motif"] == "1 3" and data[0]["gap"] == "0"
        ref = Decimal("-2.15897596273829719787004490622290806085")
        tol = Decimal("1e-32" if precision == "fp128" else "1e-14")
        assert abs(Decimal(data[0]["energy"]) - ref) < tol, data
        data = rows(run("5", "--precision", precision, "--format", "csv")[0])
        assert len(data) == 2 and {r["momentum_index"] for r in data} == {"1", "4"}
        assert all(r["s_max"] == "0.5" and r["degeneracy"] == "2" for r in data)
data = rows(run("6", "--levels", "all", "--format", "csv")[0])
assert len(data) == 13 and sum(int(r["degeneracy"]) for r in data) == 64
assert len(rows(run("6", "--levels", "3", "--format", "csv")[0])) == 3
assert rows(run("4", "--motif", "", "--format", "csv")[0])[0]["degeneracy"] == "5"
assert len(rows(run("6", "--levels", "all", "--max-motifs", "12", "--format", "csv", status=2)[0])) == 0
for budget, status in ((13, 0), (12, 2)):
    doc = json.loads(run("6", "--levels", "all", "--max-motifs", str(budget), "--format", "json", status=status)[0])
    table = doc["tables"]["levels"]
    assert doc["status"] == "complete"
    assert table["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    assert Decimal(table["summary"]["Run CPU seconds"]) >= 0
    assert Decimal(table["summary"]["Elapsed seconds"]) >= 0
selected = json.loads(run("5", "--sz", "1/2", "--format", "json")[0])["tables"]["levels"]
assert selected["metadata"]["Selected Sz"] == "0.5"
run("8", "--sz", "-1", "--format", "pretty")
for args in [("1",), ("1000001",), ("6", "--sz", "1/2"), ("6", "--motif", "1,2"),
             ("6", "--motif", "6"), ("6", "--levels", "0"), ("6", "--levels", "all", "--sz", "0"),
             ("6", "--format", "unknown"), ("6", "--motif", "1,"), ("6", "--levels")]:
    run(*args, status=1)
def records(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]


def spin_tables(args, status=0):
    doc = json.loads(run(*args, "--spin-content", "--format", "json", status=status)[0])
    assert doc["status"] == "complete"
    assert list(doc["tables"]) == ["levels", "spin_content"]
    for table in doc["tables"].values():
        assert table["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    return doc["tables"]


for precision in precisions:
    tables = spin_tables(["6", "--motif", "3", "--precision", precision])
    content = records(tables["spin_content"])
    assert [Decimal(r["spin"]) for r in content] == [0, 1, 2]
    assert all(r["multiplicity"] == "1" and r["complete"] and r["status"] == "complete" for r in content)
    assert all(r["state_id"] == "0" for r in content)
    assert int(records(tables["levels"])[0]["degeneracy"]) == 9
    odd = records(spin_tables(["5", "--precision", precision])["spin_content"])
    assert len(odd) == 2 and all(r["spin"] == 0.5 for r in odd)
    # --sz chooses motifs, not a projection of the requested decomposition.
    selected = records(spin_tables(["4", "--sz", "1", "--precision", precision])["spin_content"])
    assert [Decimal(r["spin"]) for r in selected] == [0, 1]

all_tables = spin_tables(["6", "--levels", "all"])
content = records(all_tables["spin_content"])
for level in records(all_tables["levels"]):
    dimension = sum((2*Decimal(r["spin"])+1)*int(r["multiplicity"])
                    for r in content if r["state_id"] == level["state_id"])
    assert dimension == int(level["degeneracy"])
refused = spin_tables(["6", "--levels", "all", "--max-motifs", "0"], 2)
assert not refused["levels"]["rows"] and not refused["spin_content"]["rows"]
partial = spin_tables(["6", "--levels", "all", "--max-spin-updates", "0"], 2)
assert partial["levels"]["rows"] == all_tables["levels"]["rows"]
failed = [r for r in records(partial["spin_content"]) if not r["complete"]]
assert failed and all(r["spin"] is None and r["multiplicity"] is None and r["status"] == "work_limit" for r in failed)
overflow = records(spin_tables(["300", "--motif", ",".join(map(str, range(2, 300, 3)))], 2)["spin_content"])
assert len(overflow) == 1 and overflow[0]["status"] == "count_overflow" and overflow[0]["spin"] is None
run("4", "--max-spin-updates", "0", status=1)
run("4", "--spin-content", "--max-spin-updates", "-1", status=1)

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)
    flags = ["--json", str(path/"content.json")]
    for suffix in ("csv", "tsv"):
        for name in ("levels", "spin_content"):
            flags += [f"--{suffix}-table", f"{name}={path/f'{name}.{suffix}'}"]
    run("6", "--levels", "all", "--spin-content", *flags, "--quiet", "--no-retain")
    exported = json.loads((path/"content.json").read_text())["tables"]
    for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
        for name in ("levels", "spin_content"):
            text = (path/f"{name}.{suffix}").read_text()
            actual = list(csv.DictReader((x for x in text.splitlines() if not x.startswith("#")), delimiter=delimiter))
            expected = [{k: "" if v is None else str(v).lower() if isinstance(v, bool) else str(v)
                         for k, v in r.items()} for r in records(exported[name])]
            assert actual == expected
    old = (path/"content.json").read_text()
    run("6", "--motif", "1,2", "--spin-content", "--json", str(path/"content.json"), "--force", status=1)
    assert (path/"content.json").read_text() == old
    run("4", "--motif", "2", "--spin-content", "--stream", "--no-retain", "--format", "plain")
    run("4", "--motif", "2", "--spin-content", "--format", "pretty")

print("Haldane-Shastry CLI contracts passed")
