"""Shared physical-spin table across all biquadratic solver/output paths."""
import csv
from decimal import Decimal
import json
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
    return [dict(zip([c["id"] for c in table["columns"]], r)) for r in table["rows"]]


def expected(ell):
    # Independent character polynomial in physical Sz, then highest-weight
    # subtraction c(S)=weight(S)-weight(S+1), rather than production CG fusion.
    previous, current = {}, {0: 1}
    for _ in range(ell):
        nxt = {}
        for m, count in current.items():
            for shift in (-1, 0, 1):
                nxt[m+shift] = nxt.get(m+shift, 0) + count
        for m, count in previous.items():
            nxt[m] -= count
        previous, current = current, nxt
    return {s: current.get(s, 0)-current.get(s+1, 0) for s in range(ell+1)
            if current.get(s, 0) != current.get(s+1, 0)}


def check(args, status=0):
    document = json.loads(run([*args, "--spin-content", "--format", "json"], status))
    tables = document["tables"]
    content = rows(tables["spin_content"])
    assert "independent of solver convergence" in tables["spin_content"]["metadata"]["Spin content"]
    assert tables["spin_content"]["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    states = sum((rows(tables[name]) for name in ("states", "reference", "failed") if name in tables), [])
    assert {r["state_id"] for r in content} == {r["state_id"] for r in states}
    for state in states:
        group = [r for r in content if r["state_id"] == state["state_id"]]
        ell = int(state.get("through_lines", args[0] if "--ferromagnetic" in args else 0))
        assert all(int(r["through_lines"]) == ell for r in group)
        if ell >= 46:
            assert len(group) == 1
            assert all(group[0][k] is None for k in ("spin", "multiplets", "magnetic_states"))
            assert group[0]["status"] == "total dimension exceeds uint64"
            continue
        counts = expected(ell)
        assert len(group) == len(counts)
        assert [Decimal(str(r["spin"])) for r in group] == list(counts)
        for r in group:
            s = int(Decimal(str(r["spin"])))
            assert int(r["multiplets"]) == counts[s]
            assert int(r["magnetic_states"]) == (2*s+1)*counts[s]
            assert r["status"] == "exact"
        if "multiplicity" in state:
            assert sum(int(r["magnetic_states"]) for r in group) == int(state["multiplicity"])
    streamed = json.loads(run([*args, "--spin-content", "--format", "json", "--no-retain"], status))
    assert all(streamed["tables"][name]["rows"] == table["rows"] for name, table in tables.items())
    return tables


cases = [
    ["6"], ["6", "--sectors"], ["6", "--excitations", "all", "--through-lines", "2"],
    ["6", "--singlet-excitation"], ["4", "--q-spectrum"], ["4", "--q-seed", "none"],
    ["8", "--ferromagnetic"], ["8", "--ferromagnetic", "--one-defect"],
    ["8", "--ferromagnetic", "--excitations", "all", "--through-lines", "4"],
    ["4", "--ferromagnetic", "--q-spectrum"],
    ["8", "--ferromagnetic", "--bound-pairs", "all"],
    ["8", "--ferromagnetic", "--bound-triples", "all"],
    ["10", "--ferromagnetic", "--bound-quartets", "all"],
    ["8", "--ferromagnetic", "--pair-defect", "3,2"],
    ["8", "--ferromagnetic", "--two-pairs", "1,2"],
    ["8", "--ferromagnetic", "--triple-defect", "2,1"],
]
for precision in precisions:
    for args in cases:
        check([*args, "--precision", precision])
    for args in (cases[0], cases[2], cases[3], cases[8], cases[12]):
        check([*args, "--precision", precision, "--max-iterations", "0"], 2)
    check(["4", "--q-seed", "1,2", "--precision", precision, "--max-iterations", "0"], 2)
    # Exact dimension boundary, and large-N ground/cluster overflow placeholders.
    for n in ("45", "46", "100000"):
        check([n, "--ferromagnetic", "--precision", precision])
    check(["100000", "--ferromagnetic", "--bound-quartets", "1", "--precision", precision])

assert "spin_content" not in json.loads(run(["6", "--format", "json"]))["tables"]
assert "--spin-content" in run(["--help"])
assert "Physical SU(2) content" in run(["6", "--spin-content"])
with tempfile.TemporaryDirectory() as tmp:
    base = Path(tmp)
    args = ["8", "--ferromagnetic", "--bound-pairs", "2", "--spin-content"]
    run([*args, "--no-retain", "--quiet", "--json", str(base / "spins.json"),
         "--csv-table", f"spin_content={base / 'spins.csv'}",
         "--tsv-table", f"spin_content={base / 'spins.tsv'}"])
    original = rows(json.loads((base / "spins.json").read_text())["tables"]["spin_content"])
    for suffix, delimiter in (("csv", ","), ("tsv", "\t")):
        text = (base / f"spins.{suffix}").read_text()
        assert "#" in text and "Spin content" in text
        data = list(csv.DictReader((line for line in text.splitlines() if not line.startswith("#")), delimiter=delimiter))
        assert len(data) == len(original)
        for a, b in zip(data, original):
            for key in ("state_id", "through_lines", "spin", "multiplets", "magnetic_states"):
                assert Decimal(a[key]) == Decimal(str(b[key]))
    run(["6", "--csv-table", f"spin_content={base / 'disabled.csv'}"], 1)
