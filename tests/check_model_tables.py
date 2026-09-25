"""Typed result schemas keep native scalars, auxiliary tables and failure semantics."""
import csv
from decimal import Decimal, localcontext
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

fp128, *programs = sys.argv[1:]
cases = {
    "bethe-xxx-pbc": (["6", "--roots"], {"states", "roots"}),
    "bethe-xxx-obc": (["6", "--roots"], {"states", "roots"}),
    "bethe-xxz-pbc": (["6", "--delta", "0.5", "--roots"], {"states", "roots"}),
    "bethe-xxz-obc": (["6", "--delta", "0.5", "--roots"], {"states", "roots", "boundary_roots"}),
    "bethe-lieb-liniger-pbc": (["4", "--length", "4", "--c", "1", "--roots"], {"states", "roots"}),
    "bethe-biquadratic-obc": (["4", "--roots"], {"states", "roots", "quantum_numbers"}),
    "bethe-hubbard-pbc": (["6", "--u", "4", "--roots"], {"states", "charge_roots", "spin_roots"}),
    "bethe-hubbard-obc": (["5", "--u", "4", "--roots"], {"states", "charge_roots", "spin_roots"}),
    "bethe-gaudin-yang-pbc": (["6", "--length", "6", "--c", "1", "--roots"], {"states", "charge_roots", "spin_roots"}),
    "bethe-tj-pbc": (["8", "--particles", "6", "--roots"], {"states", "first_roots", "second_roots"}),
    "bethe-sun-fermions-pbc": (["--populations", "1,0,3,1", "--length", "6", "--c", "1", "--roots"], {"states", "components", "roots"}),
    "bethe-ladder-pbc": (["6", "--rung", "0", "--roots"], {"states", "representations", "roots"}),
    "bethe-su3-pbc": (["6", "--roots"], {"states", "first_roots", "second_roots"}),
    "bethe-tb-pbc": (["4", "--roots"], {"states", "strings", "roots"}),
    "bethe-richardson": (["--levels", "0,1,2,3", "--pairs", "2", "--g", "1", "--variables"], {"states", "variables"}),
    "bethe-central-spin": (["--couplings", "1,0.7,0.3", "--field", "1", "--sz", "0", "--variables"], {"states", "variables"}),
}
env = dict(os.environ, UNI20_COLOR="never", COLUMNS="200")


def run(program, args, status=0, **kwargs):
    p = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=45, **kwargs)
    assert p.returncode == status, (program, args, p.returncode, p.stdout, p.stderr)
    return p


def records(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]


precisions = ["fp64", "long-double"]
if fp128.upper() in ("ON", "TRUE", "1"):
    precisions.append("fp128")
for executable in programs:
    program = str(Path(executable).resolve())
    name = Path(program).name
    args, expected = cases[name]
    for precision in precisions:
        common = [*args, "--precision", precision]
        doc = json.loads(run(program, [*common, "--format", "json"]).stdout)
        assert doc["status"] == "complete" and set(doc["tables"]) == expected
        main = doc["tables"]["states"]
        row = records(main)[0]
        assert row["converged"] is True
        if name == "bethe-hubbard-obc":
            assert "p" not in row and "momentum_index" not in row
        assert Decimal(row["energy"]) == Decimal(main["metadata"]["Total energy"])
        assert all(c["type"] != "string" for c in main["columns"] if c["id"] != "status")
        with localcontext() as context:
            context.prec = 70
            tolerance = Decimal("1e-29" if precision == "fp128" else "1e-12")
            if name in ("bethe-gaudin-yang-pbc", "bethe-sun-fermions-pbc"):
                if name == "bethe-gaudin-yang-pbc":
                    values = [r["k"] for r in records(doc["tables"]["charge_roots"])]
                else:
                    values = [r["rapidity"] for r in records(doc["tables"]["roots"]) if r["level"] == "0"]
                    components = records(doc["tables"]["components"])
                    assert [r["nesting_rank"] for r in components] == ["1", None, "0", "2"]
                assert abs(sum(Decimal(k)**2 for k in values) - Decimal(row["energy"])) < tolerance
            if name == "bethe-tj-pbc":
                roots = records(doc["tables"]["first_roots"])
                reconstructed = 4 - sum(1 / (Decimal(r["rapidity"])**2 + Decimal("0.25")) for r in roots)
                assert abs(reconstructed - Decimal(row["energy"])) < tolerance
            if name == "bethe-su3-pbc":
                roots = records(doc["tables"]["first_roots"])
                assert len(roots) == 4 and len(records(doc["tables"]["second_roots"])) == 2
                reconstructed = 6 - sum(1 / (Decimal(r["rapidity"]) ** 2 + Decimal("0.25")) for r in roots)
                assert abs(reconstructed - Decimal(row["energy"])) < tolerance
            if name == "bethe-tb-pbc":
                roots = records(doc["tables"]["roots"])
                assert len(roots) == 4 and len(records(doc["tables"]["strings"])) == 2
                reconstructed = Decimal(0)
                for root in roots:
                    x, y = Decimal(root["real"]), Decimal(root["imag"])
                    re, im = 1 + x*x - y*y, 2*x*y
                    reconstructed -= 4 * re / (re*re + im*im)
                assert abs(reconstructed - Decimal(row["energy"])) < tolerance
        for table_name, table in doc["tables"].items():
            assert table["summary"]["CPU time"].endswith(" s")
            assert table["summary"]["Outcome"] == "success"
            assert Decimal(table["summary"]["Run CPU seconds"]) >= 0
            assert Decimal(table["summary"]["Elapsed seconds"]) >= 0
            assert table["summary"] == main["summary"]  # One frozen run summary across all tables.
            assert table["metadata"]["Compiler"] and table["metadata"]["Platform"]
            assert table["metadata"]["Precision"] == precision
            assert all(r["state_id"] == "0" for r in records(table))
            for c in table["columns"]:
                if c["type"] == "real":
                    assert c["encoding"] == "decimal_string" and c["precision_bits"] >= 53
                    if precision == "fp128":
                        assert c["precision_bits"] == 113
        no_history = json.loads(run(program, [*common, "--format", "json", "--no-retain"]).stdout)
        assert all(no_history["tables"][key]["rows"] == doc["tables"][key]["rows"] for key in expected)
        # Final and streamed display must preserve each full-precision energy token.
        assert row["energy"] in run(program, [*common, "--format", "pretty"]).stdout
        assert row["energy"] in run(program, [*common, "--format", "plain", "--stream", "--no-retain"]).stdout
    with tempfile.TemporaryDirectory(prefix="bethe-model-tables-") as directory:
        folder = Path(directory)
        run(program, [*args, "--quiet", "--json", "all.json", "--csv", "states.csv"], cwd=folder)
        document = json.loads((folder / "all.json").read_text())
        text = "\n".join(l for l in (folder / "states.csv").read_text().splitlines() if not l.startswith("#"))
        saved = list(csv.DictReader(io.StringIO(text)))
        assert saved[0]["energy"] == records(document["tables"]["states"])[0]["energy"]
        run(program, [*args, "--table", "not_a_table", "--csv", "untouched.csv"], status=1, cwd=folder)
        assert not (folder / "untouched.csv").exists()
    if name in ("bethe-richardson", "bethe-central-spin"):
        document = json.loads(run(program, [*args, "--max-iterations", "0", "--format", "json"], status=2).stdout)
        row = records(document["tables"]["states"])[0]
        assert document["tables"]["states"]["summary"]["Outcome"] == "partial"
        assert row["converged"] is False
        if name == "bethe-richardson":
            assert Decimal(row["reached_g"]) == 0 and Decimal(row["requested_g"]) == 1
            assert Decimal(row["energy"]) == 2
        else:
            assert row["energy"] is None and row["reached_field"] is None
    if name in ("bethe-hubbard-pbc", "bethe-hubbard-obc"):
        mapped = json.loads(run(program, ["6", "--u", "-4", "--particles", "4", "--roots", "--format", "json"]).stdout)["tables"]
        assert records(mapped["states"])[0]["auxiliary_roots"] is True
        assert "auxiliary sector" in mapped["charge_roots"]["metadata"]["Roots and residuals"]
        assert len(records(mapped["charge_roots"])) == int(mapped["states"]["metadata"]["Root particles"])
        free = json.loads(run(program, ["6", "--u", "0", "--roots", "--format", "json"]).stdout)["tables"]
        assert set(free) == {"states", "free_modes"}
        assert all("quantum_number" not in r for r in records(free["free_modes"]))
    if name == "bethe-sun-fermions-pbc":
        missing = json.loads(run(program, [*args, "--max-iterations", "0", "--format", "json"], status=2).stdout)["tables"]
        assert records(missing["states"])[0]["energy"] is None and not records(missing["roots"])
        assert missing["states"]["summary"]["Outcome"] == "partial"
        free = json.loads(run(program, ["--populations", "2,0,1", "--length", "1", "--c", "0", "--roots", "--format", "json"]).stdout)["tables"]
        assert [r["component"] for r in records(free["free_modes"])] == ["0", "0", "2"]
        nested = json.loads(run(program, ["--populations", "1,1,1,1", "--length", "4", "--c", "1", "--roots", "--format", "json"]).stdout)["tables"]
        assert {r["level"] for r in records(nested["roots"])} == {"0", "1", "2", "3"}
    if name == "bethe-ladder-pbc":
        scan = json.loads(run(program, [*args, "--sectors", "--format", "json"]).stdout)["tables"]
        states = records(scan["states"])
        selected, = [r for r in states if r["selected"]]
        assert Decimal(selected["energy"]) == min(Decimal(r["energy"]) for r in states)
        assert all(r["state_id"] == selected["state_id"] for r in records(scan["roots"]))
        partial = json.loads(run(program, [*args, "--max-iterations", "0", "--format", "json"], status=2).stdout)["tables"]
        assert records(partial["states"])[0]["converged"] is False
        assert partial["states"]["summary"]["Outcome"] == "partial"
        assert records(partial["states"])[0]["energy"] is not None
        assert {r["level"] for r in records(partial["roots"])} == {"1"}
        nested = json.loads(run(program, ["8", "--rung", "0", "--roots", "--format", "json"]).stdout)["tables"]
        assert {r["level"] for r in records(nested["roots"])} == {"1", "2", "3"}
        for field, magnetization in (("0.125", "2"), ("-0.125", "-2")):
            field_tables = json.loads(run(program, ["6", "--rung", "0", "--singlets", "4",
                                                    "--field", field, "--roots", "--format", "json"]).stdout)["tables"]
            r = records(field_tables["states"])[0]
            assert r["magnetization"] == magnetization
            assert abs(Decimal(r["energy"]) - Decimal("-2.986067977499789696409173668731")) < Decimal("1e-12")
            pop = records(field_tables["representations"])
            assert int(pop[1]["population"]) - int(pop[3]["population"]) == int(magnetization)
            assert Decimal(field_tables["states"]["metadata"]["Magnetic field h"]) == Decimal(field)
        missing = json.loads(run(program, [*args, "--field", "0.5", "--max-branches", "0", "--format", "json"], status=2).stdout)["tables"]
        assert records(missing["states"])[0]["magnetization"] is None
        assert not records(missing["representations"])
        for projection in ("-2", "0", "1", "6"):
            fixed = json.loads(run(program, [*args, "--sz", projection, "--sectors", "--field", "0.5", "--format", "json"]).stdout)["tables"]
            states = records(fixed["states"])
            assert len(states) == 7 - abs(int(projection))
            assert all(r["magnetization"] == projection for r in states)
            assert Decimal(fixed["states"]["metadata"]["Requested total Sz"]) == Decimal(projection)
            selected, = [r for r in states if r["selected"]]
            assert all(r["state_id"] == selected["state_id"] for r in records(fixed["roots"]))
            pops = records(fixed["representations"])
            for r in states:
                counts = [int(p["population"]) for p in pops if p["state_id"] == r["state_id"]]
                assert counts[1] - counts[3] == int(projection)
        constrained_missing = json.loads(run(program, [*args, "--sz", "1", "--max-branches", "0", "--format", "json"], status=2).stdout)["tables"]
        assert records(constrained_missing["states"])[0]["magnetization"] is None
        assert Decimal(constrained_missing["states"]["metadata"]["Requested total Sz"]) == 1
print(f"Native model table contracts passed for {len(programs)} frontends")
