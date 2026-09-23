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
        assert Decimal(row["energy"]) == Decimal(main["metadata"]["Total energy"])
        assert all(c["type"] != "string" for c in main["columns"] if c["id"] != "status")
        with localcontext() as context:
            context.prec = 70
            tolerance = Decimal("1e-29" if precision == "fp128" else "1e-12")
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
        assert row["converged"] is False
        if name == "bethe-richardson":
            assert Decimal(row["reached_g"]) == 0 and Decimal(row["requested_g"]) == 1
            assert Decimal(row["energy"]) == 2
        else:
            assert row["energy"] is None and row["reached_field"] is None
print(f"Native model table contracts passed for {len(programs)} frontends")
