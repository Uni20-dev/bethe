"""Frontend-owned table declarations, strict selection, and independent destinations."""
import ast
import csv
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

build = Path(sys.argv[1]).resolve()
source = Path(__file__).resolve().parent.parent
env = dict(os.environ, UNI20_COLOR="never", COLUMNS="4096", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")

# Reuse the ordinary CLI smoke inputs, not a second table catalogue.
tree = ast.parse((source / "tests/check_shared_arguments.py").read_text())
base = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
            and any(isinstance(t, ast.Name) and t.id == "base_arguments" for t in n.targets))
base.update({"bethe-su3-dispersion": ["--points", "2"], "bethe-tb-dispersion": ["--points", "2"],
             "bethe-haldane-shastry-pbc": ["4"], "bethe-tb-pbc": ["4"],
             "bethe-tj-pbc": ["4"], "bethe-biquadratic-obc": ["4"]})


def run(name, args, status=0):
    p = subprocess.run([str(build / name), *args], capture_output=True, text=True, env=env, timeout=60)
    assert p.returncode == status, (name, args, p.returncode, p.stdout, p.stderr)
    return p


def document(name, args):
    return json.loads(run(name, [*args, "--format", "json"]).stdout)["tables"]


def csv_rows(path):
    return list(csv.reader(io.StringIO("\n".join(line for line in path.read_text().splitlines()
                                                if not line.startswith("#")))))


with tempfile.TemporaryDirectory(prefix="bethe-table-contracts-") as directory:
    folder = Path(directory)
    protected = folder / "protected.csv"
    other = folder / "other.csv"
    count = 0
    for file in sorted((source / "apps").glob("bethe-*.cpp")):
        name = file.stem
        help_text = run(name, ["--help"]).stdout
        if name == "bethe-hubbard-dispersion":
            assert "--csv-table" not in help_text  # Explicit single-document frontend.
            continue
        count += 1
        assert "Output tables" in help_text
        section = help_text.split("Output tables", 1)[1].split("--table", 1)[0]
        declared = set(re.findall(r"^\s*([a-z][a-z0-9_]*)(?: \[primary when available\])?: ", section, re.M))
        assert declared and "[primary when available]" in section, (name, section)
        args = base.get(name, ["4"])
        tables = document(name, args)
        assert set(tables) <= declared, (name, set(tables), declared)
        # Names rejected during parsing, even with an unrelated export and --force.
        for selector in [["--table", "removed_table"], ["--csv-table", f"removed_table={other}"],
                         ["--tsv-table", f"removed_table={other}"], ["--csv-table", f"{next(iter(declared))}="]]:
            protected.write_text("keep me")
            p = run(name, [*args, "--quiet", "--force", "--csv", str(protected), *selector], 1)
            assert protected.read_text() == "keep me" and not other.exists(), (name, selector)
            assert p.stderr

    # Each output-only family must be producible without its screen flag.
    auxiliary = [
        ("bethe-xxx-pbc", ["4"], "--roots"),
        ("bethe-xxx-obc", ["4"], "--roots"),
        ("bethe-xxz-pbc", ["4", "--delta", "0.5"], "--roots"),
        ("bethe-xxz-obc", ["4", "--delta", "0.5"], "--roots"),
        ("bethe-hubbard-pbc", ["4", "--u", "4"], "--roots"),
        ("bethe-hubbard-obc", ["4", "--u", "0"], "--roots"),
        ("bethe-su3-pbc", ["3"], "--roots"),
        ("bethe-tb-pbc", ["4"], "--roots"),
        ("bethe-gaudin-yang-pbc", ["2", "--length", "2", "--c", "1"], "--roots"),
        ("bethe-gaudin-yang-pbc", ["2", "--length", "2", "--c", "0"], "--roots"),
        ("bethe-bose-fermi-pbc", ["--bosons", "1", "--fermions", "1", "--length", "2", "--c", "1"], "--roots"),
        ("bethe-bose-fermi-pbc", ["--bosons", "1", "--fermions", "1", "--length", "2", "--c", "0"], "--roots"),
        ("bethe-sun-fermions-pbc", ["--populations", "1,1", "--length", "2", "--c", "1"], "--roots"),
        ("bethe-tj-pbc", ["4"], "--roots"),
        ("bethe-tj-pbc", ["4", "--particles", "2", "--sz", "1"], "--roots"),
        ("bethe-ladder-pbc", ["2", "--rung", "1"], "--roots"),
        ("bethe-lieb-liniger-pbc", ["2", "--length", "2", "--c", "1"], "--roots"),
        ("bethe-lieb-liniger-obc", ["2", "--length", "2", "--c", "1"], "--roots"),
        ("bethe-q-boson-pbc", ["4", "--particles", "2", "--eta", "1"], "--roots"),
        ("bethe-xyz-pbc", ["4", "--eta", "0.4", "--t", "0.7"], "--roots"),
        ("bethe-tasep-pbc", ["4", "--particles", "2"], "--roots"),
        ("bethe-asep-pbc", ["4", "--particles", "2", "--left-rate", "0.5"], "--roots"),
        ("bethe-xxz-qg-obc", ["4", "--delta", "0.5"], "--roots"),
        ("bethe-richardson", ["--levels", "0,1", "--pairs", "1", "--g", "1"], "--variables"),
        ("bethe-central-spin", ["--couplings", "1,0.7,0.3", "--field", "1", "--sz", "0"], "--variables"),
        ("bethe-sutherland-pbc", ["2", "--length", "2", "--lambda", "2"], "--pseudomomenta"),
        ("bethe-haldane-shastry-pbc", ["4"], "--spin-content"),
        ("bethe-biquadratic-obc", ["4"], "--roots"),
        ("bethe-biquadratic-obc", ["4"], "--spin-content"),
    ]
    checked = 0
    for name, args, flag in auxiliary:
        ordinary = document(name, args)
        extended = document(name, [*args, flag])
        extra = set(extended) - set(ordinary)
        assert extra, (name, args, flag)
        for table in extra:
            checked += 1
            path = folder / "selected.csv"
            json_path = folder / "ordinary.json"
            p = run(name, [*args, "--format", "json", "--force", "--csv-table", f"{table}={path}",
                           "--json", str(json_path)])
            assert set(json.loads(p.stdout)["tables"]) == set(ordinary), (name, table, p.stdout)
            assert set(json.loads(json_path.read_text())["tables"]) == set(ordinary)
            rows = csv_rows(path)
            assert rows[0] == [c["id"] for c in extended[table]["columns"]], (name, table, rows)
            assert len(rows) - 1 == len(extended[table]["rows"]), (name, table)
            # Human and streaming output must also keep file-only auxiliaries off screen.
            for live in [[], ["--stream", "--no-retain"]]:
                p = run(name, [*args, "--format", "plain", *live, "--force", "--csv-table", f"{table}={path}"])
                assert extended[table]["title"] not in p.stdout, (name, table, live, p.stdout)
            # The same table can independently go to stdout and a file.
            p = document(name, [*args, flag, "--force", "--csv-table", f"{table}={path}"])
            assert set(p) == set(extended), (name, table)

    # Available schema with no rows is valid (polarized state), not a missing table.
    path = folder / "empty.csv"
    run("bethe-xxx-pbc", ["4", "--sz", "2", "--quiet", "--csv-table", f"roots={path}"])
    assert len(csv_rows(path)) == 1

    # Outcome- and mode-dependent tables must never be synthesized to satisfy an export.
    for name, args, table in [
        ("bethe-xxx-pbc", ["4"], "spinons"),
        ("bethe-xxx-pbc", ["4", "--excitations", "1", "--spin", "1"], "failed"),
        ("bethe-hubbard-pbc", ["4", "--u", "0"], "spin_roots"),
        ("bethe-xxz-qg-obc", ["4", "--delta", "0"], "roots"),
        ("bethe-biquadratic-obc", ["4"], "string"),
    ]:
        protected.write_text("keep me")
        run(name, [*args, "--quiet", "--force", "--csv", str(protected), "--csv-table", f"{table}={other}"], 1)
        assert protected.read_text() == "keep me" and not other.exists(), (name, args, table)

print(f"Table contracts passed: {count} named frontends, {checked} auxiliary selections.")
