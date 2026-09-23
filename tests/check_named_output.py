"""Named table exports and shared CLI lifecycle; no scientific values narrowed."""
import csv
from decimal import Decimal
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program = str(Path(sys.argv[1]).resolve())
base = ["2", "--length", "4", "--lambda", "2", "--labels", "0,1", "--pseudomomenta"]
env = dict(os.environ, UNI20_COLOR="never", COLUMNS="4096")


def run(*args, status=0, **kwargs):
    p = subprocess.run([program, *args], text=True, capture_output=True, env=env, timeout=30, **kwargs)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    if status == 0:
        assert not p.stderr, p.stderr
    return p


def rows(text, delimiter=","):
    return list(csv.DictReader(io.StringIO("\n".join(x for x in text.splitlines() if not x.startswith("#"))),
                               delimiter=delimiter))


def records(table):
    return [dict(zip([c["id"] for c in table["columns"]], r)) for r in table["rows"]]


document = json.loads(run(*base, "--format", "json").stdout)
assert set(document["tables"]) == {"states", "pseudomomenta"} and document["status"] == "complete"
states, roots = (document["tables"][name] for name in ("states", "pseudomomenta"))
assert len(states["rows"]) == 1 and len(roots["rows"]) == 2
assert records(states)[0]["state_id"] == "0"
assert all(r["state_id"] == "0" for r in records(roots))
assert all(c["type"] != "string" for c in roots["columns"])
assert abs(sum(Decimal(r["k"]) ** 2 for r in records(roots)) - Decimal(records(states)[0]["energy"])) < Decimal("1e-13")
for name in document["tables"]:
    other = json.loads(run(*base, "--format", "json", "--no-retain").stdout)["tables"][name]
    assert other["rows"] == document["tables"][name]["rows"]

with tempfile.TemporaryDirectory(prefix="bethe-tables-") as directory:
    folder = Path(directory)
    p = run(*base, "--quiet", "--no-retain", "--csv", "states.csv", "--tsv-table", "pseudomomenta=roots.tsv",
            "--json", "all.json", cwd=folder)
    assert not p.stdout
    assert rows((folder / "states.csv").read_text()) == records(states)
    assert rows((folder / "roots.tsv").read_text(), "\t") == records(roots)
    assert json.loads((folder / "all.json").read_text())["tables"]["states"]["rows"] == states["rows"]
    # Reject unavailable datasets and destination aliases before creating any file.
    for extra in [("--table", "missing"), ("--csv-table", "missing=other.csv"),
                  ("--tsv", "fresh.csv"), ("--csv-table", "pseudomomenta=fresh.csv")]:
        run(*base, "--csv", "fresh.csv", *extra, status=1, cwd=folder)
        assert not (folder / "fresh.csv").exists()
    sentinel = folder / "protected.json"
    sentinel.write_text("keep me")
    for action in ["--help", "--references", "--version", "--build-info"]:
        p = run("--force", "--json", "protected.json", action, cwd=folder)
        assert sentinel.read_text() == "keep me"
        assert p.stdout and not p.stderr
    # Flag-looking file names remain values, not informational requests.
    run(*base, "--quiet", "--csv=--references", cwd=folder)
    assert rows((folder / "--references").read_text()) == records(states)
    # Literal '=' in a path is retained after the first TABLE= separator.
    run(*base, "--quiet", "--csv-table", "pseudomomenta=a=b.csv", cwd=folder)
    assert len(rows((folder / "a=b.csv").read_text())) == 2

help_text = run("--help").stdout
assert "--references" in help_text and "Used for:" not in help_text
assert run(status=1).stderr == help_text
assert "Used for:" in run("--references").stdout
for args in [("--length", "5"), ("--window", "-1"), ("--window", "1x"),
             ("--max-states", "18446744073709551616")]:
    run(*base, *args, status=1)
assert "#" not in run(*base, "--format", "csv", "--no-preamble").stdout
run(*base, "--no-retain", status=1)
run(*base, "--no-retain", "--stream", "--format", "plain")
print("Named table and shared CLI contracts passed")
