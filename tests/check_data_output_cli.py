"""Shared output contracts, exercised through the Hubbard dispersion frontend."""
import csv
from decimal import Decimal
import io
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile

program, fp128 = sys.argv[1:]
precisions = ["fp64", "long-double"]
if fp128.upper() in ("ON", "TRUE", "1"):
    precisions.append("fp128")


def run(*args, status=0, **kwargs):
    # Tests of invalid U must reach scientific validation, not duplicate-option rejection.
    interaction = [] if args[:1] == ("--u",) else ["--u", "4"]
    result = subprocess.run([program, *interaction, *map(str, args)], text=True,
                            stdout=kwargs.pop("stdout", subprocess.PIPE),
                            stderr=subprocess.PIPE, timeout=60, **kwargs)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    if status == 0:
        assert not result.stderr, result.stderr
    return result


def separated(text, delimiter=","):
    comments = dict(line[2:].split(": ", 1) for line in text.splitlines() if line.startswith("# "))
    records = csv.DictReader(io.StringIO("\n".join(line for line in text.splitlines()
                                                if not line.startswith("#"))), delimiter=delimiter)
    return comments, list(records)


def document(text):
    return json.loads(text, parse_float=Decimal)


def compare(doc, records):
    keys = [c["id"] for c in doc["columns"]]
    assert len(doc["rows"]) == len(records)
    for row, record in zip(doc["rows"], records):
        assert dict(zip(keys, ("" if v is None else str(v) for v in row))) == record


with tempfile.TemporaryDirectory(prefix="bethe-output-") as directory:
    root = Path(directory)
    for precision in precisions:
        # Apostrophes and whitespace in paths must survive provenance quoting.
        csv_path = root / f"{precision} one's table.csv"
        tsv_path = root / f"{precision}.tsv"
        json_path = root / f"{precision}.json"
        result = run("--precision", precision, "--points", 3, "--format", "plain",
                     "--csv", csv_path, "--tsv", tsv_path, "--json", json_path)
        assert "CPU time" in result.stdout and "spinon" in result.stdout
        metadata, records = separated(csv_path.read_text())
        tsv_metadata, tsv_records = separated(tsv_path.read_text(), "\t")
        doc = document(json_path.read_text())
        assert records == tsv_records and metadata == tsv_metadata
        compare(doc, records)
        for key, value in doc["metadata"].items():
            # The command has shell apostrophe escapes; CSV comments escape backslashes too.
            if key != "Command":
                assert metadata[key] == value, (key, metadata[key], value)
        for key, value in doc["summary"].items():
            assert metadata[key] == value
        assert shlex.split(doc["metadata"]["Command"]) == result.args
        assert re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ", metadata["Date"])
        assert re.fullmatch(r"[0-9a-f]{40}(-dirty)?|unavailable", metadata["Bethe revision"])
        assert re.fullmatch(r"[0-9a-f]{40}(-dirty)?|unavailable|unknown", metadata["Uni20 revision"])
        assert metadata["Program"] == "bethe-hubbard-dispersion"
        assert metadata["Precision"] == precision and metadata["Bethe version"]
        assert metadata["Status"] == "converged" and metadata["Rows"] == "9"
        assert metadata["Outcome"] == "success"
        assert Decimal(metadata["Run CPU seconds"]) >= 0
        assert Decimal(metadata["Elapsed seconds"]) >= 0
        assert metadata["Compiler"] and metadata["Platform"]
        assert re.fullmatch(r"\d+\.\d{6} s", metadata["CPU time"])
        assert csv_path.read_text().index("# CPU time:") > csv_path.read_text().index("branch,p,")
        columns = {c["id"]: c for c in doc["columns"]}
        assert columns["spin"]["encoding"] == "number" and columns["spin"]["type"] == "half_int"
        assert columns["p"]["encoding"] == "decimal_string" and columns["p"]["unit"] == "radians"
        assert columns["energy"]["nullable"] and columns["energy"]["unit"] == "t"
        assert doc["rows"][0][6] == Decimal("0.5") and doc["rows"][0][7] == "inf"
        assert doc["rows"][2][7] == "-inf"
        if precision == "fp128":
            assert columns["energy"]["precision_bits"] == 113
            assert len(records[1]["energy"]) >= 33
        elif precision == "fp64":
            assert columns["energy"]["precision_bits"] == 53
        # Retention changes storage, not values; machine stdout remains valid JSON.
        live = document(run("--precision", precision, "--points", 3, "--no-retain",
                            "--format", "json").stdout)
        assert live["rows"] == doc["rows"] and live["columns"] == doc["columns"]
        assert live["summary"]["Rows"] == "9"
        # Narrow live screens must not wrap numerical tokens into different values.
        screen = run("--precision", precision, "--points", 3, "--branch", "spinon",
                     "--stream", "--no-retain", "--format", "plain",
                     env={**os.environ, "COLUMNS": "40"}).stdout
        assert records[1]["energy"] in screen and "CPU time" in screen

    path = root / "quiet.csv"
    assert not run("--points", 2, "--quiet", "--no-retain", "--csv", path).stdout
    assert len(separated(path.read_text())[1]) == 6
    strict = root / "strict.tsv"
    text = run("--points", 2, "--format", "csv", "--no-preamble", "--tsv", strict).stdout
    assert not text.startswith("#") and "# CPU time:" not in text
    assert len(list(csv.DictReader(io.StringIO(text)))) == 6
    assert len(list(csv.DictReader(io.StringIO(strict.read_text()), delimiter="\t"))) == 6
    second_csv = root / "second.csv"
    third_csv = root / "third.csv"
    run("--points", 2, "--quiet", "--csv", second_csv, "--csv", third_csv)
    assert second_csv.read_bytes() == third_csv.read_bytes()
    assert "metadata" in document(run("--branch", "spinon", "--points", 2,
                                     "--format", "json", "--no-preamble").stdout)

    # Unconverged physics is exit 2, with matching missing cells in every sink.
    for density in ("1", "0.5"):
        path = root / f"failed-{density}.json"
        result = run("--density", density, "--branch", "holon", "--momentum", ".3",
                     "--max-iterations", 0, "--format", "csv", "--json", path, status=2)
        meta, records = separated(result.stdout)
        doc = document(path.read_text())
        compare(doc, records)
        assert doc["rows"][0][3] is None and doc["rows"][0][7] is None
        assert meta["Status"].startswith("incomplete") and "unavailable" in result.stderr
        assert doc["summary"]["Outcome"] == "partial"

    # Existing files are protected; even --force must validate all arguments first.
    existing = root / "existing.csv"
    existing.write_text("sentinel")
    assert "exists" in run("--points", 2, "--csv", existing, status=1).stderr
    for args in (("--u", 0), ("--tolerance", 0), ("--momentum", 4), ("--no-retain",),
                 ("--density", ".5", "--max-nodes", 513),
                 ("--density", ".5", "--momentum", 0)):
        new_path = root / "must-not-exist.json"
        run(*args, "--csv", existing, "--json", new_path, "--force", status=1)
        assert existing.read_text() == "sentinel" and not new_path.exists()
    run("--points", 2, "--quiet", "--csv", existing, "--force")
    assert len(separated(existing.read_text())[1]) == 6

    # Duplicate destinations, including aliases, must not truncate either file.
    existing.write_text("sentinel")
    hardlink = root / "hardlink"
    symlink = root / "symlink"
    os.link(existing, hardlink)
    symlink.symlink_to(existing)
    for alias in (existing, root / "." / existing.name, hardlink, symlink):
        run("--points", 2, "--csv", existing, "--json", alias, "--force", status=1)
        assert existing.read_text() == "sentinel"
    dangling = root / "dangling"
    dangling_target = root / "dangling-target"
    dangling.symlink_to(dangling_target)
    run("--points", 2, "--csv", dangling, "--json", dangling_target, "--force", status=1)
    assert not dangling_target.exists()
    # An output already redirected by the shell must not be overwritten internally.
    with existing.open("a") as screen:
        result = run("--points", 2, "--csv", existing, "--force", stdout=screen, status=1)
    assert "aliases stdout" in result.stderr and existing.read_text() == "sentinel"
    run("--points", 2, "--csv", root / "absent" / "out.csv", status=1)
    run("--points", 2, "--csv", root, "--force", status=1)
    run("--csv", "-", status=1)
    run("--json", status=1)
    run("--format", "bad", status=1)

    # Metadata controls stay on comment lines, even with a newline in argv.
    odd_path = root / "newline\nfile.csv"
    result = run("--points", 2, "--csv", odd_path, "--format", "json")
    assert shlex.split(document(result.stdout)["metadata"]["Command"]) == result.args
    assert len(separated(odd_path.read_text())[1]) == 6

    if Path("/dev/full").exists():
        for fmt, extra in (("csv", []), ("json", []), ("plain", []), ("pretty", []),
                           ("auto", ["--stream"]), ("plain", ["--stream"])):
            healthy = root / "healthy.json"
            with open("/dev/full", "w") as full:
                result = run("--branch", "spinon", "--points", 3, "--format", fmt,
                             "--json", healthy, "--force", *extra, stdout=full, status=1)
            assert "stdout" in result.stderr
            doc = document(healthy.read_text())
            # Failure may occur at begin, row or final flush, depending on buffering.
            assert len(doc["rows"]) <= 3
            assert int(doc["summary"]["Rows"]) == len(doc["rows"])
            assert doc["summary"]["Status"] in ("aborted", "converged")
            assert doc["summary"]["Outcome"] == (
                "failed" if doc["summary"]["Status"] == "aborted" else "success")

print("Data output CLI contracts passed")
