"""CLI contracts: native-precision tokens, conventions, metadata and failures."""
import csv
from decimal import Decimal, localcontext
import io
import subprocess
import sys

program, fp128 = sys.argv[1:]


def run(*args, status=0):
    result = subprocess.run([program, *args], capture_output=True, text=True, timeout=60)
    assert result.returncode == status, (args, result.returncode, result.stdout, result.stderr)
    if status == 0:
        assert not result.stderr, result.stderr
    return result.stdout, result.stderr


def rows(output, delimiter=","):
    assert "# CPU time:" in output and "# Energy convention:" in output, output
    assert "References:" not in output and "\x1b" not in output, output
    text = "\n".join(line for line in output.splitlines() if not line.startswith("#"))
    return list(csv.DictReader(io.StringIO(text), delimiter=delimiter))


precisions = ["fp64", "long-double"]
if fp128.upper() in ("ON", "TRUE", "1"):
    precisions.append("fp128")
else:
    _, error = run("--u", "4", "--precision", "fp128", status=1)
    assert "unavailable" in error

with localcontext() as context:
    context.prec = 60
    for precision in precisions:
        base = ("--u", "4", "--precision", precision, "--format", "csv")
        out, _ = run(*base, "--points", "3")
        data = rows(out)
        assert len(data) == 9, data
        assert {r["branch"] for r in data} == {"spinon", "holon", "antiholon"}
        assert all(r["status"] == "converged" and r["energy"] for r in data)
        assert data[0]["energy"] == data[2]["energy"] == "0"
        reference = Decimal("1.24228148941956239984591553194945497655659639")
        tolerance = Decimal("1e-31" if precision == "fp128" else "1e-13")
        assert abs(Decimal(data[1]["energy"]) - reference) < tolerance, data[1]
        for r in data:
            assert r["spin"] == ("1/2" if r["branch"] == "spinon" else "0"), r
            assert int(r["delta_n"]) == {"spinon": 0, "holon": -1, "antiholon": 1}[r["branch"]]
        symmetric = rows(run(*base, "--momentum", "1")[0])
        unshifted_out, _ = run(*base, "--momentum", "1", "--convention", "unshifted")
        assert "# Energy convention: unshifted" in unshifted_out
        for s, r in zip(symmetric, rows(unshifted_out)):
            assert s["symmetric_energy"] == r["symmetric_energy"]
            assert abs(Decimal(r["energy"]) - Decimal(s["energy"]) - 2 * int(r["delta_n"])) < tolerance
        half_fermi = rows(run(*base, "--momentum", "1", "--convention", "unshifted", "--reference", "fermi")[0])
        assert [r["energy"] for r in half_fermi] == [r["energy"] for r in symmetric]

        doped_base = (*base, "--density", "0.5", "--points", "3")
        doped_out, _ = run(*doped_base, "--reference", "fermi")
        doped = rows(doped_out)
        assert len(doped) == 9 and {r["branch"] for r in doped} == {"spinon", "holon", "charge-particle"}
        assert "# Background status: converged" in doped_out
        assert "# Fermi rapidity Q:" in doped_out and "# Mu symmetric:" in doped_out
        metadata = dict(line[2:].split(": ", 1) for line in doped_out.splitlines() if line.startswith("# "))
        mu = Decimal(metadata["Mu symmetric"])
        assert abs(Decimal(metadata["Mu unshifted"]) - mu - 2) < tolerance
        for i, r in enumerate(doped):
            assert r["status"] == "converged" and r["energy"] == r["fermi_energy"], r
            assert "energy_mesh_error" in r and "energy_quad_error" not in r
            assert r["spin"] == ("1/2" if r["branch"] == "spinon" else "0")
            assert int(r["delta_n"]) == {"spinon": 0, "holon": -1, "charge-particle": 1}[r["branch"]]
            assert abs(Decimal(r["symmetric_energy"]) - Decimal(r["fermi_energy"]) - mu * int(r["delta_n"])) < tolerance
            if i % 3 != 1:
                assert r["fermi_energy"] == "0", r
            else:
                assert Decimal(r["fermi_energy"]) > 0, r
        unshifted = rows(run(*doped_base, "--convention", "unshifted")[0])
        invariant = rows(run(*doped_base, "--convention", "unshifted", "--reference", "fermi")[0])
        for s, r, f in zip(doped, unshifted, invariant):
            assert s["fermi_energy"] == r["fermi_energy"] == f["energy"]
            assert abs(Decimal(r["energy"]) - Decimal(s["symmetric_energy"]) - 2 * int(r["delta_n"])) < tolerance

    if "fp128" in precisions:
        token = "0.5000000000000000000000000000001"
        output, _ = run("--u", "4", "--density", token, "--precision", "fp128", "--branch", "spinon",
                        "--momentum", "0", "--format", "csv")
        density_line = next(line for line in output.splitlines() if line.startswith("# Density N/L:"))
        assert abs(Decimal(density_line.split(": ", 1)[1]) - Decimal(token)) < Decimal("1e-33")

out, _ = run("--u", "4", "--points", "2", "--format", "tsv")
assert len(rows(out, "\t")) == 6
for fmt in ("plain", "pretty"):
    out, _ = run("--u", "4", "--branch", "spinon", "--points", "3", "--format", fmt)
    assert "CPU time" in out and "symmetric" in out and "spinon" in out, out
    out, _ = run("--u", "4", "--density", "0.5", "--points", "2", "--format", fmt)
    assert "CPU time" in out and "Mu symmetric" in out and "charge-particle" in out

for budget in ("--max-evaluations", "--max-levels", "--max-iterations"):
    out, error = run("--u", "4", "--branch", "holon", "--momentum", "1", "--format", "csv", budget, "0", status=2)
    data = rows(out)
    assert len(data) == 1 and not data[0]["energy"] and not data[0]["symmetric_energy"], data
    assert data[0]["status"] == ("momentum_limit" if budget == "--max-iterations" else "quadrature_limit")
    assert "unavailable" in error

for options, expected in ((["--max-nodes", "16"], "mesh_limit"),
                          (["--max-background-iterations", "0"], "density_limit"),
                          (["--max-iterations", "0"], "momentum_limit")):
    out, error = run("--u", "4", "--density", "0.5", "--branch", "holon", "--momentum", "0.3",
                     "--format", "csv", *options, status=2)
    point = rows(out)[0]
    assert point["status"] == expected and "unavailable" in error, (out, error)
    assert not any(point[column] for column in ("energy", "symmetric_energy", "fermi_energy", "parameter"))

for args in (("--u", "0"), ("--u", "-1"), ("--u", "nan"), ("--branch", "unknown"),
             ("--convention", "unknown"), ("--points", "1"), ("--points", "3", "--momentum", "1"),
             ("--momentum", "-1", "--branch", "spinon"), ("--momentum", "4"),
             ("--max-levels", "25"), ("--tolerance", "0"), ("--format", "unknown")):
    run("--u", "4", *args, status=1)
for args in (("--density", "0"), ("--density", "1.1"), ("--density", "nan"),
             ("--reference", "unknown"), ("--branch", "charge-particle"),
             ("--density", "0.5", "--branch", "antiholon"),
             ("--density", "0.5", "--max-levels", "12"),
             ("--density", "0.5", "--max-nodes", "513"),
             ("--density", "0.5", "--initial-nodes", "3"),
             ("--density", "0.5", "--branch", "spinon", "--momentum", "2"),
             ("--max-nodes", "256")):
    run("--u", "4", *args, status=1)
print("Hubbard dispersion CLI contracts passed")
