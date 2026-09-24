"""Shared two-/three-/four-string CLI, native precision and output contracts."""
from decimal import Decimal, localcontext
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, fp128, defects = sys.argv[1:]
defects = int(defects)
options = {2: "--bound-pairs", 3: "--bound-triples", 4: "--bound-quartets"}
assert defects in options
option = options[defects]
family = {2: "two", 3: "three", 4: "four"}[defects]
minimum_n = 2 * defects
scan_n = max(8, minimum_n + 2)
precisions = ["fp64", "long-double"] + (["fp128"] if fp128.upper() in ("ON", "TRUE", "1") else [])
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")


def run(args, status=0):
    p = subprocess.run([program, *args], text=True, capture_output=True, timeout=30, env=env)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def tables(args, status=0):
    doc = json.loads(run([*args, "--format", "json"], status))
    assert doc["status"] == "complete"  # Transport completion, not scientific outcome.
    for table in doc["tables"].values():
        assert table["summary"]["Outcome"] == ("success" if status == 0 else "partial")
        assert f"NOT the full {family}-defect" in table["metadata"]["Coverage"]
        assert "not a global excitation rank" in table["metadata"]["Ordering"]
        assert "exact degenerate ferro" in table["metadata"]["Gap reference"]
        assert Decimal(table["summary"]["Run CPU seconds"]) >= 0
    return doc["tables"]


def rows(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]


def tl_vacuum_characteristic(n):
    """Exact independent link-pattern matrix for sum e_i at loop weight 3.

    Faddeev-LeVerrier uses integer arithmetic; no numerical diagonalization,
    Bethe roots, or optional symbolic package is involved.
    """
    def patterns(sites):
        if not sites:
            yield ()
            return
        for k in range(1, len(sites), 2):
            for a in patterns(sites[1:k]):
                for b in patterns(sites[k+1:]):
                    yield tuple(sorted(((sites[0], sites[k]),) + a + b))
    basis = list(patterns(tuple(range(n))))
    index = {pattern: i for i, pattern in enumerate(basis)}
    dim = len(basis)
    matrix = [[0]*dim for _ in basis]
    for col, pattern in enumerate(basis):
        partner = {x: y for a, b in pattern for x, y in ((a, b), (b, a))}
        for i in range(n-1):
            if partner[i] == i+1:
                matrix[col][col] += 3
            else:
                pairs = [p for p in pattern if i not in p and i+1 not in p]
                pairs += [(i, i+1), tuple(sorted((partner[i], partner[i+1])))]
                matrix[index[tuple(sorted(pairs))]][col] += 1
    b = [[int(i == j) for j in range(dim)] for i in range(dim)]
    coefficients = [1]
    for k in range(1, dim+1):
        b = [[sum(matrix[i][l]*b[l][j] for l in range(dim)) for j in range(dim)] for i in range(dim)]
        trace = sum(b[i][i] for i in range(dim))
        assert trace % k == 0
        c = -trace // k
        coefficients.append(c)
        for i in range(dim):
            b[i][i] += c
    assert all(value == 0 for row in b for value in row)
    return coefficients


quartet_polynomial = tl_vacuum_characteristic(8) if defects == 4 else None


for precision in precisions:
    base = [str(minimum_n), "--ferromagnetic", option, "all", "--precision", precision]
    t = tables([*base, "--roots"])
    s = rows(t["states"])[0]
    with localcontext() as ctx:
        ctx.prec = 70
        tol = Decimal("1e-29" if precision == "fp128" else "1e-12")
        if defects == 2:
            exact = (15 - Decimal(17).sqrt()) / 2
        elif defects == 3:
            gap = Decimal("2.3")
            for _ in range(15):
                gap -= (((gap-17)*gap+80)*gap-106) / ((3*gap-34)*gap+80)
            exact = 5 + gap
        else:
            gap = Decimal("2.25")
            for _ in range(15):
                value, derivative = Decimal(1), Decimal(0)
                for c in quartet_polynomial[1:]:
                    derivative = derivative*gap + value
                    value = value*gap + c
                gap -= value/derivative
            exact = 7 + gap
        assert abs(Decimal(s["energy"]) - exact) < tol
        assert abs(Decimal(s["gap"]) - exact + minimum_n - 1) < tol
        assert abs(Decimal(s["tl_energy"]) - Decimal(s["gap"])) < tol
        assert abs(Decimal(s["energy"]) + 2 * Decimal(s["reference_energy"]) - Decimal(7*(minimum_n-1)) / 4) < tol
    assert s["through_lines"] == "0" and s["multiplicity"] == "1" and s["mode"] == "1"
    assert s["converged"] and len(rows(t["roots"])) == defects
    assert rows(t["reference"])[0]["multiplicity"] == {2: "55", 3: "377", 4: "2584"}[defects]
    for n in (minimum_n + 1, scan_n, 16, 65):
        count = n - (minimum_n - 1)
        args = [str(n), "--ferromagnetic", option, "all", "--precision", precision, "--roots"]
        full = tables(args)
        states = rows(full["states"])
        assert len(states) == count
        assert [int(s["mode"]) for s in states] == list(range(1, count + 1))
        assert all(s["converged"] and int(s["through_lines"]) == n - minimum_n for s in states)
        energies = [Decimal(s["energy"]) for s in states]
        assert energies == sorted(energies) and len(set(energies)) == len(energies)
        for i, string in enumerate(rows(full["string"])):
            assert int(string["state_id"]) == i and int(string["string_label"]) == count - i
            if defects == 2:
                assert int(string["deviation_sign"]) == (1 if i % 2 == 0 else -1)
                assert Decimal(string["deviation"]) * int(string["deviation_sign"]) > 0
            elif defects == 3:
                assert -3.142 < float(string["deviation_phase"]) < 3.142
                assert float(string["log_deviation"]) > 0
                assert string["deviation_real"] is not None and string["deviation_imag"] is not None
            else:
                assert int(string["inner_deviation_sign"]) == (1 if i % 2 == 0 else -1)
                assert Decimal(string["inner_deviation"])*int(string["inner_deviation_sign"]) > 0
                assert -3.142 < float(string["outer_deviation_phase"]) < 3.142
                assert float(string["inner_log_deviation"]) > 0 and float(string["outer_log_deviation"]) > 0
                assert string["outer_deviation_real"] is not None and string["outer_deviation_imag"] is not None
        roots = rows(full["roots"])
        assert len(roots) == defects * count
        for i in range(count):
            cluster = roots[defects*i:defects*(i+1)]
            assert all(r["state_id"] == str(i) for r in cluster)
            assert [int(r["index"]) for r in cluster] == list(range(defects))
            assert Decimal(cluster[-2]["u_imag"]) == Decimal(cluster[-1]["u_imag"]).copy_negate()
            assert cluster[-2]["u_real"] == cluster[-1]["u_real"]
            if defects == 3:
                assert Decimal(cluster[0]["u_real"]) == 0
            if defects == 4:
                assert cluster[0]["u_real"] == cluster[1]["u_real"]
                assert Decimal(cluster[0]["u_imag"]) == Decimal(cluster[1]["u_imag"]).copy_negate()
                assert Decimal(cluster[2]["u_real"]) > Decimal(cluster[0]["u_real"])
        first = tables([str(n), "--ferromagnetic", option, "1", "--precision", precision])
        assert rows(first["states"])[0] == states[0]
        streamed = tables([*args, "--no-retain"])
        assert all(streamed[name]["rows"] == table["rows"] for name, table in full.items())
    large = tables(["100000", "--ferromagnetic", option, "2", "--precision", precision, "--roots"])
    assert all(s["multiplicity"] is None and s["converged"] for s in rows(large["states"]))
    components = {2: ("deviation",), 3: ("deviation_real", "deviation_imag"),
                  4: ("inner_deviation", "outer_deviation_real", "outer_deviation_imag")}[defects]
    logs = ("inner_log_deviation", "outer_log_deviation") if defects == 4 else ("log_deviation",)
    assert all(all(s[key] is None for key in components) and all(float(s[key]) > 80000 for key in logs)
               for s in rows(large["string"]))
    assert abs(float(rows(large["states"])[0]["gap"]) - {2: 5/3, 3: 2, 4: 15/7}[defects]) < 4e-10
    failed = tables([str(scan_n), "--ferromagnetic", option, "2", "--max-iterations", "0",
                     "--precision", precision, "--roots"], 2)
    assert all(not s["converged"] and s["gap"] is None and s["iterations"] == "0" for s in rows(failed["states"]))
    assert rows(failed["reference"])[0]["energy"] == str(scan_n - 1)

for args in (["8", option, "1"], [str(minimum_n - 1), "--ferromagnetic", option, "1"],
             ["8", "--ferromagnetic", option, "0"],
             [str(scan_n), "--ferromagnetic", option, "all", "--max-candidates", "2"],
             ["8", "--ferromagnetic", option, "1", "--max-candidates", "0"],
             ["8", "--ferromagnetic", option, "1", "--tolerance", "nan"]):
    run(args, 1)
for selection in (["--one-defect"], ["--sectors"], ["--singlet-excitation"], ["--through-lines", "4"],
                  ["--excitations", "all"], ["--quantum-numbers", "1,2"], ["--q-spectrum"],
                  ["--q-seed", "1,2"], ["--max-attempts", "5"], ["--pair-defects", "all"],
                  ["--two-pair-states", "all"], ["--triple-defects", "all"]):
    run(["8", "--ferromagnetic", option, "1", *selection], 1)
for other in options.values():
    if other != option:
        run(["8", "--ferromagnetic", option, "1", other, "1"], 1)
assert len(rows(tables(["8", "--ferromagnetic", option, "100"])["states"])) == 9 - minimum_n
assert len(rows(tables(["100000", "--ferromagnetic", option, "1", "--max-candidates", "1"])["states"])) == 1

with tempfile.TemporaryDirectory() as tmp:
    base = Path(tmp)
    args = [str(scan_n), "--ferromagnetic", option, "2", "--roots"]
    exports = ["--json", str(base / "clusters.json")]
    for extension in ("csv", "tsv"):
        for name in ("states", "reference", "string", "roots"):
            exports += [f"--{extension}-table", f"{name}={base / f'clusters.{name}.{extension}'}"]
    run([*args, *exports])
    doc = json.loads((base / "clusters.json").read_text())
    assert len(rows(doc["tables"]["states"])) == 2
    for extension, delimiter in (("csv", ","), ("tsv", "\t")):
        for name, count in (("states", 2), ("reference", 1), ("string", 2), ("roots", 2 * defects)):
            path = base / f"clusters.{name}.{extension}"
            data = list(csv.DictReader((line for line in path.read_text().splitlines() if not line.startswith("#")), delimiter=delimiter))
            assert len(data) == count, (path, data)
assert option in run(["--help"])
assert "nachtergaele-spitzer-starr-2007" in run(["--references"])
