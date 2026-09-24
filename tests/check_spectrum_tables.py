"""Scientific contracts for named spectra, references, failed solves and root coordinates."""
from decimal import Decimal, localcontext
import json
import math
import os
from pathlib import Path
import subprocess
import sys

fp128, *programs = sys.argv[1:]
precisions = ["fp64", "long-double"] + (["fp128"] if fp128.upper() in ("ON", "TRUE", "1") else [])
cases = {
    "bethe-xxx-pbc": (["4"], [], 3),
    "bethe-xxx-obc": (["4"], [], 3),
    "bethe-xxz-pbc": (["4", "--delta", "0.5"], [], 3),
    "bethe-xxz-obc": (["4", "--delta", "0.5"], [], 2),
    "bethe-lieb-liniger-pbc": (["4", "--length", "4", "--c", "1"], ["--padding", "1"], 15),
    "bethe-biquadratic-obc": (["4"], [], 3),
}

def records(table):
    return [dict(zip([c["id"] for c in table["columns"]], row)) for row in table["rows"]]

def tables(program, args, status=0):
    p = subprocess.run([program, *args, "--format", "json"], text=True, capture_output=True,
                       env=dict(os.environ, UNI20_COLOR="never"), timeout=45)
    assert p.returncode == status, (program, args, p.returncode, p.stdout, p.stderr)
    document = json.loads(p.stdout)
    assert document["status"] == "complete"
    for table in document["tables"].values():
        assert table["summary"]["Outcome"] == ("success" if status == 0 else "partial")
    return document["tables"]

for program in programs:
    name = Path(program).name
    base, extra, count = cases[name]
    for precision in precisions:
        common = [*base, *extra, "--excitations", "all", "--roots", "--precision", precision]
        doc = tables(program, common)
        rows, reference = records(doc["states"]), records(doc["reference"])[0]
        assert len(rows) == count and all(r["converged"] for r in rows), (name, len(rows), count)
        assert [Decimal(r["energy"]) for r in rows] == sorted(Decimal(r["energy"]) for r in rows)
        assert [r["state_id"] for r in rows] == [str(i) for i in range(count)]
        assert reference["state_id"] == str(count)
        with localcontext() as context:
            context.prec = 70
            tolerance = Decimal("2e-28" if precision == "fp128" else "2e-12")
            for row in rows:
                assert abs(Decimal(row["gap"]) - Decimal(row["energy"]) + Decimal(reference["energy"])) < tolerance
            if "xxx" in name or "lieb-liniger" in name:
                for row in [*rows, reference]:
                    roots = [r for r in records(doc["roots"]) if r["state_id"] == row["state_id"]]
                    if "xxx" in name:
                        reconstructed = Decimal("1" if name.endswith("pbc") else "0.75")
                        reconstructed -= sum(2 / (1 + Decimal(r["rapidity"])**2) for r in roots)
                    else:
                        reconstructed = sum(Decimal(r["k"])**2 for r in roots)
                    assert abs(reconstructed - Decimal(row["energy"])) < tolerance, (name, reconstructed, row)
        ids = {r["state_id"] for r in [*rows, reference]}
        for table in doc.values():
            assert all(r["state_id"] in ids for r in records(table))
        if name.endswith("obc"):
            assert all("p" not in r and "momentum_index" not in r for r in rows)
        replay = tables(program, [*common, "--no-retain"])
        assert all(replay[key]["rows"] == doc[key]["rows"] for key in doc)
        # Failed candidates never enter the ranked table; estimates get distinct identities.
        failed = tables(program, ["8" if "biquadratic" in name else "6", *common[1:], "--max-iterations", "0"], status=2)
        assert all(r["converged"] for r in records(failed["states"]))
        assert not records(failed["reference"])[0]["converged"]
        assert "failed" in failed and not records(failed["failed"])[0]["converged"], (name, precision)
        assert all(r["gap"] is None for r in records(failed["states"]))
        all_states = [r for key in ("states", "reference", "failed") for r in records(failed[key])]
        assert len({r["state_id"] for r in all_states}) == len(all_states)
        # An exactly solved polarized sector can survive a failed global reference.
        if "xxx" in name or "xxz" in name or "biquadratic" in name:
            sector = ["--spin", "2"] if "xxx" in name else ["--sz", "2"] if "xxz" in name else ["--through-lines", "4"]
            vacuum = tables(program, [*common, *sector, "--max-iterations", "0"], status=2)
            assert len(records(vacuum["states"])) == 1
            assert records(vacuum["states"])[0]["converged"] and records(vacuum["states"])[0]["gap"] is None

        if "xxz" in name:
            negative = tables(program, ["6", "--delta", "-0.9", "--sectors", "--roots", "--precision", precision])
            rows = records(negative["states"])
            assert len(rows) == 7
            assert [r["sz"] for r in rows] == [-3, -2, -1, 0, 1, 2, 3]
            assert all(r["spin_reversed"] == (r["sz"] < 0) for r in rows)
            assert "rank-subtracted" in negative["roots"]["metadata"]["Residual convention"]
            for root in records(negative["roots"]):
                assert root["lambda"] is not None
                assert abs(float(root["rapidity"]) - math.sqrt(0.1/1.9)*math.tanh(float(root["lambda"]))) < 2e-15
            if name.endswith("obc"):
                assert not records(negative["boundary_roots"])
                boundary = tables(program, ["2", "--delta", "4", "--roots", "--precision", precision])
                assert not records(boundary["roots"])
                b, = records(boundary["boundary_roots"])
                assert b["kind"] == "imaginary" and b["state_id"] == "0" and b["quantum_number"] == 1
                assert abs(float(b["inverse_square"]) + 0.2) < 2e-15
                assert abs(float(b["inverse_square"]) - 0.36*math.expm1(-float(b["log_distance"]))) < 2e-15
                assert abs(Decimal(records(boundary["states"])[0]["energy"]) + Decimal("1.5")) < tolerance
        if name == "bethe-xxx-pbc":
            branch = tables(program, ["5", "--spinons", "--roots", "--precision", precision])
            assert len(records(branch["states"])) == len(records(branch["spinons"])) == 3
            for point in records(branch["spinons"]):
                k = float(point["k"])
                assert abs(float(point["epsilon_inf"]) - math.pi/2*math.sin(k)) < 2e-15
                state = records(branch["states"])[int(point["state_id"])]
                assert abs(float(point["bulk_subtracted_energy"]) - float(state["energy"]) + 5*(0.25-math.log(2))) < 2e-15
        if name == "bethe-biquadratic-obc":
            sectors = tables(program, ["4", "--sectors", "--roots", "--precision", precision])
            assert [r["through_lines"] for r in records(sectors["states"])] == ["0", "2", "4"]
            assert [r["multiplicity"] for r in records(sectors["states"])] == ["1", "8", "55"]
            overflow = tables(program, ["46", "--through-lines", "46", "--precision", precision])
            assert records(overflow["states"])[0]["multiplicity"] is None
print("Spectrum, reference and root-coordinate contracts passed")
