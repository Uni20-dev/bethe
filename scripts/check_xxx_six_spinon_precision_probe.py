#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Summarize retained opt-in precision probes; never discard failed points.

Usage: python scripts/check_xxx_six_spinon_precision_probe.py STORE_DIRECTORY
Prints compact JSON. Raw inputs, outputs, runners and receipts remain in store.
Decimal arithmetic keeps tiny densities and large inverse weights representable.
"""

import argparse
from collections import Counter
from decimal import Decimal as D, getcontext
import hashlib
import json
from pathlib import Path

getcontext().prec = 200


def relative(value, reference):
    if not reference:
        assert not value, "unexpected nonzero value at a zero reference"
        return D(0)
    return abs((value - reference) / reference)


def maximum(values):
    return format(max(values, default=D(0)), ".12E")


def group(identity):
    return "edge" if identity.startswith("stress_edge") else "deep" if identity.startswith("stress_") else "draws"


def main(root):
    def load(name):
        return json.loads((root / name).read_text())

    def rows(name):
        data = [json.loads(line) for line in (root / (name + ".jsonl")).read_text().splitlines()]
        result = {(row["id"], row.get("control", -1)): row for row in data}
        assert len(result) == len(data), "duplicate records"
        return result

    receipts = {p.name: load(p.name) for p in sorted(root.glob("*-receipt.json"))}
    runs = {}
    for receipt in receipts.values():
        for run in receipt["runs"]:
            assert run["exit_code"] == 0, run
            assert hashlib.sha256((root / (run["name"] + ".jsonl")).read_bytes()).hexdigest() == run["sha256"]
            runs[run["name"]] = run

    reference = rows("geometry-bits-512")
    assert len(reference) == 1359
    assert Counter(r["status"] for r in reference.values()) == {"inside": 53, "outside": 1306}
    geometry = {}
    for name in sorted(n for n in runs if n.startswith("geometry-")):
        data = rows(name)
        assert data.keys() == reference.keys()
        errors = {g: [] for g in ("draws", "deep", "edge")}
        beta_errors, mismatches, retries = [], [], []
        failures = []
        for key, row in data.items():
            ref = reference[key]
            if row["status"] == "error":
                failures.append(row)
                continue
            if row["status"] != ref["status"]:
                mismatches.append(row["id"])
                continue
            if row["refinements"]:
                retries.append({k: row[k] for k in ("id", "working_bits", "refinements")})
            if row["status"] == "inside":
                for value, expected in zip(row["branches"], ref["branches"], strict=True):
                    errors[group(row["id"])].append(abs(D(value["log_weight"]) - D(expected["log_weight"])))
                    beta_errors.extend(abs(D(v) - D(e)) for v, e in zip(value["beta"], expected["beta"], strict=True))
        geometry[name.removeprefix("geometry-")] = {
            "status_counts": dict(Counter(r["status"] for r in data.values())),
            "support_mismatches": mismatches,
            "failure_reasons": dict(Counter(r["reason"] for r in failures)),
            "retries": retries,
            "max_absolute_beta_error": maximum(beta_errors) if beta_errors else None,
            "max_absolute_log_weight_error": {g: maximum(v) if v else None for g, v in errors.items()},
        }

    selected = {r["id"]: r for r in load("density-selection.json")}
    assert len(selected) == 42
    density_ref = rows("density-bits-320")
    assert len(density_ref) == 84 and all(r["status"] == "value" for r in density_ref.values())
    density = {}
    for name in sorted(n for n in runs if n.startswith("density-")):
        data = rows(name)
        assert data.keys() == density_ref.keys()
        density[name.removeprefix("density-")] = {
            "status_counts": dict(Counter(r["status"] for r in data.values())),
            "max_relative_value_error": maximum(
                relative(D(r["value"]), D(density_ref[k]["value"])) for k, r in data.items() if r["status"] == "value"
            ),
            "failure_reasons": dict(Counter(r["reason"] for r in data.values() if r["status"] != "value")),
            "cpu_seconds": runs[name]["cpu_seconds"],
        }

    def screening(low, high):
        accepted = []
        for identity in sorted({key[0] for key in high}):
            budgets, values, passed = [], [], True
            for control in (0, 1):
                a, b = low[identity, control], high[identity, control]
                if a["status"] != "value" or b["status"] != "value":
                    passed = False
                    break
                av, bv = D(a["value"]), D(b["value"])
                arithmetic = (D(a["arithmetic_indicator"]) + D(b["arithmetic_indicator"]) + abs(av - bv)
                              + abs((D(a["log_scale"]) - D(b["log_scale"])).exp() - 1) * abs(bv))
                budget = D(b["contour_tail"]) + arithmetic
                passed &= arithmetic > 0 and budget <= D("1e-6") * abs(bv)
                budgets.append(budget)
                values.append(bv)
            if passed and abs(values[1] - values[0]) <= sum(budgets):
                accepted.append(identity)
        return {
            "accepted_labels": len(accepted),
            "pilot_labels": sum("_pi2_w2_" in s and s.startswith("draw_") for s in accepted),
            "stress_labels": sum(s.startswith("stress_") for s in accepted),
        }

    screens = {f"{low}/{high}": screening(rows("density-" + low), rows("density-" + high))
               for low, high in [("bits-53", "bits-64"), ("bits-64", "bits-80"), ("bits-80", "bits-96"),
                                 ("bits-96", "bits-128"), ("extended80", "fp128"), ("bits-256", "bits-320")]}

    def weighted(data, inputs):
        return {identity: D(r["value"]) * D(inputs[identity]["log_weight"]).exp()
                for (identity, control), r in data.items() if control == 1 and r["status"] == "value"}

    expected = weighted(density_ref, selected)
    pilot = {key: val for key, val in expected.items() if key.startswith("draw_") and "_pi2_w2_" in key}
    assert len(pilot) == 20
    pilot_mean = sum(pilot.values()) / 64
    native = rows("density-fp64")
    native_pilot_failures = [identity for identity in pilot if native[identity, 1]["status"] != "value"]
    density["fp64"]["pilot"] = {
        "returned_labels": len(pilot) - len(native_pilot_failures),
        "failed_labels": len(native_pilot_failures),
        "failed_reference_weight_fraction": format(
            sum(pilot[identity] for identity in native_pilot_failures) / sum(pilot.values()), ".12E"
        ),
    }
    archive = root.parent / "public-six-spinon-tutorial-20261031/pilot-09-geometry-precision/geometry-1536.json"
    archived_mean = D(json.loads(archive.read_text())["tables"]["spectrum"]["rows"][0][3]["value"])
    assert relative(pilot_mean, archived_mean) < D("1e-60")
    joint = {}
    for kind, inputs in load("joint-inputs.json").items():
        data = rows("joint-" + kind)
        points = {r["id"]: r for r in inputs["included"]}
        assert set(points) | {r["id"] for r in inputs["excluded"]} == set(selected)
        values = weighted(data, points)
        assert set(pilot) <= values.keys(), "cannot compute complete pilot mean with missing labels"
        joint[kind] = {
            "returned_labels": len(values),
            "excluded_geometry_labels": inputs["excluded"],
            "density_failure_count": sum(r["status"] != "value" for r in data.values()),
            "max_relative_weighted_error": maximum(relative(v, expected[k]) for k, v in values.items()),
            "max_relative_weighted_draw_error": maximum(relative(v, expected[k]) for k, v in values.items() if k.startswith("draw_")),
            "pilot_mean_relative_error": format(relative(sum(values[k] for k in pilot) / 64, pilot_mean), ".12E"),
        }

    strong320, strong384 = rows("stress-strong-bits-320"), rows("stress-strong-bits-384")
    assert len(strong320) == len(strong384) == 6
    strong = {
        "labels": sorted({k[0] for k in strong320}),
        "controls": [[3120, 128], [3136, 132]],
        "max_relative_320_to_384_change": maximum(relative(D(r["value"]), D(strong384[k]["value"])) for k, r in strong320.items()),
        "max_relative_change_from_original_320": maximum(relative(D(r["value"]), D(density_ref[k]["value"])) for k, r in strong384.items()),
        "screen_320_384": screening(strong320, strong384),
        "max_relative_tail_384": maximum(abs(D(r["contour_tail"]) / D(r["value"])) for r in strong384.values()),
        "max_relative_arithmetic_indicator_384": maximum(abs(D(r["arithmetic_indicator"]) / D(r["value"])) for r in strong384.values()),
    }
    strong["lower_precision"] = {}
    for kind in ("extended80", "fp128", "bits-80", "bits-96"):
        data = rows("stress-strong-" + kind)
        assert data.keys() == strong384.keys()
        strong["lower_precision"][kind] = {
            "status_counts": dict(Counter(r["status"] for r in data.values())),
            "max_relative_value_error": maximum(relative(D(r["value"]), D(strong384[k]["value"]))
                                                for k, r in data.items() if r["status"] == "value"),
            "failure_reasons": dict(Counter(r["reason"] for r in data.values() if r["status"] != "value")),
        }
    strong["lower_precision_screens"] = {
        "80/96": screening(rows("stress-strong-bits-80"), rows("stress-strong-bits-96")),
        "extended80/fp128": screening(rows("stress-strong-extended80"), rows("stress-strong-fp128")),
    }

    checked = {}
    pilot_draws = {f"draw_{r}_{s}": D(0) for r in range(4) for s in range(16)}
    for key, value in pilot.items():
        pilot_draws[key.split("_pi2_w2_")[0]] += value
    for name in sorted(n for n in runs if n.startswith("checked-")):
        data = rows(name)
        assert {k[0] for k in data} == set(pilot_draws)
        complete = all(r["status"] == "value" for r in data.values())
        checked[name] = {
            "status_counts": dict(Counter(r["status"] for r in data.values())),
            "failure_reasons": dict(Counter(r["reason"] for r in data.values() if r["status"] != "value")),
            "supported_roots": sum(r["supported_roots"] for r in data.values()),
            "max_relative_accepted_draw_error": maximum(relative(D(r["value"]), pilot_draws[k[0]]) for k, r in data.items() if r["status"] == "value"),
            "mean_relative_error": format(relative(sum(D(r["value"]) for r in data.values()) / 64, pilot_mean), ".12E") if complete else None,
            "cpu_seconds": runs[name]["cpu_seconds"],
        }

    result = {
        "schema": "bethe-six-spinon-precision-probe-v1", "date": "2026-11-01",
        "scope": "Arithmetic probes at R=1040; no claim of sampling convergence, cutoff convergence or universal minimum precision. Raw density values are not controller-accepted results.",
        "base_commit": "fc8320cc148bb3b5243b9cbbdac20bcd16f2a3f0",
        "uni20_commit": "1cb2474eb0a54e13ec7d5f2c1823752020b2425c",
        "compiler": "GNU 13.3.0 Release", "cpu": "12th Gen Intel(R) Core(TM) i7-12700KF",
        "raw_store_directory": "deep-geometry-precision-probes-20261101",
        "reference_bits": {"geometry": 512, "density": 320, "strong_contour_density": 384},
        "selection": "All 20 pilot labels, 10 prescribed stress labels, and both branches of the deepest sampled configuration at each of six other queries; selected before density evaluation.",
        "geometry_design": "4x16 exact 52-bit midpoint draws, seed 7303, seven queries, three sections, four deep configurations plus one nearly degenerate pair; 1359 solves at each format.",
        "edge_caveat": "omega=pi*sin(q/2)-2^-40 is regenerated at each format; edge errors include input rounding, which internal solver retries cannot repair.",
        "weighting_caveat": "Joint raw comparisons recombine each format's returned density and log weight using 200-digit Decimal; only full_controller_pilot exercises weighting and accumulation at the stated geometry precision.",
        "units": "bits-N is MPFR significand bits. Native extended80 has 64 significand bits; fp128 has 113. MPFR has a wider exponent range than native fp64.",
        "controls": [[1040, 32], [1056, 36]], "relative_tolerance": "1e-6",
        "geometry": geometry, "raw_density": density, "offline_recovery_screen": screens,
        "joint_raw_weighted_observations": joint, "strong_contour_reference": strong,
        "full_controller_pilot": checked, "reference_pilot_mean": format(pilot_mean, ".40E"),
        "archived_1536_geometry_256_density_pilot": {
            "sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
            "reference_mean_relative_change": format(relative(pilot_mean, archived_mean), ".12E"),
        },
        "timing_scope": "Single-process CPU seconds, two workers per batch; raw density timing excludes recovery checks and geometry. Not an end-to-end speedup benchmark.",
        "receipts": {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in receipts},
        "source_artifacts": {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                             for name in ("probe-source.cpp", "source.patch", "checked-probe-source.cpp", "checked-source.patch")},
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("store_directory", type=Path)
    main(parser.parse_args().store_directory)
