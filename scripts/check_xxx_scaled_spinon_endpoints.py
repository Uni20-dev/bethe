#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Verify the stored scaled-endpoint replay against the frozen PR #93 reference.

Usage: python scripts/check_xxx_scaled_spinon_endpoints.py STORE_DIRECTORY
Raw geometry rows stay in store; stdout is the compact validation record.
"""

import argparse
from collections import Counter
from decimal import Decimal as D, getcontext
import hashlib
import json
from pathlib import Path

getcontext().prec = 200


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    data = [json.loads(line) for line in path.read_text().splitlines()]
    result = {row["id"]: row for row in data}
    assert len(result) == len(data) == 1359
    return result


def main(root):
    frozen = root.parent / "deep-geometry-precision-probes-20261101"
    reference_path = frozen / "geometry-bits-512.jsonl"
    original_receipt = json.loads((frozen / "geometry-receipt.json").read_text())
    original_run = next(r for r in original_receipt["runs"] if r["name"] == "geometry-bits-512")
    assert digest(reference_path) == original_run["sha256"]
    reference = rows(reference_path)
    assert Counter(r["status"] for r in reference.values()) == {"inside": 53, "outside": 1306}
    receipt = json.loads((root / "receipt.json").read_text())
    results = {}
    for run in receipt["runs"]:
        path = root / (run["name"] + ".jsonl")
        assert run["exit_code"] == 0 and digest(path) == run["sha256"]
        data = rows(path)
        assert data.keys() == reference.keys()
        kind = run["name"].removeprefix("geometry-")
        differences = {"draws": [], "deep": [], "edge": []}
        beta_errors, failures, retries = [], [], []
        for identity, value in data.items():
            expected = reference[identity]
            if value["status"] == "error":
                failures.append({"id": identity, "reason": value["reason"]})
                continue
            assert value["status"] == expected["status"], (kind, identity)
            if value["refinements"]:
                retries.append({key: value[key] for key in ("id", "working_bits", "refinements")})
            if value["status"] == "inside":
                group = "edge" if identity.startswith("stress_edge") else "deep" if identity.startswith("stress_") else "draws"
                for a, b in zip(value["branches"], expected["branches"], strict=True):
                    differences[group].append(abs(D(a["log_weight"]) - D(b["log_weight"])))
                    beta_errors.extend(abs(D(x) - D(y)) for x, y in zip(a["beta"], b["beta"], strict=True))
        if kind in ("fp64", "extended80"):
            assert len(failures) == 1 and failures[0]["id"] == "stress_edge_pi2_edge_s0"
            assert not retries
        else:
            assert not failures
        results[kind] = {
            "status_counts": dict(Counter(r["status"] for r in data.values())),
            "resolved_support_mismatches": 0,
            "failures": failures, "retries": retries,
            "max_absolute_beta_error": format(max(beta_errors), ".12E"),
            "max_absolute_log_weight_error": {g: format(max(v), ".12E") if v else None for g, v in differences.items()},
            "output_sha256": run["sha256"],
        }
    result = {
        "schema": "bethe-scaled-spinon-endpoints-v1", "date": "2026-11-01",
        "scope": "Stage 1: geometry only. No density scaling, recovery-policy change, sampling convergence or cutoff-convergence claim.",
        "base_commit": receipt["base_commit"], "program_sha256": receipt["program_sha256"],
        "raw_store_directory": root.name, "reference_store_directory": frozen.name,
        "reference_sha256": original_run["sha256"], "reference_bits": 512,
        "design": "R=1040; 4x16 exact 52-bit midpoint draws, seed7303, seven queries and three sections, plus four prescribed deep configurations and one nearly degenerate pair; 1359 solves per format.",
        "representation": "Positive endpoint distance and J=1 energy share log_scale=-abs(beta). The cutoff is compared logarithmically. Ordinary projections below epsilon^2 carry explicit absolute omission bounds into support and conditioning screens.",
        "edge_caveat": "The edge energy is regenerated in each format, so differences there include input rounding. MPFR can retry locally; native formats report unresolved conditioning.",
        "results": results,
        "artifacts_sha256": {name: digest(root / name) for name in ("receipt.json", "source.patch", "xxx_spinon_endpoint.hpp", "run_geometry.py")},
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("store_directory", type=Path)
    main(parser.parse_args().store_directory)
