#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Regenerate compact two/four/six point controls, without outer integration."""
import argparse
import json
from pathlib import Path
import platform
import sys
import time

import mpmath as mp

import audit_xxx_qtm_partitions as audit
from audit_xxx_qtm_vertex import encode
import reference_xxx_qtm as qtm
import reference_xxx_qtm_vertex as bridge
from _reference_support import source_hashes


def coefficients(beta, digits):
    with mp.workdps(digits):
        _, t = bridge.shifted_points(beta, ".25")
        rows = []
        for holes, particles in bridge.assignments(len(beta)):
            result = qtm.prl([t[j]-mp.j for j in holes], [t[j] for j in particles])
            rows.append(dict(holes=holes, particles=particles,
                             reduced_coefficient=audit.reduced_coefficient(result, "prl", len(holes))))
        return rows


def fixture(name, low_digits, high_digits):
    beta = audit.FIXTURES[name]
    low, high = (coefficients(beta, digits) for digits in (low_digits, high_digits))
    with mp.workdps(high_digits):
        k0 = audit.boundary_factor(beta)
        total = k0*mp.fsum(row["reduced_coefficient"] for row in high)
        vertex = bridge.vertex_target(beta, terms=128, order=28)
        # One independent quadrature partition per sector. All partitions are
        # compared by audit_xxx_qtm_partitions.py, a deliberately larger audit.
        selected = audit.matched_partition(beta, high[0]["holes"], ".1", maxdegree=9, split=False)
        refinement = max(audit.scaled_difference(a["reduced_coefficient"], b["reduced_coefficient"])
                         for a, b in zip(low, high))
        vertex_difference = audit.scaled_difference(total, vertex["amplitude"])
        checks = dict(partition_precision=refinement < mp.power(10, -(low_digits-8)),
                      selected_jmp_prl=selected["reduced_scaled_difference"] < mp.power(10, -(high_digits-10)),
                      vertex=vertex_difference < mp.power(10, -(high_digits-10)))
        return encode(dict(fixture=name, beta=beta, spinons=len(beta),
                           digits=[low_digits, high_digits], interior_offset=".25",
                           partitions=high, boundary_sine_squared=k0, boundary_total=total,
                           vertex_target=vertex, worst_precision_difference=refinement,
                           vertex_scaled_difference=vertex_difference,
                           selected_quadrature=dict(holes=selected["holes"], offset=".1",
                               maxdegree=9, split=False,
                               jmp_reduced=selected["routes"]["jmp"]["reduced_coefficient"],
                               prl_reduced=selected["routes"]["prl"]["reduced_coefficient"],
                               scaled_difference=selected["reduced_scaled_difference"]),
                           checks=checks, passed=all(checks.values())))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--digits", type=int, default=32)
    parser.add_argument("--refined-digits", type=int, default=42)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output must not exist")
    if args.digits < 30 or args.refined_digits < args.digits+10:
        parser.error("require at least 30 decimal digits and ten extra refinement digits")
    source = source_hashes(__file__)
    start = time.perf_counter()
    results = []
    for name in ("two", "four", "six-bulk"):
        results.append(fixture(name, args.digits, args.refined_digits))
        print(name, results[-1]["checks"], flush=True)
    if source != source_hashes(__file__):
        raise RuntimeError("reference sources changed during fixture generation")
    result = dict(schema="xxx-qtm-pointwise-fixture-v1",
                  scope="All PRL partitions and one independent JMP partition per sector; no outer integral",
                  source_sha256=source, invocation=sys.argv[1:], python=platform.python_version(),
                  mpmath=mp.__version__, wall_seconds=time.perf_counter()-start,
                  fixtures=results, passed=all(row["passed"] for row in results))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as output:
        json.dump(result, output, indent=2)
        output.write("\n")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
