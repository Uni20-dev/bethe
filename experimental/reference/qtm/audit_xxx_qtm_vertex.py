#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Small pointwise QTM/vertex audit; one process, no cluster or outer integral."""
import argparse
from dataclasses import is_dataclass, asdict
import json
import platform
from pathlib import Path
import sys
import time

import mpmath as mp

from _reference_support import REFERENCE_DIR, source_hashes as public_source_hashes

import reference_xxx_qtm_vertex as bridge


FIXTURES = {
    "two": ["-1.1", ".7"],
    "four": ["-2.1", "-1.3", "-.4", ".3"],
    "six-bulk": ["-2.1", "-1.3", "-.4", ".3", "1.2", "2"],
    "six-asymmetric": ["-3", "-1.1", "-.2", ".4", "1.7", "3.2"],
    "six-clustered": ["-1.7", "-1.65", "-.4", ".3", "1.2", "2"],
}


def encode(x):
    if isinstance(x, mp.mpc):
        return dict(real=mp.nstr(x.real, mp.mp.dps), imag=mp.nstr(x.imag, mp.mp.dps))
    if isinstance(x, mp.mpf):
        return mp.nstr(x, mp.mp.dps)
    if isinstance(x, (tuple, list)):
        return [encode(v) for v in x]
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if is_dataclass(x):
        return encode(asdict(x))
    return x


def hashes():
    return public_source_hashes(__file__, REFERENCE_DIR/"test_reference_xxx_qtm_vertex.py")

def run_fixture(name, low_digits, high_digits):
    beta = FIXTURES[name]
    offsets = ("1e-4", "1e-7", "1e-10")
    raw, records, targets = [], [], []
    # All settings are declared, not stopped on the first attractive answer.
    settings = [(low_digits, e, 64, 20) for e in offsets]
    settings.append((high_digits, offsets[-1], 128, 28))
    for digits, offset, terms, order in settings:
        with mp.workdps(digits):
            start = time.perf_counter()
            def progress(done, count):
                if count > 6 and done % 5 == 0:
                    print(f"{name} dps={digits} offset={offset} partitions={done}/{count}", flush=True)
            result = bridge.partition_sum(beta, offset, progress=progress)
            target = bridge.vertex_target(beta, terms=terms, order=order)
            result.update(digits=digits, vertex=target,
                          ratio=result["total"]/target["amplitude"],
                          real_relative_error=abs(result["total"].real/target["amplitude"]-1),
                          imaginary_relative=abs(result["total"].imag/target["amplitude"]),
                          seconds=time.perf_counter()-start)
            raw.append(result)
            records.append(encode(result))
            targets.append(target)
            print(f"{name} dps={digits} offset={offset}: "
                  f"real error={mp.nstr(result['real_relative_error'],6)} "
                  f"imag={mp.nstr(result['imaginary_relative'],6)} "
                  f"condition={mp.nstr(result['cancellation_condition'],5)} "
                  f"{result['seconds']:.2f}s", flush=True)
    with mp.workdps(high_digits):
        target = targets[-1]["amplitude"]
        extrap = [bridge.richardson_real(a, b) for a, b in zip(raw[:2], raw[1:3])]
        extrap_errors = [abs(z/target-1) for z in extrap]
        precision_difference = abs(raw[-1]["total"]-raw[-2]["total"])/abs(target)
        vertex_difference = abs(targets[-1]["amplitude"]-targets[-2]["amplitude"])
        vertex_budget = sum(t["contour_tail_error"]+t["arithmetic_indicator"] for t in targets[-2:])
        real_rates = [b["real_relative_error"]/a["real_relative_error"]
                      for a, b in zip(raw[:2], raw[1:3])]
        imag_rates = [b["imaginary_relative"]/a["imaginary_relative"]
                      for a, b in zip(raw[:2], raw[1:3])]
        # These are numerical screens on these fixtures, not interval bounds.
        checks = dict(
            complete_partitions=all(len(r["pieces"]) == len(bridge.assignments(len(beta))) for r in raw),
            boundary_real=raw[-1]["real_relative_error"] < mp.mpf("1e-16"),
            boundary_imaginary=raw[-1]["imaginary_relative"] < mp.mpf("1e-7"),
            quadratic_real=all(mp.mpf("1e-7") < x < mp.mpf("1e-5") for x in real_rates),
            linear_imaginary=all(mp.mpf("1e-4") < x < mp.mpf("1e-2") for x in imag_rates),
            qtm_precision=precision_difference < mp.power(10, -(low_digits-12)),
            vertex_tail_and_precision=vertex_difference <= vertex_budget,
            extrapolated_agreement=extrap_errors[-1] < mp.mpf("1e-24"),
        )
        result = encode(dict(
            fixture=name, beta=beta, n=len(beta), records=records,
            real_boundary_extrapolations=extrap, extrapolation_relative_errors=extrap_errors,
            qtm_precision_relative_difference=precision_difference,
            vertex_refinement_difference=vertex_difference, vertex_refinement_budget=vertex_budget,
            real_error_ratios=real_rates, imaginary_ratios=imag_rates,
            checks=checks, passed=all(checks.values())))
        print(f"{name}: {'PASS' if result['passed'] else 'FAIL'} {checks}", flush=True)
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", nargs="+", choices=FIXTURES,
                        default=["two", "four", "six-bulk", "six-asymmetric"])
    parser.add_argument("--low-digits", type=int, default=35)
    parser.add_argument("--high-digits", type=int, default=50)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output must not exist")
    if args.low_digits < 30 or args.high_digits < args.low_digits+10:
        parser.error("require low-digits>=30 and at least ten extra high digits")
    source = hashes()
    wall, cpu = time.perf_counter(), time.process_time()
    fixtures = [run_fixture(name, args.low_digits, args.high_digits) for name in args.fixture]
    if source != hashes():
        raise RuntimeError("implementation changed during audit; refusing to save a mixed checkpoint")
    result = dict(
        schema="xxx-qtm-vertex-pointwise-v1",
        scope="isotropic pointwise boundary correspondence, not an outer integral or proof",
        conventions="Pauli zz QTM sum over all assignments versus converted Szz vertex momentum density",
        python=platform.python_version(), mpmath=mp.__version__, platform=platform.platform(),
        invocation=sys.argv[1:], source_sha256=source,
        wall_seconds=time.perf_counter()-wall, cpu_seconds=time.process_time()-cpu,
        processes=1, fixtures=fixtures, passed=all(f["passed"] for f in fixtures))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as f:
        json.dump(result, f, indent=2)
        f.write("\n")
    print(f"wrote {args.output}: {'PASS' if result['passed'] else 'FAIL'}", flush=True)
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
