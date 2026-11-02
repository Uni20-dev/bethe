#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Matched complex JMP/PRL partition amplitudes, with no outer quadrature.

The common imaginary translation only changes the explicit momentum factor
in the PRL expression. Compare the remaining complex coefficient B_H and
form the real-boundary amplitude analytically, without a regulator fit.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time

import mpmath as mp

import reference_xxx_qtm as qtm
import reference_xxx_qtm_vertex as bridge
from audit_xxx_qtm_vertex import FIXTURES, encode
from _reference_support import REFERENCE_DIR, REPO_ROOT, source_hashes as public_source_hashes



def decode_complex(value):
    return mp.mpc(value["real"], value["imag"])


def scaled_difference(a, b):
    """Absolute near zero, relative for magnitudes greater than one."""
    return abs(a-b)/max(1, abs(a), abs(b))


def kinematic_factor(beta, offset):
    """Common sine-squared factor for ALL assignments of the same rapidities."""
    beta, t = bridge.shifted_points(beta, offset)
    return mp.sin(mp.fsum(qtm.momentum(z) for z in t)/2)**2


def boundary_factor(beta):
    """Exactly at real beta, using p(beta/pi)=2*atan(exp(beta))."""
    beta = tuple(map(mp.mpf, beta))
    return mp.sin(mp.fsum(mp.atan(mp.exp(b)) for b in beta))**2


def reduced_coefficient(result, route, ell):
    """Remove the explicit sine factor algebraically, never by division.

    The JMP and PRL determinants/factors need not agree separately. Their
    products B_H have the same normalization and must agree pointwise.
    """
    if route == "jmp":
        return (-1)**(ell-1)*result["pair_factor"]*mp.fprod(result["determinants"])
    if route == "prl":
        return result["scalar"]*result["pair_factor"]*mp.fprod(result["determinant_slopes"])
    raise ValueError("route must be jmp or prl")


def matched_partition(beta, holes, offset, *, maxdegree=7, split=True):
    beta, t = bridge.shifted_points(beta, offset)
    holes = tuple(holes)
    assignment = next((a for a in bridge.assignments(len(beta)) if a[0] == holes), None)
    if assignment is None:
        raise ValueError("holes must be a sorted valid half-size index tuple")
    holes, particles = assignment
    u, v = tuple(t[j]-mp.j for j in holes), tuple(t[j] for j in particles)
    ell = len(holes)
    common = kinematic_factor(beta, offset)
    physical = boundary_factor(beta)
    routes = {}
    for route in ("jmp", "prl"):
        start = time.perf_counter()
        result = (qtm.jmp(u, v, maxdegree=maxdegree, split=split)
                  if route == "jmp" else qtm.prl(u, v))
        coefficient = reduced_coefficient(result, route, ell)
        routes[route] = dict(
            amplitude=result["amplitude"], reduced_coefficient=coefficient,
            boundary_amplitude=physical*coefficient,
            density=result["density"],
            determinant_coefficients=result[
                "determinants" if route == "jmp" else "determinant_slopes"],
            null_determinants=result["null_determinants"],
            null_determinants_scaled=[
                abs(d)/max(1, mp.norm(m, mp.inf)**ell)
                for m, d in zip(result["matrices"], result["null_determinants"])],
            pair_factor=result["pair_factor"], seconds=time.perf_counter()-start,
            reconstruction_scaled_error=scaled_difference(
                result["amplitude"], common*coefficient))
        if route == "jmp":
            routes[route].update(
                max_quad_error=result["max_quad_error"],
                literal_first_column_amplitude=result["literal_first_column_amplitude"],
                maxdegree=maxdegree, split=split)
        else:
            routes[route]["scalar"] = result["scalar"]
    a, b = routes["jmp"], routes["prl"]
    return dict(
        holes=holes, particles=particles, u=u, v=v, digits=mp.mp.dps,
        offset=mp.mpf(offset), sine_squared=common, boundary_sine_squared=physical,
        single_assignment_measure=qtm.measure(ell), routes=routes,
        sine_identity_scaled_error=scaled_difference(qtm.sine_squared(u, v), common),
        amplitude_scaled_difference=scaled_difference(a["amplitude"], b["amplitude"]),
        reduced_scaled_difference=scaled_difference(
            a["reduced_coefficient"], b["reduced_coefficient"]))


def verify_checkpoint(path):
    """Accept only a passed checkpoint of the same public reference sources."""
    checkpoint = json.loads(path.read_text())
    if (checkpoint.get("schema") != "xxx-qtm-vertex-pointwise-v1"
            or checkpoint.get("passed") is not True):
        raise ValueError("vertex checkpoint must have the expected schema and pass")
    expected_sources = checkpoint.get("source_sha256", {})
    required = {name for name in public_source_hashes()
                if Path(name).name.startswith("reference_")
                or Path(name).name == "_reference_support.py"}
    if not required.issubset(expected_sources):
        raise ValueError("vertex checkpoint is missing its numerical source hashes")
    for name, expected in expected_sources.items():
        source = (REPO_ROOT/name).resolve()
        if not source.is_relative_to(REPO_ROOT) or not source.is_file():
            raise ValueError(f"checkpoint source is not in this public checkout: {name}")
        if hashlib.sha256(source.read_bytes()).hexdigest() != expected:
            raise ValueError(f"checkpoint source changed: {name}")
    return checkpoint


def hashes():
    return public_source_hashes(__file__, REFERENCE_DIR/"test_reference_xxx_qtm_partitions.py")


def fresh_reference(name, digits):
    """Generate the complex boundary terms and vertex target in this checkout."""
    with mp.workdps(digits):
        beta = FIXTURES[name]
        record = bridge.partition_sum(beta, "1e-10")
        record["vertex"] = bridge.vertex_target(beta, terms=128, order=28)
        return dict(fixture=name, beta=beta, records=[encode(record)])


def run_fixture(name, previous, low_digits, high_digits, high_maxdegree=9):
    beta = FIXTURES[name]
    if beta != previous["beta"]:
        raise ValueError("fixture no longer matches the previous checkpoint")
    # Joint refinement: more arithmetic, different inner splits and degree,
    # and a smaller offset. Each setting compares the routes at IDENTICAL u,v.
    settings = [(low_digits, ".25", 7, True), (high_digits, ".1", high_maxdegree, False)]
    raw = []
    records = []
    for digits, offset, maxdegree, split in settings:
        with mp.workdps(digits):
            entries = []
            for holes, _ in bridge.assignments(len(beta)):
                entry = matched_partition(beta, holes, offset, maxdegree=maxdegree, split=split)
                entries.append(entry)
                print(f"{name} dps={digits} eta={offset} H={holes} "
                      f"JMP/PRL={mp.nstr(entry['reduced_scaled_difference'], 5)} "
                      f"{sum(r['seconds'] for r in entry['routes'].values()):.1f}s", flush=True)
            raw.append(entries)
            records.append(encode(dict(digits=digits, offset=offset, partitions=entries)))
    # Parse the earlier decimal-string evidence only inside the higher
    # precision context. Compare *each* previously saved complex term too.
    with mp.workdps(high_digits):
        earlier = previous["records"][-1]
        lookup = {tuple(p["holes"]): decode_complex(p["amplitude"]) for p in earlier["pieces"]}
        earlier_k = kinematic_factor(beta, earlier["offset"])
        now = raw[-1]
        coefficients = {p["holes"]: p["routes"]["prl"]["reduced_coefficient"] for p in now}
        physical = boundary_factor(beta)
        target = mp.mpf(earlier["vertex"]["amplitude"])
        summaries = {}
        for route in ("jmp", "prl"):
            values = [p["routes"][route]["boundary_amplitude"] for p in now]
            total = mp.fsum(values)
            summaries[route] = dict(
                boundary_total=total, vertex_complex_relative_difference=abs(total-target)/abs(target),
                sum_absolute_terms=mp.fsum(abs(z) for z in values),
                cancellation_condition=mp.fsum(abs(z) for z in values)/abs(total),
                half_sum=mp.fsum(p["routes"][route]["boundary_amplitude"]
                                 for p in now if 0 in p["holes"]))
        metrics = dict(
            worst_route_difference=max(p["reduced_scaled_difference"] for r in raw for p in r),
            worst_refined_route_difference=max(p["reduced_scaled_difference"] for p in now),
            worst_joint_refinement=max(
                scaled_difference(a["routes"][r]["reduced_coefficient"],
                                  b["routes"][r]["reduced_coefficient"])
                for a, b in zip(*raw) for r in ("jmp", "prl")),
            worst_previous_boundary_term_difference=max(
                scaled_difference(earlier_k*coefficients[h], a) for h, a in lookup.items()),
            worst_sine_identity=max(p["sine_identity_scaled_error"] for r in raw for p in r),
            worst_reconstruction=max(v["reconstruction_scaled_error"]
                                     for r in raw for p in r for v in p["routes"].values()),
            worst_complement_conjugacy=max(
                scaled_difference(c, mp.conj(coefficients[tuple(j for j in range(len(beta)) if j not in h)]))
                for h, c in coefficients.items()),
            worst_null_determinant=max(v for r in raw for p in r
                                      for route in p["routes"].values()
                                      for v in route["null_determinants_scaled"]),
            worst_jmp_reported_quad_error=max(p["routes"]["jmp"]["max_quad_error"]
                                              for r in raw for p in r),
        )
        # Numerical screens, not certified enclosures. The outer integrals
        # and the source of another implementation's discrepancy remain open.
        low_tolerance = mp.power(10, -(low_digits-8))
        high_tolerance = mp.power(10, -(high_digits-10))
        checks = dict(
            complete_partitions=all(len(r) == math.comb(len(beta), len(beta)//2) for r in raw)
                and set(coefficients) == set(lookup),
            routes_agree=metrics["worst_route_difference"] < low_tolerance,
            refined_routes_agree=metrics["worst_refined_route_difference"] < high_tolerance,
            joint_refinement=metrics["worst_joint_refinement"] < low_tolerance,
            previous_complex_terms=metrics["worst_previous_boundary_term_difference"] < high_tolerance,
            sine_identity=metrics["worst_sine_identity"] < low_tolerance,
            reconstruction=metrics["worst_reconstruction"] < low_tolerance,
            complement_conjugacy=metrics["worst_complement_conjugacy"] < high_tolerance,
            null_determinants=metrics["worst_null_determinant"] < low_tolerance,
            vertex_agreement=all(v["vertex_complex_relative_difference"] < high_tolerance
                                 for v in summaries.values()))
        result = encode(dict(
            fixture=name, beta=beta, spinons=len(beta), records=records,
            boundary_sine_squared=physical, summaries=summaries, metrics=metrics,
            vertex_target=earlier["vertex"], previous_offset=earlier["offset"],
            common_beta_measure=1/(math.factorial(len(beta))*(2*mp.pi**2)**len(beta)),
            checks=checks, passed=all(checks.values())))
        print(f"{name}: {'PASS' if result['passed'] else 'FAIL'} {checks}", flush=True)
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", nargs="+", choices=FIXTURES,
                        default=["two", "four", "six-bulk"])
    parser.add_argument("--low-digits", type=int, default=30)
    parser.add_argument("--high-digits", type=int, default=42)
    parser.add_argument("--high-maxdegree", type=int, default=9)
    parser.add_argument("--checkpoint", type=Path,
                        help="optional passed audit_xxx_qtm_vertex.py JSON from these same sources; "
                             "otherwise generate fresh boundary references")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output must not exist")
    if args.low_digits < 28 or args.high_digits < args.low_digits+10:
        parser.error("require low-digits>=28 and at least ten extra high digits")
    if args.high_maxdegree < 7:
        parser.error("high-maxdegree must be at least seven")
    source = hashes()
    wall, cpu = time.perf_counter(), time.process_time()
    checkpoint_hash = None
    if args.checkpoint is not None:
        try:
            previous = verify_checkpoint(args.checkpoint)
        except (OSError, ValueError) as error:
            parser.error(str(error))
        checkpoint_hash = hashlib.sha256(args.checkpoint.read_bytes()).hexdigest()
        fixtures = {f["fixture"]: f for f in previous["fixtures"]}
        if any(name not in fixtures for name in args.fixture):
            parser.error("each selected fixture must be present in the accepted checkpoint")
    else:
        fixtures = {name: fresh_reference(name, args.high_digits+10) for name in args.fixture}
    results = [run_fixture(name, fixtures[name], args.low_digits, args.high_digits, args.high_maxdegree)
               for name in args.fixture]
    if args.checkpoint is not None:
        verify_checkpoint(args.checkpoint)
        if checkpoint_hash != hashlib.sha256(args.checkpoint.read_bytes()).hexdigest():
            raise RuntimeError("input checkpoint changed during audit")
    if source != hashes():
        raise RuntimeError("sources changed during audit; refusing a mixed checkpoint")
    result = dict(
        schema="xxx-qtm-complex-partitions-v1",
        scope="matched complex symmetric Pauli-zz amplitudes; not an outer integral",
        index_base=0, input_rule="u_j=beta_j/pi+i*(eta-1); v_j=beta_j/pi+i*eta",
        comparison_rule="compare full complex A at matching eta, or reduced B; never abs(A) or real(A)",
        boundary_rule="A_H(0)=sin(sum_j atan(exp(beta_j)))^2 * B_H",
        python=platform.python_version(), mpmath=mp.__version__, platform=platform.platform(),
        source_sha256=source, invocation=sys.argv[1:], processes=1,
        reference_input=dict(mode="checkpoint" if args.checkpoint else "fresh",
                             checkpoint_sha256=checkpoint_hash,
                             fresh_decimal_digits=None if args.checkpoint else args.high_digits+10),
        wall_seconds=time.perf_counter()-wall, cpu_seconds=time.process_time()-cpu,
        fixtures=results, passed=all(f["passed"] for f in results))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as f:
        json.dump(result, f, indent=2)
        f.write("\n")
    print(f"wrote {args.output}: {'PASS' if result['passed'] else 'FAIL'}", flush=True)
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
