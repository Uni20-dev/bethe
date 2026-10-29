#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Independent high-precision Lee-Yang TBA reference on a composite Gauss mesh.

The production code uses a half-line trapezoid and reduced real source formulas.
This oracle uses composite Gauss-Legendre integration and the complex scattering
amplitude directly. Refining arithmetic, quadrature order and cutoff remains a
separate obligation: the requested decimal precision is not an error certificate.
Requires mpmath; it is not a production dependency.
"""
import argparse
import json
import time

import mpmath as mp


def calculate(r_text="30", dps=96, order=64, cutoff_text="4"):
    with mp.workdps(dps):
        r, cutoff = mp.mpf(r_text), mp.mpf(cutoff_text)
        if not 5 <= r <= 30 or cutoff <= 0 or order < 4 or dps < 40:
            raise ValueError("requires 5<=r<=30, positive cutoff, order>=4 and dps>=40")
        panels = int(mp.ceil(2 * cutoff))
        width = cutoff / panels
        nodes, weights = mp.gauss_quadrature(order, "legendre")
        theta = [(i + (nodes[j] + 1) / 2) * width for i in range(panels) for j in range(order)]
        weight = [weights[j] * width / 2 for _ in range(panels) for j in range(order)]
        sh, ch = [mp.sinh(x) for x in theta], [mp.cosh(x) for x in theta]
        drive = [r * x for x in ch]
        a = mp.sqrt(3) / 2
        normalization = a / mp.pi
        # Hyperbolic addition formulas avoid O(N^2) transcendental evaluations.
        operator = []
        for si, ci in zip(sh, ch):
            row = []
            for sj, cj, w in zip(sh, ch, weight):
                sm, sp = si * cj - ci * sj, si * cj + ci * sj
                cm, cp = ci * cj - si * sj, ci * cj + si * sj
                row.append(w * normalization * (cm / (sm * sm + a * a) + cp / (sp * sp + a * a)))
            operator.append(row)
        target = mp.power(10, -(dps - 15))
        counts = {"fixed_point_updates": 0, "source_queries": 0}

        def scattering(z):
            s = mp.sinh(z)
            return (s + 1j * a) / (s - 1j * a)

        def solve(source):
            correction = [mp.mpf(0)] * len(theta)
            for iteration in range(40):
                logarithm = [mp.log1p(mp.exp(-e - b - u)) for e, b, u in zip(drive, source, correction)]
                updated = [mp.fdot(row, logarithm) for row in operator]
                residual = max(abs(x - y) for x, y in zip(updated, correction))
                if residual < target:
                    return logarithm, residual
                correction = updated
                counts["fixed_point_updates"] += 1
            raise RuntimeError("reference fixed-point iteration did not converge")

        vacuum_logs, vacuum_residual = solve([mp.mpf(0)] * len(theta))
        y0 = -mp.fdot([w * e / mp.pi for w, e in zip(weight, drive)], vacuum_logs)

        def evaluate(log_d):
            counts["source_queries"] += 1
            d = mp.exp(log_d)
            beta = mp.pi / 6 + d
            source = [-2 * mp.log(abs(scattering(x + 1j * beta))) for x in theta]
            logarithm, residual = solve(source)
            continued = []
            for x, w in zip(theta, weight):
                z = 1j * beta - x
                continued.append(2 * w * normalization * (mp.cosh(z) / (mp.sinh(z)**2 + a*a)).real)
            # Positive ratio at the pole-connected root, evaluated directly
            # with high working precision rather than the production log-d form.
            log_s = mp.log((mp.sin(2*beta) + a) / (mp.sin(2*beta) - a))
            quantization = r * mp.cos(beta) - log_s + mp.fdot(continued, logarithm)
            y = 2 * r * mp.sin(beta) - mp.fdot([w * e / mp.pi for w, e in zip(weight, drive)], logarithm)
            return quantization, y, d, residual

        seed = mp.log(mp.sqrt(3)) - r * mp.sqrt(3) / 2
        root = mp.findroot(lambda x: evaluate(x)[0], (seed - 1, seed + 1),
                           tol=mp.power(10, -(dps-25)), maxsteps=40)
        quantization, y1, displacement, nonlinear_residual = evaluate(root)
        if abs(quantization) > mp.power(10, -(dps-30)):
            raise RuntimeError("reference quantization residual exceeds target")
        values = {"Y0": y0, "Y1": y1, "displacement": displacement,
                  "beta": mp.pi/6 + displacement, "gap": (y1-y0)/r,
                  "vacuum_residual": vacuum_residual, "excited_residual": nonlinear_residual,
                  "quantization_residual": abs(quantization)}
        return {"r": r_text, "dps": dps, "order": order, "cutoff": cutoff_text,
                "panels": panels, "nodes": len(theta), "mpmath_version": mp.__version__, **counts,
                **{name: mp.nstr(value, dps-10) for name, value in values.items()}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r", default="30")
    parser.add_argument("--dps", type=int, default=96)
    parser.add_argument("--order", type=int, default=64)
    parser.add_argument("--cutoff", default="4")
    args = parser.parse_args()
    start = time.monotonic()
    result = calculate(args.r, args.dps, args.order, args.cutoff)
    result["seconds"] = time.monotonic() - start
    print(json.dumps(result, indent=2))
