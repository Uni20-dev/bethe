#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Development-only Lee-Yang spin-zero one-particle TBA oracle, 5<=mL<=30.

Dorey--Tateo, hep-th/9607167, Eqs. (2.3)--(2.7). This regular infrared
branch is not a continuation through the source singularity near mL=2.53.
NumPy/SciPy are not production dependencies. Energies are bulk-subtracted.
"""
import argparse
import json
import math

import numpy as np
from scipy.optimize import brentq

from reference_lee_yang import grid, kernel, reference as vacuum_reference


def scattering(z):
    s = np.sinh(z)
    a = 1j*math.sqrt(3)/2
    return (s+a)/(s-a)


def complex_kernel(z):
    return math.sqrt(3)/(2*math.pi)*np.cosh(z)/(np.sinh(z)**2+.75)


def solve(r, cutoff, order=16, max_iterations=200):
    if not math.isfinite(r) or not 5 <= r <= 30:
        raise ValueError("excited reference oracle supports 5<=r=mL<=30")
    if not math.isfinite(cutoff) or not 0 < cutoff <= 30:
        raise ValueError("cutoff must lie in (0,30]")
    if not 2 <= order <= 64 or max_iterations < 0:
        raise ValueError("invalid quadrature or iteration controls")
    panels = math.ceil(cutoff/.75)
    theta, weight = grid(cutoff, order, panels)
    drive = r*np.cosh(theta)
    operator = (kernel(theta[:, None]-theta[None, :])
                + kernel(theta[:, None]+theta[None, :]))*weight
    total_iterations = 0

    def evaluate(log_displacement):
        nonlocal total_iterations
        displacement = math.exp(log_displacement)
        beta = math.pi/6+displacement
        source = -2*np.log(np.abs(scattering(theta+1j*beta)))
        correction = np.zeros_like(theta)
        for iteration in range(max_iterations+1):
            logarithm = np.logaddexp(0, -drive-source-correction)
            target = operator@logarithm
            residual = float(np.max(np.abs(target-correction)))
            if residual <= 2e-14:
                break
            if iteration == max_iterations:
                raise RuntimeError("reference nonlinear iteration budget exhausted")
            correction = target
        total_iterations += iteration
        # sin(2*beta)-sqrt(3)/2 without cancellation at beta -> pi/6.
        denominator = 2*math.cos(math.pi/3+displacement)*math.sin(displacement)
        log_s = math.log(denominator+math.sqrt(3))-math.log(denominator)
        convolution = float(np.dot(2*complex_kernel(1j*beta-theta).real*weight,
                                   logarithm))
        quantization = r*math.cos(beta)-log_s+convolution
        energy_over_mass = 2*math.sin(beta)-float(np.dot(weight*np.cosh(theta), logarithm))/math.pi
        return quantization, energy_over_mass, residual, beta

    # Select the pole-connected root, safely short of the singularity beta=pi/3.
    lower, upper = -r-10, math.log(math.pi/12)
    if not evaluate(lower)[0] < 0 < evaluate(upper)[0]:
        raise RuntimeError("regular one-particle root is not bracketed")
    root = brentq(lambda x: evaluate(x)[0], lower, upper, xtol=1e-13, rtol=1e-14)
    quantization, energy, residual, beta = evaluate(root)
    if abs(quantization) > 2e-12:
        raise RuntimeError("reference quantization residual too large")
    return {"r": r, "Y1": r*energy, "energy_over_mass": energy,
            "beta": beta, "pole_displacement": math.exp(root),
            "quantization_residual": abs(quantization), "nonlinear_residual": residual,
            "iterations": total_iterations, "cutoff": cutoff, "order": order,
            "panels": panels, "nodes": len(theta)}


def reference(r):
    if not math.isfinite(r) or not 5 <= r <= 30:
        raise ValueError("excited reference oracle supports 5<=r=mL<=30")
    cutoff = math.acosh(max(2, 60/r))
    a = solve(r, cutoff, order=12)
    b = solve(r, cutoff, order=24)
    c = solve(r, cutoff+2, order=24)
    mesh, tail = abs(a["Y1"]-b["Y1"]), abs(b["Y1"]-c["Y1"])
    root_mesh = abs(a["pole_displacement"]-b["pole_displacement"])
    root_tail = abs(b["pole_displacement"]-c["pole_displacement"])
    if max(mesh, tail, root_mesh, root_tail) > 2e-12:
        raise RuntimeError("reference quadrature/cutoff agreement failed")
    vacuum = vacuum_reference(r)
    c.update(Y0=vacuum["Y"], gap_over_mass=(c["Y1"]-vacuum["Y"])/r,
             mesh_difference=mesh, cutoff_difference=tail,
             source_mesh_difference=root_mesh, source_cutoff_difference=root_tail)
    return c


def check():
    x = np.array([0., .1, 1., 5.])
    assert np.max(np.abs(complex_kernel(x)-kernel(x))) < 3e-16
    beta = .6
    ratio = scattering(x-1j*beta)/scattering(x+1j*beta)
    assert np.max(np.abs(ratio.imag)) < 1e-15 and np.all(ratio.real > 0)
    assert np.max(np.abs(np.log(ratio.real)+2*np.log(np.abs(scattering(x+1j*beta))))) < 2e-15
    for r in (5, 6, 8, 10, 15, 20, 30):
        row = reference(r)
        assert math.pi/6 <= row["beta"] < math.pi/4
        assert row["gap_over_mass"] > 1
    ir = reference(30)
    expected = math.sqrt(3)*math.exp(-math.sqrt(3)*15)
    assert abs(ir["pole_displacement"]/expected-1) < 1e-8
    # The leading mu-term in E1/m is +3 exp(-sqrt(3)*r/2), not sqrt(3).
    assert abs((ir["energy_over_mass"]-1)/(3*math.exp(-math.sqrt(3)*15))-1) < .02
    try:
        solve(5, 5, max_iterations=0)
    except RuntimeError:
        pass
    else:
        raise AssertionError("zero iteration budget accepted")
    for invalid in (0, -1, 2.53, 4.99, 30.01, math.nan, math.inf):
        try:
            reference(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError("out-of-domain input accepted")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r", type=float, nargs="+", default=[5, 6, 8, 10, 15, 20, 30])
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        check()
    print(json.dumps({"source": "https://arxiv.org/html/hep-th/9607167",
                      "convention": "periodic spin-zero one-particle, bulk-subtracted Y1=L*E1_C; r=mL",
                      "scope": "regular infrared branch only, 5<=r<=30; no ultraviolet continuation",
                      "rows": [reference(r) for r in args.r]}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
