#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Independent SciPy oracle for the periodic scaling Lee-Yang ground state.

Development only, not a runtime/build dependency. Equations (133), (212),
(218) of https://arxiv.org/html/1412.8494. Outputs Y=L*E_C at r=mL;
E_C is bulk-subtracted. No excited-state or absolute bulk-energy claim.
"""
import argparse
import json
import math

import numpy as np
from scipy.special import k1, roots_legendre


def kernel(x):
    """Positive -phi/(2*pi), evaluated without large hyperbolic functions."""
    t = np.exp(-np.abs(x))
    t2 = t*t
    return math.sqrt(3)/math.pi*t*(1+t2)/(1+t2+t2*t2)


def grid(cutoff, order, panels):
    x, w = roots_legendre(order)
    width = cutoff/panels
    return ((np.arange(panels)[:, None]+(x+1)/2)*width).ravel(), np.tile(w*width/2, panels)


def solve(r, cutoff, order=16, panels=None, tolerance=2e-14, max_iterations=200):
    if not math.isfinite(r) or not 1e-6 <= r <= 50:
        raise ValueError("reference oracle supports 1e-6<=r=mL<=50")
    if not math.isfinite(cutoff) or not 0 < cutoff <= 30:
        raise ValueError("cutoff must lie in (0,30]")
    if not 2 <= order <= 64 or not 0 < tolerance < 1 or max_iterations < 0:
        raise ValueError("invalid quadrature or iteration controls")
    panels = math.ceil(cutoff/.75) if panels is None else panels
    if panels < 1 or panels*order > 2048:
        raise ValueError("reference grid must have 1..2048 nodes")
    theta, weight = grid(cutoff, order, panels)
    drive = r*np.cosh(theta)
    operator = (kernel(theta[:, None]-theta[None, :])+kernel(theta[:, None]+theta[None, :]))*weight
    # Solve for the bounded correction, avoiding subtraction of large drives.
    correction = np.zeros_like(theta)
    for iteration in range(max_iterations+1):
        logarithm = np.logaddexp(0, -drive-correction)
        target = operator@logarithm
        residual = float(np.max(np.abs(correction-target)))
        if residual <= tolerance:
            break
        if iteration == max_iterations:
            raise RuntimeError("reference nonlinear iteration budget exhausted")
        correction = target
    y = -float(np.dot(weight*drive, logarithm))/math.pi
    center = r+float(np.dot(2*kernel(theta)*weight, logarithm))
    return {"r": r, "Y": y, "effective_central_charge": -6*y/math.pi,
            "pseudoenergy_center": center, "residual": residual,
            "iterations": iteration, "cutoff": cutoff, "order": order,
            "panels": panels, "nodes": len(theta)}


def reference(r):
    if not math.isfinite(r) or not 1e-6 <= r <= 50:
        raise ValueError("reference oracle supports 1e-6<=r=mL<=50")
    cutoff = math.acosh(max(2, 60/r))
    a = solve(r, cutoff, order=12)
    b = solve(r, cutoff, order=24)
    c = solve(r, cutoff+2, order=24)
    mesh = abs(a["Y"]-b["Y"])
    tail = abs(b["Y"]-c["Y"])
    if max(mesh, tail) > 2e-12:
        raise RuntimeError("reference quadrature/cutoff agreement failed")
    c.update(mesh_difference=mesh, cutoff_difference=tail,
             leading_large_r_Y=-r*k1(r)/math.pi)
    return c


def check():
    x, w = grid(30, 24, 40)
    assert abs(float(np.dot(2*kernel(x), w))-1) < 3e-13
    assert np.all(kernel(x) > 0) and np.array_equal(kernel(x), kernel(-x))
    uv = reference(1e-5)
    assert abs(uv["effective_central_charge"]-.4) < 1e-8
    assert abs(uv["pseudoenergy_center"]-math.log((1+math.sqrt(5))/2)) < 1e-5
    ir = reference(20)
    assert abs(ir["Y"]/ir["leading_large_r_Y"]-1) < 1e-8
    for r in (.001, .01, .1, .5, 1, 2, 5, 10):
        row = reference(r)
        assert row["leading_large_r_Y"] < row["Y"] < 0
        assert 0 < row["effective_central_charge"] < .4
    try:
        solve(1, 5, max_iterations=0)
    except RuntimeError:
        pass
    else:
        raise AssertionError("zero iteration budget must fail for r=1")
    for invalid in (0, -1, math.nan, math.inf, 100):
        try:
            reference(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid oracle input accepted")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r", type=float, nargs="+", default=[1e-5, .001, .01, .1, .5, 1, 2, 5, 10, 20])
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        check()
    rows = [reference(r) for r in args.r]
    print(json.dumps({"source": "https://arxiv.org/html/1412.8494",
                      "convention": "periodic ground state, bulk-subtracted Y=L*E_C, r=mL, velocity=1",
                      "c": "-22/5", "h_min": "-1/5", "ultraviolet_c_eff": "2/5",
                      "rows": rows}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
