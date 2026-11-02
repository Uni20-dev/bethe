#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Pointwise particle-hole / vertex comparison, with no outer integration.

The QTM amplitude is symmetrized over ALL labelled hole/particle assignments.
The common contour offset tends to zero from above. The real spinon density
and its analytically derived measure conversion are independent inputs:
no fitted scale or observed integrated fraction appears here.
"""
from itertools import combinations
import math

import mpmath as mp

import _reference_support  # Locate the same checkout's public vertex controls.
import reference_xxx_qtm as qtm
from reference_xxx_six_spinon_vertex import validate
from reference_xxx_vertex_density import density_at_precision


def assignments(n):
    """Ordered (holes, particles) bipartitions; a complement is a distinct term."""
    if isinstance(n, bool) or not isinstance(n, int) or n not in (2, 4, 6):
        raise ValueError("require two, four or six spinons")
    return [(h, tuple(j for j in range(n) if j not in h))
            for h in combinations(range(n), n//2)]


def shifted_points(beta, offset):
    beta = validate(beta, 1)
    offset = mp.mpf(offset)
    if not mp.isfinite(offset) or not 0 < offset < 1:
        raise ValueError("require a real common contour offset strictly between zero and one")
    # Do not silently round a near-boundary hole onto Im(u)=-1.
    if 1-offset == 1 or offset < 64*mp.eps:
        raise ValueError("contour offset needs more working precision")
    return beta, tuple(b/mp.pi+mp.j*offset for b in beta)


def vertex_conversion(beta):
    """Multiply the Szz momentum density W_n to get the Pauli QTM amplitude.

    W_n includes 1/n!, both conservation deltas are still unintegrated.
    Szz -> Pauli gives 4, inverse Fourier measures give (2*pi)^-2,
    dp/d(beta)=sech(beta), dt=d(beta)/pi, and the QTM common measure is
    1/[n!*(2*pi)^n] AFTER summing all assignments. See the derivation note.
    """
    beta = validate(beta, 1)
    n = len(beta)
    return math.factorial(n)*(2*mp.pi**2)**n/mp.pi**2/mp.fprod(mp.cosh(b) for b in beta)


def partition_sum(beta, offset, *, route="prl", progress=None):
    beta, t = shifted_points(beta, offset)
    if route not in ("jmp", "prl"):
        raise ValueError("route must be jmp or prl")
    evaluate = qtm.prl if route == "prl" else qtm.jmp
    pieces, splits = [], assignments(len(beta))
    for holes, particles in splits:
        result = evaluate([t[j]-mp.j for j in holes], [t[j] for j in particles])
        pieces.append(dict(holes=holes, particles=particles, amplitude=result["amplitude"]))
        if progress is not None:
            progress(len(pieces), len(splits))
    total = mp.fsum(p["amplitude"] for p in pieces)
    absolute_sum = mp.fsum(abs(p["amplitude"]) for p in pieces)
    # At the real boundary complements become conjugates. Retain all terms;
    # this is only a diagnostic, never a shortcut that discards half of them.
    one_per_pair = mp.fsum(p["amplitude"] for p in pieces if 0 in p["holes"])
    lookup = {p["holes"]: p["amplitude"] for p in pieces}
    conjugacy_defect = max(abs(p["amplitude"]-mp.conj(lookup[p["particles"]]))
                          for p in pieces)/max(abs(p["amplitude"]) for p in pieces)
    density_measure = 1/(math.factorial(len(beta))*(2*mp.pi**2)**len(beta))
    return dict(total=total, pieces=pieces, offset=mp.mpf(offset),
                cancellation_condition=absolute_sum/abs(total) if total else mp.inf,
                one_per_complement_pair=one_per_pair,
                complement_conjugacy_scaled_defect=conjugacy_defect,
                pauli_rapidity_density=total*density_measure,
                density_measure=density_measure, route=route)


def vertex_target(beta, *, terms=64, order=20):
    beta = validate(beta, terms)
    result = density_at_precision(beta, terms, order)
    factor = vertex_conversion(beta)
    return dict(amplitude=factor*result.value,
                contour_tail_error=factor*result.contour_tail_error,
                arithmetic_indicator=factor*result.arithmetic_indicator,
                momentum_density=result.value, conversion_factor=factor,
                terms=terms, order=order, digits=mp.mp.dps)


def richardson_real(first, second):
    """Cancel an observed O(offset^2) real-part correction, not a rigorous bound."""
    e, f = first["offset"], second["offset"]
    if not 0 < f < e:
        raise ValueError("require two strictly decreasing positive offsets")
    return (e**2*second["total"].real-f**2*first["total"].real)/(e**2-f**2)
