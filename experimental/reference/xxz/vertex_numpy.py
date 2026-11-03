# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Optional binary64 q-product/residue reference; explicit finite truncations.

The cancelled recurrence is also used by C++; the separate mpmath reference
evaluates literal shifted residues. No outer-quadrature convergence is implied.
"""
import numpy as np


def qproduct(z, base, terms):
    out = np.ones_like(z, dtype=complex)
    for h in range(terms):
        out *= 1 - z * base**h
    return out

def double_product(z, base, terms):
    out = np.ones_like(z, dtype=complex)
    for h in range(terms):
        out *= (1 - z * base**h)**(h + 1)
    return out

def residue_ratio(u, ell, j, m, q):
    """Product-shift ratio with all negative powers of v cancelled.

    The q^2-product ratio is 1/v^2, and the theta/vacuum factor supplies
    q^2*prod(u)/v^2. The remaining powers cancel against the pole factor
    and the three spectator denominators. The limiting ratio is q^6.
    """
    r = q**4
    v = u[:, j]*r**m
    ratio = -q*q*r*np.prod(u, axis=1)/(u[:, j]*(1-r**(m+1)))
    for k in range(4):
        ratio *= 1-q**-2*v/u[:, k]
        if k != j:
            ratio /= r*v-u[:, k]
        if k != ell:
            divisor = q*q*u[:, k] if k < ell else u[:, k]
            ratio *= (1-r*v/divisor)/(1-v/divisor)
    return ratio

def vertex_amplitudes(angles, q, vacuum, poles=4, products=12, recurrence=False, diagnostics=False):
    """Vectorized direct finite-q residues, restricted to the small-|q| audit.

    angles has shape (points, spinons), xi=i*exp(-i*angle). The optional
    stable residue recurrence extends the research range to |q|<=.975.
    The product and
    residue truncations are explicit controls, not certified error bounds.
    """
    angles = np.asarray(angles, dtype=float)
    if (angles.ndim != 2 or angles.shape[1] not in (2, 4)
            or not np.all(np.isfinite(angles)) or not (-.975 if recurrence else -.2) <= q < 0
            or vacuum not in (0, 1) or poles < 1 or products < 4):
        raise ValueError('invalid vertex controls: direct |q|<=.2, recurrence |q|<=.975, vacuum 0/1, poles>=1, products>=4')
    xi = 1j * np.exp(-1j * angles)
    u = -xi * xi
    n = angles.shape[1]
    r = q**4
    qp = lambda z, base=r: qproduct(z, base, products)
    dp = lambda z: double_product(z, r, products)
    theta = lambda z: qp(q**8, q**8) * qp(z, q**8) * qp(q**8/z, q**8)
    common = np.ones(len(angles), dtype=complex)
    separation = np.full(len(angles), 1-q*q)
    common *= (-q)**(1-vacuum*(n//2)) * (1-q**-2)**(n//2-1)
    common *= qp(q*q) * qp(r)**(n-1) * (dp(r)/dp(q**6))**n
    for j in range(n):
        common /= qp(-q*u[:, j]) * qp(-q**3/u[:, j])
        for k in range(j+1, n):
            a = u[:, k]/u[:, j]
            if np.any(abs(a-1) < 1e-12):
                raise ValueError('coincident spectral nodes need a separate limiting formula')
            separation = np.minimum(separation, abs(a-1))
            common *= dp(r*a) * dp(1/a) / dp(q**6*a) / dp(q*q/a)
    if n == 2:
        monomial = np.prod(xi**np.arange(1-vacuum, 3-vacuum), axis=1)
        result = (common * monomial * theta(-q**(4*vacuum)/np.prod(u, axis=1)))[:, None]
        indicator = 128*np.finfo(float).eps*(products+poles+1)/separation[:, None]*abs(result)
        return (result, indicator) if diagnostics else result
    amplitudes = []
    indicators = []
    for ell in range(4):
        integral = np.zeros(len(angles), dtype=complex)
        absolute = np.zeros(len(angles))
        for j in range(4):
            start = int(j > ell)
            for m in range(start, min(poles, start+1) if recurrence else poles):
                v = u[:, j]*r**m
                finite_product = np.prod([1-r**h for h in range(1, m+1)])
                residue = np.full(len(angles), (-1)**m*r**(m*(m+1)//2)/(finite_product*qp(r)), dtype=complex)
                for k in range(ell):
                    residue *= q/v - 1/(q*u[:, k])
                for k in range(ell+1, 4):
                    residue *= 1/v - 1/u[:, k]
                for k in range(4):
                    if k != j:
                        residue /= qp(u[:, k]/v)
                    residue /= qp(q**-2*v/u[:, k])
                residue *= qp(-v/q, q*q) * qp(-q**3/v, q*q) * (-v)**vacuum
                residue *= theta(-q**(-2+4*vacuum)*v*v/np.prod(u, axis=1))
                integral += residue
                absolute += abs(residue)
                if recurrence:
                    for k in range(start, poles-1):
                        residue *= residue_ratio(u, ell, j, k, q)
                        integral += residue
                        absolute += abs(residue)
        powers = np.arange(1-vacuum, 5-vacuum)
        powers[ell] -= 1
        amplitudes.append(common*np.prod(xi**powers, axis=1)*integral)
        # Conservative cancellation/roundoff screen, not a certified bound or
        # an estimate of omitted products/poles. Include close-node sensitivity.
        indicators.append(abs(common)*absolute*128*np.finfo(float).eps*(products+poles+1)/separation)
    result = np.stack(amplitudes, axis=1)
    if not np.all(np.isfinite(result)):
        raise ArithmeticError('finite-q amplitude overflow or invalid arithmetic')
    return (result, np.stack(indicators, axis=1)) if diagnostics else result
