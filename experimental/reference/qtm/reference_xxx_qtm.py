#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Pointwise XXX QTM reference: independent JMP quadrature and PRL series.

Pauli zz amplitude, with the common 1/(ell!^2*(2*pi)^(2*ell)) measure exposed
separately. Caller controls mpmath precision; no outer integration is done.
This is a developer oracle, not a production solver or certified enclosure.
The published particle/hole contours use separate infinitesimal regulators
around Im(z)=+/-1/2. The wider input strips allow pointwise analytic checks
but do not certify a contour deformation.
See docs/experimental-xxx-qtm.md for conventions and limiting algebra.
"""
import math
from functools import lru_cache

import mpmath as mp


def points(holes, particles):
    """Validate the nonsingular test domain, not an outer integration contour."""
    u, v = tuple(map(mp.mpc, holes)), tuple(map(mp.mpc, particles))
    if not u or len(u) != len(v) or any(not mp.isfinite(z) for z in u+v):
        raise ValueError("require equally many finite holes and particles")
    if len(set(u)) != len(u) or len(set(v)) != len(v):
        raise ValueError("coincident same-species rapidities need a separate limit")
    if any(not -1 < z.imag < 0 for z in u) or any(not 0 < z.imag < 1 for z in v):
        raise ValueError("reference contours require holes in (-1,0), particles in (0,1)")
    if any(abs(x-y+mp.j) < 32*mp.eps for x in u for y in v):
        raise ValueError("coincident real hole/particle with separation i needs a pole limit")
    return u, v


def log_psi(z):
    """JMP (276), also the ordinary-Gamma limit of PRL (11)."""
    z = mp.mpc(z)
    a = -mp.j*z/2
    return (mp.loggamma(mp.mpf(".5")+a)+mp.loggamma(1+a)
            + 2*mp.log(mp.barnesg(1+a))-2*mp.log(mp.barnesg(mp.mpf(".5")+a)))


def momentum(z):
    return mp.pi/2-mp.j*mp.log(mp.cosh(mp.pi*(z+mp.j/2)/2)
                              /mp.cosh(mp.pi*(z-mp.j/2)/2))


def sine_squared(u, v):
    return mp.sin((mp.fsum(momentum(z) for z in v)
                   - mp.fsum(momentum(z) for z in u))/2)**2


def scattering(z):
    b = mp.j*z/2
    return mp.j*mp.exp(mp.loggamma(1+b)+mp.loggamma(mp.mpf(".5")-b)
                       -mp.loggamma(1-b)-mp.loggamma(mp.mpf(".5")+b))


def measure(ell):
    return 1/(math.factorial(ell)**2*(2*mp.pi)**(2*ell))


def jmp(holes, particles, *, maxdegree=7, split=True):
    """JMP (279)--(282), symmetrized over its distinguished first columns.

    Equation (282) is an integration-symmetry shortcut, NOT the symmetric
    pointwise amplitude. Summing over all choices of distinguished columns
    restores the amplitude convention in (243) / PRL (10).
    """
    u, v = points(holes, particles)
    ell = len(u)

    @lru_cache(maxsize=None)
    def phi_log(t, sigma):
        def term(w):
            a = sigma*mp.j*(t-w)/2
            return mp.loggamma(mp.mpf(".5")-a)+mp.loggamma(1+a)
        return mp.fsum(term(w) for w in u)-mp.fsum(term(w) for w in v)

    g = lambda a, b: mp.pi/mp.sinh(mp.pi*(b-a))
    e = lambda z: 1/z-1/(z-mp.j)
    # Split near every real part as well as zero. Increasing maxdegree and
    # disabling these splits give separate quadrature convergence controls.
    breaks = sorted(set([mp.mpf(0)]+[z.real for z in u+v])) if split else [mp.mpf(0)]
    intervals = [-mp.inf]+breaks+[mp.inf]
    errors = []

    def integrate(f):
        value, error = mp.quad(f, intervals, error=True, maxdegree=maxdegree)
        errors.append(error)
        return value

    a, b = mp.matrix(ell), mp.matrix(ell)
    amod, bmod = mp.matrix(ell), mp.matrix(ell)
    for j in range(ell):
        for k in range(ell):
            x, y = u[j], v[k]
            ratio = lambda t: mp.exp(phi_log(y, -1)-phi_log(t, 1))
            amod[j, k] = integrate(lambda t: ratio(t)*g(x, t)/(y-t))/mp.pi
            a[j, k] = g(x, y)-integrate(lambda t: ratio(t)*g(x, t)*e(t-y))/(2*mp.pi*mp.j)
            y, x = v[j], u[k]
            ratio_bar = lambda t: mp.exp(phi_log(t, -1)-phi_log(x, 1))
            bmod[j, k] = integrate(lambda t: ratio_bar(t)*g(y, t)/(x-t))/mp.pi
            b[j, k] = g(y, x)-integrate(lambda t: ratio_bar(t)*g(y, t)*e(x-t))/(2*mp.pi*mp.j)
    da, db = determinant_derivative(a, amod), determinant_derivative(b, bmod)
    signed = [(z, 1) for z in u]+[(z, -1) for z in v]
    pair = mp.exp(mp.fsum(s*t*log_psi(x-y) for x, s in signed for y, t in signed))
    scalar = (-1)**(ell-1)*sine_squared(u, v)*pair
    amplitude = scalar*da*db
    first_a, first_b = a.copy(), b.copy()
    for j in range(ell):
        first_a[j, 0], first_b[j, 0] = amod[j, 0], bmod[j, 0]
    literal = scalar*ell**2*mp.det(first_a)*mp.det(first_b)
    return dict(amplitude=amplitude, density=amplitude*measure(ell),
                matrices=(a, b), determinants=(da, db), pair_factor=pair,
                modified_columns=(amod, bmod),
                literal_first_column_amplitude=literal,
                null_determinants=(mp.det(a), mp.det(b)),
                max_quad_error=max(errors), quad_errors=errors,
                maxdegree=maxdegree, split=split, phi_cache_size=phi_log.cache_info().currsize)


def hyper_parameters(u, v, j, i=None):
    """PRL (6): q^(c)*exp(2i*gamma*d) -> ordinary parameter c/4-id/2."""
    half = mp.mpf(".5")
    if i is None:
        vv = [-mp.j*(v[j]-v[r])/2 for r in range(len(v)) if r != j]
        vu = [-mp.j*(v[j]-x)/2 for x in u]
        return [-half]+[half+d for d in vv]+vu, vv+[half+d for d in vu]
    vv = [-mp.j*(v[i]-v[r])/2 for r in range(len(v)) if r not in (i, j)]
    vu = [-mp.j*(v[i]-x)/2 for x in u]
    d = -mp.j*(v[i]-v[j])/2
    return [3*half]+[3*half+x for x in vv]+[half+d]+[1+x for x in vu], \
           [1+x for x in vv]+[2+d]+[3*half+x for x in vu]


def hyper_value_and_slope(a, b):
    """Excess-one pF(p-1)(1), plus d/dgamma of its q=e^-gamma limit.

    The first-order correction is -2*K, K=prod Gamma(b)/prod Gamma(a).
    Fixed-index q-series terms have zero first derivative: this correction
    comes from their nonuniform n~1/gamma tail. See the derivation note.
    """
    excess = mp.fsum(b)-mp.fsum(a)
    if abs(excess-1) > 100*mp.eps*max(1, mp.fsum(abs(x) for x in a+b)):
        raise ValueError("only excess-one hypergeometric series are supported")
    k = mp.gammaprod(b, a)
    return mp.hyper(a, b, 1), -2*k


def rational_r(u, v, j, i):
    d = -mp.j*(v[i]-v[j])/2
    vu = [-mp.j*(v[i]-x)/2 for x in u]
    vv = [-mp.j*(v[i]-v[r])/2 for r in range(len(v)) if r not in (i, j)]
    return (mp.fprod(x/(x+mp.mpf(".5")) for x in vu)
            * mp.fprod((x+mp.mpf(".5"))/x for x in vv)/(4*d*(1+d)))


def prl_matrix(u, v):
    """M(0), M'(0). Both r_l and S have vanishing first gamma derivative."""
    ell = len(u)
    negu, negv = tuple(-z for z in u), tuple(-z for z in v)
    ratio = [mp.fprod(scattering(y-z) for z in v)
             /mp.fprod(scattering(y-z) for z in u) for y in v]
    mat, derivative = mp.matrix(ell), mp.matrix(ell)
    for j in range(ell):
        value, slope = hyper_value_and_slope(*hyper_parameters(u, v, j))
        bar, barslope = hyper_value_and_slope(*hyper_parameters(negu, negv, j))
        mat[j, j], derivative[j, j] = bar-value*ratio[j], barslope-slope*ratio[j]
        for i in range(ell):
            if i == j:
                continue
            value, slope = hyper_value_and_slope(*hyper_parameters(u, v, j, i))
            bar, barslope = hyper_value_and_slope(*hyper_parameters(negu, negv, j, i))
            r, rb = rational_r(u, v, j, i), rational_r(negu, negv, j, i)
            mat[i, j] = -rb*bar+r*value*ratio[i]
            derivative[i, j] = -rb*barslope+r*slope*ratio[i]
    return mat, derivative


def determinant_derivative(mat, derivative):
    """Multilinearity; no inverse of the singular matrix and no finite differences."""
    pieces = []
    for k in range(mat.cols):
        changed = mat.copy()
        for j in range(mat.rows):
            changed[j, k] = derivative[j, k]
        pieces.append(mp.det(changed))
    return mp.fsum(pieces)


def prl(holes, particles):
    """Analytic isotropic limit of PRL (6)--(11), including both determinant zeros."""
    u, v = points(holes, particles)
    ell = len(u)
    sigma = mp.fsum(v)-mp.fsum(u)
    sinh_sigma = mp.sinh(mp.pi*sigma/2)
    if abs(sinh_sigma) < 32*mp.eps:
        raise ValueError("zero total-rapidity sinh denominator needs a separate removable limit")
    a, ap = prl_matrix(u, v)
    b, bp = prl_matrix(tuple(-z for z in v), tuple(-z for z in u))
    da, db = determinant_derivative(a, ap), determinant_derivative(b, bp)

    def pair(z):
        return mp.sinh(mp.pi*z/2)**2*mp.exp(log_psi(z)+log_psi(-z))

    pairs = (mp.fprod(pair(u[i]-u[j])*pair(v[i]-v[j])
                     for i in range(ell) for j in range(i+1, ell))
             /mp.fprod(pair(x-y) for x in u for y in v))
    scalar = ((mp.pi*mp.exp(log_psi(0))/2)**(2*ell)*mp.pi**2
              /sinh_sigma**2)
    amplitude = da*db*scalar*sine_squared(u, v)*pairs
    return dict(amplitude=amplitude, density=amplitude*measure(ell),
                matrices=(a, b), matrix_slopes=(ap, bp),
                determinant_slopes=(da, db), null_determinants=(mp.det(a), mp.det(b)),
                pair_factor=pairs, scalar=scalar)


def finite_q_matrix(holes, particles, gamma, *, maxterms=100000):
    """Unexpanded PRL M at gamma>0, for checking the analytic limit only.

    Direct q-series and q-Gamma products: deliberately not the production
    near-isotropic algorithm. No outer integral or finite-q Barnes factor.
    """
    u, v = points(holes, particles)
    gamma = mp.mpf(gamma)
    if not mp.isfinite(gamma) or gamma <= 0:
        raise ValueError("gamma must be positive")
    q, base = mp.exp(-gamma), mp.exp(-4*gamma)
    p, h = [mp.exp(2*mp.j*gamma*y) for y in v], [mp.exp(2*mp.j*gamma*x) for x in u]

    def s(z):
        a = mp.j*z/2
        return mp.j*mp.exp(mp.j*gamma*z)*(
            mp.qgamma(1+a, base, maxterms=maxterms)*mp.qgamma(mp.mpf(".5")-a, base, maxterms=maxterms)
            /mp.qgamma(1-a, base, maxterms=maxterms)/mp.qgamma(mp.mpf(".5")+a, base, maxterms=maxterms))

    def hyper(x, y, j, i=None):
        a, b = hyper_parameters(x, y, j, i)
        return mp.qhyper([base**z for z in a], [base**z for z in b],
                         base, base, maxterms=maxterms)

    def r(pp, hh, j, i):
        d = pp[i]/pp[j]
        return (q**2*(1-q**2)**2*d/((1-d)*(1-q**4*d))
                *mp.fprod((1-pp[i]/z)/(1-q**2*pp[i]/z) for z in hh)
                *mp.fprod((1-q**2*pp[i]/pp[k])/(1-pp[i]/pp[k])
                          for k in range(len(pp)) if k not in (i, j)))

    nu, nv = tuple(-z for z in u), tuple(-z for z in v)
    pb, hb = [1/z for z in p], [1/z for z in h]
    ratios = [mp.fprod(s(y-z) for z in v)/mp.fprod(s(y-z) for z in u) for y in v]
    result = mp.matrix(len(u))
    for j in range(len(u)):
        result[j, j] = hyper(nu, nv, j)-hyper(u, v, j)*ratios[j]
        for i in range(len(u)):
            if i != j:
                result[i, j] = (-r(pb, hb, j, i)*hyper(nu, nv, j, i)
                                +r(p, h, j, i)*hyper(u, v, j, i)*ratios[i])
    return result
