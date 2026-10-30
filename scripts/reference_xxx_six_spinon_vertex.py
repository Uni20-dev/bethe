# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Research-only XXX vertex contours for two, four and six spinons.

Finite residue sums, NOT an integrated DSF or a certified infinite series.
The six-spinon double contour separates into one-variable moments. Absolute
densities follow a fresh finite-q reduction, not a fitted sum rule or a
finite-chain extrapolation. See docs/xxx-spinon-density.md.
Requires mpmath; no production dependency.
"""
import argparse
from itertools import combinations

import mpmath as mp

from reference_xxx_four_spinon import rate


def validate(beta, terms):
    beta = tuple(map(mp.mpf, beta))
    if len(beta) not in (2, 4, 6) or any(not mp.isfinite(b) for b in beta):
        raise ValueError('require two, four or six finite real rapidities')
    if len(set(beta)) != len(beta):
        raise ValueError('coincident rapidities require a separate limiting calculation')
    if isinstance(terms, bool) or not isinstance(terms, int) or terms < 1:
        raise ValueError('residue terms must be a positive integer')
    return beta


def polynomial(beta, ell, alpha):
    """One polarization's degree n-1 polynomial in the contour integrand."""
    return (mp.fprod(alpha-beta[j]+mp.pi*1j/2 for j in range(ell))
            * mp.fprod(beta[j]-alpha+mp.pi*1j/2 for j in range(ell+1, len(beta))))


def gamma_residue(beta, j, m):
    """Direct residue of the Gamma product, before polynomial and sinh.

    Included poles are beta_j+i*pi/2-2*pi*i*m, m>=0. The opposite Gamma
    ladder is excluded. The factor 2*pi*i is the derivative of its argument,
    not an extra contour-measure factor.
    """
    x = [(beta[j]-b)/(2j*mp.pi) for k, b in enumerate(beta) if k != j]
    return (2j*mp.pi*(-1)**m * mp.gamma(m-mp.mpf('.5'))/mp.factorial(m)
            * mp.fprod(mp.gamma(a-m)*mp.gamma(m-mp.mpf('.5')-a) for a in x))


def residue_rows(beta, terms):
    """Rows indexed by polarization, then pole family j and level m."""
    beta = validate(beta, terms)
    rows = [[] for _ in beta]
    for j, b in enumerate(beta):
        x = [(b-a)/(2j*mp.pi) for k, a in enumerate(beta) if k != j]
        residue = gamma_residue(beta, j, 0)
        for m in range(terms):
            alpha = b+mp.pi*1j/2-2j*mp.pi*m
            # sinh(alpha) = i*cosh(b) on every included pole. Do not evaluate
            # periodic functions at a huge imaginary argument unnecessarily.
            common = residue*1j*mp.cosh(b)
            for ell in range(len(beta)):
                value = (mp.mpc(0) if m == 0 and j > ell
                         else common*polynomial(beta, ell, alpha))
                rows[ell].append((alpha, value))
            residue *= (-(m-mp.mpf('.5'))/(m+1)
                        * mp.fprod((m-mp.mpf('.5')-a)/(a-m-1) for a in x))
    return rows


def moments(row):
    """A+, B+, A-, B-, where A=integral f*exp(+-alpha), B=integral alpha*f*exp."""
    return tuple(mp.fsum(value*mp.exp(sign*alpha)*alpha**power for alpha, value in row)
                 for sign, power in ((1, 0), (1, 1), (-1, 0), (-1, 1)))


def separated_pair(left, right):
    ap, bp, am, bm = left
    cp, dp, cm, dm = right
    return ((bp+mp.pi*1j*ap)*cm-ap*dm-(bm+mp.pi*1j*am)*cp+am*dp)/2


def contour_amplitudes(beta, terms=256, *, direct=False):
    """Return {(positive spin positions): raw contour amplitude}.

    No averaging over polarizations, no factorial and no physical coefficient.
    `direct` evaluates the six-spinon double sum literally for small controls.
    Every result retains the specified finite pole cutoff.
    """
    beta = validate(beta, terms)
    if len(beta) == 2:
        return {(): mp.mpc(1)}
    rows = residue_rows(beta, terms)
    if len(beta) == 4:
        return {(ell,): mp.fsum(value for _, value in row) for ell, row in enumerate(rows)}
    if direct:
        return {(ell, other): mp.fsum(v*w*(a-b+mp.pi*1j)*mp.sinh(a-b)
                                      for a, v in rows[ell] for b, w in rows[other])
                for ell, other in combinations(range(6), 2)}
    values = list(map(moments, rows))
    return {(ell, other): separated_pair(values[ell], values[other])
            for ell, other in combinations(range(6), 2)}


def normalization_b():
    """Gamma(1/4)*|A_-(i*pi/2)|, independently convergent real integral."""
    integral = mp.quad(lambda t: mp.exp(t)*mp.sinh(t/2)**2
                       / (t*mp.sinh(2*t)*mp.cosh(t)), [0, 1, 4, 16, mp.inf])
    return mp.gamma(mp.mpf('.25'))*mp.exp(-integral)


def density_coefficient(n):
    """Coefficient in D_n=lim eps^n prod(cosh beta) sum |X_n|^2.

    Multiplies the full pair transition product and the raw contour norm.
    Two/four controls and the six derivation are recorded in the note.
    This fixed-rapidity limit alone does not justify a global integral.
    """
    if n == 2:
        return mp.pi**2/4
    if n == 4:
        return 1/(mp.mpf(2)**16*mp.pi**8*normalization_b()**8)
    if n == 6:
        return 1/(mp.mpf(2)**42*mp.pi**32*normalization_b()**24)
    raise ValueError('only two, four and six spinons are implemented')


def density_scale(beta):
    beta = validate(beta, 1)
    return density_coefficient(len(beta))*mp.fprod(
        rate((a-b)/(2*mp.pi)) for a, b in combinations(beta, 2))


def raw_density(beta, terms=256):
    """Pointwise D_n at finite residue cutoff; no convergence claim."""
    return density_scale(beta)*mp.fsum(abs(g)**2 for g in contour_amplitudes(beta, terms).values())


def momentum_density(n, density):
    """Szz, J=1: coefficient of product dp and both unintegrated deltas.

    Full labelled momentum cube, exactly one 1/n!, fundamental xi^2 measure.
    Not a spectral-function value. No momentum/energy deltas are evaluated.
    """
    if n not in (2, 4, 6):
        raise ValueError('only two, four and six spinons are implemented')
    return (2*mp.pi)**2/(2*mp.factorial(n)*mp.pi**(2*n))*density


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--beta', nargs='+', default=['-2.1', '-1.3', '-.4', '.3', '1.2', '2'])
    parser.add_argument('--terms', type=int, nargs='+', default=[128, 256, 512])
    parser.add_argument('--dps', type=int, default=60)
    parser.add_argument('--epsilon', nargs='*', default=[], help='also check original finite-q amplitudes')
    args = parser.parse_args()
    if args.dps < 25:
        parser.error('require dps>=25')
    mp.mp.dps = args.dps
    try:
        beta = validate(args.beta, min(args.terms))
        epsilon = list(map(mp.mpf, args.epsilon))
        if any(not mp.isfinite(e) or e <= 0 for e in epsilon):
            raise ValueError('epsilon must be finite and positive')
        if epsilon and len(beta) != 6:
            raise ValueError('--epsilon currently audits six spinons only')
    except (ValueError, TypeError) as error:
        parser.error(str(error))
    print('RESEARCH: finite contour sums; no DSF integral, certified tails or finite-N extrapolation', flush=True)
    print('beta', [mp.nstr(b, args.dps) for b in beta], 'dps', args.dps, flush=True)
    scale = density_scale(beta)
    previous = None
    for terms in args.terms:
        amplitudes = contour_amplitudes(beta, terms)
        density = scale*mp.fsum(abs(a)**2 for a in amplitudes.values())
        print('XXX terms', terms, 'D', mp.nstr(density, 30),
              'relative cutoff change', 'n/a' if previous is None else mp.nstr(abs(density/previous-1), 8), flush=True)
        previous = density
    if epsilon:
        from reference_xxz_six_spinon_vertex import SixVertex
        for e in epsilon:
            vertex = SixVertex(e)
            for vacuum in (0, 1):
                for terms in args.terms:
                    amplitudes = vertex.amplitudes(beta, vacuum, terms)
                    density = e**6*mp.fprod(mp.cosh(b) for b in beta)*mp.fsum(abs(a)**2 for a in amplitudes.values())
                    print('XXZ epsilon', mp.nstr(e, 10), 'vacuum', vacuum, 'terms', terms,
                          'D', mp.nstr(density, 30), 'ratio to last XXX cutoff', mp.nstr(density/previous, 18), flush=True)


if __name__ == '__main__':
    main()
