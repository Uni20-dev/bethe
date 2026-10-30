# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Accelerated vertex contours with propagated series-truncation bounds.

No finite-size extrapolation, no integrated DSF, no assertion of validated
floating-point interval bounds. Every report separates analytic series-tail
control from a heuristic arithmetic indicator and precision-repeat evidence.
"""
import argparse
from dataclasses import dataclass
from itertools import combinations

import mpmath as mp

from reference_hypergeometric_tail import TailDomainError, tail_multiplier
from reference_xxx_six_spinon_vertex import residue_rows, separated_pair, validate


@dataclass(frozen=True)
class ContourEstimate:
    amplitudes: dict
    errors: dict
    arithmetic_indicators: dict
    norm: object
    norm_error: object
    norm_arithmetic_indicator: object
    terms: int
    order: int
    digits: int

    def meets_tolerance(self, rtol, atol=0):
        rtol, atol = mp.mpf(rtol), mp.mpf(atol)
        if not mp.isfinite(rtol) or not mp.isfinite(atol) or rtol < 0 or atol < 0 or rtol+atol == 0:
            raise ValueError('require finite nonnegative tolerances, at least one positive')
        return self.norm_error+self.norm_arithmetic_indicator <= atol+rtol*self.norm


def ratio_parameters(beta, ell, j, *, first_moment=False, center=0):
    """Cancel exact polynomial/Gamma factors before constructing the tail.

    The uncentered residue has power n/2+1. Multiplication by alpha-center
    reduces it by one. The m=0 cancelled pole is outside the tail domain.
    """
    delta = [1j*(beta[j]-b)/(2*mp.pi) for b in beta]
    a = [d+(mp.mpf('.5') if k < ell else -mp.mpf('.5')) for k, d in enumerate(delta)]
    b = [d+(0 if k > ell else 1) for k, d in enumerate(delta)]
    power = len(beta)//2+1
    if first_moment:
        offset = 1j*(beta[j]-center)/(2*mp.pi)-mp.mpf('.25')
        a.append(offset+1)
        b.append(offset)
        power -= 1
    return a, b, power


def product_error(x, ex, y, ey):
    return abs(x)*ey+abs(y)*ex+ex*ey


def pair_error(left, right, le, re):
    # Expand the six products so that no cancellation is used in the bound.
    products = ((1, 2, 1), (0, 3, 1), (0, 2, mp.pi),
                (3, 0, 1), (2, 1, 1), (2, 0, mp.pi))
    return mp.fsum(scale*product_error(left[i], le[i], right[j], re[j]) for i, j, scale in products)/2


def pair_cross_error(left_tail, right_tail, left_round, right_round):
    products = ((1, 2, 1), (0, 3, 1), (0, 2, mp.pi),
                (3, 0, 1), (2, 1, 1), (2, 0, mp.pi))
    return mp.fsum(scale*(left_tail[i]*right_round[j]+left_round[i]*right_tail[j])
                   for i, j, scale in products)/2


def accelerated_contours(beta, terms=32, order=12):
    beta = validate(beta, terms)
    if terms < 2 or isinstance(order, bool) or not isinstance(order, int) or order < 1:
        raise ValueError('require terms>=2 and integer order>=1')
    n = len(beta)
    if n == 2:
        return ContourEstimate({(): mp.mpc(1)}, {(): mp.mpf(0)}, {(): mp.mpf(0)},
                               mp.mpf(1), mp.mpf(0), mp.mpf(0), terms, order, mp.mp.dps)
    rows = residue_rows(beta, terms+1)
    center = mp.fsum(beta)/n
    values, errors, arithmetic = [], [], []
    for ell in range(n):
        v, e, r = [], [], []
        for moment in range(1 if n == 4 else 2):
            sums, bounds, rounds = [], [], []
            for j in range(n):
                row = rows[ell][j*(terms+1):(j+1)*(terms+1)]
                samples = [value*((alpha-center) if moment else 1) for alpha, value in row]
                prefix = mp.fsum(samples[:-1])
                a, b, power = ratio_parameters(beta, ell, j, first_moment=bool(moment), center=center)
                tail = tail_multiplier(a, b, power, terms, order)
                sums.append(prefix+samples[-1]*tail.multiplier)
                bounds.append(abs(samples[-1])*tail.error_multiplier)
                rounds.append(abs(samples[-1])*tail.arithmetic_indicator
                              + 64*mp.eps*(terms+n)*(mp.fsum(abs(x) for x in samples[:-1])
                                                        + abs(samples[-1]*tail.multiplier)))
            if n == 4:
                v.append(mp.fsum(sums))
                e.append(mp.fsum(bounds))
                r.append(mp.fsum(rounds))
            else:
                # exp(+-alpha_jm)=+-i*exp(+-beta_j), exactly on the pole ladder.
                # Both B moments are centered by the same real number; it
                # cancels identically in separated_pair's differences.
                for sign in (1, -1):
                    factors = [sign*1j*mp.exp(sign*b) for b in beta]
                    v.append(mp.fsum(f*x for f, x in zip(factors, sums)))
                    e.append(mp.fsum(abs(f)*x for f, x in zip(factors, bounds)))
                    r.append(mp.fsum(abs(f)*x for f, x in zip(factors, rounds)))
        # Above: A+,A-,B+,B-. The existing pair assembler takes A+,B+,A-,B-.
        order_indices = (0,) if n == 4 else (0, 2, 1, 3)
        values.append(tuple(v[k] for k in order_indices))
        errors.append(tuple(e[k] for k in order_indices))
        arithmetic.append(tuple(r[k] for k in order_indices))
    if n == 4:
        amplitudes = {(ell,): values[ell][0] for ell in range(n)}
        error = {(ell,): errors[ell][0] for ell in range(n)}
        roundoff = {(ell,): arithmetic[ell][0] for ell in range(n)}
    else:
        amplitudes, error, roundoff = {}, {}, {}
        for ell, other in combinations(range(n), 2):
            key = ell, other
            amplitudes[key] = separated_pair(values[ell], values[other])
            error[key] = pair_error(values[ell], values[other], errors[ell], errors[other])
            roundoff[key] = (pair_error(values[ell], values[other], arithmetic[ell], arithmetic[other])
                             + pair_cross_error(errors[ell], errors[other], arithmetic[ell], arithmetic[other])
                             + 64*mp.eps*pair_error(values[ell], values[other],
                                                    tuple(map(abs, values[ell])), tuple(map(abs, values[other]))))
    norm = mp.fsum(abs(x)**2 for x in amplitudes.values())
    norm_error = mp.fsum(2*abs(amplitudes[k])*e+e*e for k, e in error.items())
    norm_arithmetic = mp.fsum(2*(abs(amplitudes[k])+error[k])*e+e*e for k, e in roundoff.items())
    return ContourEstimate(amplitudes, error, roundoff, norm, norm_error, norm_arithmetic,
                           terms, order, mp.mp.dps)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--beta', nargs='+', default=['-2.1', '-1.3', '-.4', '.3', '1.2', '2'])
    parser.add_argument('--terms', nargs='+', type=int, default=[32, 64])
    parser.add_argument('--order', nargs='+', type=int, default=[8, 12])
    parser.add_argument('--dps', nargs='+', type=int, default=[50, 80])
    parser.add_argument('--rtol', default='1e-12')
    args = parser.parse_args()
    if min(args.dps) < 25:
        parser.error('require dps>=25')
    print('RESEARCH: analytic series-tail bound evaluated in ordinary multiprecision, not interval arithmetic', flush=True)
    print('beta', args.beta, 'requested relative screen', args.rtol, flush=True)
    previous = {}
    for dps in args.dps:
        with mp.workdps(dps):
            try:
                beta = validate(args.beta, min(args.terms))
                tolerance = mp.mpf(args.rtol)
                if not mp.isfinite(tolerance) or tolerance <= 0 or min(args.order) < 1 or min(args.terms) < 2:
                    raise ValueError('require finite rtol>0, order>=1 and terms>=2')
            except (ValueError, TypeError) as error:
                parser.error(str(error))
            for terms in args.terms:
                for order in args.order:
                    try:
                        result = accelerated_contours(beta, terms, order)
                    except TailDomainError as error:
                        print('dps', dps, 'terms', terms, 'order', order, 'REJECTED:', error, flush=True)
                        continue
                    key = terms, order
                    change = ('n/a' if key not in previous else
                              mp.nstr(abs(result.norm/previous[key]-1), 8) if previous[key] else 'undefined at zero')
                    previous[key] = result.norm
                    relative = lambda x: mp.nstr(x/result.norm, 8) if result.norm else 'undefined at zero'
                    print('dps', dps, 'terms', terms, 'order', order, 'raw contour norm', mp.nstr(result.norm, 30),
                          'relative tail bound', relative(result.norm_error),
                          'roundoff indicator', relative(result.norm_arithmetic_indicator),
                          'precision change', change, 'meets requested screens', result.meets_tolerance(tolerance), flush=True)


if __name__ == '__main__':
    main()
