# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Asymptotic hypergeometric tails with a rational-defect remainder bound.

For t[m+1]/t[m]=prod(m+a)/prod(m+b), with integer power=sum(b-a)>=3.
The truncation inequality is analytic in exact arithmetic. Its evaluation
uses mpmath, NOT directed-rounding interval arithmetic; arithmetic precision
must be checked separately. No fitted coefficients or sampled stopping test.
"""
from dataclasses import dataclass
from math import comb

import mpmath as mp


class TailDomainError(ValueError):
    """The requested tail begins before the proved majorant domain."""


@dataclass(frozen=True)
class TailEstimate:
    multiplier: object
    error_multiplier: object
    arithmetic_indicator: object
    power: int
    start: int
    order: int


def multiply(a, b):
    result = [mp.mpc(0)]*(len(a)+len(b)-1)
    for j, x in enumerate(a):
        for k, y in enumerate(b):
            result[j+k] += x*y
    return result


def polynomial(parameters):
    value = [mp.mpc(1)]
    for a in parameters:
        value = multiply(value, [1, a])
    return value


def expansion(a, b, order):
    """Taylor coefficients of prod(1+a*x)/prod(1+b*x)."""
    result = [mp.mpc(0)]*(order+1)
    result[0] = 1
    for x in a:
        for k in range(order, 0, -1):
            result[k] += x*result[k-1]
    for x in b:
        for k in range(1, order+1):
            result[k] -= x*result[k-1]
    return result


def coefficients(a, b, power, order):
    """Q(m)=m*sum c[k]/m^k; cancel Q(m)-R(m)Q(m+1)-1 formally."""
    ratio = expansion(a, b, order+1)
    c = [mp.mpc(1)/(power-1)]
    for n in range(1, order+1):
        terms = []
        for k in range(n):
            choose = mp.mpf(1)
            value = mp.mpc(0)
            for j in range(n+2-k):
                value += ratio[n+1-k-j]*choose
                choose *= mp.mpf(1-k-j)/(j+1)
            terms.append(c[k]*value)
        c.append(mp.fsum(terms)/(power+n-1))
    return c


def defect_polynomial(a, b, c):
    """Exact rational numerator E(x), retaining computed low-order residues.

    delta(m)=E(1/m)/[(1/m)*B(1/m)*(1+1/m)^(K-1)].
    Mathematically E starts at degree K+2. Retaining rather than erasing the
    low-degree floating-point residues makes their effect visible.
    """
    order = len(c)-1
    ab, bb = polynomial(a), polynomial(b)
    base = multiply(bb, [comb(order-1, k) for k in range(order)])
    left = multiply(base, c)
    translated = [mp.mpc(0)]*(order+1)
    for k, value in enumerate(c):
        for j in range(order-k+1):
            translated[k+j] += value*comb(order-k, j)
    right = multiply(ab, translated)
    length = max(len(left), len(right), len(base)+1)
    result = [mp.mpc(0)]*length
    for k in range(length):
        result[k] = ((left[k] if k < len(left) else 0)
                     -(right[k] if k < len(right) else 0)
                     -(base[k-1] if 0 < k <= len(base) else 0))
    return result, left, right, base


def tail_multiplier(a, b, power, start, order=12):
    """Estimate sum(t[m], m>=start)/t[start], and bound its truncation error.

    The error includes the entire infinite remainder of the asymptotic tail,
    not just its last term. Input parameters must be finite, and start must
    exceed their moduli and one. Returns a separate roundoff *indicator*.
    """
    a, b = tuple(map(mp.mpc, a)), tuple(map(mp.mpc, b))
    if (not a or len(a) != len(b) or any(not mp.isfinite(x) for x in a+b)
            or isinstance(power, bool) or not isinstance(power, int) or power < 3
            or isinstance(start, bool) or not isinstance(start, int) or start < 2
            or isinstance(order, bool) or not isinstance(order, int) or order < 1):
        raise ValueError('require finite equal-length parameters, integer power>=3, start>=2 and order>=1')
    exponent = mp.fsum(y-x for x, y in zip(a, b))
    decay = mp.re(exponent)
    if (decay <= 2
            or abs(exponent-power) > 64*mp.eps*max(1, mp.fsum(abs(x) for x in a+b))):
        raise ValueError('power must equal sum(b-a)')
    n = mp.mpf(start)
    if n <= max(abs(x) for x in a+b):
        raise TailDomainError('tail start must exceed every parameter modulus')
    c = coefficients(a, b, power, order)
    q = n*mp.fsum(value/n**k for k, value in enumerate(c))
    defect, left, right, base = defect_polynomial(a, b, c)
    # |R(m)| <= (m/(m+1))^decay exp(C/m^2), m>=N, from the
    # convergent logarithm expansion. Sum 1/m^2 <= 1/N+1/N^2.
    # Use the parameter-derived real exponent, even when it differs from
    # the nominal integer by roundoff. The full defect covers that mismatch.
    constant = (mp.fsum(abs(x)**2/(1-abs(x)/n) for x in a+b)
                + decay/(1-1/n))/2
    envelope = mp.exp(constant*(1/n+1/n**2))
    denominator = mp.fprod(1-abs(x)/n for x in b)
    # Sum m^(-s) <= N^(-s)+N^(1-s)/(s-1), with s=decay+k-1.
    # k=0 is deliberately retained too; decay>2 ensures convergence.
    weighted = lambda values: mp.fsum(abs(value)*n**(1-k)*(1+n/(decay+k-2))
                                     for k, value in enumerate(values))
    bound = envelope/denominator*weighted(defect)
    cancellation_scale = weighted(left)+weighted(right)+weighted([0]+base)
    arithmetic = 64*mp.eps*(len(a)+order)**2*(abs(q)+envelope/denominator*cancellation_scale)
    return TailEstimate(q, bound, arithmetic, power, start, order)
