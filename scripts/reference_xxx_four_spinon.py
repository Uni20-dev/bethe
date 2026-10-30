# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Independent mpmath reference for the physical XXX four-spinon density.

Run with Python and mpmath (not a runtime dependency). Caux--Hagemans,
https://arxiv.org/abs/cond-mat/0611319, Eqs. (16)--(20). The matrix-element
series is summed as an 8F7 at unity, independently of the C++ tail recurrence.
The transition kernel uses a different analytic subtraction. The raw kernel
is converted internally to physical longitudinal W4 against labelled product(dp),
using the independently derived coefficient G(1/2)^16/3072. No sum rule is fitted.
This tests pointwise algebra and normalization, not the continuum integration.
"""
import mpmath as mp


def rate(r):
    a = 4 * abs(r)

    def f(x):
        t = mp.tanh(x)
        return ((1 + t) / mp.tanh(2*x) * (-2*mp.sin(a*x/2)**2)
                + t + t*t + 2*mp.expm1(-x)*mp.cos(a*x)) / x

    return mp.exp(-mp.log1p(1/a**2)
                  - mp.quad(f, [mp.mpf(i)/2 for i in range(201)]))


def raw_kernel(rs):
    """Unnormalized four-polarization series; arguments are rho=beta/(2*pi)."""
    r = list(map(mp.mpf, rs))
    g = [mp.mpc(0) for _ in r]
    for j in range(4):
        d = [1j*(r[j]-r[i]) for i in range(4)]
        pre = mp.cosh(2*mp.pi*r[j]) / mp.fprod(
            mp.sinh(mp.pi*(r[j]-r[i])) for i in range(4) if i != j)
        prod = mp.fprod(mp.gamma(mp.mpf('.5')+x)/mp.gamma(2+x) for x in d)
        zero = prod * mp.fprod((1+x)/(x-mp.mpf('.5')) for x in d)
        for l in range(4):
            poly = [d[i]-(mp.mpf('.5') if l > i else 0) for i in range(4) if i != l]
            if j <= l:
                g[l] += pre * zero * mp.fprod(poly)
            # Independent hypergeometric sum, not the C++ rational-tail recurrence.
            num = [mp.mpf('.5')+x for x in d] + [2+x for x in poly] + [1]
            den = [2+x for x in d] + [1+x for x in poly]
            g[l] += pre * prod * mp.fprod(1+x for x in poly) * mp.hyper(num, den, 1)
    f = mp.fprod(rate(r[j]-r[i]) for j in range(4) for i in range(j))
    return f * sum(abs(x)**2 for x in g)


def density(beta):
    """Physical Szz W4, including 1/4!, against labelled product(dp)."""
    beta = tuple(map(mp.mpf, beta))
    if len(beta) != 4 or any(not mp.isfinite(x) for x in beta):
        raise ValueError('require four finite rapidities')
    if len(set(beta)) != 4:
        raise ValueError('coincident rapidities need a limiting calculation')
    return (mp.barnesg(mp.mpf('.5'))**16 / 3072
            * raw_kernel([x/(2*mp.pi) for x in beta]))


if __name__ == '__main__':
    mp.mp.dps = 40
    for rs in [['-.6', '-.1', '.3', '.8'], ['-1', '-.4', '.2', '1.3']]:
        beta = [2*mp.pi*mp.mpf(x) for x in rs]
        print('beta', [mp.nstr(x, 40) for x in beta],
              'W4', mp.nstr(density(beta), 35), flush=True)
    # A soft pair inside the finite radial-moment audit domain.
    v = mp.mpf('.25')
    r = mp.pi/(3-2*v)*mp.exp(-12)
    lo, hi = r*(1-v), (mp.pi-r)/2
    p = [r*v, r*(1-v), lo+(hi-lo)*mp.mpf('.4')]
    p.append(mp.pi-sum(p))
    rho = [-mp.asinh(1/mp.tan(x))/(2*mp.pi) for x in p]
    print('soft-pair rho', [mp.nstr(x, 40) for x in rho], flush=True)
    print('soft-pair W4', mp.nstr(density([2*mp.pi*x for x in rho]), 40), flush=True)
