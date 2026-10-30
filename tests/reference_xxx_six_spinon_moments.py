# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Independent direct-global measure: positive PDFs and mpmath W6.

The fixture deliberately selects the background four-variable channel, but
its weight contains all overlapping one-pair and two-pair channels and all
15 choices of the native pair. No historical samples or normalization fit.
"""
from itertools import combinations
import mpmath as mp
from reference_xxx_six_spinon_spectrum import inverse, native, proposal_pdf, density_at_precision


def observation(words=(4, 32, 80, 112, 180, 220, 240), dps=85):
    with mp.workdps(dps):
        u = [(mp.mpf(w) + mp.mpf('.5')) / 256 for w in words]
        position = 15*u[0]
        pair_index = int(position)
        assert position-pair_index < mp.mpf('.5')
        pair = tuple(combinations(range(6), 2))[pair_index]
        free = [inverse(x, mp.mpf(8), mp.mpf('.5')) for x in u[1:5]]
        beta = [None]*6
        beta[pair[0]], beta[pair[1]] = [inverse(x, mp.mpf(8), mp.mpf('.5')) for x in u[5:7]]
        rest = iter(free)
        for i in range(6):
            if beta[i] is None: beta[i] = next(rest)
        q = mp.fsum(proposal_pdf([beta[k] for k in range(6) if k not in (i,j)])
                    *native(beta[i], mp.mpf(8), mp.mpf('.5'))
                    *native(beta[j], mp.mpf(8), mp.mpf('.5'))
                    for i,j in combinations(range(6), 2))/15
        measure = mp.fprod(mp.sech(b) for b in beta)/q/(4*mp.pi**2)
        energy = mp.pi/2*mp.fsum(mp.sech(b) for b in beta)
        density = density_at_precision(beta,128,24)
        return {name:mp.nstr(value,75) for name,value in
                {'I0':measure*density.value, 'I1':measure*energy*density.value,
                 'tail0':measure*density.contour_tail_error,
                 'tail1':measure*energy*density.contour_tail_error}.items()}


def fixed_observation(words=(2, 32, 80, 112, 180, 220), dps=85):
    """Independent positive-PDF assembly on a fixed momentum surface."""
    with mp.workdps(dps):
        u = [(mp.mpf(w) + mp.mpf('.5')) / 256 for w in words]
        position = 30*u[0]
        chart = int(position)
        assert position-chart < mp.mpf('.5')  # background draw, full-mixture weight
        eliminated, native_label = divmod(chart, 5)
        if native_label >= eliminated:
            native_label += 1
        beta = [None]*6
        beta[native_label] = inverse(u[5], mp.mpf(8), mp.mpf('.5'))
        free = iter(inverse(x, mp.mpf(8), mp.mpf('.5')) for x in u[1:5])
        for i in range(6):
            if i not in (eliminated, native_label): beta[i] = next(free)
        free_momentum = mp.fsum(2*mp.atan(mp.exp(b)) for b in beta if b is not None)
        lo = 2*mp.atan(mp.exp(-8))
        roots = [(s, mp.pi/2+2*mp.pi*s-free_momentum) for s in range(3)]
        roots = [(s, p) for s, p in roots if lo < p < mp.pi-lo]
        assert len(roots) == 1
        section, p = roots[0]
        beta[eliminated] = mp.log(mp.tan(p/2))
        h = []
        for i in range(6):
            for j in range(6):
                if i == j: continue
                other = [beta[a] for a in range(6) if a not in (i, j)]
                h.append(proposal_pdf(other)*native(beta[j], mp.mpf(8), mp.mpf('.5'))
                         *mp.fprod(mp.cosh(beta[a]) for a in range(6) if a != i))
        measure = 30/mp.fsum(h)/(2*mp.pi)
        energy = mp.pi/2*mp.fsum(mp.sech(b) for b in beta)
        density = density_at_precision(beta, 128, 24)
        result = {name: mp.nstr(value,75) for name,value in
                  {'m0': measure*density.value, 'm1': measure*energy*density.value,
                   'tail0': measure*density.contour_tail_error,
                   'tail1': measure*energy*density.contour_tail_error}.items()}
        result['section'] = section
        return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixed', action='store_true', help='fixed k=pi/2 rather than global fixture')
    args = parser.parse_args()
    print(fixed_observation() if args.fixed else observation(), flush=True)
