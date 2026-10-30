# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Independent coarea/proposal assembly for small physical W6 observations.

Uses direct positive PDFs and an acos energy root, with mpmath contour sums.
The exact binary midpoint inputs are shared with C++, never rounded beta files.
Run manually; mpmath is not a runtime or build dependency.
"""
from itertools import combinations
from pathlib import Path
import sys

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from reference_xxx_vertex_density import density_at_precision
from reference_xxx_spinon_proposal import inverse, native, pair_pdf, PAIRS, MATCHINGS


def proposal_pdf(beta):
    r, one = mp.mpf(8), mp.mpf(1)
    p = [native(b, r, one / 2) for b in beta]
    single = mp.fsum(pair_pdf(beta[i], beta[j], 4 * one, one / 2, r)
                    * mp.fprod(p[k] for k in range(4) if k not in (i, j))
                    for i, j in PAIRS) / 6
    double = mp.fsum(mp.fprod(pair_pdf(beta[i], beta[j], 8 * one, one, r)
                             for i, j in match) for match in MATCHINGS) / 3
    return mp.fprod(p) / 2 + single / 4 + double / 4


def observation(words, omega, dps=85):
    with mp.workdps(dps):
        free = [inverse((mp.mpf(w) + mp.mpf('.5')) / 256, mp.mpf(8), mp.mpf('.5'))
                for w in words]
        momentum = [2 * mp.atan(mp.exp(b)) for b in free]
        energy = mp.fsum(mp.pi / (2 * mp.cosh(b)) for b in free)
        lo = 2 * mp.atan(mp.exp(-8))
        total, tail, roots, sections = mp.mpf(0), mp.mpf(0), 0, []
        for section in range(3):
            half = (mp.pi / 2 + 2 * mp.pi * section - mp.fsum(momentum)) / 2
            residual = mp.mpf(omega) - energy
            if not lo < half < mp.pi - lo:
                continue
            amplitude = mp.pi * mp.sin(half)
            minimum = amplitude * mp.cos(min(half - lo, mp.pi - lo - half))
            if not minimum < residual < amplitude:
                continue
            angle = mp.acos(residual / amplitude)
            solved = [mp.log(mp.tan((half + sign * angle) / 2)) for sign in (1, -1)]
            for pair in (solved, solved[::-1]):
                beta = free + pair
                intensity = []
                for i, j in combinations(range(6), 2):
                    other = [beta[k] for k in range(6) if k not in (i, j)]
                    derivative = mp.pi / 2 * abs(mp.tanh(beta[i]) - mp.tanh(beta[j]))
                    intensity.append(proposal_pdf(other) * mp.fprod(mp.cosh(b) for b in other) * derivative)
                jacobian = 15 / mp.fsum(intensity)
                density = density_at_precision(beta, 128, 24)
                total += density.value * jacobian
                tail += density.contour_tail_error * jacobian
                roots += 1
                sections.append(section)
        return dict(value=mp.nstr(total, 75), contour=mp.nstr(tail, 75), roots=roots, sections=sections)


if __name__ == '__main__':
    cases = [((0, 1, 2, 3), '2'), ((0, 1, 254, 255), '2'),
             ((252, 253, 254, 255), '2'), ((3, 8, 247, 252), '2.24')]
    for words, omega in cases:
        print(words, omega, observation(words, omega), flush=True)
