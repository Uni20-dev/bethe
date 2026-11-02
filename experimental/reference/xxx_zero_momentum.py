#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Regenerate a regular N12 zero-momentum triplet from coarse complex seeds.

Requires mpmath>=1.2; no Bethe extension, family solver or archived data.
The total spin raising operator annihilates the singlet. This control checks
that the selected roots have momentum zero and a nonsingular Gaudin norm;
it does not fit a small numerical determinant weight to zero.
"""
import argparse
import json
from pathlib import Path
import re

import mpmath as mp
from xxx_regular_roots import equations, gaudin


def calculate(digits):
    if digits < 100:
        raise ValueError("use at least 100 decimal digits")
    with mp.workdps(digits):
        seeds = [mp.mpc('.41827'), mp.mpc('.92567'), mp.mpc('-.30676'),
                 mp.mpc('-.307074', '.999877'), mp.mpc('-.307074', '-.999877')]
        roots = list(mp.findroot(lambda *x: equations(12, x), seeds,
                                tol=mp.mpf(10)**(-digits+10), maxsteps=100))
        residual = max(abs(x) for x in equations(12, roots))
        translation = mp.fprod((x+.5j)/(x-.5j) for x in roots)
        determinant = abs(mp.det(gaudin(12, roots)))
        tolerance = mp.mpf(10)**(-digits+20)
        if residual > tolerance or abs(translation-1) > tolerance:
            raise RuntimeError("original equations or zero momentum failed")
        if not mp.isfinite(determinant) or determinant <= 1:
            raise RuntimeError("expected a resolved nonzero Gaudin determinant")
        return {'roots': [[mp.nstr(mp.re(z), digits-10),
                           mp.nstr(mp.im(z) if abs(mp.im(z)) > tolerance else 0, digits-10)]
                          for z in roots],
                'residual': mp.nstr(residual, 12),
                'translation_error': mp.nstr(abs(translation-1), 12),
                'gaudin_absolute_determinant': mp.nstr(determinant, 30)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--digits', type=int, default=120)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = calculate(args.digits)
    if args.verify:
        fixture = Path(__file__).resolve().parents[1]/'tests/xxx_zero_momentum_fixture.hpp'
        chunks = re.findall(r'(?:"[-+0-9.eE]+"\s*)+', fixture.read_text())
        strings = [''.join(re.findall(r'"([-+0-9.eE]+)"', chunk)) for chunk in chunks]
        expected = [v for pair in result['roots'] for v in pair]
        if len(strings) != len(expected):
            raise RuntimeError('unexpected fixture structure')
        with mp.workdps(args.digits):
            tolerance = mp.mpf(10)**(-min(args.digits-20, 90))
            if any(abs(mp.mpf(a)-mp.mpf(b)) > tolerance for a,b in zip(strings, expected)):
                raise RuntimeError('fixture differs from independent roots')
        print('N12 zero-momentum roots and resolved norm verified')
    else:
        print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
