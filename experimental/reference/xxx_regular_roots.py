#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Independent mpmath control for the regular N=12 XXX complex-root fixture.

Requires mpmath >= 1.2. No Bethe/Uni20 extension or calculation archive is read.
Use --verify to compare the committed C++ fixture, or --digits 130 to repeat
at higher precision. This selected-state oracle is not a general string solver.
"""
import argparse
import json
from pathlib import Path
import re

import mpmath as mp


def equations(n, roots):
    return [((x + .5j) / (x - .5j)) ** n
            * mp.fprod((x-y-1j)/(x-y+1j) for j, y in enumerate(roots) if i != j) - 1
            for i, x in enumerate(roots)]


def gaudin(n, roots):
    m = len(roots)
    result = mp.matrix(m)
    for i, x in enumerate(roots):
        for j, y in enumerate(roots):
            if i != j:
                result[i, j] = 2/(1+(x-y)**2)
        result[i, i] = n/(x*x+mp.mpf('.25')) - sum(result[i, j] for j in range(m) if j != i)
    return result


def determinant(matrix):
    # Independent direct determinant/products, with column equilibration to
    # avoid mpmath's scale-dependent pivot cutoff. No C++ logarithmic LU path.
    matrix = matrix.copy()
    scales = [max(abs(matrix[i, j]) for i in range(matrix.rows)) for j in range(matrix.cols)]
    for j, scale in enumerate(scales):
        for i in range(matrix.rows):
            matrix[i, j] /= scale
    return abs(mp.det(matrix)) * mp.fprod(scales)


def calculate(digits):
    if digits < 80:
        raise ValueError('use at least 80 decimal digits')
    with mp.workdps(digits):
        n = 12
        labels = [mp.mpf(j)-mp.mpf('2.5') for j in range(6)]
        ground = list(mp.findroot(
            lambda *x: [2*n*mp.atan(2*a)-2*mp.pi*q-sum(2*mp.atan(a-b) for b in x)
                        for a, q in zip(x, labels)],
            [mp.mpf(s) for s in ('-.65', '-.28', '-.08', '.08', '.28', '.65')],
            tol=mp.mpf(10)**(-digits+10), maxsteps=100))
        # Coarse seeds only. Solve all five complex multiplicative equations,
        # independently of the real/pair parameterization used by C++ solvers.
        seed = [mp.mpc('-.491'), mp.mpc('.098'), mp.mpc('.726'),
                mp.mpc('.422764', '.500138'), mp.mpc('.422764', '-.500138')]
        roots = list(mp.findroot(lambda *x: equations(n, x), seed,
                                 tol=mp.mpf(10)**(-digits+10), maxsteps=100))
        if max(abs(x) for x in equations(n, roots)) > mp.mpf(10)**(-digits+20):
            raise RuntimeError('reference roots failed their original equations')
        h = mp.matrix(6)
        for b, v in enumerate(roots):
            phase = ((v+.5j)/(v-.5j))**n
            minus = mp.fprod(u-v-1j for u in ground)
            plus = mp.fprod(u-v+1j for u in ground)
            for a, u in enumerate(ground):
                h[a, b] = (minus/(u-v-1j)-phase*plus/(u-v+1j))/(u-v)
        for a, u in enumerate(ground):
            h[a, 5] = 1/(u*u+mp.mpf('.25'))
        weight = n*determinant(h)**2 / (determinant(gaudin(n, ground))*determinant(gaudin(n, roots)))
        weight *= mp.fprod(abs(u*u+mp.mpf('.25')) for u in ground)
        weight /= mp.fprod(abs(v*v+mp.mpf('.25')) for v in roots)
        for family in (ground, roots):
            weight /= mp.fprod(abs(1+(x-y)**2) for i, x in enumerate(family) for y in family[:i])
        energy = n/mp.mpf(4) - sum(mp.mpf('.5')/(x*x+mp.mpf('.25')) for x in roots)
        if abs(mp.im(energy)) > mp.mpf(10)**(-digits+20):
            raise RuntimeError('nonreal reference energy')
        return {'roots': [[mp.nstr(mp.re(z), digits-10), mp.nstr(mp.im(z), digits-10)] for z in roots],
                'ground': [mp.nstr(z, digits-10) for z in ground],
                'energy': mp.nstr(mp.re(energy), digits-10), 'weight': mp.nstr(weight, digits-10)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--digits', type=int, default=110)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = calculate(args.digits)
    if args.verify:
        fixture = Path(__file__).resolve().parents[1]/'tests/xxx_complex_fixture.hpp'
        strings = re.findall(r'"([-+0-9.eE]+)"', fixture.read_text())
        reference = [v for pair in result['roots'] for v in pair]+result['ground']+[result['energy'], result['weight']]
        if len(strings) != len(reference):
            raise RuntimeError('unexpected C++ fixture structure')
        with mp.workdps(args.digits):
            tolerance = mp.mpf(10)**(-min(args.digits-20, 90))
            for actual, expected in zip(strings, reference):
                if abs(mp.mpf(actual)-mp.mpf(expected)) > tolerance*(1+abs(mp.mpf(expected))):
                    raise RuntimeError('C++ fixture differs from independent reference')
        print('N12 complex-root fixture agrees to at least 90 decimal places' if args.digits >= 110
              else f'N12 complex-root fixture agrees at {args.digits-20} decimal places')
    else:
        print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
