#!/usr/bin/env python3
"""Independent small-ring Markov-matrix reference; no Bethe roots are used.

Requires mpmath. Repeating at 100 and 130 decimal digits tests the 90-digit
L=5,N=2 reference used by the C++ precision regression.
"""
import argparse
import itertools
import json
import mpmath as mp


def markov(sites, particles):
    states = [sum(1 << j for j in occupied) for occupied in itertools.combinations(range(sites), particles)]
    index = {state: i for i, state in enumerate(states)}
    matrix = mp.matrix(len(states))
    for column, state in enumerate(states):
        for j in range(sites):
            k = (j + 1) % sites
            if state & (1 << j) and not state & (1 << k):
                following = state ^ (1 << j) ^ (1 << k)
                matrix[index[following], column] += 1
                matrix[column, column] -= 1
    assert all(sum(matrix[i, j] for i in range(len(states))) == 0 for j in range(len(states)))
    return matrix


def reference(digits):
    with mp.workdps(digits):
        eigenvalues = mp.eig(markov(5, 2), left=False, right=False)
        nonzero = [z for z in eigenvalues if abs(z) > mp.sqrt(mp.eps)]
        real = max(z.real for z in nonzero)
        leading = max((z for z in nonzero if abs(z.real - real) < mp.sqrt(mp.eps)), key=lambda z: z.imag)
        return {'decimal_digits': digits, 'real': mp.nstr(leading.real, digits-5), 'imag': mp.nstr(leading.imag, digits-5)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--digits', nargs='+', type=int, default=[100, 130])
    args = parser.parse_args()
    if min(args.digits) < 30:
        parser.error('use at least 30 decimal digits')
    values = [reference(d) for d in args.digits]
    print(json.dumps({'sites': 5, 'particles': 2, 'right_rate': 1, 'basis_size': 10, 'references': values}, indent=2))
