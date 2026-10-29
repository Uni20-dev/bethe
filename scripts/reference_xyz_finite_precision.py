#!/usr/bin/env python3
"""Independent XYZ ground energies from a spin-basis Hamiltonian and free fermions.

Uses mpmath Jacobi theta functions, not Bethe roots or Uni20. Run at two decimal
precisions and compare before updating fixtures. N=2 includes both directed bonds.
"""
import argparse
import json
import mpmath as mp


def exchange(eta, t):
    q = mp.exp(-mp.pi * t)
    return [mp.jtheta(j, mp.pi * eta, q) / mp.jtheta(j, 0, q) for j in (4, 3, 2)]


def spin_energy(n, j):
    h = mp.matrix(1 << n)
    jx, jy, jz = j
    for state in range(1 << n):
        for a in range(n):
            b = (a + 1) % n
            equal = ((state >> a) & 1) == ((state >> b) & 1)
            h[state, state] += (jz if equal else -jz) / 4
            h[state ^ (1 << a) ^ (1 << b), state] += (jx - jy if equal else jx + jy) / 4
    return mp.eigsy(h, eigvals_only=True)[0]


def xy_energy(n, j):
    jx, jy, _ = j
    # Ground-state fermion parity is (-1)^(N/2): antiperiodic for even N/2.
    offset = 0 if (n // 2) % 2 else mp.mpf('0.5')
    return -mp.fsum(mp.sqrt((jx + jy)**2 * mp.cos(2 * mp.pi * (m + offset) / n)**2
                            + (jx - jy)**2 * mp.sin(2 * mp.pi * (m + offset) / n)**2)
                    for m in range(n)) / 4


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--digits', type=int, default=140)
    args = parser.parse_args()
    mp.mp.dps = args.digits
    cases = [(4, '0.2', '2', 'spin'), (4, '0.4', '0.7', 'spin'),
             (4, '0.8', '0.3', 'spin'), (6, '0.4', '0.7', 'spin'),
             (6, '0.5', '0.7', 'xy'), (16, '0.5', '0.7', 'xy')]
    results = []
    for n, eta, t, method in cases:
        j = exchange(mp.mpf(eta), mp.mpf(t))
        energy = spin_energy(n, j) if method == 'spin' else xy_energy(n, j)
        # Independently check the parity prescription at the small XY ring.
        if method == 'xy' and n == 6:
            assert abs(energy - spin_energy(n, j)) < mp.mpf(10)**(-args.digits + 8)
        results.append(dict(sites=n, eta=eta, t=t, method=method,
                            energy=mp.nstr(energy, args.digits - 5)))
    print(json.dumps({'digits': args.digits, 'cases': results}, indent=2))


if __name__ == '__main__':
    main()
