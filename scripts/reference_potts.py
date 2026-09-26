#!/usr/bin/env python3
"""Independent critical three-state Potts clock-Hamiltonian oracle.

No Bethe equations or CFT formulas. mpmath gives arbitrary-precision energies;
--backend scipy uses sparse double-precision blocks for larger small rings.
Not a production/build dependency. Translation eigenvalue is exp(2*pi*i*k/L).
"""
import argparse
import json


def clock_block(sites, charge, momentum, real, phase):
    powers = [3**j for j in range(sites+1)]
    dimension = powers[-1]
    seen, locations, reps, periods = set(), {}, [], []
    for s in range(dimension):
        if s in seen or sum(s//powers[j] % 3 for j in range(sites)) % 3 != charge % 3:
            continue
        cycle, t = [], s
        while t not in seen:
            seen.add(t)
            cycle.append(t)
            t = (3*t) % dimension + t//powers[-2]
        if momentum*len(cycle) % sites:
            continue
        for r, t in enumerate(cycle):
            locations[t] = (len(reps), r)
        reps.append(s)
        periods.append(len(cycle))
    entries = {}

    def add(a, b, value):
        entries[a, b] = entries.get((a, b), 0)+value

    for b, s in enumerate(reps):
        digits = [s//powers[j] % 3 for j in range(sites)]
        add(b, b, sum(-2 if d == 0 else 1 for d in digits))
        for j in range(sites):
            k = (j+1) % sites
            for step in (-1, 1):
                t = s + ((digits[j]+step) % 3-digits[j])*powers[j] + ((digits[k]-step) % 3-digits[k])*powers[k]
                if t in locations:
                    a, r = locations[t]
                    add(a, b, -(real(periods[b])/real(periods[a]))**real('0.5')*phase(r))
    return len(reps), entries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sites', type=int, default=3)
    parser.add_argument('--digits', type=int, default=65)
    parser.add_argument('--backend', choices=('mpmath', 'scipy'), default='mpmath')
    args = parser.parse_args()
    if not 2 <= args.sites <= (6 if args.backend == 'mpmath' else 10) or args.digits < 20:
        parser.error('require 2<=sites<=6 (mpmath) or <=10 (scipy), and digits>=20')
    rows = []
    for q in (0, 1, -1):
        for k in range(args.sites):
            if args.backend == 'mpmath':
                import mpmath as mp
                mp.mp.dps = args.digits+15
                n, entries = clock_block(args.sites, q, k, mp.mpf,
                                         lambda r: mp.exp(2j*mp.pi*k*r/args.sites))
                h = mp.matrix(n)
                for (a, b), v in entries.items():
                    h[a, b] = v
                assert mp.norm(h-h.H) < mp.mpf(10)**(-args.digits)
                values = list(mp.eighe(h, eigvals_only=True))[:3]
                values = [mp.nstr(v, args.digits) for v in values]
            else:
                import numpy as np
                from scipy.sparse import coo_matrix
                from scipy.sparse.linalg import eigsh
                n, entries = clock_block(args.sites, q, k, float,
                                         lambda r: np.exp(2j*np.pi*k*r/args.sites))
                ij = list(entries)
                h = coo_matrix(([entries[p] for p in ij], ([a for a, b in ij], [b for a, b in ij])),
                               shape=(n, n)).tocsr()
                assert np.max(np.abs((h-h.conj().T).toarray())) < 1e-12
                values = (np.linalg.eigvalsh(h.toarray())[:3] if n < 8 else
                          np.sort(eigsh(h, k=3, which='SA', return_eigenvectors=False, tol=1e-12)))
                values = [format(v, '.16g') for v in values]
            rows.append(dict(charge=q, k=k, dimension=n, energies=values))
    print(json.dumps(dict(sites=args.sites, backend=args.backend, levels=rows), indent=2))


if __name__ == '__main__':
    main()
