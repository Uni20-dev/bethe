#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Regenerate odd-ring continuation and sector-scan fp64 references.

Optional maintainer check, not a build/test dependency. Requires NumPy/SciPy.
Run with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3
tests/reference_xxz_odd_ed.py. Peak basis size is 352716 (N=21, M=10).
Use --all-sectors --sites 13 17 --deltas -0.7 -0.99 for the sector-scan audit.
The Hamiltonian is constructed directly from spin bitstrings, without any
Bethe equations: H=sum_j (Sx_j Sx_{j+1}+Sy_j Sy_{j+1}+Delta Sz_j Sz_{j+1}).
"""

import argparse

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh


def sector_hamiltonian_parts(n, m=None):
    if m is None:
        m = n // 2
    basis = np.fromiter((b for b in range(1 << n) if b.bit_count() == m), dtype=np.int64)
    ids = np.full(1 << n, -1, dtype=np.int64)
    ids[basis] = np.arange(len(basis))
    diagonal = np.zeros(len(basis))
    rows, cols = [], []
    for j in range(n):
        k = (j + 1) % n
        anti = ((basis >> j) ^ (basis >> k)) & 1
        diagonal += (1 - 2 * anti) / 4
        source = np.flatnonzero(anti)
        target = ids[basis[source] ^ (1 << j) ^ (1 << k)]
        rows.append(source)
        cols.append(target)
    r, c = np.concatenate(rows), np.concatenate(cols)
    hopping = coo_matrix((np.full(len(r), 0.5), (r, c)), shape=(len(basis), len(basis))).tocsr()
    return hopping, diagonal


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", nargs="+", type=int, default=[13, 17, 21])
    parser.add_argument("--deltas", nargs="+", type=float, default=[-0.5, -0.97, -0.999])
    parser.add_argument("--all-sectors", action="store_true", help="check every spin-reversal-distinct sector")
    args = parser.parse_args()
    if any(n < 3 or n > 21 or n % 2 == 0 for n in args.sites):
        parser.error("this small-chain oracle requires odd N between 3 and 21")
    if any(not np.isfinite(d) or not -1 < d <= 0 for d in args.deltas):
        parser.error("references require finite -1 < Delta <= 0")
    for n in args.sites:
        sectors = range(n // 2 + 1) if args.all_sectors else [n // 2]
        for m in sectors:
            hopping, diagonal = sector_hamiltonian_parts(n, m)
            for delta in args.deltas:
                hamiltonian = hopping + diags(delta * diagonal)
                if len(diagonal) == 1:
                    e, v = np.array([delta * diagonal[0]]), np.ones((1, 1))
                else:
                    v0 = np.random.default_rng(731).normal(size=len(diagonal))
                    e, v = eigsh(hamiltonian, k=1, which="SA", tol=2e-13, ncv=min(30, len(diagonal)),
                                 maxiter=20000, v0=v0)
                residual = np.linalg.norm(hamiltonian @ v[:, 0] - e[0] * v[:, 0])
                if not np.isfinite(residual) or residual >= 1e-10:
                    raise RuntimeError(f"unconverged reference N={n}, M={m}, Delta={delta}: {residual}")
                print(f"N={n} M={m} Delta={delta} E={e[0]:.17g} residual={residual:.4g}", flush=True)


if __name__ == "__main__":
    main()
