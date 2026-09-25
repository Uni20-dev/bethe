#!/usr/bin/env python3
"""Small-chain non-Hermitian quantum-group XXZ oracle (NumPy development tool).

H=sum(SxSx+SySy+Delta SzSz)+i*sqrt(1-Delta^2)/2*(Sz_1-Sz_N).
Direct spin-basis eigensolver, independent of Bethe equations. The additional
positive-real-root iteration is exploratory, not a production solver. Never
use a Hermitian eigensolver, or interpret algebraic multiplicity as the number
of independent eigenvectors at a root of unity.
"""
import argparse
import json

import numpy as np


def hamiltonian(n, down, delta):
    if not 2 <= n <= 10 or not 0 <= down <= n or not 0 <= delta < 1:
        raise ValueError("oracle requires 2<=N<=10, 0<=down<=N, 0<=Delta<1")
    basis = [word for word in range(1 << n) if word.bit_count() == down]
    index = {word: i for i, word in enumerate(basis)}
    h = np.zeros((len(basis), len(basis)), dtype=complex)
    b = np.sqrt((1-delta)*(1+delta))/2
    for col, word in enumerate(basis):
        h[col, col] = 1j*b*(((word >> (n-1)) & 1) - (word & 1))
        for j in range(n-1):
            opposite = ((word >> j) ^ (word >> (j+1))) & 1
            h[col, col] += delta*(-.25 if opposite else .25)
            if opposite:
                h[index[word ^ (1 << j) ^ (1 << (j+1))], col] += .5
    return h


def regular_sea(n, m, delta):
    if delta == 0:
        return None  # Even-chain sea reaches infinite rapidity; treat separately.
    labels = np.arange(1, m+1)
    z = np.tan(np.pi*labels/(2*n))
    plus, minus = 1+delta, 1-delta
    for iteration in range(100000):
        phases = np.zeros(m)
        for i in range(m):
            for j in range(m):
                if i != j:
                    phases[i] += np.arctan(delta*(z[i]-z[j])/(plus-minus*z[i]*z[j]))
                    phases[i] += np.arctan(delta*(z[i]+z[j])/(plus+minus*z[i]*z[j]))
        target = (np.pi*labels+phases)/(2*n)
        residual = float(np.max(2*np.abs(np.arctan(z)-target), initial=0))
        if residual < 2e-14:
            energy = (n-1)*delta/4 - np.sum((plus-minus*z*z)/(1+z*z))
            return {"energy": float(energy), "z": z.tolist(), "residual": residual, "iterations": iteration}
        z = np.tan(target)
        if np.any(z <= 0) or np.any(minus*z*z >= plus) or not np.all(np.isfinite(z)):
            raise RuntimeError("sea iteration left the finite positive-real branch")
    raise RuntimeError("sea iteration did not converge")


def check():
    # The critical endpoint has a defective zero eigenvalue, not two eigenvectors.
    h = hamiltonian(2, 1, 0.)
    assert np.linalg.norm(h-h.conj().T) > 0
    assert np.array_equal(h@h, np.zeros((2, 2)))
    assert np.linalg.matrix_rank(h) == 1
    for delta in (.25, .6, .9):
        ev = np.linalg.eigvals(hamiltonian(2, 1, delta))
        assert np.max(np.abs(np.sort(ev.real)-[-3*delta/4, delta/4])) < 1e-12
        for n in (3, 4, 7, 8):
            # Exact one-magnon spectrum: the quantum-group descendant at E_F,
            # plus N-1 dispersive levels E_F-Delta+cos(pi*k/N).
            ef = (n-1)*delta/4
            expected = np.sort([ef, *(ef-delta+np.cos(np.pi*np.arange(1, n)/n))])
            ev = np.linalg.eigvals(hamiltonian(n, 1, delta))
            assert np.max(np.abs(np.sort(ev.real)-expected)) < 1e-10
            assert np.max(np.abs(ev.imag)) < 1e-10


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", nargs="+", type=int, default=[2, 3, 4, 6, 8])
    parser.add_argument("--delta", nargs="+", type=float, default=[.25, .6, .9])
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        check()
    rows = []
    for delta in args.delta:
        for n in args.sites:
            for down in range(n//2+1):
                h = hamiltonian(n, down, delta)
                values, vectors = np.linalg.eig(h)
                error = max(np.linalg.norm(h@vectors[:, j]-e*vectors[:, j]) for j, e in enumerate(values))
                sea = regular_sea(n, down, delta)
                if sea:
                    distance = float(np.min(np.abs(values-sea["energy"])))
                    if distance > 1e-10:
                        raise RuntimeError("regular sea energy is absent from independent spin spectrum")
                    sea["nearest_eigenvalue_distance"] = distance
                rows.append({"sites": n, "down": down, "delta": delta,
                             "eigenvalues": [[float(e.real), float(e.imag)] for e in sorted(values, key=lambda e: (e.real, e.imag))],
                             "max_right_residual": float(error), "regular_sea": sea})
    print(json.dumps({"normalization": "spin-half, J=1, imaginary opposite end fields", "rows": rows}, indent=2))


if __name__ == "__main__":
    main()
