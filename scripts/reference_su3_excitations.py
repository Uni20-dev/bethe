"""Independent momentum/Casimir-resolved color-word ED for ULS benchmarks.

Requires NumPy/SciPy; run with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1.
No Bethe equations or thermodynamic dispersion formulas enter this oracle.
"""
import argparse
import itertools
import json
import numpy as np
from scipy.sparse import coo_matrix, eye
from scipy.sparse.linalg import eigsh


def sector(length, populations, momentum, target_casimir):
    powers = 3 ** np.arange(length, dtype=np.int64)
    words = []
    for first in itertools.combinations(range(length), populations[0]):
        remaining = sorted(set(range(length)) - set(first))
        for second in itertools.combinations(remaining, populations[1]):
            digits = np.full(length, 2, dtype=np.int64)
            digits[list(first)] = 0
            digits[list(second)] = 1
            words.append(int(digits @ powers))
    words = np.array(sorted(words), dtype=np.int64)
    index = {word: i for i, word in enumerate(words)}
    digits = (words[:, None] // powers) % 3
    translated = (words * 3) % (3 ** length) + words // (3 ** (length - 1))
    translation = np.array([index[word] for word in translated])
    seen = set()
    rows, columns, values = [], [], []
    count = 0
    for start in range(len(words)):
        if start in seen:
            continue
        orbit, current = [], start
        while current not in seen:
            seen.add(current)
            orbit.append(current)
            current = translation[current]
        period = len(orbit)
        if momentum * period % length:
            continue
        rows.extend(orbit)
        columns.extend([count] * period)
        values.extend(np.exp(-2j * np.pi * momentum * np.arange(period) / length) / np.sqrt(period))
        count += 1
    basis = coo_matrix((values, (rows, columns)), shape=(len(words), count)).tocsr()

    def swaps(pairs):
        rows = []
        for i, j in pairs:
            swapped = words + (digits[:, j] - digits[:, i]) * (powers[i] - powers[j])
            rows.extend(index[word] for word in swapped)
        cols = np.tile(np.arange(len(words)), len(pairs))
        return coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(words), len(words))).tocsr()

    h = swaps([(i, (i + 1) % length) for i in range(length)])
    # The fixed-weight space also contains higher irreps. A Casimir penalty
    # selects 1 or 8, rather than mislabelling its lowest level as that irrep.
    casimir_operator = swaps(list(itertools.combinations(range(length), 2)))
    casimir_operator += (9 * length - length * length) / 6 * eye(len(words))
    selected = h + length * (casimir_operator - target_casimir * eye(len(words)))
    reduced = basis.conj().T @ selected @ basis
    assert np.max(np.abs((reduced - reduced.conj().T).data), initial=0) < 1e-12
    if count <= 16:
        energies, vectors = np.linalg.eigh(reduced.toarray())
        energy, vector = energies[0], vectors[:, 0]
    else:
        energies, vectors = eigsh(reduced, k=1, which="SA", tol=2e-12,
                                 v0=np.random.default_rng(42).normal(size=count))
        energy, vector = energies[0], vectors[:, 0]
    state = basis @ vector
    energy = float(np.vdot(state, h @ state).real)
    residual = np.linalg.norm(h @ state - energy * state)
    phase = np.vdot(state[translation], state)
    # C2=sum_{i<j} P_ij + (9L-L^2)/6, fundamental generators Tr(Ta Tb)=delta_ab/2.
    casimir_state = casimir_operator @ state
    casimir = float(np.vdot(state, casimir_state).real)
    casimir_residual = float(np.linalg.norm(casimir_state - casimir * state))
    assert residual < 1e-9 and casimir_residual < 1e-8
    assert abs(phase - np.exp(2j * np.pi * momentum / length)) < 1e-10
    return dict(energy=float(energy), casimir=casimir, residual=float(residual),
                casimir_residual=casimir_residual, dimension=count)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lengths", nargs="+", type=int, default=[6, 9, 12])
    args = parser.parse_args()
    for length in args.lengths:
        if length < 3 or length % 3:
            parser.error("lengths must be multiples of three >=3")
        n = length // 3
        ground = sector(length, (n, n, n), 0, 0)
        assert abs(ground["casimir"]) < 1e-8
        for k in sorted({1, n, length // 2}):
            state = sector(length, (n + 1, n, n - 1), k, 3)
            assert abs(state["casimir"] - 3) < 1e-8, "lowest weight-sector level is not adjoint"
            print(json.dumps(dict(length=length, momentum_index=k, ground=ground["energy"],
                                  gap=state["energy"] - ground["energy"], **state)), flush=True)


if __name__ == "__main__":
    main()
