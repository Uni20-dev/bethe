"""Shared translation projection and sparse eigensolve for independent ED audits."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


def momentum_basis(words, length, momentum, local_dimension=3):
    index = {word: i for i, word in enumerate(words)}
    translated = (words * local_dimension) % (local_dimension ** length)
    translated += words // (local_dimension ** (length - 1))
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
    return coo_matrix((values, (rows, columns)), shape=(len(words), count)).tocsr(), translation


def lowest_eigenpair(matrix):
    assert np.max(np.abs((matrix - matrix.conj().T).data), initial=0) < 1e-11
    if matrix.shape[0] <= 16:
        energies, vectors = np.linalg.eigh(matrix.toarray())
    else:
        energies, vectors = eigsh(matrix, k=1, which="SA", tol=2e-13,
                                 v0=np.random.default_rng(42).normal(size=matrix.shape[0]))
    return energies[0], vectors[:, 0]
