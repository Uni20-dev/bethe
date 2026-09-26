"""Independent spin/momentum-resolved spin-1 TB exact diagonalization.

H=sum[S.S-(S.S)^2], with no Bethe equations or spinon formulas in this oracle.
Requires NumPy/SciPy. Use OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1.
"""
import argparse
import itertools
import json
import numpy as np
from scipy.sparse import coo_matrix, eye
from spin_chain_ed import momentum_basis, lowest_eigenpair


def sector(length, spin, momentum):
    powers = 3 ** np.arange(length, dtype=np.int64)
    all_words = np.arange(3 ** length, dtype=np.int64)
    words = all_words[((all_words[:, None] // powers) % 3).sum(axis=1) == length + spin]
    index = {word: i for i, word in enumerate(words)}
    digits = (words[:, None] // powers) % 3
    basis, translation = momentum_basis(words, length, momentum)
    # Spin-1 Sz basis (-1,0,1): the allowed flip-flop matrix elements
    # in S_i.S_j equal (sqrt(2)*sqrt(2))/2=1.
    dot = np.zeros((9, 9))
    for a, b in itertools.product(range(3), repeat=2):
        col = 3 * a + b
        dot[col, col] = (a - 1) * (b - 1)
        if a < 2 and b > 0:
            dot[3 * (a + 1) + b - 1, col] = 1
        if a > 0 and b < 2:
            dot[3 * (a - 1) + b + 1, col] = 1
    bond = dot - dot @ dot
    assert np.allclose(np.linalg.eigvalsh(bond), [-6] + [-2] * 3 + [0] * 5)

    def local_sum(matrix, pairs):
        rows, cols, values = [], [], []
        for i, j in pairs:
            for row, col in zip(*np.nonzero(matrix)):
                a, b = divmod(col, 3)
                c, d = divmod(row, 3)
                selected = np.flatnonzero((digits[:, i] == a) & (digits[:, j] == b))
                target = words[selected] + (c - a) * powers[i] + (d - b) * powers[j]
                rows.extend(index[word] for word in target)
                cols.extend(selected)
                values.extend([matrix[row, col]] * len(selected))
        return coo_matrix((values, (rows, cols)), shape=(len(words), len(words))).tocsr()

    h = local_sum(bond, [(i, (i + 1) % length) for i in range(length)])
    s2 = 2 * length * eye(len(words)) + 2 * local_sum(dot, list(itertools.combinations(range(length), 2)))
    # Fixed Sz=S also contains higher total spins. Their Casimir penalty
    # exceeds H's whole [-6L,0] range, selecting total spin S exactly.
    selected = h + 4 * length * (s2 - spin * (spin + 1) * eye(len(words)))
    _, vector = lowest_eigenpair(basis.conj().T @ selected @ basis)
    state = basis @ vector
    energy = float(np.vdot(state, h @ state).real)
    residual = float(np.linalg.norm(h @ state - energy * state))
    spin_residual = float(np.linalg.norm(s2 @ state - spin * (spin + 1) * state))
    phase = np.vdot(state[translation], state)
    assert residual < 1e-8 and spin_residual < 1e-8
    assert abs(phase - np.exp(2j * np.pi * momentum / length)) < 1e-10
    return dict(energy=energy, residual=residual, spin_residual=spin_residual, dimension=basis.shape[1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lengths", nargs="+", type=int, default=[4, 6, 8, 10])
    args = parser.parse_args()
    for length in args.lengths:
        if length < 4 or length % 2:
            parser.error("lengths must be even >=4")
        ground = sector(length, 0, 0)["energy"]
        for spin, k in [(1, 1), (1, length // 2), (2, 0)]:
            state = sector(length, spin, k)
            print(json.dumps(dict(length=length, spin=spin, momentum_index=k, ground=ground,
                                  gap=state["energy"] - ground, **state)), flush=True)


if __name__ == "__main__":
    main()
