# Physical spin content of a TL level

[Model overview](biquadratic.md) · [Ferromagnetic excitations](biquadratic-ferromagnetic.md)

A through-line label $`\ell`$ is not a physical spin. Each eigenvector of the
open TL module carries a whole multiplicity space, which can contain several
different SU(2) spins and repeated copies of the same spin. This space is the
same for real-root, bound-string and Q-system levels in that module, and for
either sign of the biquadratic Hamiltonian.

The shared library now resolves that space:

```cpp
#include <bethe/temperley_lieb.hpp>
auto counts = bethe::temperley_lieb::spin_one_multiplets(4);
// counts[S] is the number of spin-S multiplets per TL eigenvector.
// Here: one singlet, three triplets, three quintets, three septets,
// and one spin-4 multiplet: 1 + 9 + 15 + 21 + 9 = 55 states.
if (counts) {
    auto triplets = (*counts)[1]; // 3, not 9 magnetic states
}
```

The argument is $`\ell`$, not the chain length (except in the ferro ground
space, where $`\ell =N`$). The vector has $`\ell +1`$ entries including zeros. It is
available for `ell=0,...,45`; larger labels return `nullopt`, consistently
with overflow of the **total** representation dimension in the existing
`spin_chain_multiplicity(3,ell)` API. This deliberately does not return a
partially representable spin distribution. Long-chain energy solvers are
unaffected: a small-ell module can still be resolved on any chain admitting it.

## Optional CLI table

```sh
bethe-biquadratic-obc 8 --ferromagnetic --bound-pairs 2 --spin-content
bethe-biquadratic-obc 6 --sectors --spin-content --json sectors.json
bethe-biquadratic-obc 8 --ferromagnetic --one-defect --spin-content \
  --csv-table spin_content=spins.csv
```

`--spin-content` works with every biquadratic state-selection mode and either
Hamiltonian sign. It adds a `spin_content` table, leaving default output unchanged.
Each row gives `state_id`, `through_lines`, physical `spin`, `multiplets`,
`magnetic_states=(2S+1)*multiplets`, and `status`. Spin uses Uni20's native
half-integer column type, exported as a decimal number. Zero-count spins are
omitted. Join by `state_id` to `states`, `reference`, or `failed`; reference
rows are included even when they repeat a state already in the main table.

For an overflowing multiplicity space there is one row with null spin and
counts, and status `total dimension exceeds uint64`. This does not change
the energy solver's outcome or exit status. For an unconverged state, exact
spin counts describe its **requested module**, not a verified eigenvalue:
check the state table's convergence status before interpreting energies.
JSON/CSV/TSV and `--no-retain` streaming use the same writer and schema.

## Why the recurrence resolves physical spins

The kernel and surjective map in
[Aufgebauer and Kluemper, Sec. 3.3, Eq. (49)](https://arxiv.org/html/1003.1932v2#S3.SS3)
give the usual TL multiplicity recurrence. In the isotropic spin-1
representation the map commutes with physical SU(2), and the removed
two-site state is a singlet. We therefore lift that recurrence from dimensions
to characters:

```math
W_0=[0],\qquad W_1=[1],\qquad W_{\ell+1}=[1]\otimes W_\ell-W_{\ell-1}.
```

Our implementation expands the tensor product using
`[1] tensor [S] = [S-1] + [S] + [S+1]` for $`S \ge 1`$, with the distinct
boundary rule `[1] tensor [0] = [1]`. All final coefficients are nonnegative
integers. The first spaces are

| ell | Multiplet counts by S=0,1,...,ell | Physical states |
| ---: | --- | ---: |
| 0 | 1 | 1 |
| 1 | 0, 1 | 3 |
| 2 | 0, 1, 1 | 8 |
| 3 | 1, 1, 2, 1 | 21 |
| 4 | 1, 3, 3, 3, 1 | 55 |

Tests reconstruct the entire energy multiset at **every physical spin**
for N=2,...,6, including odd chains. The independent oracle builds
$`-(S_{i} .S _{i +1})^{2}`$ in fixed physical magnetization spaces, then subtracts
the Sz=S+1 spectrum from Sz=S to isolate spin-S multiplets. This checks
which spins accompany each energy, not just the sum of degeneracies.
Dimension sums and overflow are also tested through the exact-integer limit.
CLI tests cover all output paths, references and failures, using an independent
magnetization-character expansion to check each spin count.

## What this does not say about excitation weights

If several TL eigenvectors share an energy, their spin counts must be added;
the API is per eigenvector, not per distinct energy. It does not select a
particular state in the degenerate ferro ground manifold. Nor does it give
matrix elements, form factors, or weights for a chosen iMPS ground state.
The presence of a spin-S multiplet is not a guarantee that a particular
operator excites it. Periodic TL representations and anisotropic physical
spin chains require separate treatment; auxiliary XXZ spin is still distinct.
