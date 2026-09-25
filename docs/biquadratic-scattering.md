# Low-energy real-root scattering windows

[Ferromagnetic overview](biquadratic-ferromagnetic.md) · [Bound pairs](biquadratic-bound-pairs.md) · [Three-defect droplets](biquadratic-bound-triples.md)

The bound branches are only part of the excitation structure. Two singlet
defects can also scatter rather than bind; three can scatter independently,
or form a pair plus a separate defect. The **positive finite-real-root
family** describes the unbound branch. A separate
[pair-plus-real-root library API](biquadratic-pair-defect.md) uses the
complex-string equations; it is not included in the real-root CLI scan.

The existing `--excitations` scan handles the full supported real family.
For a few defects on a long chain, the number of label combinations can be
large even though each Bethe solve is small. `--real-window WIDTH` restricts
that same scan to high labels, near the ferromagnetic low-energy edge:

```sh
# Two defects: ell=N-4. Scan choose(8,2)=28 candidates; retain the first ten.
build/bethe-biquadratic-obc 128 --ferromagnetic --through-lines 124 \
  --excitations 10 --real-window 8
# Odd chain, every combination in the selected window:
build/bethe-biquadratic-obc 129 --ferromagnetic --through-lines 125 \
  --excitations all --real-window 8 --roots
# Three unbound defects: ell=N-6; choose(6,3)=20 candidates.
build/bethe-biquadratic-obc 100000 --ferromagnetic --through-lines 99994 \
  --excitations all --real-window 6 --precision long-double \
  --json scattering.json --csv scattering.csv
```

The exact definition is

```math
\begin{aligned}
M&=\frac{N-\ell}{2},\\{}
I_{\min}&=N-M-\mathrm{WIDTH}+1,\qquad I_{\max}=N-M,\\{}
I_{\min}&\le I_1\lt \cdots\lt I_M\le I_{\max},\\{}
1&\le M\le\mathrm{WIDTH}\le N-M.
\end{aligned}
```

`--through-lines` must have the same parity as N. Without it, ferro scans
still default to the one-defect module ell=N-2. The window option requires
`--ferromagnetic --excitations COUNT|all`; it is not a new Hamiltonian or
solver. Without `--real-window`, the full real-family scan is unchanged.
`WIDTH=N-M` reproduces that full scan. Explicit `--quantum-numbers` and
ferro real scans now support odd and even chains; AF ground/sector helpers
and Q-system searches retain their even-chain restriction.

## What ordering and `all` mean

Every candidate in the window is solved before retaining the lowest COUNT
converged levels. `all` means **all real-root combinations in that window**,
not all physical excitations, not all real-root states of the module unless
the window spans the full label range, and not a count of bound-string modes.
The row metadata reports the first/last allowed I, width, candidate count
and convergence count. A converged window is not a certificate that the
lowest COUNT states outside it have been excluded.

Increase WIDTH and compare overlapping levels to check window coverage.
`--max-candidates` (default 10000) checks `choose(WIDTH,M)` before allocating
roots or solving a state. At fixed M and WIDTH, the number of solves does
not grow with N. Each Newton update uses O(M^2) storage and O(M^3) work.

For N=128, M=2, the first three returned gaps with WIDTH=6 are

| Integer labels | E-E0 |
| --- | --- |
| 125,126 | 2.003157340076926308... |
| 124,126 | 2.006312443850519769... |
| 124,125 | 2.008205709376788139... |

The low unbound two-defect branch approaches 2, whereas the bound-pair
edge approaches 5/3. Three unbound defects approach 3, while the
three-defect droplet approaches 2. The pair-plus-one-defect scattering
threshold is 8/3; its finite-size mixed-string states are **not** supplied
by a purely real-root scan. These comparisons refer to TL singlet insertions,
not physical spin flips. Neither these labels nor the OBC rapidities are
physical lattice momentum, and no operator spectral weights are inferred.

## Direct gaps, shared equations

The shared [quantum-group XXZ solver](../include/bethe/xxz_quantum_group.hpp)
uses the same logarithmic equations described in the [TL guide](biquadratic.md),
equivalently [Bajnok et al., Eq. (5.12)](../CITATIONS.md#bajnok-2020).
Nothing in those equations requires even N for a selected admissible real
state. The odd-chain extension is checked against the original complex
equations and independent highest-weight-module ED, not inferred solely
by removing an input check.

At Delta=3/2, in the existing angle coordinate x, evaluate directly

```math
\begin{aligned}
\delta E_{\mathrm{ref}}&=-\sum_i[(\Delta-1)+2\cos^2x_i],\\{}
\mathrm{gap}&=\mathrm{tl\_energy}=-2\delta E_{\mathrm{ref}},\\{}
E&=(N-1)+\mathrm{gap}.
\end{aligned}
```

This keeps small gaps and their ordering independent of the extensive
total-energy offset. For instance, the four lowest one-defect levels at
N=100000000 have distinct fp64 gaps even though their total energies round
to the same value. Sorting total energies would incorrectly break these
ties by labels. The common combination scanner now accepts an energy
projection; other models keep their existing default, while ferro uses
the direct shift. The generic TL adapter also uses that shift rather than
subtracting extensive XXZ energies.

The root solver, bounded-combination enumeration, heap selection, TL mapping,
and all state/root/export tables are reused. There is no separate scattering
Newton implementation. fp64, native long double and enabled fp128 all use
the same path. Roots that cannot be resolved as distinct interior points
at the selected precision are rejected, rather than reported as converged
endpoint states.

## Diagnostics and library use

Failed candidates are excluded from the ranked list, counted, and represented
by an example in `failed`; the process exits with status 2. An all-failed
window can have an empty `states` table but still a valid exact ground
`reference`. Tables and JSON/CSV/TSV exports use the existing shared
[output contract](output.md), including `state_id` links to roots and labels.

```cpp
#include <bethe/biquadratic_ferromagnetic.hpp>
auto scan = bethe::biquadratic::ferromagnetic::real_excitations_window<long double>(
    129, 125, 8, {.count=10, .max_candidates=10000});
if (scan.converged()) {
  // Sorted only within the selected window; not the full module.
  auto gap = scan.levels.front().gap;
}
```

Validation checks full/windowed agreement, direct native-precision gaps,
original complex equations, and every odd-chain real family through N=9
against module ED. For N=4,...,10, the union of the **full** two-real-root
family and the targeted bound-pair family reproduces the entire two-defect
module, eigenvalue by eigenvalue. This finite validation is not a proof of
arbitrary-N Bethe completeness. Three defects also need the
[mixed pair-plus-defect family](biquadratic-pair-defect.md); including it
reproduces the complete three-defect module in tests for N=6,...,10.
