# Non-Hermitian quantum-group XXZ chain

**Status: native positive finite-real-root library and frontend implemented
for 0<Delta<1; complex-root branches and root-of-unity representation
accounting remain follow-ups.** This is not the existing free-end XXZ model.
The separate Delta=0 construction includes the complete fixed-magnetization
spectrum and Jordan-block sizes, available through the same frontend, but
not spin-basis generalized eigenvectors.

## Command line

```sh
bethe-xxz-qg-obc 32 --delta 0.25
bethe-xxz-qg-obc 7 --delta 0.6 --through-lines 3 --roots
bethe-xxz-qg-obc 8 --delta 0.6 --numbers 1,3 --precision fp128 --json state.json
bethe-xxz-qg-obc 8 --delta 0.6 --through-lines 4 --excitations all --roots
bethe-xxz-qg-obc 4 --delta 0.25 --numbers none
bethe-xxz-qg-obc 8 --delta 0 --sz 0 --json blocks.json
bethe-xxz-qg-obc 7 --delta 0 --sz -1/2
bethe-xxz-qg-obc --references
```

For 0<Delta<1, the default sea uses ell=N mod 2 and consecutive labels I=1,...,(N-ell)/2.
`--through-lines` changes ell; `--numbers` instead supplies explicit labels
(or `none` for the polarized state). These options exclude each other. The
label ell is not an ordinary SU(2) spin or a promised root-of-unity degeneracy.
The corresponding Bethe representative has Sz=ell/2.

`--tolerance` and `--max-iterations` control the regular-root solve.
The `state` table contains energy, energy shift from the fully polarized
reference, normalized residual, iteration count and convergence status.
`--roots` adds a `roots` table with I, z, lambda and convergence status.
Failed solves leave energies missing (JSON null or empty delimited cells),
return exit status 2, and mark root coordinates as provisional last iterates.
Unrepresentable rapidities are missing rather than serialized as valid roots.
Invalid inputs exit 1 before files are opened, even with `--force`.

The [common output options](output.md) provide screen output, independent
JSON/CSV/TSV exports and streaming delivery. For separate rectangular files,
use `--csv-table state=state.csv --csv-table roots=roots.csv --roots`.
Metadata records the imaginary boundary strength, Hamiltonian normalization,
state-selection rule and numerical controls; the summary includes CPU time.

At `--delta 0`, the frontend instead lists the **complete fixed-Sz spectrum**.
`--sz` defaults to 0 for even N and +1/2 for odd N; both positive and negative
sectors are supported. The `blocks` table contains `block_id`, `energy`,
`block_size`, `zero_occupation` and comma-separated occupied `modes`.
Each row is one Jordan block, not one distinct energy: a size-two block has
one eigenvector and one generalized eigenvector. Coincident energies remain
separate rows, ordered by zero occupation and then mode labels, not energy.
Metadata includes the full sector dimension and number of size-two blocks.
No spin-basis vectors are output.

Endpoint controls `--max-blocks` (default 100000) and `--max-mode-entries`
(default 1000000) bound the complete enumeration. Budget refusal exits 1
before any output file is opened, even with `--force`. No partial spectrum is
published. The existing regular-branch `--through-lines`, `--numbers`,
`--roots`, `--tolerance`, `--max-iterations`, `--excitations` and
`--max-candidates` options are rejected at zero;
conversely `--sz` and the endpoint budget options require Delta=0. This avoids
confusing fermion-mode labels with regular Bethe labels.

## Hamiltonian and source normalization

Use spin-half operators, J=1, N>=2, and opposite imaginary end fields:

```math
\begin{aligned}
H&=\sum_{j=1}^{N-1}\left(S_j^xS_{j+1}^x+S_j^yS_{j+1}^y+\Delta S_j^zS_{j+1}^z\right)
+\frac{i\sqrt{1-\Delta^2}}2(S_1^z-S_N^z),\\{}
\Delta&=\cos\gamma,\qquad 0\lt \gamma\lt \frac\pi2.
\end{aligned}
```

The source equations are [Gainutdinov–Hao–Nepomechie–Sommese (2015),
Eqs. (1.1)–(1.3)](https://arxiv.org/html/1505.02104). Our Hamiltonian is
their Pauli-normalized Hamiltonian divided by four, with the boundary sign
reversed. Spatial reflection reverses that sign without changing energies.
The fully polarized reference is E_F=(N-1)*Delta/4. There is no translation
momentum for this open chain.

With eta=i*gamma, the multiplicative Bethe equations are

```math
\begin{aligned}
\left[\frac{\sinh(\lambda_i+i\gamma/2)}{\sinh(\lambda_i-i\gamma/2)}\right]^{2N}
&=\prod_{j\ne i}
\frac{\sinh(\lambda_i-\lambda_j+i\gamma)}{\sinh(\lambda_i-\lambda_j-i\gamma)}
\frac{\sinh(\lambda_i+\lambda_j+i\gamma)}{\sinh(\lambda_i+\lambda_j-i\gamma)},\\{}
E&=E_F-\sum_i\frac{\sin^2\gamma}{\cosh(2\lambda_i)-\cos\gamma}.
\end{aligned}
```

## Implemented regular branch

`bethe/xxz_quantum_group_critical.hpp` uses
$`z =\tanh (\lambda)/\tan (\gamma /2)`$, with ordered positive roots and
$`(1-\Delta)\,z ^{2} \lt  1+\Delta`$. Writing p=1+Delta, q=1-Delta, its logarithmic residual is

```math
F_i=4N\arctan z_i-2\pi I_i
-2\sum_{j\ne i}\left[
\arctan\!\left(\frac{\Delta(z_i-z_j)}{p-qz_iz_j}\right)
+\arctan\!\left(\frac{\Delta(z_i+z_j)}{p+qz_iz_j}\right)\right].
```

Unlike the free-end model, there is no additional boundary phase in this
equation. A simultaneous iteration solves for $`\arctan (z_{i})`$ from the remaining
terms. The residual norm is $`\max (\lvert F \rvert)/(2\,N)`$, not an energy-error estimate.
Energy is evaluated directly as
`E_F-sum((p-q*z^2)/(1+z^2))` to avoid unnecessary hyperbolic reconstruction.

`solve_real(N, Delta, labels, options)` accepts increasing positive integer
labels, M<=floor(N/2), I<=N-M and
$`I \lt  N -M +1-(N -2\,M +2)\,\gamma /\pi`$. The last condition excludes a conservative
32-epsilon relative band at the infinity threshold. These necessary bounds
do not establish completeness or guarantee convergence for every label set.
`sea_state(N, Delta, ell, options)` chooses I=1,...,(N-ell)/2; ell must have
N's parity. For an odd chain length, explicitly supply odd ell.
No unrestricted ground-state or degeneracy claim
is inferred from this label selection alone.

```cpp
#include <bethe/xxz_quantum_group_critical.hpp>
namespace qg = bethe::xxz::quantum_group::critical;
auto state = qg::sea_state<long double>(32, 0.25L, 0);
if (state.converged) {
  auto energy = *state.energy;
}
```

Existing quantum-number, compensated-summation, scalar and solver-option
facilities are reused. Native fp64, long-double and fp128 follow the same
equations; no double fallback is used. Iteration work is O(M^2), storage O(M).
Failed solves retain their last root coordinates and residual but publish no
energy. Status distinguishes iteration exhaustion from a precision/branch
limit. An unrepresentable rapidity is also a precision failure. A zero update
budget still evaluates the seed, so an analytic one-root seed or vacuum can
converge without an update.

## Independent validation

`scripts/reference_xxz_nonhermitian.py` constructs the complex spin-basis
Hamiltonian in fixed magnetization sectors and uses a general eigensolver,
not a self-adjoint solver. NumPy is only a development-oracle dependency.
`--self-test` checks the two-site spectrum, the exact one-magnon spectrum,
and the defective two-site Delta=0 endpoint. The script also compares an
exploratory real-root iteration with independent eigenvalues.

The one-down-spin spectrum is E_F together with
$`E_{F} -\Delta +\cos (\pi \,k /N)`$, k=1,...,N-1. Small-chain direct calculations through
N=8 agree with the selected sea at Delta=0.25, 0.6 and 0.9. Native tests use
two-, three- and four-site analytic energies, independent five-/six-/eight-site
matrix references and direct substitution into the complex Bethe equations
through N=32. Thus both the imaginary boundary normalization and the absence
of the free-end boundary phase are tested.

## Why root-of-unity work is separate

At Delta=0 the even-chain sea reaches infinite rapidity. Already at N=2,
the Sz=0 Hamiltonian is nonzero but H^2=0: the two algebraically repeated
zero eigenvalues have only one eigenvector. An eigenvalue list cannot stand
in for Jordan chains. Generalized eigenvectors and complete strings require
additional treatment; see [Gainutdinov–Nepomechie (2016)](https://arxiv.org/html/1603.09249).
Both papers are also recorded in [the central citation registry](../CITATIONS.md).

Do not extrapolate the massive solver's generic representation multiplicities
to these parameters. Do not replace this chain with a Hermitian matrix having
the same real eigenvalues. Even in the implemented interval, regular-root
solutions do not by themselves specify the full root-of-unity spectrum.

Additional complex-root branches remain follow-ups.
CFT fitting must identify boundary sectors and distinguish c from an effective
central charge; a numerical Casimir coefficient is not automatically c.
RSOS restrictions and periodic loop realizations are separate representations,
not a boundary-field toggle on this open spin-chain solver.

## Regular-state excitation scans

`--excitations COUNT|all` scans the regular label family at the selected
`--through-lines` value (default N mod 2). `COUNT` retains the lowest COUNT
converged states, including the sea where applicable; `all` retains every
converged candidate, **not every state of the physical sector**. Every candidate
is attempted even when only a few results are requested. `--max-candidates`
(default 10000) bounds this work before any solves or output files are opened.
It requires `--excitations`; explicit `--numbers` cannot be combined with a scan.
Budget refusal and invalid inputs exit 1, without overwriting `--force` targets.

The `levels` table is ranked by energy shift and reports `gap_from_sea`.
The `reference` table contains the selected sector's consecutive-label sea.
If a solve fails, `failed` contains the first failed candidate, unranked, with
missing energy and gap. All failed candidates are counted in the metadata,
not just that diagnostic example. Partial scans exit 2 while retaining valid
energies; gaps are missing if the sea reference failed. Successful scans exit 0
only when every candidate and the reference converged.

With `--roots`, each root row has `source` (levels/reference/failed), `state_id`
within that source, `root_id`, I, z, lambda and a convergence flag. Reference
roots can duplicate a ranked state intentionally; failed roots are provisional.
All tables support the shared independent exports, streaming and `--no-retain`.
For example, use `--csv-table levels=levels.csv --csv-table reference=sea.csv`.

For 0<Delta<1, `real_quantum_number_window(N,Delta,ell)` reports the necessary
finite-root window I=1,...,slots for M=(N-ell)/2 roots. It shares the strict
infinity-threshold test and precision margin with `solve_real`; it is not
a completeness theorem or a guarantee that every configuration converges.
`real_excitation_count(N,Delta,ell,limit)` performs allocation-free bounded
binomial counting. At very small positive Delta, finite precision can exclude
the last sea label and leave no supported configuration; that is a refusal,
not a continuation to the exact Delta=0 construction.

```cpp
namespace qg = bethe::xxz::quantum_group::critical;
auto scan = qg::real_excitations<long double>(12, 0.6L, 8,
    {.count = 10, .max_candidates = 10000});
// scan.levels: lowest converged states, including the sea if retained
// scan.candidate_count, scan.converged_count, scan.first_unconverged
```

Every label combination in the window is attempted; `count` only limits the
number retained, not the search. Use `count=real_excitation_count(...)` to
retain all candidates, subject to the same `max_candidates` budget. The
shared bounded-heap scanner orders converged levels by their energy shifts
from the polarized reference, with lexicographic labels breaking exact ties.
This avoids losing small differences to a common extensive energy offset.

The shared result field named `ground_state` is specifically the
**consecutive-label sea in the selected sector**, not a claim about the
global ground state. Each optional gap subtracts this sea's energy shift.
If the sea fails, gaps are absent; other converged energies can still be
returned. Failed candidates are excluded from the ranked list, counted, and
represented by `first_unconverged` as a diagnostic example. `converged()`
requires every candidate and the reference to converge. Even then, only this
regular-root family has been scanned: complex roots, descendants and
root-of-unity multiplicities are not supplied.

Native tests verify threshold crossings, analytic one-root energies,
all-pairs enumeration and direct substitution of every returned two-root
state into the original complex Bethe equations. Zero-budget and mixed
success/failure scans test that no failed energies or unsupported gaps are
published. The complete-spectrum mode remains specific to Delta=0.

## Exact Delta=0 spectrum and Jordan blocks

`bethe/xxz_quantum_group_free.hpp` implements the same spin-chain Hamiltonian
at Delta=0 using its number-conserving Jordan–Wigner fermions. This is not
the ordinary Hermitian open XX chain. Its one-particle matrix has hopping
1/2 and endpoint potentials -i/2 and +i/2. The characteristic polynomial is

```math
\det(xI-h)=2^{1-N}xU_{N-1}(x).
```

where U is the Chebyshev polynomial of the second kind. Thus the energies
are cos(pi*k/N), k=1,...,N-1, and an additional zero. For odd N all are
distinct. For even N, k=N/2 coincides with the extra zero; the irreducible
tridiagonal recurrence gives only one eigenvector, so this is a single
size-two Jordan block, not two ordinary zero modes. This construction is
consistent with the endpoint Bethe solutions and diagonalizability statements
in [Gainutdinov et al., Appendices C and D](https://arxiv.org/html/1505.02104).

The many-fermion Hamiltonian acts on exterior powers of this one-particle
space. Choose a subset of the nonzero modes and a zero-space occupation z:

- Odd N: z=0 or 1, each giving one size-one block.
- Even N: z=0 or 2 gives one size-one block; z=1 gives one size-two block.

Every block has energy sum_k cos(pi*k/N) over its occupied nonzero modes
and down-spin count M=(number of those modes)+z. Equal sums can occur for
different subsets: the API preserves these as separate blocks, without
floating-point degeneracy merging. The two-dimensional zero-space wedge
has zero energy and no nilpotent action; consequently no larger Jordan
blocks arise in this spin-chain Hamiltonian. These are Hamiltonian blocks,
not a classification of the entire commuting transfer-matrix family.

```cpp
#include <bethe/xxz_quantum_group_free.hpp>
namespace qg_free = bethe::xxz::quantum_group::free;
auto blocks = qg_free::sector<long double>(8, 4); // N=8, four down spins
for (auto const& b : blocks) {
  // b.energy, b.block_size (1 or 2), b.modes, b.zero_occupation, b.down
}
```

`sector(N,M)` returns the complete sector in zero-occupation/lexicographic
order, not energy order. The sum of block sizes is binomial(N,M); for even N
there are binomial(N-2,M-1) size-two blocks (zero outside the binomial range).
`SectorOptions` bounds the number of blocks (default 100000) and total stored
mode entries (default 1000000). Exceeding either throws `std::length_error`
before output allocation; no partial spectrum is returned. `block(N,modes,z)`
evaluates a selected block without a sector scan. Dispersive labels are
strictly increasing integers in 1,...,N-1 excluding N/2 for even N.
Native fp64, long-double and fp128 use the same formulas and shared compensated
summation; these mode labels are not regular real-root Bethe labels.

Tests check all sector dimensions through N=16 and native half-filled minimum
energies $`(1-\cot (\pi /(2N)))/2`$ for even N and $`(1-\csc (\pi /(2N)))/2`$ for odd N.
Independent complex spin matrices through N=8 test nullities of H-E,
(H-E)^2 and (H-E)^3 for every energy in every magnetization sector. This
checks eigenvector counts and generalized eigenspace dimensions, not merely
matching eigenvalue lists. No spin-basis Jordan vectors are returned yet.
