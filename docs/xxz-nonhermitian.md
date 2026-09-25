# Non-Hermitian quantum-group XXZ chain

**Status: native positive finite-real-root library and frontend implemented
for 0<Delta<1; complex-root branches and root-of-unity representation
accounting remain follow-ups.** This is not the existing free-end XXZ model.
The separate Delta=0 library described below includes the complete
fixed-magnetization spectrum and Jordan-block sizes, but not spin-basis
generalized eigenvectors or a command-line endpoint mode yet.

## Command line

```sh
bethe-xxz-qg-obc 32 --delta 0.25
bethe-xxz-qg-obc 7 --delta 0.6 --through-lines 3 --roots
bethe-xxz-qg-obc 8 --delta 0.6 --numbers 1,3 --precision fp128 --json state.json
bethe-xxz-qg-obc 4 --delta 0.25 --numbers none
bethe-xxz-qg-obc --references
```

The default sea uses ell=N mod 2 and consecutive labels I=1,...,(N-ell)/2.
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

## Hamiltonian and source normalization

Use spin-half operators, J=1, N>=2, and opposite imaginary end fields:

```text
H = sum_{j=1}^{N-1} (Sx_j Sx_{j+1} + Sy_j Sy_{j+1} + Delta Sz_j Sz_{j+1})
    + i*sqrt(1-Delta^2)/2 * (Sz_1-Sz_N).
Delta = cos(gamma), 0<gamma<pi/2.
```

The source equations are [Gainutdinov–Hao–Nepomechie–Sommese (2015),
Eqs. (1.1)–(1.3)](https://arxiv.org/html/1505.02104). Our Hamiltonian is
their Pauli-normalized Hamiltonian divided by four, with the boundary sign
reversed. Spatial reflection reverses that sign without changing energies.
The fully polarized reference is E_F=(N-1)*Delta/4. There is no translation
momentum for this open chain.

With eta=i*gamma, the multiplicative Bethe equations are

```text
[sinh(lambda_i+i*gamma/2)/sinh(lambda_i-i*gamma/2)]^(2N)
 = product_{j!=i} sinh(lambda_i-lambda_j+i*gamma)/sinh(lambda_i-lambda_j-i*gamma)
                 * sinh(lambda_i+lambda_j+i*gamma)/sinh(lambda_i+lambda_j-i*gamma).
E = E_F - sum_i sin(gamma)^2/(cosh(2*lambda_i)-cos(gamma)).
```

## Implemented regular branch

`bethe/xxz_quantum_group_critical.hpp` uses
`z=tanh(lambda)/tan(gamma/2)`, with ordered positive roots and
`(1-Delta)*z^2 < 1+Delta`. Writing p=1+Delta, q=1-Delta, its logarithmic residual is

```text
F_i = 4*N*atan(z_i) - 2*pi*I_i
      - 2*sum_{j!=i} [atan(Delta*(z_i-z_j)/(p-q*z_i*z_j))
                     +atan(Delta*(z_i+z_j)/(p+q*z_i*z_j))].
```

Unlike the free-end model, there is no additional boundary phase in this
equation. A simultaneous iteration solves for `atan(z_i)` from the remaining
terms. The residual norm is `max(abs(F))/(2*N)`, not an energy-error estimate.
Energy is evaluated directly as
`E_F-sum((p-q*z^2)/(1+z^2))` to avoid unnecessary hyperbolic reconstruction.

`solve_real(N, Delta, labels, options)` accepts increasing positive integer
labels, M<=floor(N/2), I<=N-M and
`I < N-M+1-(N-2*M+2)*gamma/pi`. The last condition excludes a conservative
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
`E_F-Delta+cos(pi*k/N)`, k=1,...,N-1. Small-chain direct calculations through
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

The next stages are admissible-label scans,
command-line access to the Delta=0 endpoint, and additional complex-root branches.
CFT fitting must identify boundary sectors and distinguish c from an effective
central charge; a numerical Casimir coefficient is not automatically c.
RSOS restrictions and periodic loop realizations are separate representations,
not a boundary-field toggle on this open spin-chain solver.

## Exact Delta=0 spectrum and Jordan blocks

`bethe/xxz_quantum_group_free.hpp` implements the same spin-chain Hamiltonian
at Delta=0 using its number-conserving Jordan–Wigner fermions. This is not
the ordinary Hermitian open XX chain. Its one-particle matrix has hopping
1/2 and endpoint potentials -i/2 and +i/2. The characteristic polynomial is

```text
det(x I-h) = 2^(1-N) x U_{N-1}(x),
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
energies `(1-cot(pi/(2N)))/2` for even N and `(1-csc(pi/(2N)))/2` for odd N.
Independent complex spin matrices through N=8 test nullities of H-E,
(H-E)^2 and (H-E)^3 for every energy in every magnetization sector. This
checks eigenvector counts and generalized eigenspace dimensions, not merely
matching eigenvalue lists. No spin-basis Jordan vectors are returned yet.
