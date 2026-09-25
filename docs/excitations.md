# XXX real-root excitation scans

[Back to the overview](../README.md)

An excitation scan answers a narrower question than full diagonalization:
what are the lowest energies among the supported real-root states in a chosen
total-spin family? Keeping that distinction explicit is important when
interpreting both `COUNT` and `all`. For anisotropic chains, use the separate
[XXZ excitation guide](xxz.md#real-root-excitations).

## Select a total-spin family

The periodic `bethe-xxx-pbc` and free-end `bethe-xxx-obc`
front ends enumerate supported real-root highest-weight states at a
chosen **total spin S**, with $`M =N /2-S`$ roots and $`S^z =S`$:

```sh
build/bethe-xxx-pbc 64 --excitations 10 --spin 1
build/bethe-xxx-obc 64 --excitations 10 --spin 1
build/bethe-xxx-pbc 64 --excitations all --spin 1
build/bethe-xxx-obc 64 --excitations all --spin 1
build/bethe-xxx-pbc 65 --excitations 10 --spin 1/2 --precision long-double
build/bethe-xxx-obc 32 --excitations 5 --spin 2 --roots
```

`--spin` defaults to 1 for even N and 1/2 for odd N. It must be nonnegative,
no greater than N/2, and have the same integer/half-integer parity as N/2.
`--spin` and `--max-candidates` require `--excitations COUNT|all`; `--sz` remains
the separate sector-minimum mode, not a filter for this scan.

## What is returned

The scan returns up to COUNT lowest **converged multiplets in this family**,
including its sector minimum. Each multiplet is represented once; its $`2S +1`$
SU(2) partners are not listed separately. Distinct multiplets with equal
energies are retained, including periodic reflection/momentum partners. COUNT
can cut through such degeneracies. Results are sorted by computed energy;
exact ties use lexicographic Bethe quantum numbers. Near-degenerate ordering
can change with numerical precision, and the residual is not an energy-error
bound. For odd N and S=1/2 the family includes ground-state multiplets.

Use `--excitations all` to return every converged multiplet in the selected
real-root family, without needing its size in advance. This does not include
other total-spin sectors or unsupported string states, and still respects
`--max-candidates`. For N=64, S=1 it returns 528 multiplets if all converge.

Reports give S, absolute energy, `gap=E-E0` relative to the **global ground
state of the same finite chain**, quantum numbers, and convergence diagnostics.
Periodic reports additionally give lattice momentum; open reports do not.
Gaps retain the selected arithmetic precision, including fp128. Tiny negative
gaps due to roundoff are not clamped. Root and quantum-number tables link to
levels by zero-based `state_id`; half-integers are decimal values. An empty
quantum-number set has no rows for that state. The ground reference and first
failed candidate have separate tables and distinct IDs, so neither is confused
with a ranked level. See [result tables and exports](output.md).

## What the family leaves out

This is **not a complete low-energy spectrum**, even in the requested S
sector. It enumerates only the quantum-number windows of the existing
real-root solvers. For periodic even N, S=1 gives the conventional two-spinon
triplet family; further triplets and excited singlets require additional
families. In particular, for even N and S=0 the current window contains only
the ground-state configuration. Complex/string states and infinite-root
descendants are not solved by this mode.

## Cost and candidate limits

There are $`C (N -M,M)`$ candidates for either boundary: for example, N=64, S=1
has 528. The scan solves **every** candidate before returning the requested
lowest subset, not just the first COUNT configurations. The default
`--max-candidates 10000` rejects larger families before any solve; raise it
explicitly for larger jobs, including when using `all`. A numeric COUNT and
this limit must both be positive. The limit bounds the number of configurations,
not the cost of an individual solve. A bounded heap retains only O(COUNT*M+N)
data, including the global ground reference; the underlying O(M^2) work per
iteration is unchanged.
With `all`, roots for the entire converged family are retained.

## Failure accounting

Failed candidates are excluded from the energy-ordered list, but **all**
candidate failures count toward the scan status. Reports show total and
converged candidate counts, returned multiplets, and the first failed
configuration's diagnostics. An exhausted candidate makes the ordering
incomplete and returns exit status 2, even if all displayed states converged.
If only the ground reference fails, absolute-energy ordering can still be
complete within the family, but gaps are `unavailable` and the exit status
is also 2. No unconverged ground energy is used to manufacture gaps.

## Validation

The independent small-chain tests compare the scans against exact
diagonalization through N=8, subtracting the Sz=S+1 spectrum from Sz=S to
check total-spin multiplicities, not just energy membership. They also check
momentum, exhaustive quantum-number coverage, bounded-prefix ordering,
omitted singlets, failure accounting, and an irrational gap beyond fp64.

## Library usage

```cpp
#include <bethe/heisenberg_excitations.hpp>

namespace xxx = bethe::heisenberg;
auto const spin = uni20::half_int{1};
auto const candidates = xxx::real_excitation_count(64, spin); // no solves
xxx::RealExcitationOptions selection{.count = 10, .max_candidates = 10000};
auto periodic = xxx::real_excitations<long double>(64, spin, selection);
auto open = xxx::open::real_excitations<long double>(64, spin, selection);
// Optional fourth argument: SolverOptions<Real>.
// Each scan.levels entry contains .state and .gap (std::optional<Real>).
// scan.family_converged() covers every candidate, not only retained levels.
// scan.converged() additionally requires the ground reference to converge.
```

`real_excitation_count` accepts an optional upper limit (default: maximum
`size_t`). Exceeding that limit, `max_candidates`, or representable count
throws `std::length_error`; invalid spin or scan options throw
`std::invalid_argument`. Nonfinite numerical failures propagate from the
underlying solver. Arithmetic, sorting, and gap subtraction never narrow
the selected real type to double.

Related: the [periodic quantum-number window](xxx.md#rapidity-and-quantum-number-conventions),
the [open-chain window](open-chains.md#allowed-states-and-convergence), and
[shared output and precision controls](command-line.md).
