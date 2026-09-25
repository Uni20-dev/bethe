# Periodic TASEP: relaxation gaps, not energies

**Status: native-precision library and `bethe-tasep-pbc` frontend implemented.** This is the
totally asymmetric member of the ASEP family. It is not yet a solver for
bidirectional hopping, open reservoirs, or a complete non-Hermitian spectrum.
Bidirectional hopping is available in the separate [ASEP library/frontend](asep.md).

## Physical convention

On a ring of L sites, each of N hard-core particles jumps from j to j+1 at
rate r>0 if the destination is empty. The probability vector obeys
`dP/dt=M P`: columns of M sum to zero, off-diagonal entries are transition
rates, and diagonal entries are minus total escape rates. These are inverse
time units, not quantum energies.

The stationary eigenvalue is zero. The nonzero eigenvalues have negative real
parts. We report a leading eigenvalue lambda with nonnegative imaginary part,
and the positive relaxation gap `g=-Re(lambda)`. Its conjugate is also present;
`Im(lambda)` describes oscillation, not a second contribution to the decay rate.
The relaxation time is 1/g. N=0 or L has just one configuration and no nonzero
relaxation mode, so its gap is **absent**, not a claimed zero gap.

## Command-line use

```sh
bethe-tasep-pbc 32 --particles 8 --roots
bethe-tasep-pbc 64 --particles 32 --rate 2 --precision long-double
bethe-tasep-pbc 32 --particles 8 --roots --json rates.json \
  --csv-table relaxation=rates.csv --tsv-table roots=roots.tsv
bethe-tasep-pbc --references
```

L and `--particles` are required. `--rate` defaults to one; the shared
`--precision` choices are fp64, long-double and fp128 when MPLAPACK is enabled.
The `relaxation` table has one row, with `lambda_real`, `lambda_imag`, `gap`,
`frequency`, residual, iteration counts, convergence and status. Frequency is
the nonnegative angular frequency, equal to `lambda_imag`, in inverse time—not
cycles per time. `gap=-lambda_real`; the conjugate partner is implicit. There
are no energy or momentum columns.

`has_mode=false` distinguishes the successful empty/full sector, whose mode
observables are null, from numerical failure in a nontrivial sector. The exact
stationary eigenvalue zero is recorded in metadata independently of the gap
calculation. `--roots` adds reduced-population fugacities `Z_real`, `Z_imag`;
these are not the lattice wave numbers z or necessarily the user's N particles.

Metadata includes the generator/sign convention, rate units, root reduction,
precision, all work budgets, build/invocation provenance and CPU time. All
[shared output options](output.md), including streaming and no-retain exports,
apply. Ordinary help/calculations omit the bibliography; use `--references`.

`--max-iterations`, `--max-seed-iterations` and `--max-sites` expose the library
budgets below. Numerical failures exit 2 with null JSON observables/root
coordinates (empty CSV/TSV fields). Invalid inputs exit 1 before opening files,
including existing `--force` targets. `--left-rate` and `--excitations` are not
supported: this command is a TASEP relaxation-gap calculation, not a spectrum scan.

## Library API and output contract

```cpp
#include <bethe/tasep.hpp>
auto state = bethe::tasep::relaxation_gap(32, 8, 1.0);
if (state.converged && state.gap) {
  // *state.eigenvalue is complex, *state.gap is positive.
}
```

All arithmetic stays in the selected fp64, long-double or fp128 type. Particle–hole
symmetry reduces the root count to `effective_particles=min(N,L-N)`. The retained
roots therefore describe this reduced filling, not necessarily the user's particle
coordinates. Choosing the nonnegative-frequency representative fixes a conjugate
pair convention, not an independently assigned physical momentum.

`Options<Real>` supplies a scaled equation tolerance (default 128 native epsilons),
`max_seed_iterations` (1000 all-root sweeps), `max_iterations` (100 coupled Newton
updates), and a `max_sites` work budget (256). Both iteration budgets are explicit;
there are no hidden retries. One-particle/one-hole rates are analytic and need no
iterations. `stationary_only` is a successful empty/full sector with no eigenvalue
or gap; `converged` supplies both. Failures (`seed_limit`, `iteration_limit`,
`stalled`, `precision_limit`) never publish either observable. Any retained roots
are diagnostics only. Invalid inputs and exceeding the size budget throw.

## Equations and branch selection

For unit rate and reduced population n, use the fugacity coordinate Z=2/z-1:

```text
(1-Z_j)^n (1+Z_j)^(L-n) = Y
Y = -2^L product_j (Z_j-1)/(Z_j+1)
lambda = sum_j (Z_j-1)/2.
```

These are equations (2)–(8) of
[Golinelli–Mallick (2005)](../CITATIONS.md#golinelli-mallick-2005).
Their first-relaxation branch exchanges the upper roots on either side of the
Cassini pinch relative to the n roots with largest real parts, equations
(12)–(14). We use that finite-size branch, not their asymptotic large-L gap as
an approximation. [Gwa–Spohn (1992)](../CITATIONS.md#gwa-spohn-1992) is the earlier
spectral-gap reference.

The initial Y is negative, with `log|Y|=L log(2)+n log(n/L)+(L-n) log(1-n/L)+3.5`.
All L polynomial roots are recovered by Ehrlich–Aberth iteration using the
factored logarithmic expression. Expanding this polynomial in monomials can
destroy the branch selection at larger half-filled sizes even when the final
Bethe residual looks small. Root separation and factored residuals are checked
before the literature branch is selected.
The ratio is inverted when necessary to keep its exponential bounded; this
avoids seed-iteration overflow in dilute large rings without relaxing the
residual check.

The final solve treats selected complex roots and log(Y) together, with analytic
derivatives and real/imaginary components fed to the shared recoverable pivoted
Newton solver. Imaginary residuals are reduced modulo 2*pi. The reported norm
is the largest complex logarithmic residual divided by L. Distinct roots and
the first-harmonic translation factor are checked before publishing the rate.
Residual convergence is a numerical check, not an interval certification or a
completeness proof. Larger user-selected budgets can be expensive and do not
guarantee convergence or uniform relative gap accuracy.

## Validation and next steps

Tests assemble the full physical Markov matrix independently and compare its
leading nonzero eigenvalues for every filling at L=2 through 10. A separate
90-digit ten-configuration matrix reference tests the native precision at L=5,
N=2. All enabled precisions check the original multiplicative equations,
analytic Jacobian, rate scaling, failure budgets and selected rings through
256 sites. Dense Markov diagonalization is test-only, never the solver backend.

The frontend also tests all precision modes against the independent five-site
reference, rate scaling, exact one-particle/one-hole limits, empty/full sectors,
budget/range failures, exports, citation policy and invalid-input file protection.

The [ASEP library/frontend](asep.md) uses this TASEP seed with its own bidirectional
equations and branch continuation.
