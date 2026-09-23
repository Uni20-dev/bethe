# Complex roots for the quantum-group-invariant open XXZ chain

[Biquadratic application](biquadratic.md) · [Overview](../README.md)

This is a separate solver for the **same opposite-end-field reference
Hamiltonian** used by the TL layer, not the ordinary zero-field XXZ front end.
Even N, M<=N/2 and finite Delta>1 are supported. It keeps finite complex-root
deviations; no ideal-string hypothesis is imposed. All arithmetic stays in the
selected fp64, long-double or optional binary128 precision.

## Why solve a polynomial?

Real-root logarithmic equations use ordered roots and integer phases. They
cannot simply be continued through arbitrary complex-string configurations.
Instead we encode the roots in a monic polynomial

```text
x = cosh(2u) = cos(alpha),
Q(x) = product_j (x-x_j) = x^M + sum_(k=0)^(M-1) c[k] x^k.
```

The coefficients are real; zeros can be real or form conjugate pairs.
Reflection u -> -u and periodicity u -> u+i*pi are absorbed by x. This
does not assert support for arbitrary non-self-conjugate root configurations.

[Bajnok et al., section 5](https://arxiv.org/html/1910.07805v2#S5) give a
Q-system for precisely these boundaries. We use its Wronskian relation (5.17).
With eta=acosh(Delta), s=sinh(eta), and
`x_±=Delta*x ± s*sqrt(x^2-1)`, our normalization is

```text
[P(x_+) Q(x_-) - P(x_-) Q(x_+)] / [2s sqrt(x^2-1)] = (x-1)^N.
```

P has degree N-M+1; its leading coefficient is fixed by this identity, and
we fix the P -> P+aQ freedom with P_M=0. Both sides are polynomials in x.
The implementation evaluates them without square roots: the antisymmetric
basis uses `x_+*x_-=x^2+s^2` and the recurrence

```text
D_0=0, D_1=1,
D_(k+1)=2*Delta*x*D_k-(x^2+s^2)*D_(k-1).
```

For fixed Q, descending coefficient elimination determines P using only
constant, nonzero q-number pivots. The remaining M coefficients (powers
0,...,M-2 and 2M-1) form the nonlinear system. The Jacobian differentiates
that elimination analytically. A damped Newton solve uses the raw residual
for descent, preventing normalization alone from rewarding diverging
coefficients. The stopping residual is the full Wronskian's componentwise
coefficient backward error, **not** the old logarithmic phase residual or an
energy error estimate.

The reference energy is obtained directly from Q:

```text
E_ref = (N-1)*Delta/4 - (Delta^2-1)*Q'(Delta)/Q(Delta).
```

The TL and physical biquadratic energy maps and module multiplicities are
unchanged. Root extraction is for independent validation and display, not
the energy formula.

## Numerical acceptance and limitations

A small coefficient residual is not by itself an accepted physical state.
The solver also:

- recovers separated roots using the native-precision
  [polynomial-root engine](polynomial-roots.md);
- rejects roots numerically indistinguishable from x=1, -1 or Delta
  (the excluded u=0, i*pi/2 and eta/2 solutions);
- checks the original multiplicative Bethe equations, including reflected
  scattering, using scaled products to avoid product overflow;
- propagates the root engine's estimated uncertainty through the sinh
  factors. Near a string pole, this can be much larger than machine epsilon.
  The per-equation residual must fit its own estimated bound, and the
  propagated relative uncertainty must be below 1e-4. Otherwise the state
  remains unverified, even if its polynomial residual is tiny.

The reported `bethe_residual` and `bethe_residual_bound` are maxima over
equations; acceptance is checked equation by equation. These estimates are
first-order floating-point diagnostics, **not interval certificates**. An
unresolved root or exceptional configuration is never rescued by silently
switching precision or relaxing the user's coefficient tolerance. The caller
can explicitly rerun in higher precision.

The initial implementation limits a selected solve to N<=32 to bound the
cached O(N^3) polynomial basis. Conditioning can become limiting much earlier;
this cap is not a promised working size. It is not yet a replacement for the
real solver on long chains. Continuation can reuse a converged Q at a nearby
Delta, but automated adaptive continuation and large-chain string labels are
future work.

## State discovery versus solving a state

`spectrum()` is deliberately restricted to N<=8. It samples deterministic
coefficient seeds and polynomials built from real roots and conjugate pairs.
It uses neither ED seeds nor stored eigenvalues. All attempts, including
failures and duplicate solutions, consume `max_attempts`. A zero budget
returns no discoveries. Different precisions can follow different basins.

Solutions are deduplicated by their Q coefficients, not by energy, and sorted
only after discovery. The target dimension is
`choose(N,M)-choose(N,M-1)`. Matching it is a **numerical completeness check**;
close solutions may be conservatively merged and leave a search incomplete.
An incomplete set does not establish the lowest k energies, even if sorted.

## C++ interface

```cpp
#include <bethe/biquadratic_qsystem.hpp>

std::vector<long double> seed{1.5L, -4.0L/3.0L};
auto level = bethe::biquadratic::qsystem::solve<long double>(4, seed);
if (level.reference.converged) {
  auto energy = *level.energy;
  auto copies = *level.multiplicity; // 1 for this ell=0 singlet
}
auto scan = bethe::biquadratic::qsystem::spectrum<long double>(6, 2);
// scan.complete(), expected_count, attempts, failed_attempts, states
```

The underlying `bethe::xxz::quantum_group::qsystem::{System,solve,spectrum}`
API takes Delta; `bethe::temperley_lieb::qsystem::solve` takes loop weight.
Only the biquadratic adapter fixes Delta=3/2 and attaches spin-1 multiplicities.
Selected failures preserve their coefficients and available energy estimate;
missing/pole energies are optional, not plausible-looking zeroes.

The regression suite checks the analytic Jacobian, the complex Wronskian
identity, the original Bethe equations, branch continuation, native-precision
four-site singlets, failure budgets, small TL-module ED spectra and the full
physical N=2,4,6 Hilbert-space spectra with representation weights.
Eight-site module spectra are checked against independent auxiliary XXZ ED:
with the tested extended precision and sufficient attempts all module counts
are recovered, while fp64 can leave nearly ideal strings unverified.
