# Periodic ASEP: two hopping rates

**Status: native-precision relaxation-gap library implemented; frontend pending.**
This extends [TASEP](tasep.md) to nonnegative right and left hopping rates,
including the symmetric endpoint. It does not enumerate the full spectrum or
implement open reservoirs.

## Physical convention and API

Particles hop to an empty neighbouring site on a periodic ring at rates r and s.
At least one rate must be positive; they need not sum to one. As in TASEP,
`dP/dt=M P`, with columns of M summing to zero. Eigenvalues are inverse-time
relaxation rates, not quantum energies. The stationary eigenvalue is zero.

```cpp
#include <bethe/asep.hpp>
auto state = bethe::asep::relaxation_gap(32, 8, 1.0, 0.5);
if (state.converged && state.gap) {
  // Positive decay gap: *state.gap == -state.eigenvalue->real().
  // Nonnegative angular frequency: state.eigenvalue->imag().
}
```

The arithmetic is native fp64, long-double or fp128. Reflection and particle–hole
symmetry reduce the calculation to `n=min(N,L-N)` and
`x=min(r,s)/max(r,s)`. We select the conjugate eigenvalue with nonnegative
imaginary part, without assigning it a physical momentum in the original sector.
For empty/full sectors there is only a stationary state: `stationary_only` is
successful, but both `gap` and `eigenvalue` are absent.

One-particle/one-hole results are analytic:

```text
lambda = -2(r+s) sin²(pi/L) + i |r-s| sin(2pi/L).
```

At r=s the gap is `4r sin²(pi/L)` for every nontrivial filling, with zero
frequency. These analytic cases require no iterative roots. At L=2 both
directions reach the same neighbouring site and their rates add.

## Equations and continuation

With the larger rate factored out, the Bethe equations are

```text
z_i^L = (-1)^(n-1) product_(j != i)
        [x z_i z_j -(1+x) z_i +1] / [x z_i z_j -(1+x) z_j +1]
lambda / max(r,s) = sum_j [x z_j +1/z_j -1-x].
```

These are the rate-rescaled equations (66) and (69) of
[Golinelli–Mallick (2006)](../CITATIONS.md#golinelli-mallick-2006), using
reflection to place the larger rate in the inverse-z term. The starting
branch is the existing TASEP first-relaxation solution, not a new independent
seed prescription. Adaptive continuation increases x from zero to the requested
ratio. Each step solves the coupled logarithmic equations using an analytic
Jacobian and the shared recoverable Newton solver.

Near symmetric hopping, n-1 roots approach one. Ordinary z coordinates lose
the small bias and can give a converged-looking but inaccurate frequency.
We instead store finite scaled coordinates v:

```text
delta = (max(r,s)-min(r,s))/max(r,s)
z_j = 1 + delta v_j                     (ordinary roots)
z_wave = z0 + delta v_wave              (one finite-wave root)
z0 = exp(±2 pi i/L).
```

`scaled_roots`, `wave_index`, and `wave_base=z0-1` expose this representation.
The rescaled residual divides the logarithmic equations by delta and L;
the first-harmonic driving phase is removed analytically. A shared stable
complex `log(1+z)` helper retains tiny increments. The eigenvalue is evaluated
as an analytic finite-wave contribution plus small corrections, preserving
frequencies even when the bias is only a few native epsilons. Root separation
and the translation factor are checked before publishing an observable.

The first TASEP relaxation branch is followed continuously; this is not an
exhaustive search proving that no other branch overtakes it at every size and
rate. Small-ring full-spectrum comparisons independently check the branch
selection. Numerical convergence is not an interval certification.

## Budgets and failures

`Options<Real>` supplies tolerance (128 native epsilons), `max_iterations`
(100 Newton updates per continuation attempt), `max_continuation_steps` (256
attempts, including rejected ones), and the existing `tasep::Options` as
`seed_options` (including its default 256-site work budget). Increasing the
size budget does not guarantee convergence.

`iterations` counts accepted Newton updates across all attempted stages;
`seed_iterations` and `seed_newton_iterations` report the separate seed work.
`reached_ratio` and `residual_norm` describe the last successfully reached stage,
not a claim of convergence at the requested rates when continuation fails.
Failure statuses include `seed_limit`, `iteration_limit`, `continuation_limit`,
`stalled`, and `precision_limit`. Failed solves never publish a gap or eigenvalue;
retained coordinates are diagnostic only. Invalid parameters and exceeded size
budgets throw. Overflowing physical rates are also rejected rather than reported
as a successful infinite gap.

## Validation and next step

ASEP and TASEP share an independent test-only configuration-space Markov matrix
builder. Every nontrivial filling at L=2 through 9 is compared against its full
spectrum at six hopping ratios, including both endpoints. A separate 90-digit
ten-configuration reference checks native precision at L=5, N=2, r=1, s=1/2.
All enabled precisions check original multiplicative equations, analytic
Jacobian, rate scaling, reflection/particle–hole symmetry, selected rings through
32 sites, failure budgets, and a bias only eight native epsilons from symmetry.
Dense diagonalization is never a production fallback.

Next: a frontend with both rates, shared run metadata and table exports.
`bethe-tasep-pbc` remains right-only; its options have not silently changed.
