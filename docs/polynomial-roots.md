# Recovering roots without losing the polynomial

[Coordinate wavefunctions](xxz-coordinate-wave.md) | [Mixed phantom roots](xxz-phantom.md) | [Model catalogue](models.md)

The internal `bethe::detail::recover_polynomial_roots<Real>` in
[polynomial_roots.hpp](../include/bethe/polynomial_roots.hpp) connects the
continued XXZ polynomial to explicit finite-root coordinate wavefunctions.
It accepts a monic polynomial with real coefficients, stored without its
leading 1 in ascending order, and returns complex roots in the same
coordinate. This building block is independent of the XXZ Hamiltonian.

## Iteration and numerical arithmetic

The Ehrlich-Aberth step is

```text
x_j <- x_j - Q(x_j) / [Q'(x_j) - Q(x_j)*sum_(k!=j) 1/(x_j-x_k)].
```

We use Gauss-Seidel sweeps, so later roots see earlier updates in the same
sweep. No root is explicitly divided out of Q. See
[Bini and Robol](../CITATIONS.md#bini-robol-2013), equation (3), for the
iteration and its relation to implicit deflation. Local convergence is
not a guarantee of convergence from every initial configuration.

Ordinary Horner evaluation can lose the small residual needed to refine
close roots. We therefore propagate a rounding correction alongside the
complex Horner recurrence, using real TwoSum and split-product transforms
at the selected precision. The compensated-evaluation strategy is described
by [Cameron and Graillat](../CITATIONS.md#cameron-graillat-2022).
This implementation does **not** reproduce their full running-error-bound
algorithm or claim its validated stopping guarantees.

All scalars, roots, and correction terms remain in the requested `Real`:
fp64, long double, or fp128. There is no hidden cast through double, extra
precision backend, or eigenvalue/LAPACK dependency. Compensation requires
round-to-nearest arithmetic without fast-math reassociation. The error-free
transform identities also assume no underflow or overflow; nonfinite
arithmetic is reported, but the routine is not a general underflow detector.

Before iteration, set `x=s*y` with a power-of-two s at least one, chosen
from the coefficient root bound `max_k abs(c_k)^(1/(n-k))`. This keeps the
working polynomial moderate without rounding normal coefficients during
rescaling. If rounding s upward overflows, retain the preceding power of
two. Reject scaling that makes any nonzero coefficient subnormal or zero.
Automatic seeds lie on `|y|=4`, with phase offset `1/3` of a seed spacing;
the constant and linear cases need no automatic iteration. Optional caller
seeds are supplied in the original x coordinate.

## What a resolved result means

`status == PolynomialRootStatus::resolved` requires all of the following:

- Every root's `abs(Q)/sum_k abs(c_k)*abs(y)^k` is at most the requested
  tolerance, including the monic term. An exact zero numerator and zero
  magnitude sum have residual zero.
- Reconstructing `product_j (y-y_j)` agrees with every scaled coefficient
  to the requested tolerance, using denominator `max(1,abs(c_k_scaled))`.
  This catches duplicated solutions that individually satisfy Q=0 but
  omit another root. Factor ordering greedily maximizes the product of
  distances to already selected roots (Leja-style); multiplication also
  carries rounding corrections.
- Root estimates are separated relative to a numerical uncertainty model.
  With `B=sum abs(c_k)*abs(y)^k`, its derivative-magnitude polynomial B',
  and `u=16*n*epsilon`, use
  `rho=(abs(Q)+u*B)/(abs(Q')-u*B')`. A nonpositive denominator gives infinite
  uncertainty. Every pair must have distance greater than
  `8*(rho_i+rho_j)`, and the largest radius in original coordinates must
  be finite.

The last check deliberately includes a coefficient/roundoff scale even
when the computed residual is zero. The rho values and separation ratio
are **diagnostics, not rigorous inclusion disks**. In particular, a small
backward error need not mean a tolerance-sized root error. Coefficient
conditioning, input uncertainty, and near-multiple roots remain important.
Roots are not snapped onto the real axis or paired by averaging.

Other statuses distinguish exhausted iteration budgets, unresolved clusters,
stalled updates, and nonfinite/unrepresentable arithmetic. Results do not
silently retry with looser tolerances or another precision. Diagnostics
are reset to unavailable (infinity, or zero separation) if a failed partial
sweep leaves no completed assessment of the returned roots.

## Budget and usage

```cpp
using Real = long double;
std::array<Real, 2> coefficients{-Real{2}, Real{0}}; // Q(x)=x^2-2
auto result = bethe::detail::recover_polynomial_roots<Real>(coefficients);
if (result.status != bethe::detail::PolynomialRootStatus::resolved) {
    // Report unresolved root recovery; do not accept a physical state.
}
```

Defaults are `tolerance=128*epsilon`, `max_iterations=200`, `max_degree=64`.
The degree limit is checked before allocating working vectors. A zero
iteration budget assesses only the seeds (automatic or supplied); exact
resolved seeds may succeed without a sweep. The count records attempted
full-root sweeps, including a failed partial sweep. Per-sweep work is O(n²)
and scratch memory O(n). Root order is unspecified.

For the XXZ continuation's affine polynomial, convert each recovered root
with `z=center+coordinate_scale*x`, then `v=-(1-i*z)/(1+i*z)` and pass the
momentum factors to [CoordinateBetheWave](xxz-coordinate-wave.md). Check
the original Bethe equations and any singular factors separately; recovering
Q's roots is not proof of a physical Bethe vector or a sector minimum.

## Independent checks and current coverage

Native-precision tests include exact mixed real/complex roots, roots of
unity through degree 24, rescaling, explicit seeds and budgets, repeated
roots, close distinguishable roots, duplicated rootwise solutions, and
unrepresentable scaling. The close-root and degree-24 cases detect the
failures of uncorrected evaluation and naive factor multiplication.

The full connection is also tested, not just polynomial residuals:

- Recover continued odd-ring polynomials at N=5,7,9,11 in every sector
  with 2 through floor(N/2) roots, at `Delta=-0.2,-0.7,-0.9,-0.999`;
  build nonzero coordinate vectors and apply the periodic spin Hamiltonian.
- Reduce mixed phantom polynomials at `(N,r,p)=(9,3,1),(11,4,1),(11,3,2),
  (13,4,2),(13,3,3)`. Recovered finite vectors satisfy the twisted Hamiltonian;
  their nonzero dressed vectors satisfy the periodic Hamiltonian at the
  same energy. The test caches amplitudes only to avoid repeated evaluation.

These are numerical validations of the selected branches, not a theorem
that every regular reduced state has a nonzero phantom lift. General
mixed/singular-state acceptance and sector-minimum tracking remain open;
public solver domains and automatic scan acceptance are unchanged.
