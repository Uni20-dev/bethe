# A regular-state test without extracting Bethe roots

[Odd-ring continuation and scans](xxz-odd-continuation.md) | [Wronskian diagnostic](xxz-wronskian.md)

The internal `check_regular_polynomial<Real>` in
[xxz_regularity.hpp](../include/bethe/xxz_regularity.hpp) tests a sufficient
**regular-case** nonzero-vector criterion using the solver's polynomial
coefficients. It supplements the quantum Wronskian and does not require its
generic-anisotropy assumption. Singular limiting states are a different
problem; they are not declared unphysical when this test cannot resolve them.

## The norm formula and our deduction

[Korepin's historical account](../CITATIONS.md#korepin-2009), section IV,
equation (23), expresses the scalar product of an on-shell XXZ Bethe vector
and its algebraic dual as a product of one-particle and scattering factors
times the determinant of the logarithmic Bethe-equation Jacobian. The roots
are taken distinct. A nonzero scalar product implies a nonzero ket even
without identifying the algebraic dual with its Hermitian conjugate.

Our deduction is that it is enough to verify, on an **exact solution**:

1. Conventional rapidities are finite and distinct modulo their period.
2. The one-particle factors and both pair-scattering factors are nonzero.
3. The full logarithmic Bethe Jacobian is nonsingular.

Then the regular algebraic Bethe construction produces a nonzero
eigenvector. This is a state-existence test, not a lowest-energy theorem.
Its numerical implementation is not an interval-arithmetic certificate:
finite residuals and rank indicators only test numerical versions of the
hypotheses, at the caller's precision and requested residual tolerance.

## Expressing those conditions in polynomial coordinates

The input is the monic polynomial F(x), with physical rapidity `z=o+s*x`.
Let `a=sqrt((1+Delta)/(1-Delta))`. The conventional multiplicative rapidity is

```text
exp(2*lambda) = (a+z)/(a-z).
```

Roots at z=+-a therefore require infinite-rapidity limits. The one-particle
singularities `lambda=+-i*gamma/2`, where `Delta=cos(gamma)`, map to z=+-i.
Homogeneous Horner evaluation tests F at these points without forming
large coordinate ratios. Each `*_margin` for a value is its magnitude
divided by the coefficient-magnitude bound evaluated at
`max(1,abs(x))`. The unit-radius floor also detects small endpoint values
without cancellation, such as `F(x)=x^M` near x=0; a cancellation-only
ratio would misleadingly return one. Real coefficients make the
two driving values at +-i conjugates, so only one needs evaluation.

The remaining root conditions are checked in the quotient algebra C[x]/F:

- Multiplication by F' is invertible exactly when F has distinct roots.
- Multiplication by the self-scattering-removed K_- is invertible exactly
  when it has no common root with F. At a root x_i,

```text
K_-(x_i) = product_(j!=i)
    [1+Delta-(1-Delta)*z_i*z_j-i*Delta*(z_i-z_j)].
```

Thus the second matrix detects vanishing pair-scattering factors. The
opposite-sign factor is covered by conjugation for a real polynomial.
`PolynomialBetheSystem::scattering_remainder` exposes K_- modulo F before
the driving factor is multiplied in; it uses the same divided-Horner
identity that removes self-scattering without a numerical division.

For the Jacobian test, use **all M coefficient equations**, not the reduced
M-1 equations and imposed momentum used in the continuation driver.
At an exact solution, evaluating their first variations at the roots is a
Vandermonde change of basis. The root-to-coefficient map is nonsingular for
distinct roots. Dividing the rootwise equations by their nonzero driving
and scattering factors, and changing from x to lambda, then gives the
logarithmic Bethe Jacobian. All these changes are invertible under the
preceding conditions. Therefore nonsingularity of the full coefficient
Jacobian is equivalent to the determinant condition needed here.

This deduction deliberately excludes collisions, infinite rapidities,
singular driving roots, and exact strings; a limiting construction is
required when a prefactor or change of variables vanishes. A root-of-unity
anisotropy by itself is **not** a reason to reject a regular finite-root
state. This differs from applying the generic-q Wronskian theorem outside
its stated hypotheses.

## Numerical rank, not a terminal linear solve

All arithmetic stays in fp64, native long double, or fp128. Complex
multiplication matrices use the real block representation

```text
[ Re(A)  -Im(A) ]
[ Im(A)   Re(A) ].
```

After normalizing the largest matrix entry to one, complete-pivot
elimination returns

```text
smallest_absolute_pivot / (matrix_order * max(1, element_growth)).
```

This is a rank-resolution indicator, **not a reciprocal condition number**
or a rigorous error bound. Zero/nonfinite pivots or arithmetic give an
unresolved result. Every value/rank margin must exceed `64*epsilon`, and
the full polynomial residual must satisfy the caller's tolerance (default
`32*epsilon`). The routine does not loosen tolerances to obtain a pass.

An exact singular matrix is an expected input here. Uni20's ordinary dense
solve uses a terminal singular-matrix error policy, which aborts under its
default configuration. This rank calculation does not invoke that solve
or mutate Uni20's process-wide error policy. The singular-matrix regressions
run under the unmodified default policy in all supported scalar types.

The result distinguishes `regular_on_shell`, `off_shell`,
`exceptional_or_unresolved`, `ill_conditioned`, and `nonfinite`. Invalid
arguments throw. A small endpoint or resultant margin is not proof of an
exact singularity: numerical resolution can also be inadequate. Uncomputed
rank margins stay zero and an unavailable residual stays infinite.

## Use in the sector scan and tests

An optional final `rotation` argument applies the same criterion to the
[constant-twist finite-root equations](xxz-phantom.md). It defaults to one;
the sector scan still uses the periodic, zero-twist check. This extension
tests the finite reduced state, not its singular phantom dressing.

Each equation-converged entry in the internal all-sector scan now carries
an optional `regularity` result alongside its Wronskian diagnostic.
`regular_states_complete` says whether every entry passed this regular
criterion. Neither diagnostic changes the candidate minimum or substitutes
for sector-minimum tracking. Failed continuations have neither diagnostic.

Tests cover all folded sectors through N=9 at five couplings; selected
N=13,17,21 branches at Delta=-0.97 and -0.999; and explicit repeated,
infinite, singular-driving, and exact-string roots. Root-of-unity tests
separate regular vacua/one-magnon states from infinite-root states.

An independent native-precision five-site check extracts the two roots
from an analytic quadratic only in the **test**, constructs the coordinate
Bethe wavefunction, verifies it is nonzero, and applies the spin Hamiltonian
directly. This covers both real and conjugate-pair roots. A separate affine
coordinate check compares K_- against the explicit pair factor above.

An [explicit spin-helix construction](xxz-spin-helix.md) now handles the
all-phantom collision separately. Public odd-ring integration still requires
a policy for the other exceptional states this regular criterion excludes,
together with sector-minimum tracking. A singular regular-root norm formula
alone is never a reason to discard a state.
