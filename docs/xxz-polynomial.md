# Polynomial XXZ equations across root collisions

[Negative-anisotropy work](xxz-negative.md) | [Model catalogue](models.md)

The internal building block
`bethe::xxz::detail::PolynomialBetheSystem<Real>` in
[xxz_polynomial.hpp](../include/bethe/xxz_polynomial.hpp) expresses the
periodic Bethe equations in coefficients rather than individually labelled
roots. It includes analytic coefficient derivatives, energy, and momentum.
It accepts `-1<Delta<=0`, including the free-fermion point for initialization.

This is **not yet a production ground-state solver**. The polynomial can
represent real roots and conjugate pairs, and the tests follow small-ring
sector minima through collisions. General continuation, numerical conditioning,
physical-state admissibility, and global sector selection still need work.
The public APIs and frontends remain unchanged.

## Starting equations and parity

Use the coordinate Bethe ansatz in
[Caux's notes](../CITATIONS.md#caux-xxz-coordinate), Eqs. `xxz.be` and `xxz.e`,
at zero twist, J=1. We restore the ferromagnetic energy `N*Delta/4` to the
notes' shifted Hamiltonian. Set

```text
exp(i*k_j) = -(1-i*z_j)/(1+i*z_j),
D(z,w) = 1+Delta-(1-Delta)*z*w,
f_±(z,w) = D(z,w) ± i*Delta*(z-w),
sigma = (-1)^(N-M-1).
```

Clearing the scattering denominators gives, for each root,

```text
(1+i*z_j)^N * product_(l!=j) f_-(z_j,z_l)
  - sigma*(1-i*z_j)^N * product_(l!=j) f_+(z_j,z_l) = 0.
```

The `(-1)^N` in the momentum parametrization is essential on odd rings.
This equation comes from the periodic coordinate ansatz; it does not rely
on an even-ring ground-state identification theorem. Clearing denominators
alone is not a proof of admissibility at singular configurations.

## Removing the self-scattering factor before reduction

Write the monic root polynomial with **real** coefficients in ascending order:

```text
Q(z) = product_j (z-z_j) = z^M + sum_(r=0)^(M-1) c_r*z^r.
```

Real coefficients allow conjugate pairs. They do not represent arbitrary
non-self-conjugate complex root configurations. Define

```text
A_±(z) = 1+Delta ± i*Delta*z,
B_±(z) = (1-Delta)*z ± i*Delta,
D_self(z) = A_±(z)-B_±(z)*z = 1+Delta-(1-Delta)*z^2.
```

The product over **all** roots is `B_±^M Q(A_±/B_±)`. At a root it includes
`D_self`, whereas the Bethe equation excludes self-scattering. Retaining that
factor would introduce false roots at
`z=±sqrt((1+Delta)/(1-Delta))`.

Instead form the polynomial divided difference

```text
K_±(z) = [B_±(z)^M*Q(A_±(z)/B_±(z)) - B_±(z)^M*Q(z)] / D_self(z).
```

This is a polynomial identity, not a numerical division by a possibly zero
quantity. For a simple root `z_j`, it evaluates to the required product
over `l!=j`. The equations become divisibility by Q of

```text
F(z) = (1+i*z)^N*K_-(z) - sigma*(1-i*z)^N*K_+(z).
```

At multiple roots, divisibility also imposes derivative conditions. These
retain information lost by merely evaluating the same root equation twice.
They describe the continuous coefficient branch being followed, not a general
classification of all singular Bethe solutions.

## Evaluation and derivatives

All computations are in the quotient algebra modulo Q. For either sign,
initialize `K=0`, `Z=1`, `P=1`. For `r=0,...,M-1`, using the previous values
on each right-hand side, update

```text
K <- A*K + P*Z,
Z <- z*Z + c_(M-r-1),
P <- B*P.
```

Reduce each result modulo Q. This divided-Horner recurrence constructs K
without forming high-degree polynomials or dividing by `D_self`. Repeated
multiplication by `1+i*z` gives
`G=(1+i*z)^N*K_- mod Q`. Since Q is real, the other term is its
coefficient-wise complex conjugate. Thus the M raw residuals are

```text
Im(G_r), sigma=+1;
Re(G_r), sigma=-1.
```

The factor of two (or 2i) does not affect their zeros. Directional forward
differentiation supplies an analytic real M-by-M Jacobian. It differentiates
the reduction modulo Q too: treating that reduction as constant would omit
part of the derivative. No finite differences are used by the module.

The diagnostic `norm` is the largest absolute raw residual divided by
`max_r |G_r|`. A vanishing denominator raises an error; it is not interpreted
as convergence. The Jacobian differentiates the **raw residual**, not this
normalization. This coefficient diagnostic is neither a rootwise phase error
nor an energy error estimate.

## Observables without root extraction

The same coefficients give

```text
E = N*Delta/4 + M*(1-Delta) + 2*Im[Q'(i)/Q(i)],
exp(i*P_total) = Q(-i)/Q(i) = conjugate(Q(i))/Q(i).
```

Evaluating Q and Q' uses compensated sums of the coefficients with exact
powers of i. Complex division is scaled to avoid squaring large components.
`Q(i)=0` is rejected as an observable pole. Root extraction is unnecessary
for these quantities and would be ill-conditioned at a collision.

`momentum_defect` returns `|exp(i*N*P_total)-1|`, a necessary translation
check. It is **not sufficient** to prove a nonzero physical Bethe vector or
to identify a sector minimum. A recorded 21-site failed continuation has
raw coefficient residuals below 1e-7 but a momentum defect above 0.1;
the regression explicitly retains this counterexample.

## Independent collision benchmark

For N=5, M=2, momentum `2*pi/5`, write `C=cos(pi/5)`. In a relative-separation
basis `r=1,2`, direct reduction of the spin Hamiltonian gives

```text
H - E_ferro = [[-Delta, C], [C, -2*Delta-C]],
e = (-3*Delta-C-sqrt((Delta+C)^2+4*C^2))/2,
u = 4*C^2/(2*C^2-2*Delta-e),
Q(z) = z^2 + tan(pi/5)*u*z + 1-u.
```

At `Delta=-C`, the two roots coincide at `-tan(pi/10)`. They are real above
that value and a conjugate pair below it. The tests check the discriminant,
energy, polynomial residual, quantized momentum, and nonsingular coefficient
Jacobian on both sides and at the collision, in all three precisions.

Other tests compare the reduced equations with direct scattering products at
arbitrary real/conjugate roots, check analytic Jacobians against finite
differences, and follow every `M<=N/2` sector for N=3,5,7,9 through 41 couplings
from 0 to -0.999. Exact spin-basis diagonalization supplies energies and joint
energy/translation checks. The continuation loop is test-only: it is not an
API or a claimed large-chain algorithm.

The implementation checkpoint passes the full 412-test GCC 13 Debug suite
with fp128 and the 285-test Clang 20 Release suite without MPLAPACK. Of these,
19 and 13 cases respectively test this building block and its audit examples.

## Remaining implementation work

1. Add adaptive continuation with a shared iteration budget, checks of the
   intended momentum branch, and explicit diagnostics when conditioning or
   admissibility prevents certification. Investigate a better-conditioned
   polynomial basis or factored variables for larger chains; do not accept
   a small raw coefficient residual as a physical-state certificate.
2. Establish sector-minimum tracking through singular/root-of-unity cases
   and compare sectors for the global ground state. The existing
   smallest-|Sz| rule is not valid for all negative-coupling odd rings.
3. Integrate the real and polynomial representations into the existing
   ground-state API/frontends, reporting the representation and residual
   convention explicitly. Only then widen the documented input range.
