# XYZ: even periodic ground branch

**Status: library ground-branch solver; no frontend or excited-spectrum coverage.**
Native fp64, long-double and fp128 implementations cover the symmetric regular-root
branch of even periodic chains with real `0<eta<1` and rectangular `tau=i*t`,
`t>0`. Independent small-chain diagonalization validates the ground-state selection
on the tested parameter grid; convergence is not a general completeness proof.

## Convention and initial target

Use spin operators `S=sigma/2` and

```text
H = sum_j (Jx Sx_j Sx_(j+1) + Jy Sy_j Sy_(j+1) + Jz Sz_j Sz_(j+1)).
Jx=theta4(eta)/theta4(0), Jy=theta3(eta)/theta3(0), Jz=theta2(eta)/theta2(0)
theta_j(u) = vartheta_j(pi*u | i*t), t>0.
```

The parametrization is equation (2) of
[Zhang–Klümper–Popkov (2024)](../CITATIONS.md#zhang-klumper-popkov-2024).
Their Hamiltonian uses Pauli matrices, so their energies are four times ours.
L=2 includes both periodic bonds. Real eta and rectangular `tau=i*t` give the
initial Hermitian family. Generic XYZ does not conserve total Sz; the future
frontend must not label roots as a fixed-spin-flip sector.

The first solver target is the even periodic chain, starting with regular
ground-state roots and independent small-chain Hamiltonian checks. The
broader catalogue target remains finite-size energies and excitations; neither
two-site identities nor the trigonometric limit alone meet that target.

## Implemented numerical layer

`bethe/detail/elliptic_theta.hpp` evaluates all four theta functions and their
first two derivatives at complex u. Derivatives are with respect to **u**, not
pi*u. Its `ThetaJet` stores a common real `log_scale`:

```text
actual derivative[r] = exp(log_scale) * derivative[r],  r=0,1,2.
```

Keep this representation for quotients: `theta_ratio` reconstructs only the
ratio, including derivative ratios. Computing unscaled theta values first can
underflow even when the required ratio is perfectly representable.

The evaluator uses Fourier series for t>=1 and a Poisson-transformed Gaussian
series for t<1. Argument reduction restores quasi-periodic factors and their
derivatives. Paired sine/sinh expressions and a local even-function expansion
preserve tiny arguments; compensated sums retain native fp64, long-double,
and fp128 arithmetic. No intermediate conversion to double is used.
Definitions and identities are documented by
[DLMF §20.2](https://dlmf.nist.gov/20.2) and
[§20.7](https://dlmf.nist.gov/20.7).

This internal interface requires finite arguments and finite t>0. It rejects
period reductions beyond its scalar/integer resolution and unrepresentable
scaled results or quotients. Components smaller than the representable range
of the common scale can round to zero; the representation is not arbitrary
precision and does not promise relative accuracy for such components.
Very large log scales also limit absolute reconstruction accuracy. General
complex tau is not supported.

## XYZ normalization hooks

`bethe::xyz::couplings(eta,t)` returns the three real coefficients.
`candidate_energy<Real>(N,roots,eta,t)` evaluates the regular expression for
even N>=2, N/2 supplied roots, and real 0<eta<1:

```text
g(u) = theta1(eta)/theta1'(0) * theta1'(u)/theta1(u)
E = N*g(eta)/4 + sum_j [g(lambda_j-eta/2)-g(lambda_j+eta/2)]/2.
```

This follows equations (8), (22), (26), and section VI of the cited paper with
the spin normalization above. **It is not a solve operation:** no Bethe
equation, root sum, or physical-state condition is checked. It keeps the complex
result, rejects poles, and does not regularize singular bound pairs.

Tests use independent 90-digit theta values, dual representations, differentiated
period identities, tiny arguments, zeros, and large/small t. Couplings are
checked against high-precision values and XXX/XXZ/XY limits. All four two-site
energy expressions agree with the Bell-state eigenvalues of the physical
Hamiltonian in the declared normalization.

## Remaining checkpoints

1. Extend robustness at extreme eta/t and larger sizes. The present direct seed
   avoids a separate continuation solver; an XXZ continuation path may help
   difficult parameters. Roots shifted by periods must be handled together
   with their phase parameter, not independently.
2. Add the ground-branch frontend using the shared CLI and output facilities.
3. Handle singular/bound-pair solutions separately before claiming spectrum
   completeness. The modern chiral construction explicitly warns about missing
   states; use [Baxter's formulation](../CITATIONS.md#baxter-1973) as well.
4. Expose only validated branches through a dedicated frontend with the shared
   CLI, output tables, precision controls, and citation registry. Follow with
   excited families, coupling inversion/axis mappings, and odd rings as
   separately validated capabilities.

## Ground-state library API

```cpp
#include <bethe/xyz.hpp>
auto state = bethe::xyz::ground_state(16, 0.4, 0.7);
if (state.converged) { /* *state.energy is the total energy */ }
```

The solver fixes `M=N/2`, `xi=0`, and pairs imaginary roots `lambda=i*x`
with their negatives (including zero for odd M). Thus their sum is exactly
zero. Positive roots remain ordered in `0<x<t/2`; no independent period
folding is performed. With centered integer/half-integer labels I, it solves

```text
N phi_(eta/2)(x_j) - sum_(k!=j) phi_eta(x_j-x_k) = 2*pi*I_j
phi_a(x) = 2 arg theta1(a+i*x).
```

These are the regular equations (46)–(48) of the cited paper on this branch.
An analytic reflection-reduced Jacobian feeds the same physical-domain damped
Newton driver used by the continuum models. The reported residual is
`max(abs(F))/(N*(1-eta))`: the extra scale prevents accepting a small raw
residual simply because the equations flatten near eta=1. Cancellation there
can still cause stagnation; use higher precision or inspect a failure, rather
than interpreting retained diagnostic roots as a solved state.

`energy` is absent on failure. Status distinguishes iteration limit, stalled
iteration and numerical range/precision failure. Invalid physical inputs throw.
Default tolerance is 32 native epsilons; zero iterations still allows the
exact two-site solution. Momentum is pi for odd M and zero for even M, from
the paired-root translation factors. Sz is not a quantum number of this solver.

Regression tests compare against independently assembled full spin Hamiltonians
for N=2,4,6,8 at 60 combinations of eta and t, verify the original multiplicative
complex equations and Jacobian, and compare the t=100 trigonometric limit with
the existing XXZ solver through N=32. Both signs of Jz and the XY point are
included, with equation/limit tests in every enabled native precision.
