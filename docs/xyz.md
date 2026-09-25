# XYZ: numerical foundation and implementation plan

**Status: not yet a spectrum solver or frontend.** The first checkpoint adds
native-precision elliptic functions, coupling normalization, and a regular-root
energy expression. It does not select or certify eigenstates, nor claim
ground-state or full-spectrum coverage.

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

1. Implement regular elliptic Bethe residuals and analytic Jacobians, with
   explicit quasi-period/root-sum conditions. Roots shifted by periods must
   be handled together with the associated phase parameter, not independently.
2. Establish the even-chain ground branch, including a verified mapping from
   the existing XXZ rapidity convention for continuation. Compare energies
   against independent physical-Hamiltonian diagonalization for several sizes
   and anisotropies; residual convergence alone is insufficient.
3. Handle singular/bound-pair solutions separately before claiming spectrum
   completeness. The modern chiral construction explicitly warns about missing
   states; use [Baxter's formulation](../CITATIONS.md#baxter-1973) as well.
4. Expose only validated branches through a dedicated frontend with the shared
   CLI, output tables, precision controls, and citation registry. Follow with
   excited families, coupling inversion/axis mappings, and odd rings as
   separately validated capabilities.
