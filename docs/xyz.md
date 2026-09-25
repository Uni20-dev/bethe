# XYZ: even periodic ground branch

**Status: ground-branch library and `bethe-xyz-pbc` frontend; no excited-spectrum coverage.**
Native fp64, long-double and fp128 implementations cover the symmetric regular-root
branch of even periodic chains with real $`0\lt \eta \lt 1`$ and rectangular `tau=i*t`,
$`t \gt 0`$. Independent small-chain diagonalization validates the ground-state selection
on the tested parameter grid; convergence is not a general completeness proof.

## Convention and initial target

Use spin operators $`S =\sigma /2`$ and

```math
\begin{aligned}
H&=\sum_j\left(J_xS_j^xS_{j+1}^x+J_yS_j^yS_{j+1}^y+J_zS_j^zS_{j+1}^z\right),\\{}
J_x&=\frac{\theta_4(\eta)}{\theta_4(0)},\qquad
J_y=\frac{\theta_3(\eta)}{\theta_3(0)},\qquad
J_z=\frac{\theta_2(\eta)}{\theta_2(0)},\\{}
\theta_j(u)&=\vartheta_j(\pi u\mid it),\qquad t>0.
\end{aligned}
```

The parametrization is equation (2) of
[Zhang–Klümper–Popkov (2024)](../CITATIONS.md#zhang-klumper-popkov-2024).
Their Hamiltonian uses Pauli matrices, so their energies are four times ours.
L=2 includes both periodic bonds. Real eta and rectangular `tau=i*t` give the
initial Hermitian family. Generic XYZ does not conserve total Sz; roots are
not labeled as a fixed-spin-flip sector.

The first solver target is the even periodic chain, starting with regular
ground-state roots and independent small-chain Hamiltonian checks. The
broader catalogue target remains finite-size energies and excitations; neither
two-site identities nor the trigonometric limit alone meet that target.

## Command-line use

```sh
bethe-xyz-pbc 16 --eta 0.4 --t 0.7 --roots
bethe-xyz-pbc 32 --eta 0.7 --t 2 --precision long-double
bethe-xyz-pbc 16 --eta 0.4 --t 0.7 --roots --json xyz.json \
  --csv-table states=xyz.csv --tsv-table roots=roots.tsv
bethe-xyz-pbc --references
```

Both `--eta` and `--t` are required; this first interface uses elliptic
parameters, not an inversion from arbitrary Jx/Jy/Jz. Eta below/above 1/2 gives
positive/negative Jz; eta=1/2 is the anisotropic XY chain. Exactly eta=0 or 1
is excluded: use the XXX/XXZ tools for those endpoint limits.

The `states` table reports total/per-site energy, momentum, residual, iterations,
convergence and status. `--roots` adds a `roots` table with decimal labels I and
separate `lambda_real`/`lambda_imag` columns. Metadata records the elliptic
parameters, three exchange constants, spin convention, precision, solver
controls, provenance and CPU time via the common [output system](output.md).
References appear only with `--references`, not during ordinary calculations.

All enabled precisions (`fp64`, `long-double`, and MPLAPACK-backed `fp128`) parse
and solve natively. `--max-iterations 0` checks only the seed; it is exact for
two sites, but generally fails for larger rings. A failed solve exits 2 and
exports missing energies, momenta and root coordinates (JSON null, CSV/TSV
empty). Numerical diagnostics and labels remain. Invalid input exits 1 before
any output file is opened, including existing targets supplied with `--force`.
Odd rings, `--sz` sectors and `--excitations` are not supported.

## Implemented numerical layer

`bethe/detail/elliptic_theta.hpp` evaluates all four theta functions and their
first two derivatives at complex u. Derivatives are with respect to **u**, not
pi*u. Its `ThetaJet` stores a common real `log_scale`:

```math
\mathrm{actual\ derivative}[r]=e^{\mathrm{log\_scale}}\,\mathrm{derivative}[r],
\qquad r=0,1,2.
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

```math
\begin{aligned}
g(u)&=\frac{\theta_1(\eta)}{\theta_1'(0)}\frac{\theta_1'(u)}{\theta_1(u)},\\{}
E&=\frac N4g(\eta)+\frac12\sum_j[g(\lambda_j-\eta/2)-g(\lambda_j+\eta/2)].
\end{aligned}
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
2. Add coupling inversion/axis mappings, with explicit conventions and tests.
3. Handle singular/bound-pair solutions separately before claiming spectrum
   completeness. The modern chiral construction explicitly warns about missing
   states; use [Baxter's formulation](../CITATIONS.md#baxter-1973) as well.
4. Extend the frontend only with separately validated excited families and
   odd-ring capabilities; the existing frontend deliberately exposes ground
   states only.

## Ground-state library API

```cpp
#include <bethe/xyz.hpp>
auto state = bethe::xyz::ground_state(16, 0.4, 0.7);
if (state.converged) { /* *state.energy is the total energy */ }
```

The solver fixes $`M =N /2`$, `xi=0`, and pairs imaginary roots $`\lambda =i \,x`$
with their negatives (including zero for odd M). Thus their sum is exactly
zero. Positive roots remain ordered in $`0\lt x \lt t /2`$; no independent period
folding is performed. With centered integer/half-integer labels I, it solves

```math
\begin{aligned}
N\phi_{\eta/2}(x_j)-\sum_{k\ne j}\phi_\eta(x_j-x_k)&=2\pi I_j,\\{}
\phi_a(x)&=2\arg\theta_1(a+ix).
\end{aligned}
```

These are the regular equations (46)–(48) of the cited paper on this branch.
An analytic reflection-reduced Jacobian feeds the same physical-domain damped
Newton driver used by the continuum models. The reported residual is
$`\max (\lvert F \rvert)/(N \,(1-\eta))`$: the extra scale prevents accepting a small raw
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
