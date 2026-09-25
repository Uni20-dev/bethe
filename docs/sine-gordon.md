# Sine-Gordon: bulk-subtracted finite-volume vacuum

**Status: native-precision kernel, vacuum-energy library and `bethe-sine-gordon-vacuum` frontend implemented.** This is a field
theory calculation, not a classical sine-Gordon PDE evolution or a finite-site
spin-chain diagonalization.

## Running the vacuum solver

```sh
bethe-sine-gordon-vacuum --length 1 --p 2
bethe-sine-gordon-vacuum --mass 2 --length 0.5 --p 0.5 --json vacuum.json
bethe-sine-gordon-vacuum --length 1 --p 1 --precision fp128 --csv vacuum.csv
```

`--length` and `--p` are required; `--mass` defaults to one. All real inputs
are parsed directly in the chosen precision: fp64 (default), long-double or
fp128 (with MPLAPACK). The last example is the free Dirac check, not an
interacting calculation. Interacting fp128 runs can take minutes.

The single `vacuum` table reports `casimir_energy` (E_C), `scaling_function`
(Y), `effective_central_charge`, separate convergence diagnostics and work
counters. These are not a total absolute energy or an energy per lattice site.
The metadata records physical conventions, both numerical contour shifts,
input controls, provenance and CPU time. See [output and exports](output.md)
for JSON/CSV/TSV, streaming and overwrite protection. Literature appears with
`--references`, not on each ordinary run.

`--tolerance` is an absolute tolerance in Y, default $`262144\,\epsilon`$ in the
selected precision. `--contour-shift` and `--initial-cutoff` override the
automatic contour and rapidity cutoff choices described below.
`--initial-intervals` (64) and `--max-intervals` (2048) control rapidity
resolution; `--max-iterations` (10000) is the total nonlinear-update budget
across all meshes and both contours. `--max-cutoffs` (3) limits rapidity
cutoff attempts per contour. Independently, `--max-kernel-evaluations`
(100000), `--max-kernel-levels` (16) and `--max-fourier-cutoffs` (64) limit
each Fourier kernel table construction. Kernel accuracy is assigned from the
vacuum error budget, not a separate user tolerance.

Exit status is 0 on convergence, 1 for invalid input/output errors, and 2 for
an incomplete solve. On status 2 all three physical observables are missing
(JSON null, empty CSV/TSV cells); diagnostics remain available. Unknown error
estimates may be infinite. Numerical validation precedes opening output files,
including when `--force` is given.

## Physical convention

Use the canonically normalized kinetic term `(1/2)(partial phi)^2`, coupling
$`0 \lt  \beta ^{2} \lt  8\,\pi`$, and the parameter

```math
p=\frac{\beta^2}{8\pi-\beta^2}.
```

Then p>1 is repulsive, 0<p<1 attractive, and p=1 the free massive Dirac point.
In Rutkevich's normalization, his beta squared is our $`\beta ^{2}/(8\,\pi)`$,
his xi is p, and his gamma is `pi/(1+p)`. State these conventions explicitly;
the different meanings of beta must not be mixed.

Take the soliton mass M as the physical input scale, set velocity and hbar to
one, and use $`u =M \,L`$. The observable is the **bulk-subtracted** vacuum
energy $`E_{C} =E_{0} -L \,e_{\mathrm{bulk}}`$, or the dimensionless scaling function $`Y =L \,E_{C}`$.
No nonuniversal absolute bulk energy is inferred. In the attractive regime M
still denotes the soliton mass, not necessarily the lightest breather mass.

The starting point is [Destri–de Vega](../CITATIONS.md#destri-de-vega-1992).
For explicit shifted-contour formulas and checks on normalization we use
[Rutkevich (2020)](../CITATIONS.md#rutkevich-2020), equations (22), (68), (69).

## Implemented kernel

```cpp
#include <bethe/sine_gordon.hpp>
auto g = bethe::sine_gordon::scattering_kernel(std::complex<double>(0.8, 0.4), 2.7);
if (g.converged) {
  // *g.value is the complex kernel, not an energy.
}
```

The kernel is

```math
G_p(z)=\int_0^\infty\frac{\cos(kz)}{2\pi}
\frac{\sinh[(p-1)\pi k/2]}{\sinh(p\pi k/2)\cosh(\pi k/2)}\,dk.
```

Its Fourier representation is used only inside
`|Im(z)| < pi*min(1,p)`. Invalid couplings, nonfinite arguments, points outside
this strip and nonpositive tolerances throw. We enforce the same strip at p=1
even though the identically zero kernel can itself be continued further.

`KernelOptions<Real>` provides an absolute component tolerance (1024 native
epsilons), up to 100000 integrand evaluations, 16 quadrature levels, and 64
Fourier-cutoff attempts. All arithmetic stays in the chosen fp64, long-double
or fp128 type. The implementation reuses the common compensated tanh-sinh
quadrature; it does not introduce a double-only special-functions dependency.

The hyperbolic ratio is evaluated with `expm1` and decaying exponentials.
For complex arguments, growth from $`\cos (k \,z)`$ is combined with the Fourier
decay *before* either factor is evaluated. This avoids an overflowing cosh
multiplied by an underflowing ratio on otherwise legal contours. The same
formula retains the small multiplier for p only a few native epsilons from one.

Let `d=pi*min(1,p)-|Im(z)|>0`. Beyond a cutoff K, an absolute tail bound is

```math
\frac{e^{-dK}}{\pi d(1-e^{-p\pi K})}.
```

The returned `tail_bound` and `quadrature_error` are separate. The latter is a
successive-mesh estimate, not a rigorous enclosure. A kernel value is optional
and is published only when their sum meets the requested tolerance. Failure
statuses distinguish cutoff, quadrature and finite-range limitations. Unknown
error diagnostics remain infinite rather than appearing to be zero.

Independent checks include the exact transforms

```math
\begin{aligned}
G_{1/2}(z)&=-\frac1{2\pi\cosh z},\\{}
G_2(z)&=\frac{z}{2\pi^2\sinh z},\qquad G_2(0)=\frac1{2\pi^2},\\{}
G_1(z)&=0,\qquad \left.\frac{\partial G_p(0)}{\partial p}\right|_{p=1}=\frac18.
\end{aligned}
```

The first two are test references, not implementation shortcuts. Tests also
cover evenness, complex conjugation, a separate 90-digit quadrature reference,
overflow-prone Fourier factors, native near-free couplings, and failure budgets.

## Vacuum-energy API

```cpp
#include <bethe/sine_gordon_vacuum.hpp>
auto state = bethe::sine_gordon::vacuum_energy(1.0, 2.0, 0.7); // M, L, p
if (state.converged) {
  // *state.casimir_energy = E0-L*e_bulk
  // *state.scaling_function = L*state.casimir_energy
  // *state.effective_central_charge = -6*Y/pi
}
```

The library accepts finite positive M, L and p. It follows the untwisted vacuum
branch in both regimes; difficult parameters can exhaust a numerical budget.
It does not return the bulk energy density or an excited-state spectrum.

`VacuumOptions<Real>` controls:

- `tolerance`: absolute accuracy target in Y, default 262144 native epsilons.
  This is **not** a relative guarantee for exponentially small infrared energies.
- `initial_intervals=64`, `max_intervals=2048`: an even, uniform rapidity mesh,
  doubled until two successive energy differences are below tolerance/4.
  The explicit work ceiling is 16384 intervals; it is not a convergence promise.
- `max_iterations=10000`: total damped nonlinear updates across all meshes,
  cutoffs and both contours, rather than a hidden budget per retry.
- `max_cutoffs=3`: cutoff attempts per contour, each separately mesh-resolved.
  Consecutive cutoffs differ by one rapidity unit and must agree within tolerance/4.
- Optional `initial_cutoff` and `contour_shift`; defaults estimate a safe tail
  scale and use $`\eta =\pi \,\min (1,p)/4`$. The second independently resolved contour is
  $`3\,\eta /4`$, and must agree in Y within tolerance/2.
- `kernel`: Fourier evaluation/level/cutoff budgets per table. Its tolerance is
  assigned internally from the vacuum tolerance and rapidity cutoff.

`nonlinear_residual`, `kernel_error`, `mesh_error`, `cutoff_error` and
`contour_error` are reported separately. These mesh/cutoff comparisons are
numerical estimates, not rigorous interval bounds. `intervals` and `cutoff`
report the larger final resolved sizes of the two contours; work counters
include unsuccessful earlier attempts. A small nonlinear residual cannot by
itself produce `converged=true`.

All three observables are optional and remain absent on `kernel_limit`,
`iteration_limit`, `mesh_limit`, `cutoff_limit`, `contour_limit` or
`precision_limit`. Invalid parameters throw before iteration; unrepresentable
scaled lengths or energies fail without publishing values. A branch guard keeps
the logarithm argument in its open right half-plane. A transient pseudoenergy
with negative real part alone is not rejected: this occurs during perfectly
regular attractive ultraviolet iterations.

## Shifted-contour vacuum NLIE

For $`0 \lt  \eta \lt  \pi \,\min (1,p)/2`$, use

```math
\begin{aligned}
A(\theta)&=\log(1+e^{-\epsilon(\theta)}),\\{}
\epsilon(\theta)&=-iu\sinh(\theta+i\eta)
-\int G_p(\theta-t)A(t)\,dt
+\int G_p(\theta-t+2i\eta)\overline{A(t)}\,dt,\\{}
Y&=-\frac u\pi\operatorname{Im}\int\sinh(\theta+i\eta)A(\theta)\,d\theta.
\end{aligned}
```

The finite contour shift regularizes the counting-function logarithms without
crossing kernel poles. Do not move directly to eta=pi/2 in the attractive
regime. The chosen eta is numerical, not a physical model parameter: agreement
at two independently resolved shifts is an important check.

The implementation tabulates real and shifted kernels on a difference grid,
using shared native Gauss-Legendre quadrature and compensated sums. The Fourier
multiplier is computed once per quadrature point for all displacements;
short phase recurrences reduce repeated trigonometric calls. The standalone
tanh-sinh kernel provides an independent check of this tabulation. Its Fourier
tail bound is unchanged. Two successive Fourier-mesh agreements are required.

The unknown is the correction to the analytically known driving term, updated
with half damping. Vacuum parity `epsilon(-theta)=conj(epsilon(theta))` reduces
the convolution cost. The complex logarithm utility is shared with ASEP.
Convolution is currently quadratic in mesh size: a demanding fp128 solve can
take minutes in a Debug build. No FFT or lower-precision fallback is used.

## Validation and remaining scope

All precisions check the free-Dirac integral against an independent high-precision
reference, and the interacting p=2 result against an independent, real-valued
D3 TBA calculation with two separately refined Gauss meshes. The latter uses
the equations summarized in [Hegedus (2026)](../CITATIONS.md#hegedus-2026),
(2.6)-(2.8), with the magnon plateau integrated analytically.

Further tests cover attractive/repulsive couplings, mass-length scaling, altered
contours, the ultraviolet $`c_{\mathrm{eff}} \to 1`$ limit, and the leading soliton plus
breather wrapping correction at p=1/2. They also distinguish nonlinear,
kernel, rapidity-mesh and cutoff failures, and verify that failed states have
no published observables.

Frontend regressions additionally check native parsing/output, mass-length
scaling, budget forwarding, missing failed observables and shared export contracts.

This first target is the untwisted zero-topological-charge vacuum. Excited
states require additional source terms and branch/quantization bookkeeping;
large-volume Bethe–Yang equations alone are not an exact finite-volume solver.
