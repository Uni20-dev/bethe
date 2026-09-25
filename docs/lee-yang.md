# Scaling Lee–Yang model: periodic finite-volume energies

**Status: native fp64/long-double/fp128 periodic ground-state TBA library and
frontend implemented and checked against an independent oracle.** This is a continuum
field-theory calculation, not a finite spin-chain or RSOS Hamiltonian solver.
The regular spin-zero one-particle library and frontend are also available
for $`5\le mL\le30`$; continuation to small volume is not yet implemented.

## Command line

```sh
bethe-lee-yang-vacuum --length 1
bethe-lee-yang-vacuum --mass 2 --length 0.5 --precision fp128 --json vacuum.json
bethe-lee-yang-vacuum --length 0.001 --csv vacuum.csv --tsv vacuum.tsv
bethe-lee-yang-vacuum --references
```

The particle mass defaults to 1; circumference is required. The numerical
options below are available as `--tolerance`, `--initial-intervals`,
`--max-intervals`, `--max-iterations`, `--max-cutoffs`,
`--max-kernel-products`, and `--initial-cutoff`.

The `vacuum` table contains the bulk-subtracted `casimir_energy`,
`scaling_function` Y, and computed `effective_central_charge`, followed by
nonlinear residual/error, mesh and cutoff errors, the direct tail bound,
cutoff, intervals, iteration/cutoff/work counts and convergence status.
Metadata labels $`c=-22/5`$, $`h_{\min}=-1/5`$ and $`c_{\mathrm{eff}}^{\mathrm{UV}}=2/5`$ as **theory data**;
these are not inferred from a finite-r calculation. CPU time and numerical
controls accompany all shared screen/file outputs.

Successful solves exit 0. Numerical failures exit 2 and publish no physical
observables; unavailable diagnostics are also missing (JSON null or empty
delimited cells), not infinity strings. Invalid input exits 1 before opening
files, including existing `--force` targets. The [common output options](output.md)
provide independent JSON/CSV/TSV exports, streaming and optional table retention.

### One-particle level and gap

```sh
bethe-lee-yang-excited --length 5
bethe-lee-yang-excited --mass 2 --length 2.5 --precision fp128 --json levels.json
bethe-lee-yang-excited --length 10 --csv-table gap=gap.csv --tsv-table source=source.tsv
bethe-lee-yang-excited --references
```

This tool computes one zero-momentum excited level, not an arbitrary
excitation scan. Both mass and length must be positive and $`5\le mL\le30`$.
The same numerical controls apply to each state separately; their tolerance
defaults to $`65536\,\epsilon`$. `--max-root-iterations` (default 256 per grid)
controls the excited source quantization. The default work budget is
1000000000 kernel products **per state**, not a combined allowance.

Three tables are always present:

- `levels`: vacuum and one-particle bulk-subtracted energies and scaling
  functions, with independent convergence diagnostics and work counters.
- `source`: $`\beta`$, its displacement from $`\pi/6`$, quantization residual,
  propagated source-error estimate, root-iteration count and status.
- `gap`: $`E_{1,C}-E_{0,C}`$, `scaled_gap=L*gap`, estimated `gap_error` in energy
  units, and convergence status.

The energy gap is **not** the bulk-subtracted one-particle energy. The tool
subtracts the scaling functions before dividing by L. Its gap-error estimate
sums both verified cutoff errors (which already include nonlinear and mesh
errors) and subtraction roundoff, then divides by L. It is not a rigorous
error certificate. The excited table does not label its energy as an
effective central charge.

If either state fails, the tool exits 2 and leaves the gap/error fields
missing. A successfully computed level remains available. For example,
`--max-root-iterations 0` can leave a converged vacuum alongside an
unavailable excited level. Invalid input, including unsupported mL, exits 1
before opening output files even with `--force`. Shared CPU-time metadata
covers both solves; references appear only with `--references`.

## Physics and normalization

The non-unitary minimal model $`M(2,5)`$ has $`c=-22/5`$ and lowest conformal weight
$`h_{\min}=-1/5`$. Its massive integrable perturbation has one particle of mass $`m`$.
Our first target is the periodic circle of circumference L, with velocity 1,
positive m and L, and no defects. The source is
[Bajnok–el Deeb–Pearce (2015)](../CITATIONS.md#bajnok-el-deeb-pearce-2015),
especially Eqs. (1)–(3), (133), (151)–(152), (212), and (218).

Writing $`r=mL`$, the source-free ground-state TBA and energy are

```math
\begin{aligned}
\epsilon(\theta)&=r\cosh\theta+
\int K(\theta-\theta')\log\!\left(1+e^{-\epsilon(\theta')}\right)\,d\theta',\\{}
K(x)&=-\frac{\phi(x)}{2\pi}
=\frac{\sqrt3}{2\pi}\frac{\cosh x}{\sinh^2x+3/4},\\{}
Y(r)&=LE_C(L)=-\frac r{2\pi}\int\cosh\theta\,
\log\!\left(1+e^{-\epsilon(\theta)}\right)\,d\theta.
\end{aligned}
```

All integrals run over the real line. The sign is important: $`K`$ is positive.
$`E_C`$ is the bulk-subtracted finite-size energy, not an absolute extensive
vacuum energy. No bulk constant is added in this first implementation.
The ultraviolet limit is $`Y\to-\pi/15`$ and therefore
$`c_{\mathrm{eff}} (r)=-6Y /\pi \to 2/5`$. This is **not** the central charge: the relation is
$`c=c_{\mathrm{eff}}^{\mathrm{UV}}+24h_{\min}=-22/5`$. The lowest state belongs to the nontrivial
conformal sector, not the identity sector. Excited-state TBA needs additional
sources and their quantization conditions; changing particle number in this
source-free equation is not sufficient.

## Numerical formulation

The implementation solves for
the bounded correction $`u=\epsilon-r\cosh\theta`$, using even parity to work on
the positive half-line. The folded kernel is $`K(\theta-y)+K(\theta+y)`$.
For large arguments use $`t=e^{-\lvert x\rvert}`$ and

```math
K(x)=\frac{\sqrt3}{\pi}\frac{t(1+t^2)}{1+t^2+t^4}.
```

This avoids overflow from squared hyperbolic functions. Evaluate the
logarithm with a stable softplus routine. Since K integrates to 1 and the
solution has $`\epsilon\ge r\gt 0`$, the continuum fixed-point map on nonnegative $`u`$
has sup-norm Lipschitz constant at most $`1/(1+e^r)\lt 1/2`$. A discretization
must still check its own quadrature and roundoff; the nonlinear residual
alone is not an energy error estimate.

The oracle uses composite Gauss–Legendre rules supplied by SciPy. The native
solver uses a uniform trapezoidal rapidity mesh and cached difference/sum
kernels, giving an independent discretization. It reuses stable thermal
factors, scalar math and compensated summation from the existing solvers,
without sine-Gordon contour parameters: this is a real scalar TBA. Nonlinear,
mesh and cutoff diagnostics are separate. Two successive mesh agreements and
an enlarged cutoff verification are required, with explicit budgets and
missing energies on failure. The tolerance targets absolute error in Y,
avoiding an implicit change of tolerance when m and L are rescaled.

The direct omitted energy tail has a useful bound: for cutoff $`B\gt 0`$,
$`\lvert Y_{\mathrm{tail}}\rvert\le\coth(B)e^{-r\cosh B}/\pi`$. This bounds the omitted integral
using $`\epsilon\ge r\cosh\theta`$; it does **not** alone bound the response of
the retained solution to the omitted tail. The cutoff verification must
account for that response too.

## Reproducible independent oracle

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 scripts/reference_lee_yang.py --self-test
```

NumPy and SciPy are development-only dependencies. The oracle covers
$`10^{-6}\le r\le50`$ and checks orders 12 versus 24 and a cutoff enlarged by 2.
Its self-test checks kernel normalization, positivity and parity, the UV
golden-ratio plateau, the UV effective charge, the leading large-volume
term $`Y\sim-rK_1(r)/\pi`$, invalid inputs and iteration-budget exhaustion.
$`K_1`$ here is the modified Bessel function, not the scattering kernel.
Mesh/cutoff differences are diagnostics, not certified error bounds.

Representative independently computed fp64 values (use about 1e-12 absolute
tolerance when comparing Y, not all printed digits as exact references):

| $`r`$ | $`Y=LE_C`$ |
| ---: | ---: |
| 0.001 | -0.20943937151285794 |
| 0.01 | -0.20942648592833846 |
| 0.1 | -0.20835015785667577 |
| 0.5 | -0.19017376405211875 |
| 1 | -0.15320688011013006 |
| 2 | -0.0812547352032044 |
| 5 | -0.006408441082412908 |
| 10 | -0.000059359269102363974 |

## Native library

```cpp
#include <bethe/lee_yang.hpp>
auto state = bethe::lee_yang::ground_state<long double>(1.0L, 1.0L);
if (state.converged) {
  auto Y = *state.scaling_function;
  auto E_C = *state.casimir_energy;
  auto c_eff = *state.effective_central_charge;
}
```

Inputs are mass m and length L; both and their product must be finite,
representable and positive. There is no double fallback. All three physical
outputs are optional and absent on failure. `kernel(x)` also exposes the
positive scattering kernel for finite real x.

`Options<Real>` defaults and meanings:

| Option | Default | Meaning |
| --- | ---: | --- |
| tolerance | 8192 epsilon | Absolute target in Y, required in (0,1) |
| initial_intervals | 32 | Starting intervals on [0,B] |
| max_intervals | 2048 | Maximum intervals per grid; hard cap 8192 |
| max_iterations | 1000 | Fixed-point updates per grid; zero still tests the seed |
| max_cutoffs | 3 | Total cutoff trials, including the initial one |
| max_kernel_products | 200000000 | Total folded kernel-times-logarithm terms across all grids |
| initial_cutoff | automatic | Positive B; verification increases it by 1 |

The automatic cutoff is at least 2 and makes the free driving energy at the
boundary larger than a tolerance-dependent logarithmic target. The driving
term is evaluated through logarithms to avoid overflow from an intermediate
cosh when r is tiny. Grids double their interval count; at least three
converged grids are needed for two successive agreements at each cutoff.
A single cutoff trial can never publish an energy.

The nonlinear error estimate uses the discrete kernel row-sum bound and the
energy's sensitivity to the correction u, with a native roundoff floor.
`nonlinear_residual` itself remains a sup-norm equation residual, not an
energy error. `mesh_error` includes adjacent-grid changes and their nonlinear
errors; `cutoff_error` includes the change between independently mesh-refined
cutoffs, numerical errors and the direct tail bound. These estimates are not
rigorous interval-arithmetic certificates.

`Status` distinguishes `converged`, `iteration_limit`, `mesh_limit`,
`cutoff_limit`, `work_limit`, and `precision_limit`. Invalid inputs throw
`std::invalid_argument`. Failure states retain input mass/length, diagnostics and work
counters, not a publishable approximate energy. A diagnostic not yet evaluated
is infinity. `iterations` sums updates over all grids; `intervals` and `cutoff`
describe the last attempted grid. Storage is O(N), work O(N^2) per fixed-point
evaluation, with the work budget checked before each evaluation.

Tests cover the independent values above, UV/IR limits, exact rescaling at
fixed mL, native round trips and agreement between different initial meshes
and cutoffs at the selected precision. Every failure status has a regression,
including a too-small cutoff, exhausted work and unrepresentable mL.

## First excited state: independent infrared oracle

**Native regular one-particle library and frontend implemented.**
`scripts/reference_lee_yang_excited.py` solves the spin-zero one-particle
branch for $`5\le r=mL\le30`$, independently of the production C++ vacuum solver.
It reuses the Python reference quadrature and vacuum calculation for gaps.

The source is [Dorey–Tateo (1996)](../CITATIONS.md#dorey-tateo-1996),
Eqs. (2.3)–(2.7). With $`\theta_{0} =i \,\beta`$, define

```math
\begin{aligned}
S(z)&=\frac{\sinh z+i\sqrt3/2}{\sinh z-i\sqrt3/2},\\{}
s(\theta,\beta)&=-2\log\lvert S(\theta+i\beta)\rvert,\\{}
\epsilon(\theta)&=r\cosh\theta+s(\theta,\beta)+\int K(\theta-y)L(y)\,dy,\\{}
L(\theta)&=\log\!\left(1+e^{-\epsilon(\theta)}\right),\\{}
0&=r\cos\beta-\log S(2i\beta)+\int K(i\beta-y)L(y)\,dy,\\{}
Y_1&=LE_{1,C}=2r\sin\beta-\frac r{2\pi}\int\cosh\theta\,L(\theta)\,d\theta,\\{}
\frac{\mathrm{gap}}m&=\frac{Y_1-Y_0}{r}.
\end{aligned}
```

Here the integrals are over the full real line, and K uses its analytic
continuation for the quantization equation. Conjugate symmetry makes that
integral real. $`E_{1,C}`$ and $`E_{0,C}`$ share the vacuum's bulk subtraction; $`E_{1,C}`$ alone
is not the excitation gap. The regular source pair collides with scattering
singularities near $`r=2.53`$; these equations are not a UV continuation algorithm.

Our oracle brackets $`\beta`$ between $`\pi/6`$ and $`\pi/4`$, away from that collision.
For $`r\ge5`$ the bare driving term plus source stays positive on this bracket;
this gives a contractive inner real-axis fixed-point problem. Brent's method
solves the outer quantization condition in $`\log (\beta -\pi /6)`$. In particular,
the denominator of $`S(2i\beta)`$ is evaluated as
$`2\,\cos (\pi /3+d)\,\sin (d)`$, $`d=\beta-\pi/6`$, rather than subtracting nearly equal
sines. The upper r limit keeps the fp64 energy correction resolvable.

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 scripts/reference_lee_yang_excited.py --self-test
```

Checks compare quadrature orders 12 and 24 and an enlarged cutoff, for both
energy and source displacement. They also check source reality, the complex
kernel's real-axis reduction, nonlinear/quantization residuals, exhausted
iterations and domain rejection. Large-r checks use
$`\beta-\pi/6\sim\sqrt3\,e^{-\sqrt3r/2}`$ and
$`E_{1,C}/m-1\sim3e^{-\sqrt3r/2}`$; the latter still has subleading integral
corrections. Agreement diagnostics are not certified error bounds.

| $`r`$ | $`Y_1=LE_{1,C}`$ | $`(E_1-E_0)/m`$ |
| ---: | ---: | ---: |
| 5 | 5.146781165267865 | 1.0306379212700556 |
| 6 | 6.076715733314169 | 1.0132130508527144 |
| 8 | 8.019437400300701 | 1.0024791202717358 |
| 10 | 10.004545725250384 | 1.0004605084519487 |
| 20 | 20.00000175700884 | 1.0000000880377056 |

Use about 2e-12 absolute tolerance in Y1 for comparisons, not every displayed
digit as an exact reference. This oracle is a development dependency only.

## Native one-particle library

```cpp
#include <bethe/lee_yang_excited.hpp>
auto excited = bethe::lee_yang::one_particle<long double>(1.0L, 5.0L);
auto vacuum = bethe::lee_yang::ground_state<long double>(1.0L, 5.0L);
if (excited.converged && vacuum.converged) {
  auto gap = *excited.casimir_energy - *vacuum.casimir_energy;
}
```

The library has the oracle's explicit $`5\le mL\le30`$ domain. It computes native
fp64, long-double or fp128 values without a double fallback. The optional
outputs are `scaling_function` $`Y_1`$, `casimir_energy` $`E_{1,C}`$, $`\beta`$ and
`pole_displacement` $`\beta-\pi/6`$. All are absent on failure. There is deliberately
no excited-state `effective_central_charge` field: this is an energy level,
not the vacuum scaling function.

`OneParticleOptions<Real>` inherits the vacuum controls, with default
tolerance $`65536\,\epsilon`$, work budget 1000000000 kernel products, and an
additional `max_root_iterations=256` per grid. The absolute tolerance still
targets Y1. Root iterations sum over grids; inner iterations and kernel
products include every source-root evaluation. The work budget is checked
before each integral-operator evaluation. An exhausted root or inner
iteration budget returns `iteration_limit`; other statuses match the vacuum.

The shared real TBA iteration now includes a source-dependent contraction
estimate and sensitivity of the continued-kernel integral. The outer solve
uses safeguarded secant updates in log displacement, falling back within its
sign bracket. It retains the last resolved slope when changes in the
quantization residual become comparable to inner-solve roundoff, rather
than differentiating numerical noise. Acceptance checks both residual and
successive source positions. `quantization_residual` is dimensionless;
`source_error` estimates its propagated absolute Y1 error, including the
change of source position, using a conservative response factor and the
local secant slope. These are numerical estimates, not certified bounds.

The shared mesh/cutoff controller checks the displacement as well as energy.
The direct energy-tail bound is multiplied by the maximum of exp(-source)
on the allowed beta bracket; cutoff enlargement also checks the response
of the coupled source root. `nonlinear_error` includes source-location and
energy-roundoff contributions on an accepted grid. Failed solves retain
diagnostics, not provisional physical outputs. A gap calculation must check
both states and propagate both error estimates.

Tests cover independent oracle energies, analytic source/kernel identities,
default-tolerance agreement between different meshes and cutoffs in all
three precisions, mass/length rescaling, every failure status, and a 51-point
fp64 scan of the supported interval. The vacuum regressions exercise the
same refactored iteration and refinement code.

## Next checkpoints

Continuation through the source collision toward the UV identity sector,
moving particles, additional particles, boundaries and defects remain
separate extensions. None should silently reuse the source-free equation.
