# Scaling Lee–Yang model: periodic finite-volume ground state

**Status: native fp64/long-double/fp128 periodic ground-state TBA library
implemented and checked against an independent oracle; frontend is next.** This is a continuum
field-theory calculation, not a finite spin-chain or RSOS Hamiltonian solver.

## Physics and normalization

The non-unitary minimal model M(2,5) has c=-22/5 and lowest conformal weight
h_min=-1/5. Its massive integrable perturbation has one particle of mass m.
Our first target is the periodic circle of circumference L, with velocity 1,
positive m and L, and no defects. The source is
[Bajnok–el Deeb–Pearce (2015)](../CITATIONS.md#bajnok-el-deeb-pearce-2015),
especially Eqs. (1)–(3), (133), (151)–(152), (212), and (218).

Writing r=mL, the source-free ground-state TBA and energy are

```text
epsilon(theta) = r cosh(theta) + integral K(theta-theta') log(1+exp(-epsilon(theta'))) dtheta'
K(x) = -phi(x)/(2 pi) = sqrt(3)/(2 pi) cosh(x)/(sinh(x)^2+3/4)
Y(r) = L E_C(L) = -r/(2 pi) integral cosh(theta) log(1+exp(-epsilon(theta))) dtheta
```

All integrals run over the real line. The sign is important: K is positive.
E_C is the bulk-subtracted finite-size energy, not an absolute extensive
vacuum energy. No bulk constant is added in this first implementation.
The ultraviolet limit is Y -> -pi/15 and therefore
`c_eff(r)=-6Y/pi -> 2/5`. This is **not** the central charge: the relation is
`c=c_eff(UV)+24*h_min=-22/5`. The lowest state belongs to the nontrivial
conformal sector, not the identity sector. Excited-state TBA needs additional
sources and their quantization conditions; changing particle number in this
source-free equation is not sufficient.

## Numerical formulation

The implementation solves for
the bounded correction u=epsilon-r*cosh(theta), using even parity to work on
the positive half-line. The folded kernel is K(theta-y)+K(theta+y).
For large arguments use t=exp(-abs(x)) and

```text
K(x) = sqrt(3)/pi * t*(1+t^2)/(1+t^2+t^4).
```

This avoids overflow from squared hyperbolic functions. Evaluate the
logarithm with a stable softplus routine. Since K integrates to 1 and the
solution has epsilon>=r>0, the continuum fixed-point map on nonnegative u
has sup-norm Lipschitz constant at most 1/(1+exp(r))<1/2. A discretization
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

The direct omitted energy tail has a useful bound: for cutoff B>0,
`abs(Y_tail) <= coth(B)*exp(-r*cosh(B))/pi`. This bounds the omitted integral
using epsilon>=r*cosh(theta); it does **not** alone bound the response of
the retained solution to the omitted tail. The cutoff verification must
account for that response too.

## Reproducible independent oracle

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 scripts/reference_lee_yang.py --self-test
```

NumPy and SciPy are development-only dependencies. The oracle covers
1e-6<=r<=50 and checks orders 12 versus 24 and a cutoff enlarged by 2.
Its self-test checks kernel normalization, positivity and parity, the UV
golden-ratio plateau, the UV effective charge, the leading large-volume
term `Y ~ -r*K_1(r)/pi`, invalid inputs and iteration-budget exhaustion.
K_1 here is the modified Bessel function, not the scattering kernel.
Mesh/cutoff differences are diagnostics, not certified error bounds.

Representative independently computed fp64 values (use about 1e-12 absolute
tolerance when comparing Y, not all printed digits as exact references):

| r | Y=L E_C |
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

## Next checkpoints

`bethe-lee-yang-vacuum` will follow the existing sine-Gordon finite-volume
naming and shared CLI/table/run-context facilities. Known CFT constants must
be labelled theory data, not fitted results. The identity-sector state and
other excited-state source terms, then boundaries or defects, are separate
extensions; they must not silently reuse the source-free ground-state equation.
