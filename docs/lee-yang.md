# Scaling Lee–Yang model: periodic finite-volume ground state

**Status: source conventions and independent numerical oracle established;
native C++ solver and frontend are the next checkpoint.** This is a continuum
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

The following choices and checks are our implementation design. Solve for
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
solver should use a uniform rapidity mesh and cached difference/sum kernels,
giving an independent discretization. Reuse scalar math and compensated
summation from the existing solvers, but do not inherit sine-Gordon contour
parameters: this is a real scalar TBA. Keep nonlinear, mesh and cutoff
diagnostics separate. Require successive mesh agreements and an enlarged
cutoff verification, with explicit budgets and missing energies on failure.
Convergence checks should target absolute error in dimensionless Y, avoiding
an implicit change of tolerance when m and L are rescaled.

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

## Planned native/public contract

- `bethe/lee_yang.hpp`: scalar kernel and periodic ground-state TBA, native
  fp64/long-double/fp128 throughout, with no double fallback.
- Input m and L; record r=mL. Return optional Y, E_C=Y/L and c_eff(r), plus
  separate numerical diagnostics and work counters. Known CFT constants
  must be labelled theory data, not fitted numerical results.
- `bethe-lee-yang-vacuum`: one frontend, following the existing sine-Gordon
  finite-volume naming convention and shared CLI/table/run-context facilities.
- Validate independent values above, UV/IR limits, mass/length rescaling,
  native-precision round trips, and individual numerical-budget failures.
- Later checkpoints: the identity-sector state and other excited-state
  source terms, followed by boundaries or defects when needed. They must
  not silently reuse the source-free ground-state equations.
