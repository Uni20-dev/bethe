# Equal-coupling Bose–Fermi mixture

The [worked tutorial](tutorials/bose-fermi.md) compares repulsion and selected
compositions, with plots, downloadable roots and Lieb–Liniger checks.

The library `bethe/bose_fermi.hpp` implements periodic ground states of an
equal-mass gas containing scalar bosons and spinless fermions. Bose–Bose and
Bose–Fermi contact strengths must be equal and repulsive. This is a distinct
graded model, not another color of the SU(kappa) Fermi gas. The frontend is
`bethe-bose-fermi-pbc`.

## Command-line use

```sh
bethe-bose-fermi-pbc --bosons 2 --fermions 3 --length 5 --c 1 --roots
bethe-bose-fermi-pbc --bosons 4 --fermions 0 --length 4 --c 2 --json bosons.json
bethe-bose-fermi-pbc --bosons 1 --fermions 1 --length 1 --c 1e-40 --precision fp128
bethe-bose-fermi-pbc --bosons 2 --fermions 3 --length 5 --c 1 --precision digits-100 --roots
bethe-bose-fermi-pbc --references
```

Both populations, circumference, and coupling are required. The defaults are
fp64 and ground-state calculation; there is no excitation mode yet.
`--max-iterations` and `--tolerance` expose the library controls described below.
The common [output options](output.md) support screen presentation, CPU timing,
metadata, and JSON/CSV/TSV files, including streaming and `--no-retain`.
Literature is printed only with `--references`, not ordinary help.

Alongside fp64, platform-native long double and optional fp128, the frontend
accepts explicit MPFR precision (`bits-N` or `digits-N`) when
[MPFR scalars](building.md#enable-arbitrary-precision) are
configured. This covers the existing mixed, pure-species, free and vacuum
paths. More arithmetic bits do not extend the supported shells or guarantee
convergence within the iteration budget.

`states` contains total energy, momentum/index, residual, iterations, and
convergence status. With `--roots`, interacting calculations add `charge_roots`
and `auxiliary_roots`; the pure-boson reduction has an empty auxiliary table.
Free calculations instead add `free_modes`, with species, integer mode, and
physical k. Bosons each occupy mode zero. Half-integer labels are decimal in
exports. Use `--json result.json` for all tables, or named destinations such as
`--csv-table charge_roots=charge.csv`.

Failed solves leave energy, physical momentum, and numerical roots missing
(JSON null, CSV/TSV empty); requested Bethe labels remain available. Exit codes
are 0 for convergence, 2 for numerical failure, and 1 for invalid input or
output errors. Parameters and supported shells are checked before opening
output files, including with `--force`.

## Hamiltonian and supported sectors

In units $`\hbar ^{2}/(2m)=1`$, circumference $`\ell`$, the first-quantized Hamiltonian is

```math
H=-\sum_j\partial_j^2+2c\sum_{i\lt j}\delta(x_i-x_j),\qquad c\ge0.
```

Antisymmetry removes the contact interaction between identical fermions.
For $`c \gt 0`$ and both species present, the current branch requires **odd N_f**;
$`N_{b}`$ may be odd or even. Pure bosons, pure fermions, vacuum, and $`c =0`$ accept
arbitrary counts. Unequal masses/couplings, attraction, other interacting
periodic shells, excitations, boundaries, and thermodynamics are not supported.

The equations and normalization follow
[Imambekov–Demler, Annals of Physics (2006)](../CITATIONS.md#imambekov-demler-2006-applications),
equations (3), (28)–(34); the original short report is
[Phys. Rev. A (2006)](../CITATIONS.md#imambekov-demler-2006).
For N=N_b+N_f charge roots k and M=N_b auxiliary roots lambda:

```math
\begin{aligned}
\ell k_j+\sum_a2\arctan\!\left(\frac{2(k_j-\lambda_a)}c\right)&=2\pi I_j,\\{}
\sum_j2\arctan\!\left(\frac{2(\lambda_a-k_j)}c\right)&=2\pi J_a,\\{}
E&=\sum_jk_j^2.
\end{aligned}
```

Our $`J`$ is the negative of the paper's auxiliary logarithmic label. There is
**no auxiliary–auxiliary scattering term**.

Our parity audit gives `2I = M (mod 2)` and `2J = N (mod 2)` from the rational
equations. The centered choices $`I_{j} =j -(N -1)/2`$, $`J_{a} =a -(M -1)/2`$ satisfy both
precisely when N_f is odd. We reject even interacting mixed shells rather
than treating a centered solution with incompatible parity as periodic.
This branch has zero total momentum. At c=0, even fermion populations have
two lowest filled shells; we choose the positive-momentum representative,
`momentum_index=N_f/2`, with `P=2*pi*momentum_index/ell`.

## C++ API

```cpp
#include <bethe/bose_fermi.hpp>

auto state = bethe::bose_fermi::ground_state(2, 3, 5.0, 1.0);
// N_b, N_f, ell, c; double, long double, or enabled Uni20 fp128.
if (state.converged) {
    auto energy = *state.energy;
    auto const& k = state.momenta;
    auto const& lambda = state.auxiliary;
}
```

`SolverOptions<Real>` is the shared solver control: default residual tolerance
32 epsilon and 10000 accepted Newton updates. Zero updates only checks the
seed. Failures return `iteration_limit`, `stalled`, or `precision_limit`, with
no energy. Root arrays on failure, when available, are diagnostics only.
Invalid parameters and unsupported shells throw before solving.

Runtime precision requires options and physical inputs at the same finite
precision; exact-only, unset and mixed-precision inputs are rejected:

```cpp
auto p = uni20::Precision::bits(256);
bethe::bose_fermi::SolverOptions<uni20::mpreal> options(p);
auto state = bethe::bose_fermi::ground_state(
    2, 3, uni20::mpreal("5", p), uni20::mpreal("1", p), options);
```

The shared continuum Newton driver uses normal Uni20 dispatch for MPFR
updates, with CPU LU available when MPLAPACK is disabled. No global precision or exponent-range setting is changed.

The pure-boson branch delegates to Lieb–Liniger. Free states store boson zero
momenta followed by the occupied fermion modes in `momenta`; exact integer
fermion labels are in `free_fermion_modes`. These are not interacting charge
labels. The mixed branch returns ordered real charge and auxiliary roots and
their separate half-integer labels. Pure-boson labels are Lieb–Liniger labels;
no auxiliary roots are introduced for this reduction.

For independent compositions or couplings, capture each request and its options
by value in the existing [bounded Async grid](async-execution.md#observable-grids).
Each solve owns its Newton scratch and runs its dependent iterations sequentially.
Wait on the final grid result; failed solves retain their diagnostics and missing
energy. The same grid can contain different explicit precisions. No additional
scheduler or model-specific execution wrapper is needed.

## Numerical method and verification

Our mixed solver uses dimensionless roots $`q =\ell \,k`$, $`g =c \,\ell`$, and auxiliary
coordinates scaled by `max(1,g)`. Reflection reduction enforces exact paired
roots and explicit central zeros. The analytic Jacobian couples only the two
root families. It uses the shared rational scattering kernel, compensated
sums, normal Uni20 LU dispatch, and physical-domain damped Newton driver. Small-g
rank phases are canceled analytically before summing complementary angles;
charge residuals are divided by `max(sqrt(g),abs(q))` there, and otherwise
equations are scaled by N. Returned physical roots are checked again after
dimensional conversion. Unresolved subnormal or overflowing energies are
not published.

Tests include original complex rational equations, numerical differentiation
of the Jacobian, scaling, free-shell energies, weak-coupling slopes, strong
fermionization, and native-precision two-particle energies. One fermion agrees
with the Lieb–Liniger ground energy; one boson agrees with the corresponding
Gaudin–Yang impurity. A separate momentum-space contact Hamiltonian for two
bosons and three fermions supplies decreasing variational upper bounds as
the cutoff grows. Those finite-cutoff energies are not claimed to be exact.

MPFR tests exercise 128, 256 and 333 bits, including original complex equations,
pure-boson and one-fermion Lieb–Liniger reductions, free energies, failed solves,
finite-precision central zeros and serial/concurrent agreement. Independent
512-bit two-body bisection checks accuracy growth beyond binary128. At
$`c=10^{-60}`$ and length one, 333-bit arithmetic resolves
$`E/(2c)=1-c/12+O(c^2)`$; binary128 cannot resolve that correction. The mixed
$`N_b=2,N_f=3`$ calculation also resolves the first-order energy shift $`14c`$
above its nonzero free Fermi sea at that coupling. Input, root and state exports
retain the selected precision, including free boson zeros.
