# Rational Gaudin central spin

[Overview](../README.md) · [Model catalogue](models.md) · [Bibliography](../CITATIONS.md)

`bethe-central-spin` finds the lowest energy at a **specified total spin
projection** for

```math
H=BS_0^z+\sum_{j=1}^{N_b}A_j\,\mathbf S_0\cdot\mathbf S_j.
```

Every spin, including the central spin 0, is spin 1/2. There are no bath
fields, bath–bath interactions, or energy shifts. Couplings may have either
sign but must be finite, nonzero and mutually distinct. The field B can be
positive, negative, or exactly zero. This is a star of interactions, not a
chain: PBC/OBC, lattice momentum and an energy density are not part of the
interface.

## First calculations

```sh
build/bethe-central-spin --couplings 1,0.7,0.3 --field 1 --sz 0
build/bethe-central-spin --couplings 1,-2,0.3,-0.1 --field 0 --sz 1/2
build/bethe-central-spin --couplings 1 --field 1 --sz 0 --precision fp128 --variables
```

`--sz` includes the central spin, and is required. For Nb bath spins its
allowed values run from `-(Nb+1)/2` to `(Nb+1)/2` in unit steps. A sector
minimum need not be the global ground state; compare sectors if that is
what the application needs. In particular, an antiferromagnetic central
spin with several bath spins is not generally a singlet at zero field.

The one-bath Sz=0 answer is

```math
E=-\frac A4-\frac{\sqrt{A^2+B^2}}2.
```

Thus A=B=1 gives E=-0.9571067811865475244..., and A=1, B=0 gives E=-3/4.
The polarized sectors instead have `E=sum A_j/4 +/- B/2`. An empty list,
`--couplings ''`, represents an isolated central spin and requires Sz=+/-1/2.

The report uses Uni20 presentation, retains the selected fp64/long-double/fp128
precision, and includes CPU time. `--references` displays literature from
the shared citation registry separately from help. There is no excitation scan,
rapidity output, wavefunction, form factor or dynamics API yet.

## Equations and conventions

The Hamiltonian, all-down-reference Bethe equations and energy are fixed by
[Faribault–Schuricht, Eqs. (2)–(6)](https://arxiv.org/html/1306.2541v2), setting
the bath field to zero. For M up spins and B>=0,

```math
\begin{aligned}
\epsilon_0&=0,\qquad\epsilon_j=-\frac1{A_j},\\
-2B+\sum_i\frac1{\lambda_\alpha-\epsilon_i}
-2\sum_{\beta\ne\alpha}\frac1{\lambda_\alpha-\lambda_\beta}&=0,\\
E&=\frac12\sum_\alpha\frac1{\lambda_\alpha}-\frac B2+\frac14\sum_j A_j.
\end{aligned}
```

Spin reversal maps `(B,Sz)` to $(-B,-S^z)$, without changing energy. The
implementation therefore solves at h=|B|, using `M=(Nb+1)/2+Sz` for B>=0
and the spin-reversed M for B<0.

Direct rapidities can collide or escape to infinity. Following the
quadratic eigenvalue-variable approach of
[Faribault et al.](https://arxiv.org/html/1103.0472v2), introduce

```math
\begin{aligned}
A_\star&=\max_j|A_j|,\quad q_0=0,\quad q_j=A_\star/A_j,\quad z_\alpha=-A_\star\lambda_\alpha,\\
t&=\frac{A_\star}{2h+A_\star},\quad p=1-t,\qquad
v_i=t\sum_\alpha\frac1{q_i-z_\alpha}.
\end{aligned}
```

Algebra applied to the rational equations in the convention above gives

```math
\begin{aligned}
v_i^2-pv_i-t\sum_{j\ne i}\frac{v_i-v_j}{q_i-q_j}&=0,\\
\sum_i v_i&=Mp,\\
E&=hv_0+\frac{A_\star}{2}v_0-\frac h2+\frac14\sum_j A_j.
\end{aligned}
```

These v variables are **not occupations or local magnetizations**, and can
be negative. They remain finite at exactly B=0, where t=1 and p=0. No small
field is substituted for a requested zero field. Numerically p is computed
separately from the field ratio, retaining a tiny nonzero field even when
t rounds to one. Unrepresentable ratios/pole separations are rejected.

The number constraint is enforced by eliminating an initially occupied
bath variable. Keeping the central variable independent preserves its tiny
high-field value. All quadratic equations enter a rectangular Newton
correction using native-precision QR; no equation is discarded, no normal
equations are formed, and long-double does not require a LAPACK backend.

## Which state is followed?

The starting point t=0 is infinite positive field. In any nonpolarized
sector the central spin is down, and the M largest A_j have bath spins up.
“Largest” is algebraic, not absolute value. Distinct couplings make that
product-state minimum unambiguous. It is **not** generally the occupation
seed which minimizes the Richardson pairing Hamiltonian, even though the
models share Gaudin eigenstates.

There is a useful state-selection argument specific to the star geometry.
A local z-axis phase on each bath spin makes every off-diagonal exchange
matrix element negative, for either sign of its A_j. For nonzero couplings,
the spin-configuration graph in each nonpolarized fixed-M sector is
connected. Perron–Frobenius therefore gives a unique lowest state in that
sector at every finite field: the ground branch cannot cross another
state in the same sector. Zero-field SU(2) degeneracies between different
magnetization sectors do not invalidate that argument.

Adaptive predictor/corrector continuation follows this branch from t=0 to
the requested t. It limits the distance from the predicted variables and
checks product-state variational bounds, an interaction-norm lower bound,
and the ground energy's field concavity. These are numerical safeguards,
not a claim that a residual alone certifies an eigenstate's ordering.
Polarized and isolated-spin cases are evaluated analytically.

## Budgets, failures and library use

```cpp
#include <bethe/central_spin.hpp>
std::vector<long double> a{1.L, 0.7L, 0.3L};
auto s = bethe::central_spin::sector_ground_state<long double>(
    a, 1.L, uni20::half_int{0});
if (s.converged) {
    auto energy = *s.energy;
}
```

`SolverOptions<Real>` has `residual_tolerance` (default 32 machine epsilons),
`max_iterations` (default 10000 attempted Newton corrections, including
retries) and `max_stages` (default 10000 attempted stages). The tolerance
is a polynomial backward residual, **not an energy-error bound**. Close
couplings can amplify energy error; use higher precision and compare
results when a problem is ill-conditioned.

Zero budget suppresses continuation. The seed then has no finite reached
field and no energy: `State::reached_field` and `State::energy` are empty
optionals. If a later stage fails, these contain the last accepted field
and its energy, never an energy silently attributed to the requested B.
`converged` requires reaching the target. CLI exit codes are 0 for success,
2 for an incomplete solve and 1 for invalid input/numerical exceptions.
Exact analytic cases still succeed with zero budgets.

`State` retains the requested Sz/field/couplings, the spin-reversal flag,
attempt counters, reached residual and number-constraint error. Its v
array is central first, then in the **original input coupling order**;
for negative B it describes the spin-reversed frame. The array is empty
for analytic polarized cases, which do not need continuation variables.

## Validation and next scope

Tests construct the spin Hamiltonian independently and compare every small
magnetization sector for positive, negative, mixed-sign and clustered
couplings, including zero and negative fields. Native-precision checks
cover the two-spin formula, the original one- and two-root rational equations,
the analytic Jacobian/tangent, spin reversal, input permutation, overall
energy scaling, the high-field seed and zero-field multiplet degeneracy.
A 23-spin bath exercises the solver beyond those small ED fixtures.
CLI tests cover precision-preserving output, citations, CPU time and
honest incomplete-solve reporting.

Repeated bath couplings need confluent/derivative Gaudin equations or
grouped bath-spin representations, with their multiplicities and sectors
handled explicitly. A zero coupling adds a decoupled spin and requires
minimization over its allowed projections. Neither case is silently
perturbed into the distinct-coupling problem. Higher local spins, general
linear combinations of Gaudin charges, excited states and observable
matrix elements remain separate future work.
