# Open XXZ chains with free ends

[Back to the overview](../README.md)

The `bethe-xxz-obc` front end combines anisotropy with free ends. It uses
the same scalar precision and reporting controls as the other programs, but
its equations include reflected scattering and a boundary phase. Simply
doubling the periodic chain length would miss that phase.

## Model and first calculations

```math
\begin{aligned}
H&=\sum_{i=0}^{N-2}\left(S_i^xS_{i+1}^x+S_i^yS_{i+1}^y+\Delta S_i^zS_{i+1}^z\right),\\{}
J&=1,\quad h=0,\quad \Delta\gt -1\ \text{(ground states)},\quad N\ge2,
\quad\text{no boundary fields}.
\end{aligned}
```

```sh
build/bethe-xxz-obc 64 --delta 0.5
build/bethe-xxz-obc 65 --delta 0.5 --sz -1/2 --roots
build/bethe-xxz-obc 16 --delta 0.75 --sectors --precision long-double
build/bethe-xxz-obc 16 --delta 3 --roots
build/bethe-xxz-obc 33 --delta -0.9 --sz -1/2 --roots
# In a binary128-enabled build:
build/bethe-xxz-obc 64 --delta 0.999999999999999999999999 --precision fp128
```

The default is the ground state, represented by Sz=0 for even N and Sz=1/2
for odd N. `--sz` selects a sector minimum; `--sectors` returns one minimum
per magnetization sector. Negative sectors use spin reversal. There is no
lattice momentum for free ends, and generic XXZ states are not classified
by total spin S or SU(2) multiplets. There is no `--spin`, `--spinons`, or
boundary-selection option.

For $`-1\lt \Delta \lt 0`$, the solver uses positive hyperbolic lambda coordinates
and Newton iteration. `--roots` prints both lambda and
$`z =s \,\tanh (\lambda)`$, with $`s =\sqrt{(1+\Delta)/(1-\Delta)}`$. The reported residual
uses **rank-subtracted equations divided by N*s**, not the nonnegative
solver's $`\max \lvert F \rvert /(2\,N)`$. This scaling keeps the stopping test meaningful as
Delta approaches -1, where the z roots shrink. See the
[negative-anisotropy guide](xxz-negative.md) for the equations and derivation.
Exactly Delta=-1 and lower anisotropies are not included by this extension.

For $`\Delta \gt 1`$, even zero-magnetization ground states need a distinguished
boundary root. `--roots` prints its inverse square $`y =1/z_{B} ^{2}`$ and logarithmic
distance w separately from the bulk real roots. Negative y means that z_B is
imaginary; it is not a real rapidity or a failed solve. See the
[massive boundary-root guide](xxz-open-massive.md) for the coordinate crossing,
finite-size deviation, residual convention, and continuation diagnostics.

The default arithmetic is fp64; long-double and optional fp128 retain native
precision in parameters, roots, energies, gaps, and output. See the shared
[CLI controls](command-line.md) for tolerances, CPU time, and presentation.

## Real-root excitations

Both excitation scans and specified quantum-number lists still require
$`0\le \Delta \le 1`$; negative and massive ground-state support do not extend that family.

```sh
build/bethe-xxz-obc 16 --delta 0.5 --excitations 10 --sz 1
build/bethe-xxz-obc 8 --delta 0.5 --excitations all --sz 2 --roots
build/bethe-xxz-obc 8 --delta 0.75 --quantum-numbers 1,3 --roots
```

`--excitations COUNT|all` scans the supported finite-real-root family at fixed
Sz, including its sector minimum. Without `--sz`, the scan uses Sz=1 for
even N and Sz=1/2 for odd N. This is different from the ground-state default.
Gaps are relative to the **global ground energy**, not the sector minimum.

For M=N/2-|Sz|, the supported labels are distinct increasing positive integers:

```math
\begin{aligned}
1&\le I\le N-M,\qquad I\lt I_\infty,\\{}
I_\infty&=N-M+1-\frac{(N-2M+1)\gamma}{\pi},\\{}
\gamma&=\arccos\Delta.
\end{aligned}
```

The strict second bound comes from sending one rapidity to infinity with the
others finite. The implementation uses the equivalent expression
$`I_{\mathrm{infinity}} =(N +1)/2+(N -2\,M +1)\,\arcsin (\Delta)/\pi`$. For 0<Delta<1 it also excludes
a conservative $`32\,\epsilon \,\max (1,I_{\mathrm{infinity}})`$ band below this threshold.
At Delta=0, the exact integer condition is $`2\,I \lt N +1`$; at Delta=1, the existing
open XXX window `1..N-M` is used exactly. The vacuum M=0 is a single empty
configuration.

This is a **restricted family, not a complete Sz spectrum**. Complex strings,
infinite rapidities, and other root branches are excluded, including the
XXX infinite-root descendants at Delta=1. Even at the free-fermion point
Delta=0 the scan is only the finite-real-rapidity subset. For even N at Sz=0
it contains only the ground configuration. `all` means all configurations
in this family, not all physical eigenstates.

The family size changes with anisotropy: N=8, Sz=2 has 6 candidates at Delta=0,
10 at Delta=1/2, and 15 at Delta=1. At N=8, Sz=1, there are 4 candidates at
Delta=1/2 and 10 just above it, as additional roots move in from infinity.
The output records the actual window and candidate count.

Every configuration is solved before retaining the lowest COUNT energies;
`--max-candidates` (default 10000) limits this exhaustive work, also with `all`.
Preflight counting is allocation-free and overflow-checked. Unconverged
candidates are excluded and the ordering is marked incomplete. A failed ground
reference leaves valid candidate energies available but gaps unavailable.
Either failure gives exit status 2.

`--quantum-numbers` solves one explicit label list in the same window, inferring
M and Sz=N/2-M; use `none` for the polarized vacuum. It is mutually exclusive
with `--sz`, `--sectors`, and `--excitations`. The latter cannot be combined
with `--sectors` either.

## Library usage

```cpp
#include <bethe/xxz_open.hpp>
#include <bethe/xxz_excitations.hpp>

namespace obc = bethe::xxz::open;
auto ground = obc::ground_state<long double>(64, 0.5L);
auto sector = obc::sector_ground_state<long double>(65, 0.5L,
                                                   uni20::half_int::parse("-1/2"));
auto sectors = obc::sector_ground_states<long double>(16, 0.75L);
auto massive = obc::ground_state<long double>(16, 3.0L);
auto negative = obc::ground_state<long double>(33, -0.9L);
auto scan = obc::real_excitations<long double>(
    16, 0.5L, uni20::half_int{1}, {.count = 10, .max_candidates = 10000});
auto count = obc::real_excitation_count(8, 0.5L, uni20::half_int{2});
auto window = obc::real_quantum_number_window(8, 0.5L, uni20::half_int{2});
obc::QuantumNumbers numbers{uni20::half_int{1}, uni20::half_int{3}};
auto state = obc::solve_real<long double>(8, 0.75L, numbers);
```

The three ground-state functions return `obc::GroundState<Real>` (or a vector
of it). Alongside delta, bulk scaled roots and labels, Sz, spin reversal,
energy, and convergence diagnostics, this includes `log_rapidities`, a
`GroundResidualConvention` tag, an optional `boundary_root`,
`root_delta`, and a `GroundSolveStatus` distinguishing convergence, budget
exhaustion, and a stalled line search. Bulk arrays exclude the distinguished
root and its label. `log_rapidities` stores the solver's lambda coordinates
only for negative Delta; do not reconstruct them from rounded z values.
The convention tag distinguishes `logarithmic_phase` (0<=Delta<=1),
`negative_rank_scaled` (-1<Delta<0), and `massive_regularized` (Delta>1).
At $`\Delta \le 1`$, `boundary_root` is empty. At $`0\le \Delta \le 1`$, all numerical
results retain the original real-root solver's values. Always check
`converged`: on failure, the reported energy is an unconverged estimate,
with the residual evaluated at the returned coordinates and requested Delta.

`solve_real` and excitation scans retain `obc::RealState<Real>` with no boundary
root. Ground states no longer implicitly convert to `RealState`: code with an
explicit old return type should use `GroundState` or `auto`, and account for
`boundary_root` when counting roots. Neither type has **momentum members**.
`solve_real` also accepts shared `bethe::SolverOptions<Real>` and an optional
initial-root span. Enumeration, bounded ranking, and convergence bookkeeping
are shared with the other models. Numerical code has no presentation dependency.

## Boundary equations and limiting cases

This section describes the all-real solver for $`0\le \Delta \le 1`$; the
[massive equations and regularization](xxz-open-massive.md) are separate.
The [negative-Delta scaled equations](xxz-negative.md) are also separate;
the excitation labels below must not be extrapolated to that regime.

Use the same scaled coordinate as periodic XXZ:
$`z =\tanh (\lambda)/\tan (\gamma /2)`$, not the conventional rapidity lambda. Physical
roots are positive and satisfy $`(1-\Delta)\,z ^{2}\lt 1+\Delta`$. Writing
$`p =1+\Delta`$, $`q =1-\Delta`$, and $`c =q /p`$, the equations and energy are

```math
\begin{aligned}
A_{ij}&=\frac{\Delta(z_i-z_j)}{p-qz_iz_j},&
B_{ij}&=\frac{\Delta(z_i+z_j)}{p+qz_iz_j},\\{}
F_i&=4N\arctan z_i+4\arctan(cz_i)-2\pi I_i
-2\sum_{j\ne i}(\arctan A_{ij}+\arctan B_{ij}),\\{}
E&=\frac{(N-1)\Delta}{4}-\sum_i\frac{p-qz_i^2}{1+z_i^2}.
\end{aligned}
```

Both direct and reflected self-scattering terms are excluded. The boundary
term `4*atan(c*z)` is essential even with zero boundary fields. At Delta=0
scattering vanishes but the boundary term does not: the roots become
$`z_{i} =\tan (\pi \,I_{i} /[2\,(N +1)])`$, with energy `-sum_i cos(pi*I_i/(N+1))`.
These are the open-chain standing waves, with denominator **N+1**, not N.

The sector minimum fills I=1,...,M. Simultaneous updates start at zero unless
finite nonnegative in-branch guesses are supplied. The implementation moves
the boundary correction to the right using
$`\arctan (z)-\arctan (c \,z)=\arctan (2\,\Delta \,z /(p +q \,z ^{2}))`$; the resulting driving denominator
is $`2\,(N +1)`$. This keeps updates on the supported positive finite branch and
solves Delta=0 in at most one update. At **exactly** Delta=1 it delegates to
the open XXX solver, preserving its roots $`z =2\,\lambda_{\mathrm{XXX}}`$ and diagnostics.
Nearby anisotropies are not snapped to the endpoint.

Convergence uses $`\max \lvert F \rvert /(2\,N)`$, with default tolerance 32 times the selected
type's epsilon. It is not an energy-error bound. Both residual and energy
describe the returned iterate, even on budget exhaustion. Work per update
is O(M^2), or O(M) at Delta=0, with O(M) state storage. Retaining all sectors
costs O(N^2) storage. Excitation scans retain
O(min(COUNT,candidates)*M+N), including their ground reference.

## Validation

Independent small-chain exact diagonalization checks sector minima and every
enumerated excitation, including multiplicities, through N=9. Analytic checks
include $`E_0 (N =2)=-1/2-\Delta /4`$, its root
$`z =\sqrt{(1+\Delta)/(3-\Delta)}`$, and
$`E_0 (N =3)=-(\Delta +\sqrt{\Delta ^{2}+8})/4`$. Irrational-energy and gap tests discriminate
against narrowing higher precision to double. Further checks cover independent
hyperbolic residuals, the XX and XXX limits, both sides of an infinity threshold,
larger chains, exhausted budgets, and full-precision plain/pretty output.

The negative-Delta public API is checked separately against all spin-basis
sector minima through N=9 at four couplings. Tests also check spin reversal,
the global ground-state selection, the residual-convention tags, independent
scaled residuals at successful and failed iterates, and native-precision
two-/three-site energies. Endpoint tests reach $`\Delta =-1+128\,\epsilon`$ on
even and odd open chains. CLI tests retain all printed z/lambda digits in
both wide and narrow terminal reports, and keep negative excitations rejected.

Related: [periodic XXZ](xxz.md), [free-end XXX](open-chains.md), and
[references and provenance](../CITATIONS.md).
