# Periodic XXZ chains

[Back to the overview](../README.md)

The anisotropy changes the equations and which quantum numbers belong to
the supported real-root family. This guide starts with ground states, then
builds up to excitation scans. Unlike XXX multiplet scans, XXZ states are
labelled by magnetization Sz, not total spin S.

## Model and supported regime

`bethe-xxz-pbc` is a separate periodic-chain front end for

```math
\begin{aligned}
H&=\sum_i\left(S_i^xS_{i+1}^x+S_i^yS_{i+1}^y+\Delta S_i^zS_{i+1}^z\right),
\qquad J=1,\quad h=0,\\{}
\text{ground states/sectors:}\quad&
\begin{cases}\Delta\gt -1,&N\text{ even},\\ \Delta\ge0,&N\text{ odd}.\end{cases}
\end{aligned}
```

The anisotropy is required. This periodic-chain implementation supports
ground states and magnetization-sector minima for every finite Delta>=0.
Even rings additionally support $`-1\lt \Delta \lt 0`$ using rank-subtracted,
scaled equations. Negative-Delta odd rings are
[deferred until needed](xxz-negative.md#deferred-odd-ring-work-restart-checklist);
their internal candidate solver is not exposed through this API or CLI.
They are explicitly rejected rather than assigned the wrong ground sector.
For $`0\le \Delta \le 1`$ it also supports a restricted finite-real-root excitation
family, not the complete spectrum or strings. For open boundaries,
use the separate [`bethe-xxz-obc` front end](xxz-open.md).
That free-end solver supports $`\Delta \gt -1`$ ground states for either parity, including the
massive boundary root; its real-root excitation modes require $`0\le \Delta \le 1`$.
Exactly Delta=-1 and lower anisotropies are not included. The
existing analytic thermodynamic `bethe::xxz::spinon_energy` retains its wider
$`-1 \lt  \Delta \le 1`$ domain; it is independent of this finite-size solver.

A separate [spin-helix library API](xxz-spin-helix.md) constructs explicit
eigenstates at commensurate couplings on either parity, including the
isotropic endpoints. These are not general ground-state or excitation scans.

## Ground states and sectors

```sh
build/bethe-xxz-pbc 64 --delta 0.5
build/bethe-xxz-pbc 65 --delta 0.5 --sz -1/2 --roots
build/bethe-xxz-pbc 16 --delta 0.75 --sectors --precision long-double
build/bethe-xxz-pbc 4 --delta 0 --format pretty
build/bethe-xxz-pbc 64 --delta 2 --roots
build/bethe-xxz-pbc 64 --delta -0.9 --roots
build/bethe-xxz-pbc 15 --delta 10 --sectors --precision long-double
# In a binary128-enabled build:
build/bethe-xxz-pbc 64 --delta 0.999999999999999999999999 --precision fp128
```

For ground-state and sector-minimum runs, `--sz` and `--sectors` work as in
the XXX programs. [Precision, tolerance, update budgets, CPU time, roots,
and plain/pretty formatting](command-line.md) use the shared controls.
Default precision is fp64. `--sz` and `--sectors` are mutually exclusive.
With neither ground-state selector, even N uses
Sz=0 and odd N returns one Sz=1/2 ground-state representative. Spin reversal
supplies negative sectors; the all-sector scan reuses those partners. At
generic anisotropy the states are labeled by Sz, not total spin S, and no
SU(2) multiplet degeneracies are inferred. `--spin` and `--spinons` are not
options of this front end; use `--sz` for excitation scans as described below.
Those scans have their own default sector, distinct from the ground-state default.

For negative Delta, the CLI reports the scaled residual convention explicitly
and `--roots` prints both z and the native hyperbolic lambda. For even rings,
the ground-sector labels are symmetric and total momentum is pi for odd M,
zero for even M, with $`M =N /2-\lvert S^z \rvert`$. The smallest-|Sz| ground-state rule is
**not** applied to negative-Delta odd rings: there the global minimum can
instead be polarized, and conjugate root pairs require a different representation.

## Library usage

Library usage stays in the same `bethe::bethe` target:

```cpp
#include <bethe/xxz.hpp>

auto gs = bethe::xxz::ground_state<long double>(64, 0.5L);
auto sector = bethe::xxz::sector_ground_state<long double>(65, 0.5L,
                                                          uni20::half_int::parse("-1/2"));
auto sectors = bethe::xxz::sector_ground_states<long double>(16, 0.75L);
bethe::xxz::SolverOptions<long double> options{.max_iterations = 20000};
auto precise = bethe::xxz::ground_state<long double>(64, 0.99L, options);
auto negative = bethe::xxz::ground_state<long double>(64, -0.9L, options);
```

Ground-state functions now return `xxz::GroundState<Real>` (or a vector of
it), distinct from `xxz::RealState<Real>` used by `solve_real` and excitation
scans. Both retain $`\delta`$, scaled roots, Bethe quantum numbers, Sz,
spin-reversal flag, energy, momentum/index, and convergence diagnostics.
Ground results additionally retain `log_rapidities` for negative Delta,
`GroundSolveStatus` (converged, iteration limit, or stalled), and a
`GroundResidualConvention` tag (`logarithmic_phase` or `negative_rank_scaled`).
Do not reconstruct lambda from rounded z values. Nonnegative ground-state
results preserve all numerical values of the previous solver; their lambda
array is empty. There is **no implicit conversion** between ground results
and real excitation states. Code naming the old ground return type should
use `GroundState` or `auto` and inspect the residual convention.

The common iteration
options live in `<bethe/solver.hpp>` as `bethe::SolverOptions<Real>`; the
model-namespace aliases remain available. Numerical code has no CLI or
presentation dependency.

## Scaled rapidities and Bethe equations

For $`0 \le \Delta \lt  1`$, define $`\gamma =\arccos (\Delta)`$ and store
$`z =\tanh (\lambda)/\tan (\gamma /2)`$, **not the conventional rapidity lambda**.
The output labels this coordinate as scaled rapidity z. Its real branch has
$`(1-\Delta)\,z ^{2} \lt  1+\Delta`$. In this coordinate the equations and energy are

```math
\begin{aligned}
A_{ij}&=\frac{\Delta(z_i-z_j)}{1+\Delta-(1-\Delta)z_iz_j},\\{}
F_i&=2N\arctan z_i-2\pi I_i-\sum_{j\ne i}2\arctan A_{ij},\\{}
E&=\frac{N\Delta}{4}-\sum_i\frac{1+\Delta-(1-\Delta)z_i^2}{1+z_i^2},\\{}
P&=\pi M-\frac{2\pi}{N}\sum_i I_i\pmod{2\pi}.
\end{aligned}
```

The sector minimum has $`M =N /2-\lvert S^z \rvert`$ roots and consecutive labels
`I_i=i-(M-1)/2-(N mod 2)/2`, with zero-based i. For odd N this selects the
reflection-related minimum with labels centered at -1/2. Momentum indices
are calculated with integer arithmetic. These ground-state labels do not
constitute an admissibility rule for XXZ excited states; the excitation window
is defined below.

This rational form is algebraically equivalent to the logarithmic hyperbolic
Bethe equations, but avoids cancellation in expressions such as
$`\cosh (2\,\lambda)-\Delta`$ near Delta=1. At Delta=0 scattering vanishes and the
free-fermion solution takes at most one update. At **exactly** Delta=1 the
solver delegates to the existing XXX sector solver and returns its coordinate
$`z =2\,\lambda_{\mathrm{XXX}}`$; nearby anisotropies are not snapped to the endpoint. All
calculations and parameter parsing use the selected real precision.

### Easy-axis ground states: Delta>1

The massive regime has a different real-rapidity contour. Set
$`\Delta =\cosh (\eta)`$ and store $`z =\tan (\lambda)/\tanh (\eta /2)`$, where
$`-\pi /2\lt \lambda \lt \pi /2`$. This z again tends to $`2\,\lambda_{\mathrm{XXX}}`$ as Delta approaches
one. Unlike the massless contour, it has no finite bound on z. The CLI
continues to label the stored variable **scaled rapidity z**, not lambda.

The trigonometric Bethe equations are Eq. (1.2) of
[Dugave et al.](../CITATIONS.md#dugave-2015) at zero twist. Their Hamiltonian
uses Pauli matrices; our spin-1/2 energy is one quarter of their J=1 energy
at zero field. Algebraically, define the continuous real scattering phase

```math
\begin{aligned}
\phi_{ij}&=\mathrm{atan2}\!\left(z_i-z_j,1+\frac1\Delta+\left(1-\frac1\Delta\right)z_iz_j\right),\\{}
F_i&=2N\arctan z_i-2\pi I_i-2\sum_{j\ne i}\phi_{ij},\\{}
E&=\left(\frac N4-M\right)\Delta+\sum_i\frac{z_i^2-1}{z_i^2+1}.
\end{aligned}
```

The same consecutive ground-state labels and momentum formula above apply.
In conventional lambda coordinates, $`2\,\phi_{\mathrm{ij}}`$ is the continuous lift of
$`2\,\arctan (\tan (\lambda_{i} -\lambda_{j})/\tanh (\eta))`$. Keeping the atan2 quadrant matters:
root differences can exceed pi/2, so replacing it by a principal atan yields
the wrong logarithmic branch. Dividing both atan2 arguments by Delta and
collecting the explicit anisotropy term in E avoids avoidable large-Delta
overflow. Nonrepresentable final energies remain errors.

All `ground_state`, `sector_ground_state`, and `sector_ground_states` entry
points support this regime, including odd N and spin-reversed sectors.
General `solve_real`, `real_quantum_number_window`, and excitation APIs still
reject Delta>1: validating a sector minimum does not classify excited labels.
The same restriction applies to `--quantum-numbers` and `--excitations`.
No strings, thermodynamic mode, or spontaneous-symmetry-broken state is implied.

## Convergence and validation

For nonnegative Delta, convergence uses $`\max \lvert F \rvert /N`$ with the same default 32-epsilon tolerance as XXX,
not an energy-error bound. Both residual and energy describe the returned
iterate, including on budget exhaustion. The solver starts at zero roots
and performs simultaneous updates. It uses O(M^2) work per sweep (O(M) at
Delta=0) and O(M) state storage; retaining every sector uses O(N^2) storage.
Invalid inputs throw `std::invalid_argument`; nonfinite arithmetic or leaving
the finite real-root branch throws `std::runtime_error`. CLI exit statuses
remain 0 for convergence, 2 for an unconverged result, and 1 for errors.

For $`-1\lt \Delta \lt 0`$ on even rings, Newton solves the rank-subtracted equations
in lambda, scaled by $`s =\sqrt{(1+\Delta)/(1-\Delta)}`$ to remain discriminating
near Delta=-1. The residual tag is `negative_rank_scaled`, meaning the
equations in the [negative-anisotropy guide](xxz-negative.md), not $`\max \lvert F \rvert /N`$.
That guide derives the equations from [Kozlowski](../CITATIONS.md#kozlowski-2017),
whose periodic ground-state identification explicitly assumes even length.
This path starts from a free-fermion hyperbolic seed, uses O(M^2) workspace
and O(M^3) work per Newton update, and counts accepted updates against the
shared iteration budget. A stalled line search is distinguished from budget
exhaustion in the result and both CLI formats. Energy and residual always
refer to the returned coordinates at the requested Delta, even on failure.

Normalization checks include $`E_0 (N =2)=-1-\Delta /2`$ (the periodic bond is
counted twice), $`E_0 (N =3)=-1/2-\Delta /4`$, and
$`E_0 (N =4)=-(\Delta +\sqrt{\Delta ^{2}+8})/2`$. Tests compare every sector minimum and
its momentum against independent bit-basis exact diagonalization through
N=9, verify the XX free-fermion and exact XXX limits, independently evaluate
the hyperbolic equations, and exercise native long-double/fp128 precision,
near-endpoint anisotropies, spin reversal, larger chains, and CLI formatting.
For Delta>1 the suite additionally checks the original complex trigonometric
equations (independent of atan2), including a state whose scattering phase
crosses the principal-atan boundary. Sector ED includes Delta=1.01,2,10;
larger systems through N=128 and native-precision approach to XXX are tested.
As Delta tends to infinity, the roots approach
`lambda_j=pi*(I_j-sum(I)/N)/(N-M)` and $`E /\Delta \to N /4-M`$.
Tests check this limit in odd/even sectors and exercise finite answers near
the largest representable anisotropy, as well as explicit overflow errors.

Negative-Delta API tests compare every even-ring sector through N=8 with
independent spin-basis energies and joint energy/translation spectra. Native
two-/four-site formulas, both endpoints, spin reversal, failed iterates, and
chains through N=64 are checked in all supported precisions. CLI tests preserve
z/lambda digits in wide and narrow reports. Separate tests ensure the new
ground result preserves nonnegative solver values exactly and cannot be
implicitly converted to an excitation state. Negative excitation labels
and negative odd-ring ground calls remain rejected, including vacua/scans.

## Real-root excitations

This section applies only to $`0\le \Delta \le 1`$.

```sh
build/bethe-xxz-pbc 64 --delta 0.5 --excitations 10
build/bethe-xxz-pbc 16 --delta 0.75 --sz -1 --excitations all --roots
build/bethe-xxz-pbc 8 --delta 0.5 --sz 2 --excitations all --max-candidates 100
build/bethe-xxz-pbc 8 --delta 0.5 --quantum-numbers -3/2,1/2
```

`--excitations COUNT|all` scans the entire supported quantum-number family,
retaining the lowest COUNT converged states or all of them. The sector minimum
is included. Without `--sz`, the excitation default is Sz=1 for even N and
Sz=1/2 for odd N. Negative Sz uses spin reversal with unchanged energy and
momentum. The reported gap is **E-E0 relative to the global ground state**,
not the minimum of the selected sector. Distinct degenerate states are kept.

`--max-candidates` defaults to 10000 and is checked before either a ground
solve or candidate allocation. `all` does not bypass it. Counting is
overflow-checked; limiting the returned count reduces retained storage, not
the number of solves. Failed candidates are excluded from ranking but counted
and reported, and make the ordering incomplete. If the ground solve fails,
valid candidate energies remain available but their gaps are unavailable.
Either failure gives exit status 2.

### The anisotropy-dependent window

For M=N/2-|Sz|, use distinct sorted quantum numbers of parity
`2*I == N-M-1 (mod 2)`. The supported family is the intersection

```math
\begin{aligned}
|2I|&\le N-M-1&&\text{(conventional XXX window)},\\{}
|I|&\lt I_\infty,\\{}
I_\infty&=\frac12\left[N-M+1-\frac{(N-2M+2)\gamma}{\pi}\right],\\{}
\gamma&=\arccos\Delta.
\end{aligned}
```

The second bound follows by sending one conventional rapidity to infinity
while the others remain finite in the logarithmic equations. The code uses
the equivalent $`2\,I_{\mathrm{infinity}} = N /2 + (N -2\,M +2)\,\arcsin (\Delta)/\pi`$, evaluated in
the selected precision. For 0<Delta<1, labels within a conservative
$`32\,\epsilon \,\max (1,2\,I_{\mathrm{infinity}})`$ band below this doubled bound are also excluded
to keep marginal infinite roots out of the finite branch. Delta=0 uses exact
integer comparisons; Delta=1 uses the existing XXX window exactly.

This is an explicitly restricted family, **not a complete Sz spectrum or a
classification of every real solution**. It excludes complex/string roots,
infinite rapidities, and at Delta=1 the SU(2) descendants. At Delta=0 it still
enumerates only the finite-real-rapidity subset, not the full free-fermion
spectrum. Family sizes can change with anisotropy as roots reach infinity.
For example, N=8, Sz=2 has 6 candidates at Delta=1/2 and 15 at Delta=1.
Even N at Sz=0 has only the single ground configuration in this family;
requesting `all` cannot produce excited states there.

### Reports and specified states

The report displays the selected window, candidate/convergence counts,
Sz, momentum, energy, gap, and quantum numbers. `--quantum-numbers` solves one
explicit configuration in the same window, with M inferred from the list and
Sz=N/2-M; use `none` for the vacuum. It cannot be combined with `--sz`,
`--sectors`, or `--excitations`. Excitation scans likewise cannot be combined
with `--sectors`.

### Excitation API

```cpp
#include <bethe/xxz_excitations.hpp>

auto scan = bethe::xxz::real_excitations<long double>(
    64, 0.5L, uni20::half_int{1}, {.count = 10, .max_candidates = 10000});
auto count = bethe::xxz::real_excitation_count(64, 0.5L, uni20::half_int{1});
auto window = bethe::xxz::real_quantum_number_window(64, 0.5L, uni20::half_int{1});
bethe::xxz::QuantumNumbers numbers{uni20::half_int::parse("-3/2"),
                                  uni20::half_int::parse("1/2")};
auto state = bethe::xxz::solve_real<long double>(8, 0.5L, numbers);
```

`solve_real` also accepts solver options and an optional initial-root span in
the same scaled coordinate. Excitation enumeration, bounded energy ranking,
and convergence bookkeeping are shared with XXX; the anisotropy-dependent
window and state interpretation remain in the XXZ module. Retained scan
storage is O(min(COUNT,candidates)*M+N), including the ground reference.
Tests check every returned energy and momentum with ED multiplicities through
N=9, independently enumerate the window at Delta=0,1/2,1, test both sides of
an infinity threshold, and verify native-precision gaps and exact XXX matching.

Related: [XXX conventions](xxx.md), [analytic spinon dispersions](spinons.md#thermodynamic-dispersion),
and [References and provenance](../CITATIONS.md).
