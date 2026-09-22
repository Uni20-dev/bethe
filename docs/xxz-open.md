# Open XXZ chains with free ends

[Back to the overview](../README.md)

The `bethe-xxz-obc` front end combines anisotropy with free ends. It uses
the same scalar precision and reporting controls as the other programs, but
its equations include reflected scattering and a boundary phase. Simply
doubling the periodic chain length would miss that phase.

## Model and first calculations

```text
H = sum_(i=0)^(N-2) [Sx_i Sx_(i+1) + Sy_i Sy_(i+1) + Delta Sz_i Sz_(i+1)],
J=1, h=0, 0 <= Delta <= 1, no boundary fields, N >= 2.
```

```sh
build/bethe-xxz-obc 64 --delta 0.5
build/bethe-xxz-obc 65 --delta 0.5 --sz -1/2 --roots
build/bethe-xxz-obc 16 --delta 0.75 --sectors --precision long-double
# In a binary128-enabled build:
build/bethe-xxz-obc 64 --delta 0.999999999999999999999999 --precision fp128
```

The default is the ground state, represented by Sz=0 for even N and Sz=1/2
for odd N. `--sz` selects a sector minimum; `--sectors` returns one minimum
per magnetization sector. Negative sectors use spin reversal. There is no
lattice momentum for free ends, and generic XXZ states are not classified
by total spin S or SU(2) multiplets. There is no `--spin`, `--spinons`, or
boundary-selection option.

The default arithmetic is fp64; long-double and optional fp128 retain native
precision in parameters, roots, energies, gaps, and output. See the shared
[CLI controls](command-line.md) for tolerances, CPU time, and presentation.

## Real-root excitations

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

```text
1 <= I <= N-M,
I < I_infinity,
I_infinity = N-M+1 - (N-2*M+1)*gamma/pi,
gamma = acos(Delta).
```

The strict second bound comes from sending one rapidity to infinity with the
others finite. The implementation uses the equivalent expression
`I_infinity=(N+1)/2+(N-2*M+1)*asin(Delta)/pi`. For 0<Delta<1 it also excludes
a conservative `32*epsilon*max(1,I_infinity)` band below this threshold.
At Delta=0, the exact integer condition is `2*I<N+1`; at Delta=1, the existing
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
auto scan = obc::real_excitations<long double>(
    16, 0.5L, uni20::half_int{1}, {.count = 10, .max_candidates = 10000});
auto count = obc::real_excitation_count(8, 0.5L, uni20::half_int{2});
auto window = obc::real_quantum_number_window(8, 0.5L, uni20::half_int{2});
obc::QuantumNumbers numbers{uni20::half_int{1}, uni20::half_int{3}};
auto state = obc::solve_real<long double>(8, 0.75L, numbers);
```

`obc::RealState<Real>` holds delta, scaled roots, labels, Sz, spin-reversal
flag, energy, and convergence diagnostics, with **no momentum members**.
`solve_real` also accepts shared `bethe::SolverOptions<Real>` and an optional
initial-root span. Enumeration, bounded ranking, and convergence bookkeeping
are shared with the other models. Numerical code has no presentation dependency.

## Boundary equations and limiting cases

Use the same scaled coordinate as periodic XXZ:
`z=tanh(lambda)/tan(gamma/2)`, not the conventional rapidity lambda. Physical
roots are positive and satisfy `(1-Delta)*z^2<1+Delta`. Writing
`p=1+Delta`, `q=1-Delta`, and `c=q/p`, the equations and energy are

```text
A_ij = Delta*(z_i-z_j)/(p-q*z_i*z_j),
B_ij = Delta*(z_i+z_j)/(p+q*z_i*z_j),
F_i  = 4*N*atan(z_i) + 4*atan(c*z_i) - 2*pi*I_i
       - 2*sum_(j != i) [atan(A_ij) + atan(B_ij)],
E    = (N-1)*Delta/4 - sum_i (p-q*z_i^2)/(1+z_i^2).
```

Both direct and reflected self-scattering terms are excluded. The boundary
term `4*atan(c*z)` is essential even with zero boundary fields. At Delta=0
scattering vanishes but the boundary term does not: the roots become
`z_i=tan(pi*I_i/[2*(N+1)])`, with energy `-sum_i cos(pi*I_i/(N+1))`.
These are the open-chain standing waves, with denominator **N+1**, not N.

The sector minimum fills I=1,...,M. Simultaneous updates start at zero unless
finite nonnegative in-branch guesses are supplied. The implementation moves
the boundary correction to the right using
`atan(z)-atan(c*z)=atan(2*Delta*z/(p+q*z^2))`; the resulting driving denominator
is `2*(N+1)`. This keeps updates on the supported positive finite branch and
solves Delta=0 in at most one update. At **exactly** Delta=1 it delegates to
the open XXX solver, preserving its roots `z=2*lambda_XXX` and diagnostics.
Nearby anisotropies are not snapped to the endpoint.

Convergence uses `max|F|/(2*N)`, with default tolerance 32 times the selected
type's epsilon. It is not an energy-error bound. Both residual and energy
describe the returned iterate, even on budget exhaustion. Work per update
is O(M^2), or O(M) at Delta=0, with O(M) state storage. Retaining all sectors
costs O(N^2) storage. Excitation scans retain
O(min(COUNT,candidates)*M+N), including their ground reference.

The separate [massive ground-state module](xxz-open-massive.md) now handles
`Delta>1` in the library, including the continued boundary root. Its CLI and
unified ground-state API integration are pending; the commands and all-real
API described above still require `0<=Delta<=1`.

## Validation

Independent small-chain exact diagonalization checks sector minima and every
enumerated excitation, including multiplicities, through N=9. Analytic checks
include `E0(N=2)=-1/2-Delta/4`, its root
`z=sqrt((1+Delta)/(3-Delta))`, and
`E0(N=3)=-(Delta+sqrt(Delta^2+8))/4`. Irrational-energy and gap tests discriminate
against narrowing higher precision to double. Further checks cover independent
hyperbolic residuals, the XX and XXX limits, both sides of an infinity threshold,
larger chains, exhausted budgets, and full-precision plain/pretty output.

Related: [periodic XXZ](xxz.md), [free-end XXX](open-chains.md), and
[references and provenance](../CITATIONS.md).
