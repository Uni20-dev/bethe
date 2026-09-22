# Massive free-end XXZ: ground-state boundary root

[Open XXZ guide](xxz-open.md) | [Model catalogue](models.md)

The library module `bethe/xxz_open_massive.hpp` implements sector minima at
`Delta>1`, for either chain parity and either sign of Sz. It uses the same
spin-1/2 Hamiltonian, J=1 and zero boundary fields as the gapless open solver.
This is an implementation checkpoint: the existing `bethe-xxz-obc` executable
and `xxz::open::ground_state` still use the `0<=Delta<=1` path. Integrating the
new state representation into their reports is the next step, not a new
executable or a different physical model.

```cpp
#include <bethe/xxz_open_massive.hpp>

auto state = bethe::xxz::open::massive::ground_state<long double>(16, 3.0L);
auto sector = bethe::xxz::open::massive::sector_ground_state(
    17, 3.0L, uni20::half_int::parse("-1/2"));
```

All arithmetic, including the boundary-root equation and Newton matrix, uses
the selected real type: fp64, long-double, or optional fp128. Check
`converged`, `status`, and `residual_norm` before using an energy. This module
does not implement excited states, boundary fields, or general complex strings.

## Why extending the real interval is insufficient

Write `M=N/2-|Sz|` and `Delta=cosh(eta)`. Bulk roots use
`z=tan(lambda)/tanh(eta/2)`. When `M<N/2`, the ground state has M positive
bulk roots. For even N with Sz=0, its last root must be tracked separately:
as anisotropy grows it passes through infinity in z and becomes purely
imaginary. This is not a singularity of the energy.

The two-site example makes the issue explicit:

```text
z_B^2 = (1+Delta)/(3-Delta),
E0    = -1/2-Delta/4.
```

The root is real below Delta=3, infinite at Delta=3, and imaginary above it.
Changing a range check on an all-real solver cannot represent all three cases.
The transition anisotropy depends on chain length; 3 is not a general cutoff.

The reflection equations follow [Mei and Bolech](../CITATIONS.md#mei-2017),
Eqs. (11)-(12). The massive boundary-root interpretation and exponentially
small finite-size deviations are discussed by
[Grijalva, De Nardis, and Terras](../CITATIONS.md#grijalva-2019), Sec. 4.3.2.
Their Pauli-matrix Hamiltonian at zero boundary fields is divided by four
and shifted by `(N-1)*Delta/4` to obtain ours. The regularization below is
an algebraic implementation of the finite-size equations, not an ideal-string
or thermodynamic approximation.

## Regularized coordinate

For the distinguished root, define

```text
r = (Delta-1)/(Delta+1),
y = 1/z_B^2,
w = -log(1+y/r^2),    y = r^2*expm1(-w).
```

The ground-state branch has `y>-r^2`. Negative w describes a real z_B; w=0
is its crossing through infinity; positive w describes an imaginary z_B.
Large positive w retains the deviation from the boundary pole even when
`y` rounds to `-r^2`. Do not reconstruct w from that rounded value.

`State::rapidities` and `quantum_numbers` contain only bulk roots and their
labels `1,...,M-1` when a boundary root is present, or `1,...,M` otherwise.
The optional `boundary_root` stores w as `log_distance`, y as
`inverse_square`, and the last label M. The label tracks continuation from
the real branch; it does not assert that the continued root is real.

## Coupled equations

Set `P=1+1/Delta`, `Q=1/Delta-1`. For bulk roots, the half scattering phases
are evaluated with `atan2`, retaining their winding. A bulk root z sees the
distinguished root and its reflection through their combined phase

```text
S_B(z,y) = atan2(2*z*(P*y-Q), (P^2-z^2)*y+1-Q^2*z^2).
```

Near the pole at strong anisotropy, both arguments in this expression suffer
cancellation. The implementation expands them in `1/Delta` and
`r^2*exp(-w)` before evaluation; it does not evaluate the displayed differences
naively.

The bulk residual used in the solver is

```text
R_i = [2*N*atan(z_i) - 2*atan(r*z_i) - pi*I_i
       - sum_(j != i, bulk) {atan2(z_i-z_j,P-Q*z_i*z_j)
                           +atan2(z_i+z_j,P+Q*z_i*z_j)}
       - S_B(z_i,y)] / N.
```

Omit `S_B` if there is no distinguished root. Direct and reflected
self-scattering are both excluded.

For even N=2M, the unregularized equation for the last root vanishes at
`1/z_B=0` regardless of Delta. Divide out that vanishing factor **before**
taking the limit, or infinity becomes a spurious solution. Define the analytic
continuation

```text
A(y,a) = atan(sqrt(y)*a)/sqrt(y),       y>0,
         a,                           y=0,
         atanh(sqrt(-y)*a)/sqrt(-y),   y<0,
a_j    = (P+i*z_j)/(1-i*Q*z_j).

R_B = [-N*A(y,1) + A(y,1/r) + sum_bulk Re A(y,a_j)] / N.
```

The implementation evaluates the real parts without complex transcendental
functions. In particular, for w>0 let `t=sqrt(-expm1(-w))`; then
`r*A(y,1/r)=(w/2+log1p(t))/t`. This retains the boundary-pole logarithm
when t rounds to one. Small-argument formulas are used near the crossing to
avoid manufacturing a tiny residual through cancellation.

The energy is

```text
E = [(N-1)/4-M]*Delta + sum_bulk (z_j^2-1)/(z_j^2+1)
    + (1-y)/(1+y),
```

with the last term omitted if there is no boundary root. Evaluate its
denominator as `1+y=4/Delta/(1+1/Delta)^2+r^2*exp(-w)` near the pole.

## Iteration and diagnostics

Damped Newton uses bulk angles `atan(z)` and w, with a central-difference
Jacobian in native precision. Continuation starts at `min(Delta,1.1)` and
increases `Delta-1` by a factor of 1.25 per converged stage. Ordered bulk
roots and, while real, a largest distinguished root are enforced during
backtracking. After the crossing, continuation carries w relative to the
moving pole rather than rounding away the small deviation.

Each accepted Newton update consumes the single global iteration budget;
there is no hidden unbudgeted seed solve. On exhaustion, energy and residual
are still evaluated for the requested Delta at the returned roots.
`root_delta` records the continuation stage reached. `status` distinguishes
convergence, budget exhaustion, and a stalled line search.

The residual is `max(|R_i|,|R_B|)`, with the boundary entry absent when
appropriate. Its boundary component is **regularized**, unlike the old
all-real logarithmic residual. The default tolerance is 32 native epsilons.
It is not an energy-error bound, or a bound on the exponentially small
ground/first-excited splitting; this solver returns only the ground branch.
Work is O(M^3) per Newton update, with O(M^2) dense workspace.

## Validation and remaining integration

Typed tests check every magnetization sector against independent small-chain
exact diagonalization through N=9, original complex reflection equations,
analytic two- and three-site energies in native precision, the root crossing,
spin reversal, exhausted-budget consistency, and longer chains whose boundary
deviation is smaller than machine epsilon. Further tests check continuity to
XXX without snapping Delta to one, and strong-coupling root limits through
Delta=10^9. The near-crossing residual test
specifically detects cancellation that can be hidden by double-precision ED.

Next integrate the explicit boundary coordinate and status into the existing
open-XXZ frontend and ground-state API, retaining the old all-real API for
its restricted excitation family. Do not advertise massive excitation scans
or infer the first-excited splitting from this ground-only solver.
