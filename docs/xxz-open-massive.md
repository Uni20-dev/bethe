# Massive free-end XXZ: ground-state boundary root

[Open XXZ guide](xxz-open.md) | [Model catalogue](models.md)

The library module `bethe/xxz_open_massive.hpp` implements sector minima at
$`\Delta \gt 1`$, for either chain parity and either sign of Sz. It uses the same
spin-1/2 Hamiltonian, J=1 and zero boundary fields as the gapless open solver.
The existing `bethe-xxz-obc` executable and the three ground-state functions in
`bethe::xxz::open` select this module for Delta>1. The gapless and exact XXX
paths are unchanged; no new executable is needed.

```cpp
#include <bethe/xxz_open.hpp>

auto state = bethe::xxz::open::ground_state<long double>(16, 3.0L);
auto sector = bethe::xxz::open::sector_ground_state(
    17, 3.0L, uni20::half_int::parse("-1/2"));
```

All arithmetic, including the boundary-root equation and Newton matrix, uses
the selected real type: fp64, long-double, or optional fp128. Check
`converged`, `status`, and `residual_norm` before using an energy. This module
does not implement excited states, boundary fields, or general complex strings.

## Why extending the real interval is insufficient

Write $`M =N /2-\lvert S^z \rvert`$ and $`\Delta =\cosh (\eta)`$. Bulk roots use
$`z =\tan (\lambda)/\tanh (\eta /2)`$. When $`M \lt N /2`$, the ground state has M positive
bulk roots. For even N with Sz=0, its last root must be tracked separately:
as anisotropy grows it passes through infinity in z and becomes purely
imaginary. This is not a singularity of the energy.

The two-site example makes the issue explicit:

```math
z_B^2=\frac{1+\Delta}{3-\Delta},\qquad E_0=-\frac12-\frac\Delta4.
```

The root is real below Delta=3, infinite at Delta=3, and imaginary above it.
Changing a range check on an all-real solver cannot represent all three cases.
The transition anisotropy depends on chain length; 3 is not a general cutoff.

The reflection equations follow [Mei and Bolech](../CITATIONS.md#mei-2017),
Eqs. (11)-(12). The massive boundary-root interpretation and exponentially
small finite-size deviations are discussed by
[Grijalva, De Nardis, and Terras](../CITATIONS.md#grijalva-2019), Sec. 4.3.2.
Their Pauli-matrix Hamiltonian at zero boundary fields is divided by four
and shifted by $`(N -1)\,\Delta /4`$ to obtain ours. The regularization below is
an algebraic implementation of the finite-size equations, not an ideal-string
or thermodynamic approximation.

## Regularized coordinate

For the distinguished root, define

```math
r=\frac{\Delta-1}{\Delta+1},\qquad y=\frac1{z_B^2},\qquad
w=-\log(1+y/r^2),\qquad y=r^2\operatorname{expm1}(-w).
```

The ground-state branch has $`y \gt -r ^{2}`$. Negative w describes a real z_B; w=0
is its crossing through infinity; positive w describes an imaginary z_B.
Large positive w retains the deviation from the boundary pole even when
$`y`$ rounds to $`-r ^{2}`$. Do not reconstruct w from that rounded value.

`GroundState::rapidities` and `quantum_numbers` contain only bulk roots and their
labels `1,...,M-1` when a boundary root is present, or `1,...,M` otherwise.
The optional `boundary_root` stores w as `log_distance`, y as
`inverse_square`, and the last label M. The label tracks continuation from
the real branch; it does not assert that the continued root is real.

## Coupled equations

Set $`P =1+1/\Delta`$, $`Q =1/\Delta -1`$. For bulk roots, the half scattering phases
are evaluated with `atan2`, retaining their winding. A bulk root z sees the
distinguished root and its reflection through their combined phase

```math
S_B(z,y)=\operatorname{atan2}\!\left(2z(Py-Q),(P^2-z^2)y+1-Q^2z^2\right).
```

Near the pole at strong anisotropy, both arguments in this expression suffer
cancellation. The implementation expands them in $`1/\Delta`$ and
$`r ^{2}\,\exp (-w)`$ before evaluation; it does not evaluate the displayed differences
naively.

The bulk residual used in the solver is

```math
\begin{aligned}
R_i=\frac1N\biggl[&2N\arctan z_i-2\arctan(rz_i)-\pi I_i\\{}
&-\sum_{\substack{j\ne i\\j\ \mathrm{bulk}}}
\left\{\operatorname{atan2}(z_i-z_j,P-Qz_iz_j)
+\operatorname{atan2}(z_i+z_j,P+Qz_iz_j)\right\}-S_B(z_i,y)\biggr].
\end{aligned}
```

Omit $`S_{B}`$ if there is no distinguished root. Direct and reflected
self-scattering are both excluded.

For even N=2M, the unregularized equation for the last root vanishes at
$`1/z_{B} =0`$ regardless of Delta. Divide out that vanishing factor **before**
taking the limit, or infinity becomes a spurious solution. Define the analytic
continuation

```math
\begin{aligned}
A(y,a)&=\begin{cases}
\arctan(\sqrt y\,a)/\sqrt y,&y>0,\\{}
a,&y=0,\\{}
\operatorname{artanh}(\sqrt{-y}\,a)/\sqrt{-y},&y<0,
\end{cases}\\{}
a_j&=\frac{P+iz_j}{1-iQz_j},\\{}
R_B&=\frac1N\left[-NA(y,1)+A(y,1/r)+\sum_{\mathrm{bulk}}\operatorname{Re}A(y,a_j)\right].
\end{aligned}
```

The implementation evaluates the real parts without complex transcendental
functions. In particular, for w>0 let `t=sqrt(-expm1(-w))`; then
`r*A(y,1/r)=(w/2+log1p(t))/t`. This retains the boundary-pole logarithm
when t rounds to one. Small-argument formulas are used near the crossing to
avoid manufacturing a tiny residual through cancellation.

The energy is

```math
E=\left(\frac{N-1}{4}-M\right)\Delta
+\sum_{\mathrm{bulk}}\frac{z_j^2-1}{z_j^2+1}+\frac{1-y}{1+y}.
```

with the last term omitted if there is no boundary root. Evaluate its
denominator as $`1+y =4/\Delta /(1+1/\Delta)^{2}+r ^{2}\,\exp (-w)`$ near the pole.

## Iteration and diagnostics

Damped Newton uses bulk angles `atan(z)` and w, with a central-difference
Jacobian in native precision. Continuation starts at `min(Delta,1.1)` and
increases $`\Delta -1`$ by a factor of 1.25 per converged stage. Ordered bulk
roots and, while real, a largest distinguished root are enforced during
backtracking. After the crossing, continuation carries w relative to the
moving pole rather than rounding away the small deviation.

Each accepted Newton update consumes the single global iteration budget;
there is no hidden unbudgeted seed solve. On exhaustion, energy and residual
are still evaluated for the requested Delta at the returned roots.
`root_delta` records the continuation stage reached. `status` distinguishes
convergence, budget exhaustion, and a stalled line search.

The residual is $`\max (\lvert R_{i} \rvert,\lvert R_{B} \rvert)`$, with the boundary entry absent when
appropriate. Its boundary component is **regularized**, unlike the old
all-real logarithmic residual. The default tolerance is 32 native epsilons.
It is not an energy-error bound, or a bound on the exponentially small
ground/first-excited splitting; this solver returns only the ground branch.
Work is O(M^3) per Newton update, with O(M^2) dense workspace.

## Command-line output

```sh
build/bethe-xxz-obc 16 --delta 3 --roots --precision long-double
build/bethe-xxz-obc 5 --delta 2 --sectors --roots --format plain
```

`--roots` writes bulk `roots` and separate `boundary_roots` tables, linked to
`states` by `state_id`. The boundary schema includes `quantum_number`,
`inverse_square`, `log_distance`, and `kind`.
Its `kind` describes the returned coordinate (real, imaginary, or infinity),
not an exact symbolic determination of the finite-chain crossing. Pretty
output uses separate bulk and boundary tables and retains full-precision
numeric tokens even on narrow terminals. A boundary-only state is not labeled
fully polarized.

Massive reports name the regularized residual convention and show `Root Delta`.
The typed state table includes `root_delta`, `converged`, and descriptive
`status` columns, distinguishing iteration limits from stalled line searches.
CSV/TSV and JSON exports preserve these diagnostics and native-precision
coordinates; see [output schemas](output.md). A nonconverged calculation exits with status 2;
invalid requests exit with status 1. Neither output mode invents momentum.

## Validation and remaining scope

Typed tests check every magnetization sector against independent small-chain
exact diagonalization through N=9, original complex reflection equations,
analytic two- and three-site energies in native precision, the root crossing,
spin reversal, exhausted-budget consistency, and longer chains whose boundary
deviation is smaller than machine epsilon. Further tests check continuity to
XXX without snapping Delta to one, and strong-coupling root limits through
Delta=10^9. The near-crossing residual test
specifically detects cancellation that can be hidden by double-precision ED.

API tests verify dispatch, spin reversal, and exact preservation of the
gapless solver's values. CLI tests exercise all precision modes, roots,
sector tables, narrow-terminal formatting, CPU time, and exhausted budgets.
The old all-real API retains its restricted excitation family. Massive
excitation scans and the first-excited splitting remain unimplemented.
