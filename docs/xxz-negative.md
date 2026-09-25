# Negative XXZ anisotropy: ground-root engine and remaining work

[Model catalogue](models.md) | [Periodic XXZ](xxz.md) | [Open XXZ](xxz-open.md)

The [free-end ground-state API and frontend](xxz-open.md) now support
$-1<\Delta <0$, for either parity of N and every physical Sz sector. The
[periodic API and frontend](xxz.md) support the same interval on **even rings**.
The `bethe::xxz::detail::negative_ground_roots` engine in
[xxz_negative.hpp](../include/bethe/xxz_negative.hpp) solves ground-state
sectors for $-1<\Delta <0$ on even periodic rings and free-end chains of either
parity. It works in fp64, long-double, and optional fp128. Negative-Delta odd
rings are still not public. Their remaining physical-state and lowest-state
selection work is **deferred until a concrete need arises**; see the
[restart checklist](#deferred-odd-ring-work-restart-checklist). The internal
implementation and its tests are retained. Excitations and $\Delta \le -1$
require separate work.

The Hamiltonian is unchanged:

```math
\begin{aligned}
H&=\sum_{\mathrm{bonds}}\left(S_i^xS_j^x+S_i^yS_j^y+\Delta S_i^zS_j^z\right),
\qquad S=\frac12,\\
M&=\frac N2-|S^z|.
\end{aligned}
```

## Sources and geometry

The periodic equations follow [Kozlowski](../CITATIONS.md#kozlowski-2017),
Eqs. (0.4), (0.7), with the Pauli Hamiltonian divided by four. That paper
explicitly assumes **even length**; its ground-state identification cannot
simply be applied to odd rings. The open equations follow the reflection
equations of [Mei and Bolech](../CITATIONS.md#mei-2017), Eqs. (11)-(12),
with the normalization already derived in the [open guide](xxz-open.md).
The rank subtraction and residual scaling below are our algebraic/numerical
implementation of these equations.

## Why subtract the phase ranks?

For negative Delta the ground roots approach a shrinking interval as
Delta approaches -1. A small residual in the unscaled equations can then
mean only that the interval is small, rather than that the roots are accurate.
Define

```math
s=\sqrt{\frac{1+\Delta}{1-\Delta}},\qquad
b=\sqrt{1+\Delta}\sqrt{1-\Delta},\qquad a=-\Delta,\qquad
z_i=s\tanh\lambda_i.
```

The engine works with finite ordered lambda, positive for free ends. On this
branch the scattering denominators are positive. For distinct roots define
the complementary half-scattering phase

```math
\begin{aligned}
p(x)&=\frac{\arctan(s\tanh x)}s,\\
h(x)&=\frac{\operatorname{sgn}(x)}s\operatorname{atan2}\!\left(b,a|\tanh x|\right).
\end{aligned}
```

The original half-scattering phase is
`-sign(x)*pi/2 + s*h(x)`. Its rank-dependent constant cancels the consecutive
ground labels analytically. For PBC those labels are
$I_{i} =i -(M -1)/2$, with zero-based i. For OBC they are $I_{i} =i +1$;
the positive reflected roots and the boundary phase contribute the remaining
constants. This yields the following equations, all divided by N and s
relative to the original logarithmic equations:

```math
\begin{aligned}
R_i^{\mathrm{PBC}}&=\frac1N\left[Np(\lambda_i)-\sum_{j\ne i}h(\lambda_i-\lambda_j)\right],\\
R_i^{\mathrm{OBC}}&=\frac1N\left[2Np(\lambda_i)-\frac2s\arctan\!\left(\frac s{\tanh\lambda_i}\right)
-\sum_{j\ne i}\{h(\lambda_i-\lambda_j)+h(\lambda_i+\lambda_j)\}\right].
\end{aligned}
```

Self-scattering, including the reflected self, is excluded. These formulas
apply to the ordered ground branch, not arbitrary excitation labels.
The analytic Jacobian uses

```math
h'(x)=-\frac{(1-\Delta)a\,\operatorname{sech}^2x}{b^2+a^2\tanh^2x}.
```

The periodic diagonal is $N \,\operatorname{sech} (\lambda_{i})^{2}/(1+z_{i} ^{2})$ minus the sum of
direct h' terms. The open diagonal doubles that driving term, adds
$2\,\operatorname{sech} (\lambda_{i})^{2}/(\tanh (\lambda_{i})^{2}+s ^{2})$, and subtracts direct and reflected
h' terms. Off-diagonal entries are direct h', minus reflected h' for OBC.
Every Jacobian entry is divided by N, matching R.

The energy is evaluated without subtracting $1-\tanh (\lambda)^{2}$:

```math
E=\frac{\mathrm{bonds}\,\Delta}{4}
-\sum_i\frac{(1+\Delta)\operatorname{sech}^2\lambda_i}{1+z_i^2},
\qquad
\mathrm{bonds}=\begin{cases}N,&\mathrm{PBC},\\N-1,&\mathrm{OBC}.\end{cases}
```

## Numerical contract

The initial roots are the analytic free-fermion roots at Delta=0, expressed
in lambda. Newton iteration uses a dense Uni20 solve and a residual-decreasing
line search that preserves root order and the physical domain. Each accepted
step consumes one iteration. There is no hidden continuation budget.

The result retains both z (`rapidities`) and lambda (`log_rapidities`), labels,
energy, accepted iteration count, and status. `residual_norm` is $\max \lvert R_{i} \rvert$
in the **scaled equations above**, evaluated at the returned roots. It is not
interchangeable with the old nonnegative solver's $\max \lvert F_{i} \rvert /N$ diagnostic.
Always check `converged`; an exhausted budget or stalled line search returns
the last iterate, not a claimed solution. Odd periodic input is rejected.

The public free-end result retains both coordinate arrays and carries
`GroundResidualConvention::negative_rank_scaled`, including on failure and
for the polarized vacuum. Its `root_delta` equals the requested Delta; this
engine uses a direct solve, not coupling continuation. The CLI reports the
residual convention and appends lambda to the negative-Delta root table.
Its default ground-state sector remains Sz=0 (even N) or Sz=1/2 (odd N);
the odd-*periodic* caveat below does not apply to an open chain.

The periodic result is now a distinct `xxz::GroundState<Real>`, with the
same two coordinate arrays and its own residual/status tags; real excitation
states remain `xxz::RealState<Real>`. There is no implicit conversion between
them. On an even ring, the symmetric labels give momentum zero for even M
and pi for odd M. All periodic negative-Delta entry points consistently
reject odd N, including polarized sectors and scans, pending the odd-ring work.

Tests independently check all sectors through N=9 for OBC and even PBC at
four negative couplings, using spin-basis exact diagonalization. Periodic
momentum is checked against a joint energy/translation spectrum. Additional
tests cover the original unscaled equations, finite-difference Jacobians,
spin reversal, analytic N=2 and N=3 energies, invalid input, and iteration
budgets. Native-precision checks reach $\Delta =-1+128\,\epsilon$ through N=64,
and ensure that the zero-step seed does not falsely pass the scaled residual.

An independent first-order check uses the staggered rotation to the isotropic
ferromagnet at Delta=-1. In its symmetric fixed-Sz state,
$<S^z_{i} \,S^z_{j} \ge [(N -2M)^{2}-N]/[4\,N \,(N -1)]$ for distinct sites, fixing the
energy slope as Delta increases. This distinguishes the sector minimum from
merely returning the polarized limiting energy.

The initial internal-engine checkpoint passed all 393 tests in the GCC 13 Debug/fp128 build and all
272 in the Clang 20 Release build without MPLAPACK. Sixteen and eleven of
those cases, respectively, belong to this new engine and odd-ring audit.

Free-end public integration adds typed tests of every returned sector, global
selection, native-precision energies near both Delta=-1 and Delta=0, residual
reconstruction, exhausted budgets, and the separation from excitation scans.
The free-end checkpoint passed 445 GCC/fp128 cases and 307 Clang cases. The
CLI checks also pass in an applications-only build using the published Uni20
pin, without a sibling checkout or GoogleTest.

The even-periodic integration checkpoint passes 457 GCC 13 Debug/fp128 tests
and 315 Clang 20 Release tests, including its additional typed API/translation
checks and full-precision CLI cases. The periodic CLI also passes against
the pinned Uni20 build. Shared XXZ definitions now live below the solvers in
`xxz_common.hpp`, allowing both real and polynomial engines to remain
independent of the public ground-state dispatcher.

## Why odd rings are not public

An odd ring is not bipartite. Two issues need explicit handling:

1. The real-root continuation can encounter colliding roots and require
   complex conjugate pairs. A real-coordinate range change is insufficient;
   regularized complex-root variables or a polynomial formulation are needed.
   The [polynomial building block](xxz-polynomial.md) now supplies coefficient
   equations, analytic derivatives, and observables. An internal
   [adaptive odd-ring driver](xxz-odd-continuation.md) now follows these
   coefficients with momentum, conditioning, and energy-concavity checks;
   selected sector energies agree with independent spin-basis calculations
   through 21 sites. An independent [Wronskian diagnostic](xxz-wronskian.md)
   and explicit phantom-point eigenvectors now support the admissibility
   work, but do not yet provide a general acceptance policy.
   It is not yet a production odd-ring solver.
2. The global minimum need not lie in the smallest-|Sz| sector. An independent
   ED regression for N=5, Delta=-0.9 finds the fully polarized energy -1.125
   below the |Sz|=1/2 sector minimum (approximately -0.9898034892).
   The internal [all-sector candidate scan](xxz-odd-continuation.md#comparing-sectors-without-assuming-the-answer)
   now compares every folded branch and does not select an incomplete
   minimum. Its result is separate from physical-state validation.

The second observation alone rules out routing all negative couplings through
the existing `ground_state` wrapper's sector choice. The internal work is
retained for a future implementation, but no complete odd-ring phase diagram
is claimed by the small-system audit.

## Deferred odd-ring work: restart checklist

**Decision:** leave negative-Delta odd periodic ground states unsupported in
the public API and CLI, and revisit only when those results are needed.
This does not defer negative-Delta even rings or free ends, which are already
supported. It does not remove the separate explicit spin-helix API. Other
models can proceed independently; there is no commitment to an experimental
odd-ring frontend or to completing every singular-state family first.

### Preserved checkpoint

Through commit `d7fc028`, the internal implementation has polynomial
continuation, all-sector candidate comparison, regular-state and Wronskian
diagnostics, explicit helix states, mixed-phantom dressing and bounded
nonzero-amplitude checks, native-precision root recovery, and independent
free-sea/helix variational upper bounds. Follow the
[continuation guide](xxz-odd-continuation.md) for equations, test coverage,
counterexamples, and links to the individual modules. Existing regression
tests remain part of the normal suite; deferral is not permission to break
these components.

### Issues to resolve before public ground-state support

1. **Lowest-state selection within a sector.** The driver follows one
   free-sea branch at fixed momentum. A physical, converged eigenstate need
   not be the sector minimum. Continuity, concavity, and variational bounds
   reject some wrong branches but do not establish that no lower branch or
   momentum sector was missed. The documented N=17 branch jump is a concrete
   regression, not evidence of a complete classification.
2. **A single physical-state acceptance policy.** Regularity, helix matching,
   and mixed-phantom witnesses cover different cases. Generic-q Wronskian
   assumptions cannot be applied indiscriminately at roots of unity;
   repeated roots, exact strings, and other singular limits can remain
   unresolved. Witness budgets and conditioning limits must remain visible,
   and an unresolved check must not be interpreted as a zero or absent state.
3. **Recoverable numerical failures and practical costs.** The continuation
   and Wronskian code still call Uni20's ordinary dense solve, whose default
   singular-matrix policy can abort. An expected singularity needs a safe
   failure path before public exposure, without changing process-global
   policy. Establish a useful size/precision envelope: polynomial conditioning
   and exponential subset work in the optional wavefunction checks matter.
4. **Ground-state versus candidate API semantics.** Keep polynomial
   coefficients and affine coordinates authoritative; complex-root recovery
   is optional and may be unresolved. Extend result/status and residual tags
   rather than squeezing this representation into real-root arrays. Compare
   every required magnetization sector for a global ground state, retaining
   failed or unresolved sectors instead of selecting from an incomplete set.
   The current internal `lowest_index` compares converged branch energies;
   it is not by itself a public ground-state acceptance flag.

### Where to resume

First specify the actual requested lengths, magnetization sectors, couplings,
and accuracy. Then audit branch selection over that domain against independent
spin-basis spectra: dense coupling sweeps on small rings, selected larger
sparse calculations, momentum competitors, and neighborhoods of known
collisions. Use the existing `tests/reference_xxz_odd_ed.py` and typed tests
as the starting point. Cross-check precision and continuation step sensitivity.
This investigation should determine whether the present branch is adequate or
whether competing branches or a more systematic Q-system treatment are needed;
adding another necessary acceptance inequality alone does not settle it.

A future numerical release need not wait for a new general completeness
theorem. It does need a defensible state-selection method, discriminating
validation over its stated scope, and honest unresolved/failure behaviour.
Only after that decision should API/CLI integration and any experimental
opt-in be designed. Negative-Delta excitation scans remain separate work.
