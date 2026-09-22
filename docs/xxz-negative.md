# Negative XXZ anisotropy: ground-root engine and remaining work

[Model catalogue](models.md) | [Periodic XXZ](xxz.md) | [Open XXZ](xxz-open.md)

This is an implementation checkpoint, **not a newly supported frontend mode**.
The internal `bethe::xxz::detail::negative_ground_roots` engine in
[xxz_negative.hpp](../include/bethe/xxz_negative.hpp) solves ground-state
sectors for `-1<Delta<0` on even periodic rings and free-end chains of either
parity. It works in fp64, long-double, and optional fp128. The existing public
ground-state functions and executables still require nonnegative Delta.
Public integration, odd periodic sectors, and global sector selection remain
part of this extension; excitations and `Delta<=-1` require separate work.

The Hamiltonian is unchanged:

```text
H = sum_bonds [Sx_i*Sx_j + Sy_i*Sy_j + Delta*Sz_i*Sz_j],  S=1/2.
M = N/2 - |Sz|.
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

```text
s = sqrt((1+Delta)/(1-Delta)),
b = sqrt(1+Delta)*sqrt(1-Delta),  a = -Delta,
z_i = s*tanh(lambda_i).
```

The engine works with finite ordered lambda, positive for free ends. On this
branch the scattering denominators are positive. For distinct roots define
the complementary half-scattering phase

```text
p(x) = atan(s*tanh(x))/s,
h(x) = sign(x)*atan2(b, a*abs(tanh(x)))/s.
```

The original half-scattering phase is
`-sign(x)*pi/2 + s*h(x)`. Its rank-dependent constant cancels the consecutive
ground labels analytically. For PBC those labels are
`I_i=i-(M-1)/2`, with zero-based i. For OBC they are `I_i=i+1`;
the positive reflected roots and the boundary phase contribute the remaining
constants. This yields the following equations, all divided by N and s
relative to the original logarithmic equations:

```text
PBC: R_i = [N*p(lambda_i) - sum_(j!=i) h(lambda_i-lambda_j)]/N,

OBC: R_i = [2*N*p(lambda_i) - 2*atan(s/tanh(lambda_i))/s
            - sum_(j!=i) {h(lambda_i-lambda_j)+h(lambda_i+lambda_j)}]/N.
```

Self-scattering, including the reflected self, is excluded. These formulas
apply to the ordered ground branch, not arbitrary excitation labels.
The analytic Jacobian uses

```text
h'(x) = -(1-Delta)*a*sech(x)^2 / [b^2+a^2*tanh(x)^2].
```

The periodic diagonal is `N*sech(lambda_i)^2/(1+z_i^2)` minus the sum of
direct h' terms. The open diagonal doubles that driving term, adds
`2*sech(lambda_i)^2/(tanh(lambda_i)^2+s^2)`, and subtracts direct and reflected
h' terms. Off-diagonal entries are direct h', minus reflected h' for OBC.
Every Jacobian entry is divided by N, matching R.

The energy is evaluated without subtracting `1-tanh(lambda)^2`:

```text
E = bonds*Delta/4 - sum_i (1+Delta)*sech(lambda_i)^2/(1+z_i^2),
bonds = N (PBC), N-1 (OBC).
```

## Numerical contract

The initial roots are the analytic free-fermion roots at Delta=0, expressed
in lambda. Newton iteration uses a dense Uni20 solve and a residual-decreasing
line search that preserves root order and the physical domain. Each accepted
step consumes one iteration. There is no hidden continuation budget.

The result retains both z (`rapidities`) and lambda (`log_rapidities`), labels,
energy, accepted iteration count, and status. `residual_norm` is `max|R_i|`
in the **scaled equations above**, evaluated at the returned roots. It is not
interchangeable with the old nonnegative solver's `max|F_i|/N` diagnostic.
Always check `converged`; an exhausted budget or stalled line search returns
the last iterate, not a claimed solution. Odd periodic input is rejected.

Tests independently check all sectors through N=9 for OBC and even PBC at
four negative couplings, using spin-basis exact diagonalization. Periodic
momentum is checked against a joint energy/translation spectrum. Additional
tests cover the original unscaled equations, finite-difference Jacobians,
spin reversal, analytic N=2 and N=3 energies, invalid input, and iteration
budgets. Native-precision checks reach `Delta=-1+128*epsilon` through N=64,
and ensure that the zero-step seed does not falsely pass the scaled residual.

An independent first-order check uses the staggered rotation to the isotropic
ferromagnet at Delta=-1. In its symmetric fixed-Sz state,
`<Sz_i*Sz_j>=[(N-2M)^2-N]/[4*N*(N-1)]` for distinct sites, fixing the
energy slope as Delta increases. This distinguishes the sector minimum from
merely returning the polarized limiting energy.

This checkpoint passes all 393 tests in the GCC 13 Debug/fp128 build and all
272 in the Clang 20 Release build without MPLAPACK. Sixteen and eleven of
those cases, respectively, belong to this new engine and odd-ring audit.

## Why odd rings are still work in progress

An odd ring is not bipartite. Two issues need explicit handling:

1. The real-root continuation can encounter colliding roots and require
   complex conjugate pairs. A real-coordinate range change is insufficient;
   regularized complex-root variables or a polynomial formulation are needed.
   The [polynomial building block](xxz-polynomial.md) now supplies coefficient
   equations, analytic derivatives, and observables. An internal
   [adaptive odd-ring driver](xxz-odd-continuation.md) now follows these
   coefficients with momentum, conditioning, and energy-concavity checks;
   selected sector energies agree with independent spin-basis calculations
   through 21 sites. It is not yet a production odd-ring solver.
2. The global minimum need not lie in the smallest-|Sz| sector. An independent
   ED regression for N=5, Delta=-0.9 finds the fully polarized energy -1.125
   below the |Sz|=1/2 sector minimum (approximately -0.9898034892).

The second observation alone rules out routing all negative couplings through
the existing `ground_state` wrapper's sector choice. These are not grounds
for silently dropping odd rings from the project: the next work is their
root/state classification, followed by public integration with honest root
representation and residual diagnostics. No complete odd-ring phase diagram
is claimed by the small-system audit.
