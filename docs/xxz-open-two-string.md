# A targeted complex-root singlet on long open chains

[Biquadratic application](biquadratic.md) · [Q-system search](xxz-open-qsystem.md)

This page describes the AF singlet branch. The same regularized solver also
supports [ferromagnetic empty-sea bound pairs](biquadratic-bound-pairs.md),
including signed deviations and other string labels on odd/even chains,
and [a bound pair plus one real defect](biquadratic-pair-defect.md).

The Q-system is useful for discovering small spectra, but finding every
polynomial is an expensive way to obtain one low-lying level on a long chain.
The separate two-string solver targets an ell=0 branch with **one conjugate
pair and N/2-2 positive real roots**, for even N>=4. In the spin-1 biquadratic
chain this is a physical singlet with multiplicity one.

```sh
build/bethe-biquadratic-obc 128 --singlet-excitation
build/bethe-biquadratic-obc 512 --singlet-excitation --precision long-double --roots
build/bethe-biquadratic-obc 128 --singlet-excitation --json singlet.json
```

This branch is the lowest excited singlet in the small-chain ED checks.
Following it to longer chains is a **targeted branch calculation**, not a
proof that no other singlet lies below it, and not a global first-excitation
claim across TL modules. It is separate from `--excitations`, whose real-root
family does not contain this state, and from the multistart `--q-spectrum`.

## Coordinates that retain the string deviation

We use the same opposite-end-field XXZ reference, Delta=cosh(eta)>1, as the
TL layer. These are not the ordinary zero-boundary-field XXZ equations.
Starting from [Bajnok et al., equation (5.12)](https://arxiv.org/html/1910.07805v2#S5),
write

```text
u_j = i*alpha_j/2,          0<alpha_1<...<alpha_r<pi, r=N/2-2,
u_± = (eta+d)/2 ± i*a/2,   0<a<pi, 0<d<eta,
L = -log(d).
```

The real labels are I_j=j and the two-string label is 1. The positive-deviation
branch is fixed for this singlet API; more strings are not implemented.
The separate bound-pair API allows negative deviations and other string labels.
The equations below are our regularization of the
published equations, not a claim that the paper proves its energy ordering.

For a nearly ideal string, forming `u_+ + u_- - eta` can give zero by rounding,
although the true deviation is positive. We therefore solve for **L**, and
evaluate the singular factor through
`log(sinh(d)) = -L + log(sinh(d)/d)`. The latter correction has a regular
small-d limit. Even if `exp(-L)` underflows, L remains finite and meaningful.
We never set d=0 as a physical approximation to obtain convergence.

Define

```text
Theta(beta;w) = 2 atan2(sin(beta/2), tanh(w)*cos(beta/2)),
G(beta;w) = log|sinh(w+i*beta/2)|,
w_+ = (3 eta+d)/2,  w_- = (eta-d)/2,
S(beta) = Theta(beta;w_+) + Theta(beta;w_-),
C(a;d) = 2 atan2(tanh(d/2)*cos(a/2), sin(a/2)).
```

The real-root equations are

```text
2N Theta(alpha_j;eta/2)
 - sum_(k!=j) [Theta(alpha_j-alpha_k;eta)+Theta(alpha_j+alpha_k;eta)]
 - S(alpha_j-a) - S(alpha_j+a) = 2 pi j.
```

The pair's phase and log-modulus equations are

```text
2N [Theta(a;eta+d/2)+C(a;d)] - 2 Theta(2a;eta)
 - sum_j [S(a-alpha_j)+S(a+alpha_j)] = 2 pi,

2N [G(a;eta+d/2)-G(a;d/2)] - log(sinh(2eta+d)) - L
 + log(sinh(d)/d)
 - sum_(j,sign=±) [G(a+sign*alpha_j;w_+)-G(a+sign*alpha_j;w_-)] = 0.
```

These retain the finite deviation. The reflected self-scattering is essential:
dropping it changes both equations. The real roots need not all lie below a;
the continuous phase convention also covers their crossing of the center.

## Algorithm, precision and output

A damped Newton solver uses an analytic Jacobian in the selected native
fp64, long-double or optional fp128 type. An ideal-string solve initializes
the coordinates; acceptance always uses the finite-deviation equations.
Both stages share the user's `--max-iterations` budget, with no hidden extra
updates. A zero budget returns the unconverged seed. Residuals are divided
by 2N; the default tolerance is 32 machine epsilon, not an energy-error bound.

There is no site cutoff. Storage is O(N^2) and a dense Newton update costs
O(N^3); long chains are not cost-free. Conditioning, initialization or limited
precision can still prevent convergence. Failures retain consistent estimates
and diagnostics, return exit 2, and never emit a verified gap. No fallback to
another precision or relaxation of the requested tolerance occurs.

Named tables are `states`, `reference`, `string`, and optionally `roots`.
The `string` table always carries `center=a` and `log_deviation=L`.
`deviation=d` is null on underflow, not a reported exact zero. Displayed u
coordinates are rounded; **do not reconstruct d from their difference**.
The state table reports physical energy, TL and auxiliary XXZ energies,
E-E0, multiplicity, and separate phase/modulus diagnostics. Gaps are emitted
only if both the excitation and global ground reference converge.

The energy and multiplicity maps are unchanged:
`E_biquadratic = 2 E_ref - 7(N-1)/4`, ell=0, multiplicity=1.

```cpp
#include <bethe/biquadratic_two_string.hpp>
auto state = bethe::biquadratic::two_string::singlet<long double>(128);
if (state.reference.converged) {
  auto energy = state.energy;
  auto L = state.reference.log_deviation;
}
```

The lower APIs are `bethe::temperley_lieb::two_string::singlet(N,loop_weight)`
and `bethe::xxz::quantum_group::two_string::singlet(N,Delta)`, with loop weight
>2 and Delta>1 respectively. The public biquadratic adapter fixes Delta=3/2.

Tests cover the exact four-site singlet at native precision, analytic
Jacobians, independent small-chain ED through N=10, original complex Bethe
equations, longer chains through N=512, logarithmic coordinates below the
underflow threshold, iteration budgets, failure output and typed exports.
This does not yet provide arbitrary complex-root states, multiple strings,
singlet excitation ranks, or odd-chain spinons.
