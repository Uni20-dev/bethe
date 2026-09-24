# Four-defect bound droplets

[Ferromagnetic overview](biquadratic-ferromagnetic.md) · [Triples](biquadratic-bound-triples.md) · [Triple plus defect](biquadratic-triple-defect.md)

Four TL singlet insertions can bind into one droplet, or **quartet**. This
is the last of the five four-defect string topologies, alongside four real
roots, one pair plus two real roots, two pairs, and one triple plus a real
root. The selected-state library API accepts odd/even N>=8 in fp64, native
long double, or enabled fp128:

```cpp
#include <bethe/biquadratic_ferromagnetic.hpp>
namespace ferro = bethe::biquadratic::ferromagnetic;
auto state = ferro::bound_quartet<long double>(128, 1);
if (state.reference.converged) {
    auto gap = state.tl_energy; // E-(N-1), computed without extensive subtraction
    auto inner_L = state.reference.inner_log_deviation;
    auto outer_L = state.reference.outer_log_deviation;
    auto outer_phi = state.reference.outer_deviation_phase;
}
```

`mode=1,...,N-7` selects the four-string branch, not a rank in the full
module or excitation spectrum. The module is `ell=N-8`; these are TL
singlet insertions, not four physical spin flips. CLI integration is a
separate next step. The lower API is
`bethe::xxz::quantum_group::four_string::bound_quartet(N,Delta,mode,options)`
for finite Delta>1, with no blanket convergence or branch-existence guarantee.

## Two deviations, four real unknowns

The reference is quantum-group XXZ with opposite end fields, at Delta=3/2
for the biquadratic model. Define eta=acosh(Delta), and use canonical roots

```text
u = (eta+d)/2 + i*a/2,       conjugate(u),
w = u+eta+z,                conjugate(w),
d = sigma*exp(-L_inner),    z = exp(-L_outer+i*phi),
J = N-6-mode,               sigma = (-1)^(N-J-1) = (-1)^(mode+1).
```

The inner pair has a signed real deviation, while the outer pair needs a
complex relative deviation. The solver retains both logarithms, the sign,
and phi even when the deviations underflow or round away in the roots.
It requires `0<a<pi`, `0<a+2 Im(z)<pi`, and `|d|,|z|<eta/4`.

Our regularization of
[Bajnok et al., Eq. (5.12)](../CITATIONS.md#bajnok-2020) cancels the singular
internal scattering in the sum of the u and w equations. Write

```text
D0 = log sinh(u+eta/2) - log sinh(u-eta/2),
D1 = log sinh(w+eta/2) - log sinh(w-eta/2),       D = D0+D1,
A  = log sinh(2eta+z),
B  = log sinh(3eta+d+i*a+z) - log sinh(eta+d+i*a+z),
K  = log sinh(2eta+i*a+z) - log sinh(i*a+z),
R  = log sinh(3eta+d+z) - log sinh(eta+d+z),
q  = log sinh(2eta+d),
Q  = log sinh(4eta+d+2 Re(z)) - log sinh(2eta+d+2 Re(z)),
c_d = log(sinh(d)/d),       c_z = log(sinh(z)/z),
t0 = Theta(2a;eta),         t1 = Theta(2a+4 Im(z);eta).
```

Theta and wrap are as in the pair/triple guides. The four residuals are

```text
f0 = 2 Im(D)+pi - [t0+t1+2 Im(B)+2 Im(K)+pi*(J+1)]/N,
f1 = Re(D)  - [q+L_inner-Re(c_d)+Q+2 Re(B)+2 Re(R)]/(2N),
f2 = Re(D1) - [Re(A)+L_outer-Re(c_z)+Re(B)+Re(K)+Re(R)+Q]/(2N),
f3 = wrap(2N Im(D1)-Im(A)+phi+Im(c_z)-Im(B)-Im(K)-Im(R)-t1+pi)/(2N).
```

The unsquared equations fix sigma; the product phase alone cannot recover
it. The analytic Jacobian differentiates these regularized expressions.
The ideal-string initializer reduces the center equation to

```text
N*Theta(a;2eta) - sum_{k=1,2,3} Theta(2a;k*eta) = pi*J.
```

Initialization and the full finite-deviation solve share one Newton budget;
only the latter can certify convergence. Work and storage per step do not
grow with N. This reuses the common string Newton driver and correction
functions, plus the stable inner-pair energy. The outer conjugate pair adds
its real energy contribution directly. Failed results retain consistent
estimates and residuals, not verified energies.

## Threshold and coverage

The thermodynamic droplet formula
`g_M=2*sinh(eta)*tanh(M*eta/2)` gives **g4=15/7** at Delta=3/2; see
[Nachtergaele, Spitzer and Starr, Theorem 2.1](../CITATIONS.md#nachtergaele-spitzer-starr-2007).
The lowest quartet branch approaches this value from above. It lies below
the separated triple-plus-single threshold 3 and the two-pair threshold
10/3, but the global positive gap remains the one-defect edge, approaching 1.

All five families together reproduce independent four-defect TL-module ED
**as a multiset**, leaving no unmatched levels at these tested sizes:

| N | Four real | Pair + two real | Two pairs | Triple + real | Quartet | Module dimension |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 1 | 6 | 1 | 5 | 1 | 14 |
| 9 | 5 | 20 | 3 | 12 | 2 | 42 |
| 10 | 15 | 45 | 6 | 21 | 3 | 90 |

This is finite-size numerical evidence, not a general completeness proof
or a full spin-chain spectrum. Individual quartet branches match distinct
ED levels at Delta=1.25,1.5,2,3; mode 1 is the module minimum in the ferro
mapping at these sizes. Tests also verify the original inner/outer complex
equations with both inner signs, every Jacobian column, the threshold and
both logarithmic deviations through N=100000, and failure diagnostics.
