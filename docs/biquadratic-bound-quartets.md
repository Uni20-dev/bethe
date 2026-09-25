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
module or excitation spectrum. The module is $\ell =N -8$; these are TL
singlet insertions, not four physical spin flips. The lower API is
`bethe::xxz::quantum_group::four_string::bound_quartet(N,Delta,mode,options)`
for finite Delta>1, with no blanket convergence or branch-existence guarantee.

## Command line and output

```sh
bethe-biquadratic-obc 128 --ferromagnetic --bound-quartets 4
bethe-biquadratic-obc 9 --ferromagnetic --bound-quartets all --roots
bethe-biquadratic-obc 128 --ferromagnetic --bound-quartets 4 --precision fp128 --json quartets.json
```

`COUNT` selects the first modes, capped at the available $N -7$; `all`
means all quartet labels, **not all four-defect excitations**. The usual
`--max-candidates` limit applies. This selection cannot be combined with
other state families or an explicit through-line sector.

The shared bound-cluster frontend writes `states`, the exact ground-space
`reference`, and `string` tables; `--roots` adds the four reconstructed roots
per state. The string table preserves both logarithmic deviations, the inner
sign and the outer phase. Underflowed deviation components are null, not a
claim of an exact ideal string. JSON/CSV/TSV exports and `--no-retain`
streaming use the same tables as the pair and triple frontends. For example,
`--csv-table string=quartets.string.csv` exports the string parameters.
Failed solves return exit status 2, keep diagnostic estimates, and leave
the verified `gap` null.

## Two deviations, four real unknowns

The reference is quantum-group XXZ with opposite end fields, at Delta=3/2
for the biquadratic model. Define eta=acosh(Delta), and use canonical roots

```math
\begin{aligned}
u&=\frac{\eta+d}{2}+\frac{ia}{2},\qquad \overline u,\\
w&=u+\eta+z,\qquad \overline w,\\
d&=\sigma e^{-L_{\mathrm{inner}}},\qquad z=e^{-L_{\mathrm{outer}}+i\phi},\\
J&=N-6-\mathrm{mode},\qquad \sigma=(-1)^{N-J-1}=(-1)^{\mathrm{mode}+1}.
\end{aligned}
```

The inner pair has a signed real deviation, while the outer pair needs a
complex relative deviation. The solver retains both logarithms, the sign,
and phi even when the deviations underflow or round away in the roots.
It requires $0<a <\pi$, `0<a+2 Im(z)<pi`, and $\lvert d \rvert,\lvert z \rvert <\eta /4$.

Our regularization of
[Bajnok et al., Eq. (5.12)](../CITATIONS.md#bajnok-2020) cancels the singular
internal scattering in the sum of the u and w equations. Write

```math
\begin{aligned}
D_0&=\log\sinh(u+\eta/2)-\log\sinh(u-\eta/2),\\
D_1&=\log\sinh(w+\eta/2)-\log\sinh(w-\eta/2),\qquad D=D_0+D_1,\\
A&=\log\sinh(2\eta+z),\\
B&=\log\sinh(3\eta+d+ia+z)-\log\sinh(\eta+d+ia+z),\\
K&=\log\sinh(2\eta+ia+z)-\log\sinh(ia+z),\\
R&=\log\sinh(3\eta+d+z)-\log\sinh(\eta+d+z),\\
q&=\log\sinh(2\eta+d),\\
Q&=\log\sinh(4\eta+d+2\operatorname{Re}z)-\log\sinh(2\eta+d+2\operatorname{Re}z),\\
c_d&=\log\frac{\sinh d}{d},\qquad c_z=\log\frac{\sinh z}{z},\\
t_0&=\Theta(2a;\eta),\qquad t_1=\Theta(2a+4\operatorname{Im}z;\eta).
\end{aligned}
```

Theta and wrap are as in the pair/triple guides. The four residuals are

```math
\begin{aligned}
f_0&=2\operatorname{Im}D+\pi-\frac{t_0+t_1+2\operatorname{Im}B+2\operatorname{Im}K+\pi(J+1)}N,\\
f_1&=\operatorname{Re}D-\frac{q+L_{\mathrm{inner}}-\operatorname{Re}c_d+Q+2\operatorname{Re}B+2\operatorname{Re}R}{2N},\\
f_2&=\operatorname{Re}D_1-\frac{\operatorname{Re}A+L_{\mathrm{outer}}-\operatorname{Re}c_z+\operatorname{Re}B+\operatorname{Re}K+\operatorname{Re}R+Q}{2N},\\
f_3&=\frac{\operatorname{wrap}(2N\operatorname{Im}D_1-\operatorname{Im}A+\phi+\operatorname{Im}c_z-\operatorname{Im}B-\operatorname{Im}K-\operatorname{Im}R-t_1+\pi)}{2N}.
\end{aligned}
```

The unsquared equations fix sigma; the product phase alone cannot recover
it. The analytic Jacobian differentiates these regularized expressions.
The ideal-string initializer reduces the center equation to

```math
N\Theta(a;2\eta)-\sum_{k=1}^{3}\Theta(2a;k\eta)=\pi J.
```

Initialization and the full finite-deviation solve share one Newton budget;
only the latter can certify convergence. Work and storage per step do not
grow with N. This reuses the common string Newton driver and correction
functions, plus the stable inner-pair energy. The outer conjugate pair adds
its real energy contribution directly. Failed results retain consistent
estimates and residuals, not verified energies.

## Threshold and coverage

The thermodynamic droplet formula
$g_{M} =2\,\sinh (\eta)\,\tanh (M \,\eta /2)$ gives **g4=15/7** at Delta=3/2; see
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
