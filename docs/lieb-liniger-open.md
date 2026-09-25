# Lieb–Liniger bosons in a hard-wall box

[Periodic gas](lieb-liniger.md) · [Model catalogue](models.md)

The library supports repulsive bosons with Dirichlet walls at x=0 and x=L,
in the same units as the ring: `H=-sum d_j^2 + 2c sum_(i<j) delta(x_i-x_j)`.
Both L and c must be finite and positive. Ground states, specified real-root
states and finite-window excitation scans work in fp64, native long double
and enabled fp128. The executable is `bethe-lieb-liniger-obc`.

## Start with a calculation

```sh
build/bethe-lieb-liniger-obc 4 --length 4 --c 1 --roots
build/bethe-lieb-liniger-obc 4 --length 4 --c 1 --excitations all --padding 1
build/bethe-lieb-liniger-obc 4 --length 4 --c 1 --quantum-numbers 1,2,4,5 --precision fp128
build/bethe-lieb-liniger-obc 4 --length 4 --c 1 --excitations 2 --padding 1 \
  --roots --json box.json --csv-table roots=box.roots.csv
```

Both `--length` and `--c` are required. The second example solves all five
configurations of four occupied labels in `1,...,5`, including the ground
state. `all` refers to this window, not the infinite spectrum. The same command
on the ring frontend instead has 15 candidates because ring padding is two-sided.

The shared frontend supplies CPU timing, `--references`, precision selection,
JSON/CSV/TSV, and optional `--no-retain` streaming. Main rows include energy,
gap (for ranked scans), residual, iteration count and convergence. Unlike the
ring, there are no total-momentum columns. `--roots` adds positive wave numbers
and integer labels. Reference and failed rows have distinct `state_id` values
for joins with roots. Unconverged energies are estimates; a partial solve exits
with status 2 and scan failures are not ranked. Invalid input exits with status 1.

## Library

```cpp
#include <bethe/lieb_liniger_open.hpp>
namespace box = bethe::lieb_liniger::open;
auto ground = box::ground_state<long double>(4, 4.L, 1.L);
auto excited = box::solve_real<long double>(4.L, 1.L,
    {uni20::half_int{1}, uni20::half_int{2}, uni20::half_int{4}, uni20::half_int{5}});
auto scan = box::real_excitations<long double>(4, 4.L, 1.L, 2, {.count=10});
// Inspect ground.converged, excited.converged and scan.converged().
```

## Labels, reflection and energy

The hard-wall Bethe equations originate with
[Gaudin (1971)](../CITATIONS.md#gaudin-1971). In dimensionless variables
q=kL and g=cL we use

```math
q_j+\sum_{l\ne j}\left[\arctan\!\left(\frac{q_j-q_l}{g}\right)
+\arctan\!\left(\frac{q_j+q_l}{g}\right)\right]=\pi I_j.
```

The second phase is reflected scattering. The sum excludes j=l in **both**
terms: there is no self-image scattering at 2q_j. This agrees with the
ground-state equation (8) of
[Reichert et al. (2019)](../CITATIONS.md#reichert-2019), where the self-image
term in the all-index sum is explicitly cancelled.

The labels are strictly increasing positive integers for either particle-number
parity, with ground labels 1,...,N. The positive roots satisfy
`0<k_1<...<k_N` and `E=sum k_j^2`. They are standing-wave Bethe coordinates:
there is no conserved total momentum, so the box State intentionally has no
`momentum` or `momentum_index` field. Reflecting the ring's roots and doubling
its size would introduce a different self-scattering equation.

## Excitations and numerical contract

Padding P means choose N labels from `1,...,N+P`, giving `choose(N+P,N)`
candidates. This is a **single upper edge**, unlike the ring's two-edge padding.
The ground state is included. COUNT retains the lowest converged results after
solving the entire window; it does not reduce the solve count. The common scan
budget rejects oversized windows before solving. No finite window is the full
continuum spectrum. Vacuum N=0 is allowed and has one empty-label state.

The two boundary geometries share the physical-domain damped Newton driver,
stable rational scattering kernel, parameter validation, compensated sums and
finite-window enumerator, as well as the CLI/reporting implementation. The hard-wall model owns its reflection phases,
Jacobian, labels, seed and observables. For g<1 it removes the exact rank phase
analytically, using target $`\pi \,(I_{j} -j)`$ with zero-based j and complementary
arctangents. The reported residual is
`max_j |F_j|/max(1,|target_j|,|q_j|)`, not an energy-error bound.

`max_iterations=0` evaluates the seed only. Initial wave numbers may be supplied
for restart. After conversion back to physical k, the residual is recomputed
from those returned coordinates. Distinct positive roots must remain representable;
at extremely weak coupling, roots near the same free-boson box mode can become
indistinguishable in finite precision. Such a seed is rejected, rather than
silently treated as a valid coincident-root solution. Increase precision when
needed. Unconverged results remain estimates and scans keep failure diagnostics.
Attraction, exactly c=0/infinity, general Robin walls and thermodynamics are not
implemented by this solver.

## Independent checks

For N=1, k=pi I/L exactly. For N=2, A=q_2+q_1 and D=q_2-q_1 independently solve
$`x +2 \arctan (x /g)=\pi K`$, with K=I_1+I_2 and I_2-I_1 respectively. Tests bisect these
scalar equations without the many-body Newton system. They also check the
original complex exponential equations, every analytic Jacobian column, finite
windows and ordering, restart/budget behavior and length scaling.

The ground energy tends to $`N \,\pi ^{2}/L ^{2}`$ at weak repulsion (with first-order
shift $`3c \,N \,(N -1)/(2L)`$), and to `pi^2*sum(j^2)/L^2` in the impenetrable limit.
Ground-state tests cover up to N=32 and couplings from 1e-8 to 1e8 in all
enabled precisions. These checks do not claim that every extreme parameter
combination or arbitrarily large excited label is numerically resolved.
