# Lieb–Liniger at zero temperature

`bethe-lieb-liniger-dispersion` gives the repulsive bulk ground state and
Lieb's type-I and type-II excitation curves at fixed density. The underlying
C++ API supports fp64, native long double, and enabled fp128 throughout.

Use the same units as the [finite ring](lieb-liniger.md):
`H = -sum d²/dx² + 2c sum delta(x_i-x_j)`, with `c>0` and density `n>0`.
The bulk energy is extensive: the result is **energy per length**, not energy
per particle. Ring and hard-wall systems have the same bulk limit; boundary
energies are not included here.

## Command-line examples

```sh
bethe-lieb-liniger-dispersion --c 4 --points 65 --csv curves.csv
bethe-lieb-liniger-dispersion --c 1 --density 0.5 --branch type-i --p-max 8
bethe-lieb-liniger-dispersion --c 4 --branch type-ii --momentum 1 --precision fp128
bethe-lieb-liniger-dispersion --c 4 --points 1001 --json curves.json --quiet --no-retain
```

By default both curves have 33 equally spaced physical momenta from 0 to
`2*pi*n`. Type I actually extends to infinity; `--p-max` changes its grid
endpoint without changing the type-II interval. `--momentum` instead evaluates
one requested p on each selected branch. It cannot be combined with `--points`
or `--p-max`. Negative-momentum curves follow by parity; this interface uses p>=0.

The `dispersion` table includes `p`, `p_over_pi_n`, `energy`, the particle/hole
`rapidity`, its positive `edge_distance`, error estimates, inversion iterations,
and status. Metadata contains Q, ground energy per length, chemical potential,
coupling, density, precision, solver controls, and provenance; the final summary
includes CPU time. All [shared output options](output.md) apply. File exports
and machine stdout receive each row as it is calculated; `--no-retain` does not
accumulate an extra result vector. Human output uses the Uni20 presentation layer.

Unconverged energies and rapidities are null in JSON and empty in CSV/TSV,
with exit status 2. A completed output document does **not** mean every point
converged: inspect its scientific outcome and each row's status.

## Equations and implementation

Occupied rapidities fill `[-Q,Q]`. The root density satisfies

```text
rho(k) - integral[-Q,Q] c rho(k') / (pi (c²+(k-k')²)) dk' = 1/(2 pi)
n = integral[-Q,Q] rho(k) dk
E/L = integral[-Q,Q] k² rho(k) dk.
```

The [original solution](../CITATIONS.md#lieb-liniger-1963) and
[Caux's derivation](../CITATIONS.md#caux-lieb-liniger) explain why the
normalization determines the Fermi rapidity Q. Q is a rapidity endpoint,
not the physical excitation momentum.

The dressed energy obeys the same integral operator with source `k²-mu`.
We solve separately for the dressed `k²` and unit sources and impose
`epsilon(Q)=0` to obtain the chemical potential. See also equations
(2.87)–(2.88) of [Franchini's lecture notes](../CITATIONS.md#franchini-2011).
The unit source is simply `2 pi rho`, so it needs no additional linear solve.

For our positive scattering phase `theta(x)=2 atan(x/c)`, define

```text
p_dressed(k) = k + integral[-Q,Q] theta(k-k') rho(k') dk'
p_dressed(Q) = pi*n.
```

Type I moves a particle from the right Fermi edge to `k>=Q`:
`p=p_dressed(k)-pi*n`, `E=epsilon(k)`, with `p>=0`.
Type II leaves a hole at `-Q<=k<=Q` and puts a particle at the right edge:
`p=pi*n-p_dressed(k)`, `E=-epsilon(k)`, with `0<=p<=2*pi*n`.
These are **fixed-N excitation gaps**, not energies for adding/removing a bare
particle: no chemical-potential offset should be added to the exported energy.
See [Lieb II](../CITATIONS.md#lieb-1963-excitations), and Caux's
[type-I](https://integrability.org/g_l_e_I.html) and
[type-II](https://integrability.org/g_l_e_II.html) derivations.

Both branches are gapless at p=0. Type II is symmetric about p=pi*n and also
vanishes at p=2*pi*n. In the Tonks limit they become `p*(2*pi*n+p)` and
`p*(2*pi*n-p)`, respectively. For very small positive p, the code solves for
the distance from Q and factors energy/phase differences by that distance.
Thus the gap can remain nonzero even when the reported absolute rapidity
rounds to Q. `edge_distance` retains that information. Type-II momenta above
pi*n are evaluated by reflection about the branch midpoint.

The two converged background meshes are separately inverted at the **same
physical p**. Point energies have an additional relative mesh-error test;
background convergence alone does not certify an entire excitation curve.
The reported `energy_error` combines mesh differences, momentum inversion
residuals propagated with the local velocity, and a floating-point allowance.
`momentum_error` is the larger inversion residual in physical units, not an
independent rigorous bound on the integral equation. Exact branch endpoints
are analytic zero-energy points once the background is converged.

All calculations use `gamma=c/n` and rapidities measured in units of n.
Thus `Q=n Qbar`, `E/L=n³ ebar`, and `mu=n² mubar`. The solver shares native-real
Gauss–Legendre quadrature with the doped Hubbard implementation, the pivoting
linear solver with other Bethe systems, and the stable rational scattering
kernel with finite continuum solvers.

## C++ API

```cpp
#include <bethe/lieb_liniger_thermo.hpp>

auto result = bethe::lieb_liniger::thermo::ground_state(4.0, 1.0);
if (result.converged) {
    auto q = *result.fermi_rapidity;
    auto energy_per_length = *result.energy_per_length;
    auto mu = *result.chemical_potential;
}

namespace ll = bethe::lieb_liniger::thermo;
ll::Solver<double> solver(4.0, 1.0); // cache both converged background meshes
auto particle = solver.at_momentum(ll::Branch::type_i, 1.0);
auto hole = solver.at_momentum(ll::Branch::type_ii, 1.0);
// Check point.converged before dereferencing point.energy.
```

The scalar type can also be `long double` or Uni20's enabled fp128 type.
`Options<Real>` controls `tolerance`, `initial_nodes`, `max_nodes`, and the
total density-iteration budget `max_iterations` across all meshes.
Defaults are `4096*epsilon`, 16, 256, and 128. Mesh sizes double; the final
size need not equal `max_nodes` if it is not on that doubling sequence.
The supported node budget is 4–512.
`max_momentum_iterations` (default 160) limits inversion updates for both
meshes combined, separately for each point. At the CLI these budgets are
`--max-background-iterations` and `--max-iterations`, respectively.

Convergence requires both normalization and agreement of successive resolved
meshes for Q, energy, and chemical potential. `mesh_error` is eight times
their largest relative change; `density_error` is the relative normalization
residual. These are numerical diagnostics, not rigorous error bounds.
A mesh must also resolve the scattering-kernel width before it is eligible
for the comparison. A single mesh can never report success.

Failure statuses distinguish `mesh_limit`, `density_limit`, `momentum_limit`, and
`precision_limit`. Under-resolved linear systems, including loss of positive
root density, trigger mesh refinement; exhaustion returns `mesh_limit`.
No physical observable is populated on failure. Invalid parameters throw.
Very weak coupling may exceed the mesh budget; no weak-coupling approximation
is silently substituted. Overflow or underflow in physical-unit conversion
also fails explicitly. Exactly zero/infinite c, attraction, finite temperature,
and boundary corrections are outside this API.
For T>0 at fixed chemical potential, use the separate
[Yang–Yang equilibrium library](lieb-liniger-thermal.md).

Tests cover native polynomial quadrature, density and dressed-energy equations,
finite-ring extrapolation, scale covariance, the thermodynamic derivative
`mu = d(E/L)/dn` at fixed c, weak/strong limits, and failure budgets.
The excitation tests additionally compare an independent backflow linear
system, increasing finite rings with explicitly changed Bethe labels, Tonks
curves, density scaling, reflection, and tiny positive gaps in every precision.
The frontend tests exercise all three export formats, metadata, streaming,
native input/output, invalid arguments, and numerical failures.

These elementary curves do not supply form factors, spectral weights, a full
multiparticle enumeration, or a momentum-resolved hard-wall spectrum.
