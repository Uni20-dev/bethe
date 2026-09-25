# Lieb–Liniger at zero temperature

The library now solves the repulsive bulk ground state at fixed density.
This is the background needed for Lieb's type-I and type-II excitation curves;
the curves and their frontend are a subsequent checkpoint, not yet implemented.

Use the same units as the [finite ring](lieb-liniger.md):
`H = -sum d²/dx² + 2c sum delta(x_i-x_j)`, with `c>0` and density `n>0`.
The bulk energy is extensive: the result is **energy per length**, not energy
per particle. Ring and hard-wall systems have the same bulk limit; boundary
energies are not included here.

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
```

The scalar type can also be `long double` or Uni20's enabled fp128 type.
`Options<Real>` controls `tolerance`, `initial_nodes`, `max_nodes`, and the
total density-iteration budget `max_iterations` across all meshes.
Defaults are `4096*epsilon`, 16, 256, and 128. Mesh sizes double; the final
size need not equal `max_nodes` if it is not on that doubling sequence.
The supported node budget is 4–512.

Convergence requires both normalization and agreement of successive resolved
meshes for Q, energy, and chemical potential. `mesh_error` is eight times
their largest relative change; `density_error` is the relative normalization
residual. These are numerical diagnostics, not rigorous error bounds.
A mesh must also resolve the scattering-kernel width before it is eligible
for the comparison. A single mesh can never report success.

Failure statuses distinguish `mesh_limit`, `density_limit`, and
`precision_limit`. Under-resolved linear systems, including loss of positive
root density, trigger mesh refinement; exhaustion returns `mesh_limit`.
No physical observable is populated on failure. Invalid parameters throw.
Very weak coupling may exceed the mesh budget; no weak-coupling approximation
is silently substituted. Overflow or underflow in physical-unit conversion
also fails explicitly. Exactly zero/infinite c, attraction, finite temperature,
and boundary corrections are outside this API.

Tests cover native polynomial quadrature, density and dressed-energy equations,
finite-ring extrapolation, scale covariance, the thermodynamic derivative
`mu = d(E/L)/dn` at fixed c, weak/strong limits, and failure budgets.
