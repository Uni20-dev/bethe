# Finite-temperature Lieb–Liniger equilibrium

The library implements the repulsive Yang–Yang equation in the **grand-canonical
ensemble**: give c>0, T>0, and chemical potential mu; obtain the density rather
than prescribing it. A fixed-density wrapper instead accepts c,T,n and solves
for mu. `bethe-lieb-liniger-thermal` exposes both ensembles, including temperature
scans. Chemical-potential scans are not yet implemented.

The Hamiltonian is `H=-sum d_i²+2c sum delta(x_i-x_j)`, with `hbar=2m=k_B=1`.
Pressure and energy are per length, density is particles per length, and entropy
is also per length. The grand potential per length is **minus the pressure**.
The uniform thermodynamic limit is not a finite-ring spectrum or a trapped gas.

## Command-line examples

```sh
# Grand-canonical equilibrium, with density determined by mu.
build/bethe-lieb-liniger-thermal --c 4 --temperature 1 --mu -1

# Fixed-density temperature scan; mu is solved independently at each point.
build/bethe-lieb-liniger-thermal --c 4 --temperature 0.5 \
  --temperature-end 2 --points 9 --density 1 --max-nodes 512 --csv thermal.csv

# Extended precision and simultaneous structured exports.
build/bethe-lieb-liniger-thermal --c 4 --temperature 1 --density 0.1 \
  --precision fp128 --json thermal.json --tsv thermal.tsv

build/bethe-lieb-liniger-thermal --references
```

Choose exactly one of `--mu` and `--density`. Without `--temperature-end`, the
tool evaluates one point. With it, `--points` defaults to 17 samples, includes
both endpoints, and permits ascending or descending positive temperatures.
The specified c and ensemble parameter stay fixed throughout the scan.

The named table `thermodynamics` reports temperature, chemical potential,
density, pressure, energy/length, entropy/length, quadrature nodes and cutoff,
error estimates, iteration/evaluation counts, and status. Mesh/cutoff diagnostics
are from the final inner solve, not totals; Newton updates and equilibrium
evaluations are totals for that row. The density mismatch applies only to
fixed-density runs. After a failed canonical solve, inner diagnostics and mu
are unavailable; after a failed grand-canonical solve, the supplied mu remains
known. Neither publishes failed physical observables.

The numerical controls mirror the C++ options below: `--tolerance`,
`--initial-nodes`, `--max-nodes`, `--initial-cutoff`, `--max-cutoffs`, and
`--max-iterations`. Fixed-density runs additionally accept `--density-tolerance`
and `--max-evaluations`. Native defaults depend on the selected precision.

The shared [output options](output.md) provide presentation-layer screen output,
CSV/TSV/JSON files, provenance metadata and CPU timing. Rows stream to attached
exports during the scan; `--no-retain` avoids keeping a second copy in memory.
Failure at a point does not stop subsequent temperatures: its physical values
are null in JSON or empty in CSV/TSV, and the final exit code is 2. Invalid input
returns 1, before opening output files. References appear only with `--references`.

## Equations

With `C(x)=c/(pi*(c²+x²))`, the Yang–Yang equation is

```text
epsilon(k) = k²-mu - integral C(k-q) T log(1+exp(-epsilon(q)/T)) dq
filling(k) = 1/(1+exp(epsilon(k)/T))
rho_total(k) = 1/(2*pi) + integral C(k-q) filling(q) rho_total(q) dq.
```

All integrals cover the real axis. The filling has a Fermi-like form because
it counts occupied Bethe quantum numbers; the physical particles remain bosons.
Pressure integrates `T log(1+exp(-epsilon/T))/(2*pi)`. Particle and energy
densities integrate `filling*rho_total` and `k²*filling*rho_total` respectively.
Entropy integrates the binary entropy of the filling times `rho_total`.

Sources: [Yang–Yang (1969)](../CITATIONS.md#yang-yang-1969),
[Guan–Batchelor (2011)](../CITATIONS.md#guan-batchelor-2011), equations (2)–(5),
and [Franchini's notes](../CITATIONS.md#franchini-2011), section 2.11.
We solve the nonlinear equation numerically; no finite-c polylog approximation
is substituted. The ideal-Fermi fugacity series supplies an independent test
in the Tonks limit.

## C++ interface

```cpp
#include <bethe/lieb_liniger_thermal.hpp>

namespace yy = bethe::lieb_liniger::thermal;
auto state = yy::equilibrium(4.0, 1.0, -1.0); // c, T, mu
if (state.converged) {
    auto pressure = *state.pressure;
    auto density = *state.density;
    auto energy_per_length = *state.energy_per_length;
    auto entropy_per_length = *state.entropy_per_length;
}
```

Use `long double` or enabled Uni20 fp128 scalars for extended arithmetic.
There is no narrowing to fp64 inside the solver or its quadrature rules.
The default `Options<Real>` has relative tolerance `16384*epsilon`, 16 initial
nodes, at most 256 nodes, 512 accepted Newton updates across the whole solve,
and at most four rapidity cutoffs. Node budgets may be increased to 512;
cutoff budgets to 32. Zero iteration/cutoff budgets are allowed for diagnostics.

`initial_cutoff`, if specified, is a positive physical rapidity. Otherwise the
initial domain is chosen from T, mu, and the requested tolerance. Successive
cutoffs grow by 5/4. Mesh sizes double at each fixed cutoff, starting again at
`initial_nodes`; a non-doubling `max_nodes` need not itself be sampled.

## Three distinct convergence checks

1. A damped Newton solve converges the discretized pseudoenergy equation.
   Its Jacobian is also the dressed linear operator for `rho_total`.
2. Successive Gauss–Legendre meshes must agree for **all four** thermodynamic
   observables at a fixed cutoff.
3. Two separately mesh-converged cutoffs must agree for all four observables.

`nodes` counts the positive-half nodes: evenness supplies the negative half
through `C(k-q)+C(k+q)`. `cutoff` is physical; `cutoffs` counts attempted domains.
`mesh_error` and `cutoff_error` are eight times the largest relative observable
change. They are estimates, not rigorous bounds on omitted tails or quadrature.
The library reports no physical observable until both checks pass. Small c,
low T, or a large user-selected initial cutoff may need more nodes than allowed.

Calculations scale energies by `S=max(T,abs(mu))` and rapidities by `sqrt(S)`.
`nonlinear_residual` is the largest scaled-equation residual divided by
`max(1,abs(epsilon),abs(k²-mu))` in these dimensionless variables. Its target is
`tolerance*(T/S)/64`, so the equation's precision requirement tracks temperature.

Stable shared thermal factors avoid exponentially overflowing expressions.
Entropy is calculated directly, not by subtracting `p+e-mu*n` at low T.
Native Gauss–Legendre rules, compensated sums, rational scattering kernels, and
recoverable pivoting solves are shared with existing continuum/Hubbard solvers.

`iteration_limit`, `mesh_limit`, `cutoff_limit`, and `precision_limit` distinguish
failure paths; all four optional observables remain empty on failure. A failed
coarse mesh is retried at higher order; if the last mesh cannot resolve the
arithmetic, the result is `precision_limit`. Invalid parameters throw.
No asymptotic zero-density or T=0 answer is silently substituted on underflow.
Use the separate [zero-temperature API](lieb-liniger-thermo.md) for T=0.

Validation includes the original integral equations, scale covariance,
`p+e-mu*n=T*s`, pressure derivatives at fixed c, a native-precision fugacity
series in the Tonks limit, decreasing-temperature comparison with the ground
state, and independently exhausted nonlinear, mesh, and cutoff budgets.

## Fixed density

```cpp
auto canonical = yy::at_density(4.0, 1.0, 1.0); // c, T, n
if (canonical.converged) {
    auto mu = canonical.state->chemical_potential;
    auto free_energy_per_length = mu * *canonical.state->density
                               - *canonical.state->pressure;
}
```

`DensityOptions<Real>` contains the inner `equilibrium` options, an outer
relative density `tolerance` (default `65536*epsilon`), and `max_evaluations`
(default 128). Each evaluation is a complete grand-canonical solve with its
own inner budgets; the inner tolerance is tightened to at most one eighth of
the density tolerance. `evaluations` counts those solves and `iterations`
sums their accepted Newton updates.

The chemical potential starts at the classical-gas estimate, but every density
used is calculated from the full Yang–Yang equation. Exponentially expanding
steps bracket the target; Illinois regula falsi then interpolates inside that
bracket, falling back to its midpoint when rounding places a trial at an end.
The stopping test requires a relative density mismatch at most half the outer
tolerance, leaving room for the inner numerical error estimate. Small bracket
width alone is not accepted as convergence.

`density_error` is the last successfully evaluated relative density mismatch.
`density_limit` means the outer evaluation budget was exhausted. Inner failures
retain their original status; no failed density is treated as zero and no
partially converged thermodynamics are published. The optional `state` is set
only on success, and its density is the calculated value, not overwritten by
the requested density. Extreme parameters can still encounter the inner mesh
or precision limits described above. Tests cover native-precision round trips,
scale covariance, dilute and degenerate gases, and independent outer/inner
budget exhaustion.
