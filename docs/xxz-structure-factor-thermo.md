# Massive XXZ two-spinon spectral intensity

[XXZ family](xxz-family.md) · [Spinon dispersion](xxz-dispersion.md) ·
[Output and files](output.md)

`bethe-xxz-structure-factor-thermo` calculates the zero-temperature, zero-field
**transverse two-spinon contribution** for the infinite antiferromagnetic XXZ
chain with $`\Delta\gt 1`$. It adds operator-dependent intensity to the continuum
kinematics supplied by `bethe-xxz-dispersion`.

```math
H=J\sum_j\left(S_j^xS_{j+1}^x+S_j^yS_{j+1}^y+\Delta S_j^zS_{j+1}^z\right),
\qquad S^\alpha=\sigma^\alpha/2,\quad J\gt 0.
```

The result is unbroadened and keeps its physical normalization. Four and higher
spinons are omitted; the two-spinon result is never rescaled to exhaust a sum
rule. The longitudinal response is a different calculation at finite anisotropy.

## Calculate a cut or grid

```sh
bethe-xxz-structure-factor-thermo --delta 1.5 --momentum 1.5 --omega 2
bethe-xxz-structure-factor-thermo --delta 1.5 --momentum 1.5 --points 501 --csv cut.csv
bethe-xxz-structure-factor-thermo --delta 2 --momentum-points 129 --points 501 --threads 4 --json grid.json
bethe-xxz-structure-factor-thermo --delta 1.5 --momentum 1.5 --points 501 --tsv-table continuum=bounds.tsv
bethe-xxz-structure-factor-thermo --references
```

The first command gives approximately 4.15974928658 for the spectral density.
Momenta are in radians per site, in [0,2π]. `--momentum` requests a single
momentum; otherwise `--momentum-points` samples the full interval, including both
periodic endpoints. Count that seam once when integrating or plotting a periodic
grid. Frequencies are excitation energies above the ground state.

`--omega` requests one frequency. Otherwise `--points`, `--omega-min` and
`--omega-max` define a uniform frequency grid; the default upper endpoint is
twice the maximum single-spinon energy. A rectangular grid has at most one
million points. Uniform grids are useful for plots, but endpoint-aware quadrature
is preferable for integrated sum rules.

The frontend declares two tables:

| Table | Columns and purpose |
| --- | --- |
| `spectrum` (primary) | `k`, `omega`, `density`, `estimated_error`, `evaluations`, `status` |
| `continuum` | `k`, `sheet`, `present`, `lower`, `upper`; the two spectral sheets separately |

`--csv-table continuum=FILE` writes those bounds independently of the screen
selection. Unknown table names are errors. JSON/CSV/TSV, metadata, streaming,
`--no-retain`, overwrite protection and screen selection follow the
[common output contract](output.md). `--threads` uses the shared scheduler for
bounded, ordered blocks of independent points; it does not change their weights.
The caller collects rows in momentum-major order. CPU time covers the spectral
grid and ordered row output.

## Momentum, channel and normalization

The response is $`S^{-+}_2(k,\omega)`$, with the raising operator on the ground
state. In finite-volume notation before taking the thermodynamic limit,

```math
S_k^+=\frac{1}{\sqrt{N}}\sum_j e^{-ikj}S_j^+,
\qquad
S^{-+}_2(k,\omega)=2\pi\sum_{n\in\text{two spinons}}
\left|\langle n|S_k^+|0\rangle\right|^2\delta(\omega-E_n+E_0).
```

The density column is **S itself**, so frequency integrals use
$`d\omega/(2\pi)`$. The complete transverse response, averaged over a full
Brillouin zone, has zeroth moment 1/2. The two-spinon part is smaller:

```math
\int_0^{2\pi}\frac{dk}{2\pi}
\int_0^\infty\frac{d\omega}{2\pi}S^{-+}_2(k,\omega)\leq\frac12.
```

The massive-chain expression uses the translation-symmetric average of the two
Néel vacua. Its physical momentum is $`k=\pi-Q`$ relative to the convention in
[Caux–Mossel–Pérez Castillo](../CITATIONS.md#caux-mossel-perez-castillo-2008).
Both spectral sheets contribute with **unequal weights**. Folding their energy
supports into a two-site Brillouin zone does not make the spectral intensity
π-periodic. The implementation retains the full physical momentum convention.

Only in the isotropic limit does SU(2) relate this transverse response to the
[XXX longitudinal two-spinon result](xxx-structure-factor-thermo.md) by a factor
of two. Do not apply that conversion at general Δ. Exchange scaling is
$`S_J(k,\omega)=S_1(k,\omega/J)/J`$.

## Numerical status and precision

| Status | Meaning |
| --- | --- |
| `converged` | Interior density and an estimated numerical error are available |
| `outside_continuum` | Exactly zero contribution from both sheets at the supplied point |
| `threshold` | The supplied frequency equals a sheet endpoint at working precision; density and error are absent |
| `numerical_failure` | Budget or arithmetic resolution is insufficient; density and error are absent |

An exact endpoint is not assigned a finite value, even when a one-sided limit
exists. An internal sheet endpoint also receives `threshold`: the other sheet
may be regular there, but the complete point value is deliberately left absent.
Thresholds and outside-support points are valid output statuses. Numerical
failures retain their rows and make the frontend exit with status 2.

Native `fp64`, `long-double` and configured `fp128` are supported. With MPFR,
`--precision bits-N` or `digits-N` constructs all continuous inputs at the
requested precision; no ambient precision is changed. MPLAPACK is not needed
for this scalar kernel. For example:

```sh
bethe-xxz-structure-factor-thermo --delta 1.001 --momentum 1.5 --omega 1.8 --precision bits-160
```

`--tolerance` controls the **log-kernel series tail**, with a default of 4096
machine epsilons and a minimum of 64 epsilons. `--max-terms` bounds that series
separately for each contributing sheet. Neither is a promise of a relative
error in the final density. The exported error combines its tail estimate with
an arithmetic-conditioning estimate; it is not a rigorous enclosure, an input
uncertainty estimate, or an estimate of omitted spinon sectors.

Near an edge, subtraction in the inverse dispersion amplifies roundoff. Near
Δ=1, the mass becomes exponentially small. Use more precision when those scales
matter. An unrepresentable mass or bandwidth is rejected during preparation;
unresolved interior arithmetic is not silently turned into zero density.

## C++ interface

```cpp
#include <bethe/xxz_thermodynamic_structure_factor.hpp>

bethe::xxz::ThermodynamicTwoSpinonStructureFactor<double> spectrum(1.5);
auto point = spectrum(1.5, 2.0);
auto sheets = spectrum.sheets(1.5);
if (point.value) {
    // Use *point.value, with *point.error and point.status.
}
```

The immutable preparation can be shared among independent queries. For MPFR,
construct Δ, finite exchange and query coordinates at the same explicit
precision. The default exact exchange 1 is accepted. Optional controls use
`TwoSpinonStructureFactorOptions<uni20::mpreal>(precision)` and a tolerance at
that precision. `converged()` on a returned point means there was no numerical
failure; inspect `value` as well, because a threshold deliberately has none.
The [observable-grid API](parallel-grids.md) composes queries with Async and
waits on the final result; no scheduler-wide drain is required.

## Implementation and independent checks

The spectral sheets and weights follow published Eqs. (44), (51), (56) of
[Caux–Mossel–Pérez Castillo (2008)](../CITATIONS.md#caux-mossel-perez-castillo-2008).
The implementation rationalizes cancelling Jacobian factors, inverts the
elliptic map with a restricted Carlson RF helper, and uses Uni20's scaled theta
constants/jets. It sums the nondecaying harmonic part of the transition kernel
analytically and bounds the remaining geometrically decaying series tail.

At a sheet's closing momentum its width vanishes quadratically. Independently
rounding the upper and lower bounds can invert that narrow interval. On the
branch with lower angular bound $`\ell=(1+k')sc`$, the code instead uses
$`u=\sqrt{\ell^2+(s^2-k'c^2)^2}`$, where $`s^2+c^2=1`$ and $`k'`$ is the
single-spinon gap divided by its maximum energy. Both bounds are then scaled
by twice that maximum energy. This keeps the interval ordered; a width below
the chosen precision may round to zero. Increase precision to resolve it.

The standalone `scripts/reference_xxz_two_spinon.py` provides a direct mpmath
elliptic/theta implementation and a separate original q-product vertex route
from [Abada–Bougourzi–Si-Lakhal](../CITATIONS.md#abada-bougourzi-si-lakhal-1997),
Eqs. (2.6)–(2.7) restricted to two spinons. The latter searches momentum roots
and computes the Jacobian directly. Its mesh is an interior-point control,
not a general threshold root finder. A NumPy integral over the complete
spectral-angle domain checks the two-vacuum measure before resolving k and ω.
It is restricted to Δ≥1.2 and has explicit product/grid controls.

With mpmath and NumPy installed, run from the source directory:

```sh
python3 scripts/reference_xxz_two_spinon.py --delta 1.5 --momentum 1.5 --omega 2 --digits 120
python3 -m unittest discover -s scripts -p test_reference_xxz_two_spinon.py -v
```

Regression tests include independent point values, near-edge references down
to relative offsets of 10^-24 in wider arithmetic, both sheets, exchange and
momentum symmetries, explicit failures, and MPFR precision ownership. Integrated
checks compare the resolved density against the original vertex measure,
approach the independent XXX result from Δ>1, and test saturation towards the
Ising limit. These controls do not certify arbitrary points or replace the
reported numerical diagnostics. [Validation provenance](validation/xxz-two-spinon-spectral-weight.json)
records the tested configurations and source correspondence.
