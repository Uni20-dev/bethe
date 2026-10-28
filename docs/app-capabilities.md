# App execution and arbitrary-precision coverage

This inventory separates concurrent execution from arbitrary-precision support.
All numerical apps offer native precision choices, subject to the build; run
`bethe-info` for configured formats and providers. `mpreal` means explicit
`bits-N`/`digits-N` precision, not a promise that a chosen mesh or iteration
budget reaches that accuracy. See [precision controls](command-line.md#select-arithmetic-precision)
and [Async ownership](async-execution.md).

The table records implemented paths and the next concrete audit. An app can
remain serial when its work is too cheap to benefit from scheduling; independent
parameter sweeps can still use the generic library grid. Native-only apps below
do not yet expose MPFR at their frontend, even if some underlying helpers do.
Names omit the `bethe-` executable prefix. Linked pending PRs record separate,
validated checkpoints; they are not capabilities of an installed build until
merged and rebuilt. Each PR records its precision/provider configurations and
independent numerical or serial/concurrent validation.

| Apps | Existing execution path | MPFR frontend support | Next check / dependency |
| --- | --- | --- | --- |
| `xxx-pbc`, `xxx-obc` | Continuous real-root scans; owning Async scan API | Yes, MPFR scalars | Existing concurrent-precision scan tests pass; retain owned solver options through final consumers |
| `xxx-structure-factor` | Async preparation, fused root/form-factor runners, ordered reduction | Yes, MPLAPACK MPFR required | Existing mixed-precision fused/prepared comparisons pass, including failed references; no extra wrapper needed |
| `lieb-liniger-pbc`, `lieb-liniger-obc` | Shared continuous real-root scanner; pending [PR #17](https://github.com/Uni20-dev/bethe/pull/17) exposes owning Async scans | Yes, MPLAPACK MPFR required | Integrate reviewed ownership, downstream-consumer and concurrent-precision checks |
| `richardson`, `central-spin` | Individual solver calls; pending [PR #17](https://github.com/Uni20-dev/bethe/pull/17) supplies an owning grid-sweep example and full-state composition checks | Yes, MPFR scalars | Integrate reviewed example/tests; keep each solve's continuation sequential |
| `xxx-structure-factor-thermo` | Bounded ordered grid, direct table collection | Native only | Precision-aware transition-rate integrals; no new state solver needed |
| `lieb-liniger-dispersion` | Prepared solver and bounded grids; shared Async backgrounds for library clients | Native only on merged baseline; [PR #9](https://github.com/Uni20-dev/bethe/pull/9) adds MPFR | Integrate reviewed thermodynamic precision support |
| `lieb-liniger-thermal` | Serial grid; pending [PR #11](https://github.com/Uni20-dev/bethe/pull/11) adds bounded independent temperature solves | Native only | Explicit-precision quadrature and stopping rules after the thermodynamic precision checkpoint |
| `hubbard-dispersion` | Bounded independent queries; pending [PR #14](https://github.com/Uni20-dev/bethe/pull/14) adds shared Async preparation and one branch-spanning window | Native only | Precision-aware kernels and stable gamma-ratio phases |
| `hubbard-continuum` | Serial searches; pending [PR #12](https://github.com/Uni20-dev/bethe/pull/12) adds bounded independent extrema queries | Native only | Integrate reviewed concurrency; same thermodynamic precision dependencies |
| `xxz-pbc`, `xxz-obc` | Shared continuous real-root scans; owning parameter sweeps validated in pending PR #15 | Native only on merged baseline; [PR #15](https://github.com/Uni20-dev/bethe/pull/15) adds MPFR with MPLAPACK | Integrate reviewed precision support for existing ground/excitation families |
| `q-boson-pbc` | Shared real-root excitation scans | Native only | Runtime precision and complete scan diagnostics |
| `su3-pbc` | Individual nested solves; owning mixed-precision sweeps validated in pending PR #13 | Native only on merged baseline; [PR #13](https://github.com/Uni20-dev/bethe/pull/13) adds MPFR with MPLAPACK | Integrate reviewed balanced-singlet precision support |
| `tb-pbc` | Individual ground-state solve | Native only | Explicit precision, finite string deviations and complex scalar representation; shared Newton helper awaits PR #9 |
| `hubbard-pbc`, `hubbard-obc` | Individual finite-size solves; owning parameter sweeps validated in pending PR #16 | Native only on merged baseline; [PR #16](https://github.com/Uni20-dev/bethe/pull/16) adds MPFR with MPLAPACK | Integrate reviewed precision for existing sectors, both interaction signs and symmetry mappings |
| `gaudin-yang-pbc` | Individual nested solves; owning mixed-precision sector sweeps validated in pending PR #18 | Native only on merged baseline; [PR #18](https://github.com/Uni20-dev/bethe/pull/18) adds MPFR with MPLAPACK | Integrate reviewed precision for existing repulsive shells, free/polarized limits and requested-coupling failure diagnostics |
| `sun-fermions-pbc`, `tj-pbc`, `ladder-pbc` | Composition audit pending | Native only | Shared explicit-precision Newton helper awaits PR #9; audit nested continuation and scratch ownership |
| `bose-fermi-pbc` | Individual nested solves; owning mixed-precision requests validated in pending PR #19 | Native only on merged baseline; [PR #19](https://github.com/Uni20-dev/bethe/pull/19) adds MPFR with MPLAPACK | Integrate reviewed mixed/pure/free-state precision and shared continuum-provider guard |
| `su3-dispersion` | Serial analytic CLI; owning mixed-precision Async queries validated in pending PR #20 | Native only on merged baseline; [PR #20](https://github.com/Uni20-dev/bethe/pull/20) adds MPFR scalars, no MPLAPACK | Integrate reviewed band/bound precision; individual formulas remain too cheap to justify per-point scheduling |
| `haldane-shastry-pbc`, `sutherland-pbc` | Serial analytic CLIs; owning mixed-precision window/parameter scans validated in pending PR #21 | Native only on merged baseline; [PR #21](https://github.com/Uni20-dev/bethe/pull/21) adds MPFR scalars, no MPLAPACK | Integrate reviewed precision and exact-integer gap/order handling; cheap formulas stay serial |
| `tb-dispersion` | Serial analytic CLI; owning mixed-precision Async queries validated in pending PR #22 | Native only on merged baseline; [PR #22](https://github.com/Uni20-dev/bethe/pull/22) adds MPFR scalars, no MPLAPACK | Integrate reviewed shared-band precision; XXZ/XYZ elliptic models still need separate audits |
| `sine-gordon-dispersion` | Serial analytic CLI; owning mixed-precision Async queries validated in pending PR #23 | Native only on merged baseline; [PR #23](https://github.com/Uni20-dev/bethe/pull/23) adds MPFR scalars, no MPLAPACK | Integrate reviewed masses, strict breather selection and pair thresholds; cheap formulas stay serial |
| `xyz-dispersion` | Serial CLI; pending [PR #24](https://github.com/Uni20-dev/bethe/pull/24) reuses owning bound branches and validates shared Async preparation / parameter curves | Native only | Integrate reviewed reuse and measured curve-level concurrency; explicit-precision elliptic/theta functions remain necessary |
| `xxz-dispersion` | Immutable band with cheap queries; bulk preparation remains serial | Native only | Explicit-precision theta constants and bulk quadrature; a per-momentum task is unlikely to help after preparation |
| `sine-gordon-bethe-yang`, `biquadratic-obc`, `xxz-qg-obc` | Sector-specific solver/composition audit pending | Native only | Treat real, complex and singular sectors separately |
| `lee-yang-vacuum`, `lee-yang-excited`, `sine-gordon-vacuum`, `sine-gordon-excited`, `kondo-response` | Integral-equation composition audit pending | Native only | Preserve coupled iterations and continuation; audit discretization and complex functions |
| `asep-pbc`, `tasep-pbc`, `potts-pbc`, `xyz-pbc` | Complex/nonlinear composition audit pending | Native only | Branch tracking and scalar-function coverage before advertising MPFR |

`bethe-info` reports capabilities and has no numerical parallelism to add.
The inventory covers 43 numerical frontends; it does not extend any model's
physical state family. Generic scalar functions and numerical primitives belong
in Uni20; model equations and physical normalization remain in Bethe.

## Validation records

- [Async execution](async-execution.md) describes scan, preparation and grid
  lifetime/failure contracts, including completion of all admitted work.
- [Parallel grids](parallel-grids.md) records independent scalar comparisons,
  precision coverage and benchmark conditions.
- [MPFR builds](building.md#enable-arbitrary-precision) describes provider
  requirements. Native, wide and sanitizer CI exercise different configurations.
- Thermodynamic XXX CLI tests compare complete serial, threaded and non-retained
  rows, including grids larger than the admission window, mixed successful and
  failed points, threshold counts, native precisions and output preservation.

The pending MPFR checkpoints #9, #13, #15, #16, #18, #19, #20, #21, #22 and #23 verify precision growth
against independent references beyond binary128, typed endpoint/failure values,
CLI export/provider contracts and ownership in composed Async calculations.
The pending scheduling checkpoints #11, #12 and #14 compare complete results
against serial/baseline calculations and document measured wall/CPU and memory
behavior. Finite Hubbard's scalar/LU work is independent of the complex
Gamma-ratio functions still required by thermodynamic Hubbard.

The composition audit in pending #17 adds finite Lieb–Liniger owning scan APIs
and a compiled Richardson/central-spin sweep. Deferred input destruction,
downstream consumers, one/multiple workers, failed solves, cancellation and
independent MPFR precisions are covered. Existing XXX scan and structure-factor
composition tests were rerun; their APIs need no new wrapper. No new precision
support or parallel speedup is claimed by that audit.

Gaudin–Yang's pending #18 verifies the unreduced complex Bethe equations and
two-body accuracy growth at 128/256/333 bits. Free states and zero-update
diagnostics also work with MPFR scalars alone at the library level; Newton
updates and the MPFR frontend require MPLAPACK. Mixed interacting populations
remain odd, and continuation within each independent solve remains sequential.

Bose–Fermi's pending #19 validates mixed-species equations, Lieb–Liniger
reductions, free states and missing-energy failures. The shared continuum
Newton guard keeps zero-update diagnostics usable in scalar-only MPFR builds;
actual updates and the frontend require MPLAPACK. Tests resolve a weak-coupling
shift above a nonzero Fermi sea beyond binary128 resolution and compare complete
states under owning, mixed-precision grid composition. Mixed interacting shells
still require odd fermion number.

Haldane–Shastry and Sutherland's pending #21 inherit #20's scalar CLI dispatcher.
The spin model keeps exact motif ordering; the gas keeps exact MPFR integer gap
coefficients until physical scaling. Full MPFR tests run without MPLAPACK,
including differential-Hamiltonian references, large-exponent small gaps,
precision growth, budget refusals and independent owning Async spectra. Native
ED/wavefunction checks and 561 unchanged exported rows validate the existing
formats. Individual analytic CLI calculations remain serial.

TB dispersion's pending #22 uses the same scalar dispatcher from #20 and
extends the common sine/massive-band helper rather than duplicating kinematics.
Independent dressed-sea and continuum references, scalar-only MPFR physics,
positive momenta at the MPFR exponent floor and owning mixed-precision queries
pass. Existing TB/XXZ/XYZ checks and 1,710 identical native output rows validate
shared-band consumers. This does not enable XXZ/XYZ's elliptic models in MPFR.

Sine-Gordon particle dispersion's pending #23 also uses #20's scalar dispatcher.
Independent algebraic breather masses, precision-sensitive stable-species
selection, unequal-mass threshold minimization and full owning Async queries
pass at 128/256/333 bits. Full MPFR physics runs in scalar-only sanitizers;
1,095 native/binary128 output rows match the baseline. The scattering phase,
Bethe–Yang and finite-volume integral equations remain separate precision work.

XYZ's pending #24 prepares the bound coefficient once per branch and reuses it
for both parity copies. Independent elliptic/limiting references, owning Async
readers and failed preparations/queries pass in native/wide/sanitizer builds.
End-to-end binary128 curve time falls from 5.42 to 0.27 seconds with identical
output; an independent 2,097,664-energy parameter sweep falls from 4.480 to
0.585 seconds with eight workers. These are distinct reuse and scheduling
measurements; the CLI stays serial. Commands and limits are in the PR's
`docs/parallel-grids.md`. It does not extend the theta helper to MPFR.

Update this inventory with each reviewed app checkpoint. A proposed capability
or a pending audit is not evidence of implemented support.
