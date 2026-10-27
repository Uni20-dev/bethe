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
Names omit the `bethe-` executable prefix.

| Apps | Existing execution path | MPFR frontend support | Next check / dependency |
| --- | --- | --- | --- |
| `xxx-pbc`, `xxx-obc` | Continuous real-root scans; owning Async scan API | Yes, MPFR scalars | Independent sector/sweep composition and configuration coverage |
| `xxx-structure-factor` | Async preparation, fused root/form-factor runners, ordered reduction | Yes, MPLAPACK MPFR required | Preserve determinant precision and failure diagnostics in composed sweeps |
| `lieb-liniger-pbc`, `lieb-liniger-obc` | Shared continuous real-root scanner | Yes, MPLAPACK MPFR required | Public composition of independent solves/scans |
| `richardson`, `central-spin` | Individual solver calls; generic grid can compose independent solves | Yes, MPFR scalars | Parameter sweeps while retaining sequential continuation inside each solve |
| `xxx-structure-factor-thermo` | Bounded ordered grid, direct table collection | Native only | Precision-aware transition-rate integrals; no new state solver needed |
| `lieb-liniger-dispersion` | Prepared solver and bounded grids; shared Async backgrounds for library clients | Native only on merged baseline; [PR #9](https://github.com/Uni20-dev/bethe/pull/9) adds MPFR | Integrate reviewed thermodynamic precision support |
| `lieb-liniger-thermal` | Serial temperature grid | Native only | Independent temperature solves and explicit-precision quadrature/stopping rules |
| `hubbard-dispersion` | Bounded independent momentum queries | Native only | Shared Async preparation; precision-aware kernels and stable gamma-ratio phases |
| `hubbard-continuum` | Serial momentum searches | Native only | Independent extrema queries; same thermodynamic precision dependencies |
| `xxz-pbc`, `xxz-obc` | Shared continuous real-root scans | Native only | Explicit-precision root solving and owning composition |
| `q-boson-pbc` | Shared real-root excitation scans | Native only | Runtime precision and complete scan diagnostics |
| `su3-pbc`, `tb-pbc` | Individual nested ground-state solves | Native only | Precision controls and independent solves; supported sectors remain unchanged |
| `hubbard-pbc`, `hubbard-obc` | Individual finite-size solves | Native only | Runtime-precision finite equations; independent of the thermodynamic gamma-ratio work |
| `sun-fermions-pbc`, `gaudin-yang-pbc`, `bose-fermi-pbc`, `tj-pbc`, `ladder-pbc` | Composition audit pending | Native only | Related real/nested solvers, continuation and scratch ownership |
| `su3-dispersion`, `tb-dispersion`, `haldane-shastry-pbc`, `sutherland-pbc` | Analytic/formula-based paths | Native only | Low-cost precision extensions; measure before parallelizing individual points |
| `xxz-dispersion`, `xyz-dispersion`, `sine-gordon-dispersion` | Independent queries; prepared-context audit pending | Native only | Elliptic/theta/integral accuracy and immutable scratch ownership |
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

Update this inventory with each reviewed app checkpoint. A proposed capability
or a pending audit is not evidence of implemented support.
