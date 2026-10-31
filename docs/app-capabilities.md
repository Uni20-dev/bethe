# App execution and arbitrary-precision coverage

This inventory separates concurrent execution from arbitrary-precision support.
All numerical apps offer native precision choices, subject to the build; run
`bethe-info` for configured formats and providers. `mpreal` means explicit
`bits-N`/`digits-N` precision, not a promise that a chosen mesh or iteration
budget reaches that accuracy. See [precision controls](command-line.md#select-arithmetic-precision)
and [Async ownership](async-execution.md).

The table records implemented execution paths and numerical limits. An app can
remain serial when its work is too cheap to benefit from scheduling; independent
parameter sweeps can still use the generic library grid. All 44 numerical
frontends below expose explicit MPFR when that scalar family is configured.
Names omit the `bethe-` executable prefix. This inventory describes the current
source tree. Rebuild an older
installation to obtain these capabilities. Linked PRs preserve the original
precision/provider configurations and independent numerical or serial/concurrent
validation; their measurements are historical records, not new integration runs.

| Apps | Existing execution path | MPFR frontend support | Scope and numerical limits |
| --- | --- | --- | --- |
| `xxx-pbc`, `xxx-obc` | Continuous real-root scans; owning Async scan API | Yes, MPFR scalars | Owning options survive through final consumers; real-root families only |
| `xxx-structure-factor` | Async preparation, fused root/form-factor runners, ordered reduction | Yes, MPFR scalars; CPU LU available | Two-spinon family; failed references and factors remain explicit |
| `lieb-liniger-pbc`, `lieb-liniger-obc` | Shared continuous real-root scanner; owning Async scans ([#17](https://github.com/Uni20-dev/bethe/pull/17)) | Yes, MPFR scalars; CPU LU available | Finite label windows; independent requests can use different precisions |
| `richardson`, `central-spin` | Individual solver calls; owning grid-sweep example and full-state composition checks ([#17](https://github.com/Uni20-dev/bethe/pull/17)) | Yes, MPFR scalars | Continuation within each solve remains sequential |
| `xxx-structure-factor-thermo` | Bounded ordered two-spinon grid; shared Async sector, node and replicate work | Yes, scalar MPFR | Two spinons by default; selected 2/4/6 partial sums with separate quadrature, sampling, contour, arithmetic and cutoff diagnostics |
| `xxx-moments-thermo` | Shared Async sector queries; bounded node and replicate work | Yes, scalar MPFR | Direct fixed-momentum/global M0/M1; physical channel/exchange conversion and optional six-spinon covariance; finite domains require refinement |
| `lieb-liniger-dispersion` | Prepared solver and bounded grids; shared Async backgrounds | Yes, MPFR scalars; CPU LU available ([#9](https://github.com/Uni20-dev/bethe/pull/9)) | Arithmetic precision does not replace mesh, cutoff or inversion convergence |
| `lieb-liniger-thermal` | Bounded independent CLI temperature solves; owning mixed-precision library grids | MPFR scalars, grand-canonical and fixed-density equilibrium | Mesh, rapidity cutoff and density inversion budgets remain independent of arithmetic precision |
| `hubbard-dispersion` | Shared Async preparation and one branch-spanning query window ([#14](https://github.com/Uni20-dev/bethe/pull/14)) | Yes, MPFR scalars; doped backgrounds use CPU/provider LU | Explicit quadrature, stable gamma-ratio phase and shared immutable kernel |
| `hubbard-continuum` | Bounded independent extrema queries ([#12](https://github.com/Uni20-dev/bethe/pull/12)) | Yes, no MPLAPACK required | Explicit constituent and extrema precision; independent work budgets |
| `xxz-pbc`, `xxz-obc` | Continuous real-root scans; owning parameter sweeps ([#15](https://github.com/Uni20-dev/bethe/pull/15)) | Yes, MPFR scalars; CPU LU available | Ground/excitation domains are unchanged by precision or concurrency |
| `q-boson-pbc` | Continuous real-root scans; owning mixed-precision queries ([#25](https://github.com/Uni20-dev/bethe/pull/25)) | Yes, MPFR scalars; CPU LU available | Ground and canonical states retain full zero-update/failure diagnostics |
| `su3-pbc` | Individual nested solves; owning mixed-precision sweeps ([#13](https://github.com/Uni20-dev/bethe/pull/13)) | Yes, MPFR scalars; CPU LU available | Balanced singlets retain their existing physical scope |
| `tb-pbc` | Single-request CLI stays serial; owning mixed-precision Async size sweeps compare complete states ([#37](https://github.com/Uni20-dev/bethe/pull/37)) | MPFR scalars, including finite string deviations and owning complex rapidities | Even-ring singlet ground-state branch; the residual tolerance is not an energy-error bound |
| `hubbard-pbc`, `hubbard-obc` | Individual finite-size solves; owning parameter sweeps ([#16](https://github.com/Uni20-dev/bethe/pull/16)) | Yes, MPFR scalars; CPU LU available | Supported sectors, both interaction signs and symmetry mappings |
| `gaudin-yang-pbc` | Individual nested solves; owning mixed-precision sector sweeps ([#18](https://github.com/Uni20-dev/bethe/pull/18)) | Yes, MPFR scalars; CPU LU available | Restricted repulsive shells; failed continuation can stop below requested coupling |
| `sun-fermions-pbc` | Single-request CLI stays serial; owning mixed-precision population/coupling grids preserve complete states ([#38](https://github.com/Uni20-dev/bethe/pull/38)) | Yes, MPFR scalars; CPU LU available | Existing odd-population interacting branch, free limits and partial accepted stages; independent two-body and boson-reduction checks |
| `tj-pbc` | Single-request CLI stays serial; owning one/two-worker mixed-precision grids preserve complete states | Native and explicit MPFR; CPU LU available | Existing doped parity restrictions remain; higher precision does not resolve the known large dilute-sector seed/Newton limitation |
| `ladder-pbc` | Owning prepared branch catalogue with dependent Async field/rung queries; complete mixed-precision results checked against direct scans | Native and explicit MPFR; CPU LU available | Fixed-Sz/singlet scans retain their direct paths, physical branch scope and whole-scan budgets; one precision per preparation |
| `bose-fermi-pbc` | Individual nested solves; owning mixed-precision requests ([#19](https://github.com/Uni20-dev/bethe/pull/19)) | Yes, MPFR scalars; CPU LU available | Restricted mixed shells; zero-update diagnostics remain explicit |
| `su3-dispersion` | Serial analytic CLI; owning mixed-precision Async queries ([#20](https://github.com/Uni20-dev/bethe/pull/20)) | Yes, MPFR scalars | Individual formulas remain too cheap to justify per-point scheduling |
| `haldane-shastry-pbc`, `sutherland-pbc` | Serial analytic CLIs; owning mixed-precision window/parameter scans ([#21](https://github.com/Uni20-dev/bethe/pull/21)) | Yes, MPFR scalars | Gap/order labels remain exact integers; cheap formulas stay serial |
| `tb-dispersion` | Serial analytic CLI; owning mixed-precision Async queries ([#22](https://github.com/Uni20-dev/bethe/pull/22)) | Yes, MPFR scalars | Analytic spinon kinematics; no operator-dependent spectral weights |
| `sine-gordon-dispersion` | Serial analytic CLI; owning mixed-precision Async queries ([#23](https://github.com/Uni20-dev/bethe/pull/23)) | Yes, MPFR scalars | Particle kinematics are distinct from scattering and finite-volume numerical convergence |
| `xyz-dispersion` | Serial CLI; owning bound branches and shared Async preparation / parameter curves ([#24](https://github.com/Uni20-dev/bethe/pull/24)) | Explicit MPFR scalars | Shared Uni20 complex theta jets; owning mixed-precision bound readers; no MPLAPACK required |
| `xxz-dispersion` | Immutable band with cheap queries; bulk preparation remains serial | Explicit scalar MPFR | Scaled theta constants, finite quadrature and analytic limits; no MPLAPACK required |
| `xxz-qg-obc` | Continuous regular-root scans; owning regular/free-sector queries ([#26](https://github.com/Uni20-dev/bethe/pull/26)) | Yes, MPFR scalars | Critical regular-root/free-endpoint scope; cheap free blocks stay serial |
| `sine-gordon-bethe-yang` | Bounded independent level scans ([#27](https://github.com/Uni20-dev/bethe/pull/27)) | Explicit MPFR scalars | No MPLAPACK; Fourier, iteration and omitted-wrapping limits remain separate |
| `biquadratic-obc` | Continuous real-root and owned bound-cluster scans; independent AF sector minima ([#35](https://github.com/Uni20-dev/bethe/pull/35)) | Yes, all existing modes with MPFR scalars; CPU LU available | Restricted candidate families, ranking and explicit failure diagnostics; selected states and Q-system/singlet paths stay direct |
| `lee-yang-vacuum`, `lee-yang-excited` | Owning mixed-precision parameter sweeps; optional vacuum/excited overlap | Yes, MPFR scalars | Nonlinear, mesh, cutoff and work convergence are separate from arithmetic bits |
| `sine-gordon-vacuum` | Serial single-scale CLI; owning mixed-precision Async parameter sweeps | Explicit MPFR/MPC; no MPLAPACK required | Independent kernel, nonlinear, mesh, cutoff and contour controls; more bits can require larger work budgets |
| `sine-gordon-excited` | Owning mixed-precision parameter sweeps; optional independent vacuum/excited overlap | Explicit MPFR/MPC; no MPLAPACK required | Selected p>=1, I=1/2 or 3/2 family; independent Fourier, hole, nonlinear and refinement budgets; coupled iterations remain sequential |
| `kondo-response` | Bounded field scans and owning Async queries ([#29](https://github.com/Uni20-dev/bethe/pull/29)); single-point/default one-worker requests stay direct | Yes, MPFR scalars; no MPLAPACK required | Precision-aware log-Gamma amplitude, Euler-series tails and endpoint-aware quadrature; independent work budgets and explicit failure diagnostics |
| `potts-pbc` | Bounded charged-momentum scans and owning mixed-precision levels after one vacuum solve | Yes, MPFR scalars; no MPLAPACK required | Explicit constants, twists, auxiliary roots and complete failure diagnostics; vacuum and charged one-hole family only |
| `asep-pbc` | Serial CLI; owning mixed-precision shared TASEP seeds and dependent Async rate grids | Explicit MPFR/MPC; CPU/provider real LU | Existing continuation and work budgets; tiny-bias drift remains explicit |
| `tasep-pbc` | Serial CLI; owning mixed-precision size/filling requests | Explicit MPFR/MPC; CPU/provider real LU | Existing gap branch and work budgets; no full-spectrum claim |
| `xyz-pbc` | Serial CLI; owning mixed-precision size/eta/t sweeps | Explicit MPFR; CPU/provider LU | Shared complex theta jets, precision-aware roots, Newton controls and failure results |

`bethe-info` reports capabilities and has no numerical parallelism to add.
The inventory covers 44 numerical frontends; it does not extend any model's
physical state family. Generic scalar functions and numerical primitives belong
in Uni20; model equations and physical normalization remain in Bethe.

## Validation records

[Historical precision and execution checks](validation-records.md) retain the
feature-checkpoint comparisons, provider configurations and test counts.
See [testing and CI](testing.md) for the current validation workflow and
[benchmarks](benchmarks.md) for measured workloads.

<a id="betheyang-precision"></a>
[Bethe–Yang precision checkpoint](validation-records.md#betheyang-precision).
