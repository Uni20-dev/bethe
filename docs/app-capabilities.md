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
Names omit the `bethe-` executable prefix. This inventory describes the current
source tree, including the integrated checkpoints #9–#38. Rebuild an older
installation to obtain these capabilities. Linked PRs preserve the original
precision/provider configurations and independent numerical or serial/concurrent
validation; their measurements are historical records, not new integration runs.

| Apps | Existing execution path | MPFR frontend support | Next check / dependency |
| --- | --- | --- | --- |
| `xxx-pbc`, `xxx-obc` | Continuous real-root scans; owning Async scan API | Yes, MPFR scalars | Preserve owned solver options through final consumers |
| `xxx-structure-factor` | Async preparation, fused root/form-factor runners, ordered reduction | Yes, MPLAPACK MPFR required | Preserve mixed-precision and failed-reference coverage |
| `lieb-liniger-pbc`, `lieb-liniger-obc` | Shared continuous real-root scanner; owning Async scans ([#17](https://github.com/Uni20-dev/bethe/pull/17)) | Yes, MPLAPACK MPFR required | Preserve downstream-consumer and concurrent-precision checks |
| `richardson`, `central-spin` | Individual solver calls; owning grid-sweep example and full-state composition checks ([#17](https://github.com/Uni20-dev/bethe/pull/17)) | Yes, MPFR scalars | Keep continuation within each solve sequential |
| `xxx-structure-factor-thermo` | Bounded ordered grid, direct table collection | Native only | Precision-aware transition-rate integrals |
| `lieb-liniger-dispersion` | Prepared solver and bounded grids; shared Async backgrounds | Yes, MPLAPACK MPFR required ([#9](https://github.com/Uni20-dev/bethe/pull/9)) | Preserve explicit-precision thermodynamic validation |
| `lieb-liniger-thermal` | Bounded independent CLI temperature solves; owning mixed-precision library grids | MPFR + MPLAPACK, grand-canonical and fixed-density equilibrium | Mesh, rapidity cutoff and density inversion budgets remain independent of arithmetic precision |
| `hubbard-dispersion` | Shared Async preparation and one branch-spanning query window ([#14](https://github.com/Uni20-dev/bethe/pull/14)) | Native only | Precision-aware kernels and stable gamma-ratio phases |
| `hubbard-continuum` | Bounded independent extrema queries ([#12](https://github.com/Uni20-dev/bethe/pull/12)) | Native only | Thermodynamic kernel and gamma-ratio precision |
| `xxz-pbc`, `xxz-obc` | Continuous real-root scans; owning parameter sweeps ([#15](https://github.com/Uni20-dev/bethe/pull/15)) | Yes, MPLAPACK MPFR required | Preserve existing ground/excitation restrictions and low-precision failure diagnostics |
| `q-boson-pbc` | Continuous real-root scans; owning mixed-precision queries ([#25](https://github.com/Uni20-dev/bethe/pull/25)) | Yes, MPLAPACK MPFR required | Preserve ground/canonical-state precision and scalar-only zero-update diagnostics |
| `su3-pbc` | Individual nested solves; owning mixed-precision sweeps ([#13](https://github.com/Uni20-dev/bethe/pull/13)) | Yes, MPLAPACK MPFR required | Balanced singlets retain their existing physical scope |
| `tb-pbc` | Single-request CLI stays serial; owning mixed-precision Async size sweeps compare complete states ([#37](https://github.com/Uni20-dev/bethe/pull/37)) | MPFR + MPLAPACK, including finite string deviations and owning complex rapidities | Preserve the even-ring singlet ground-state branch; the residual tolerance is not an energy-error bound |
| `hubbard-pbc`, `hubbard-obc` | Individual finite-size solves; owning parameter sweeps ([#16](https://github.com/Uni20-dev/bethe/pull/16)) | Yes, MPLAPACK MPFR required | Preserve existing sectors, both interaction signs and symmetry mappings |
| `gaudin-yang-pbc` | Individual nested solves; owning mixed-precision sector sweeps ([#18](https://github.com/Uni20-dev/bethe/pull/18)) | Yes, MPLAPACK MPFR required | Preserve repulsive-shell restrictions and requested-coupling failure diagnostics |
| `sun-fermions-pbc` | Single-request CLI stays serial; owning mixed-precision population/coupling grids preserve complete states ([#38](https://github.com/Uni20-dev/bethe/pull/38)) | Yes, MPLAPACK MPFR required | Existing odd-population interacting branch, free limits and partial accepted stages; independent two-body and boson-reduction checks |
| `tj-pbc` | Single-request CLI stays serial; owning one/two-worker mixed-precision grids preserve complete states | Native and explicit MPFR with MPLAPACK; library free/XXX/zero-update paths also work without that provider | Existing doped parity restrictions remain; higher precision does not resolve the known large dilute-sector seed/Newton limitation |
| `ladder-pbc` | Owning prepared branch catalogue with dependent Async field/rung queries; complete mixed-precision results checked against direct scans | Native and explicit MPFR with MPLAPACK; library product/zero-branch paths work without the provider | Fixed-Sz/singlet scans retain their direct paths, physical branch scope and whole-scan budgets; one precision per preparation |
| `bose-fermi-pbc` | Individual nested solves; owning mixed-precision requests ([#19](https://github.com/Uni20-dev/bethe/pull/19)) | Yes, MPLAPACK MPFR required | Preserve mixed-shell restrictions and scalar-only zero-update diagnostics |
| `su3-dispersion` | Serial analytic CLI; owning mixed-precision Async queries ([#20](https://github.com/Uni20-dev/bethe/pull/20)) | Yes, MPFR scalars | Individual formulas remain too cheap to justify per-point scheduling |
| `haldane-shastry-pbc`, `sutherland-pbc` | Serial analytic CLIs; owning mixed-precision window/parameter scans ([#21](https://github.com/Uni20-dev/bethe/pull/21)) | Yes, MPFR scalars | Preserve exact-integer gap/order handling; cheap formulas stay serial |
| `tb-dispersion` | Serial analytic CLI; owning mixed-precision Async queries ([#22](https://github.com/Uni20-dev/bethe/pull/22)) | Yes, MPFR scalars | XXZ/XYZ elliptic models still need separate precision audits |
| `sine-gordon-dispersion` | Serial analytic CLI; owning mixed-precision Async queries ([#23](https://github.com/Uni20-dev/bethe/pull/23)) | Yes, MPFR scalars | Scattering phases and finite-volume equations remain separate precision work |
| `xyz-dispersion` | Serial CLI; owning bound branches and shared Async preparation / parameter curves ([#24](https://github.com/Uni20-dev/bethe/pull/24)) | Native only | Explicit-precision elliptic/theta functions |
| `xxz-dispersion` | Immutable band with cheap queries; bulk preparation remains serial | Native only | Model adoption of theta constants and finite quadrature; see upstream prerequisites below |
| `xxz-qg-obc` | Continuous regular-root scans; owning regular/free-sector queries ([#26](https://github.com/Uni20-dev/bethe/pull/26)) | Yes, MPFR scalars | Preserve critical regular-root/free-endpoint scope; cheap free blocks stay serial |
| `sine-gordon-bethe-yang` | Bounded independent level scans ([#27](https://github.com/Uni20-dev/bethe/pull/27)) | Native only | MPFR phase integration and a model-level precision audit |
| `biquadratic-obc` | Continuous real-root and owned bound-cluster scans; independent AF sector minima ([#35](https://github.com/Uni20-dev/bethe/pull/35)) | Native only | Model-level explicit precision; selected states and dependent Q-system/singlet paths stay direct |
| `lee-yang-vacuum`, `lee-yang-excited` | Owning mixed-precision parameter sweeps; optional vacuum/excited overlap | Yes, MPFR scalars | Keep nonlinear, mesh, cutoff and work convergence separate from arithmetic bits |
| `sine-gordon-vacuum`, `sine-gordon-excited` | Owning parameter sweeps; optional vacuum/excited overlap ([#28](https://github.com/Uni20-dev/bethe/pull/28)) | Native only | Explicit precision through controls/constants/reductions and complex kernels/quadrature; coupled iterations remain sequential |
| `kondo-response` | Bounded field scans and owning Async queries ([#29](https://github.com/Uni20-dev/bethe/pull/29)); single-point/default one-worker requests stay direct | Native only | Explicit-precision scaled positive-real Gamma, Euler-series tails and quadrature |
| `potts-pbc` | Bounded charged-momentum scans and owning mixed-precision levels after one vacuum solve | Yes, MPFR scalars; no MPLAPACK required | Explicit constants, twists, auxiliary roots and complete failure diagnostics; vacuum and charged one-hole family only |
| `asep-pbc` | Serial CLI; owning shared TASEP seed and dependent Async rate grids ([#31](https://github.com/Uni20-dev/bethe/pull/31)) | Native only | Explicit complex arithmetic, seed/continuation controls and diagnostics |
| `tasep-pbc`, `xyz-pbc` | Serial CLIs; owning size/filling/eta/t sweeps with a compiled example ([#32](https://github.com/Uni20-dev/bethe/pull/32)) | Native only | Explicit complex arithmetic, controls and diagnostics; finite XYZ also needs complex theta arguments/derivatives |

`bethe-info` reports capabilities and has no numerical parallelism to add.
The inventory covers 43 numerical frontends; it does not extend any model's
physical state family. Generic scalar functions and numerical primitives belong
in Uni20; model equations and physical normalization remain in Bethe.

## Validation records

The records below summarize validation at each original feature checkpoint.
The stacked integration additionally runs the routine native, wide-precision,
sanitizer and documentation CI suites on the combined trees.

- [Async execution](async-execution.md) describes scan, preparation and grid
  lifetime/failure contracts, including completion of all admitted work.
- [Parallel grids](parallel-grids.md) records independent scalar comparisons,
  precision coverage and benchmark conditions.
- [MPFR builds](building.md#enable-arbitrary-precision) describes provider
  requirements. Native, wide and sanitizer CI exercise different configurations.
- Thermodynamic XXX CLI tests compare complete serial, threaded and non-retained
  rows, including grids larger than the admission window, mixed successful and
  failed points, threshold counts, native precisions and output preservation.

The MPFR checkpoints #9, #13, #15, #16, #18, #19, #20, #21, #22, #23, #25 and #26 verify precision growth
against independent references beyond binary128, typed endpoint/failure values,
CLI export/provider contracts and ownership in composed Async calculations.
The scheduling checkpoints #11, #12, #14 and #27 compare complete results
against serial/baseline calculations and document measured wall/CPU and memory
behavior. Finite Hubbard's scalar/LU work is independent of the complex
Gamma-ratio functions still required by thermodynamic Hubbard.

The composition audit in #17 adds finite Lieb–Liniger owning scan APIs
and a compiled Richardson/central-spin sweep. Deferred input destruction,
downstream consumers, one/multiple workers, failed solves, cancellation and
independent MPFR precisions are covered. Existing XXX scan and structure-factor
composition tests were rerun; their APIs need no new wrapper. No new precision
support or parallel speedup is claimed by that audit.

Gaudin–Yang's #18 verifies the unreduced complex Bethe equations and
two-body accuracy growth at 128/256/333 bits. Free states and zero-update
diagnostics also work with MPFR scalars alone at the library level; Newton
updates and the MPFR frontend require MPLAPACK. Mixed interacting populations
remain odd, and continuation within each independent solve remains sequential.

Bose–Fermi's #19 validates mixed-species equations, Lieb–Liniger
reductions, free states and missing-energy failures. The shared continuum
Newton guard keeps zero-update diagnostics usable in scalar-only MPFR builds;
actual updates and the frontend require MPLAPACK. Tests resolve a weak-coupling
shift above a nonzero Fermi sea beyond binary128 resolution and compare complete
states under owning, mixed-precision grid composition. Mixed interacting shells
still require odd fermion number.

Haldane–Shastry and Sutherland's #21 inherit #20's scalar CLI dispatcher.
The spin model keeps exact motif ordering; the gas keeps exact MPFR integer gap
coefficients until physical scaling. Full MPFR tests run without MPLAPACK,
including differential-Hamiltonian references, large-exponent small gaps,
precision growth, budget refusals and independent owning Async spectra. Native
ED/wavefunction checks and 561 unchanged exported rows validate the existing
formats. Individual analytic CLI calculations remain serial.

TB dispersion's #22 uses the same scalar dispatcher from #20 and
extends the common sine/massive-band helper rather than duplicating kinematics.
Independent dressed-sea and continuum references, scalar-only MPFR physics,
positive momenta at the MPFR exponent floor and owning mixed-precision queries
pass. Existing TB/XXZ/XYZ checks and 1,710 identical native output rows validate
shared-band consumers. This does not enable XXZ/XYZ's elliptic models in MPFR.

Sine-Gordon particle dispersion's #23 also uses #20's scalar dispatcher.
Independent algebraic breather masses, precision-sensitive stable-species
selection, unequal-mass threshold minimization and full owning Async queries
pass at 128/256/333 bits. Full MPFR physics runs in scalar-only sanitizers;
1,095 native/binary128 output rows match the baseline. The scattering phase,
Bethe–Yang and finite-volume integral equations remain separate precision work.

XYZ's #24 prepares the bound coefficient once per branch and reuses it
for both parity copies. Independent elliptic/limiting references, owning Async
readers and failed preparations/queries pass in native/wide/sanitizer builds.
End-to-end binary128 curve time falls from 5.42 to 0.27 seconds with identical
output; an independent 2,097,664-energy parameter sweep falls from 4.480 to
0.585 seconds with eight workers. These are distinct reuse and scheduling
measurements; the CLI stays serial. Commands and limits are in the PR's
`docs/parallel-grids.md`. It does not extend the theta helper to MPFR.

The theta-constant prerequisite is proposed separately in
[Uni20 PR #103](https://github.com/Uni20-dev/uni20/pull/103). It provides scaled
zero-argument rectangular constants at native and explicit MPFR precision,
with independent product/Gamma references and accuracy growth through 1,024
bits. Local native/wide/sanitizer checks and upstream CI pass; independent
review and merge remain pending. Bethe has not changed its Uni20 pin or enabled MPFR for XXZ/XYZ.
XXZ bulk-energy integration and XYZ complex theta arguments/derivatives remain
separate work.

The shared finite-integration prerequisite is proposed in
[Uni20 PR #104](https://github.com/Uni20-dev/uni20/pull/104). Array-valued tanh-sinh
quadrature carries explicit precision through endpoint distances, accumulation
and recoverable error diagnostics. Local native/wide/sanitizer checks pass,
including reflected/translated singularities, cutoff sensitivity and accuracy
growth through 1,024 bits. Upstream CI passes; independent review and merge remain pending.
Its estimates are not certified bounds. Bethe still needs model-level adoption,
kernel/tolerance audits and native-output comparisons for sine-Gordon, XXZ bulk,
Hubbard and Kondo; no new app capability is implied by this upstream proposal.

Update this inventory with each integrated app checkpoint. A proposed capability
or a pending audit is not evidence of implemented support.

Kondo's [#29](https://github.com/Uni20-dev/bethe/pull/29) adds inclusive
signed field scans and streams complete rows through the existing bounded grid.
Every query owns its series/quadrature scratch; owning deferred/threaded
composition and failed responses retain all diagnostics. Native/wide/sanitizer
checks pass 15/21/15; 30 single-point cases preserve merged-main outputs.
Substantial high-field curves scale 7.6–7.7x at eight workers, including a
65-point binary128 curve taking 74.02 s direct versus 9.64 s parallel. Cheap
low-field curves saturate earlier; timings and CPU costs are in that PR's
`docs/parallel-grids.md`. Kondo's Gamma amplitude is positive-real, not complex.
MPFR requires precision-aware asymptotic order/shift, series-tail and quadrature
error control, as well as explicit constants and reductions; scheduling does
not remove that numerical work.

Q-boson's #25 retains the existing continuous scanner and compares
complete states, gaps and failed references in serial/threaded execution.
Independent occupation-Hamiltonian energies and the unreduced Bethe equations
validate MPFR, including a weak-deformation correction beyond binary128 and
positive energies at eta=1e-7000. The scalar-only library supports analytic
limits and zero-update diagnostics through #19's shared guard; actual Newton
updates and the CLI require MPLAPACK. No new scheduling mechanism is needed.

Quantum-group XXZ's #26 validates both the regular critical real-root
family and complete free-endpoint sectors with MPFR scalars alone. Independent
spin-Hamiltonian eigenvalues and complex reflection equations verify precision
growth at 128/256/333 bits. All 45 states of a 16-site family, including partial
and failed scans, agree under serial/threaded execution. Owning regular/free
queries also pass deferred and threaded composition tests. Native output
comparisons retain 681 identical rows; no massive or new complex-root family
is introduced.

Sine-Gordon Bethe–Yang's #27 uses the existing bounded grid for
independent interacting levels. One worker, a single level and the free point
stay direct and serial. Native/wide/sanitizer CLI checks preserve all fields,
partially converged scans and ordered output. Matched baseline comparisons
retain 495 selected rows and 2,066 benchmark rows. Eight-worker wall times
improve from 6.33 to 0.85 seconds for 2,049 binary64 levels and from 30.82 to
5.71 seconds for 17 binary128 levels. CPU/memory costs, a roughly 4% one-worker
binary128 overhead and measurement limits are recorded in the PR's
`docs/parallel-grids.md`. Numerical kernels and precision support are unchanged.

Potts's #30 retains one vacuum solve and reuses each charged-momentum
result for both conjugate charge rows. The existing bounded grid streams
ordered rows, including missing-reference and iteration failures; default
one-worker and selected single-state requests stay serial. Complete owning
states match under deferred/one/four-worker execution. Native/wide/sanitizer
checks, 33 matched-baseline cases and tutorial exports pass. All 8,100 timed
rows match baseline output; eight-worker whole-curve speedups are 7.4x in
binary64/binary80 and 6.7x in binary128. CPU/memory costs and shared-workstation
limits are documented in the PR's `docs/parallel-grids.md`. The underlying
real-root XXZ equations and precision support are unchanged.

ASEP's #31 extracts the existing rate solver around an owning
`PreparedGap`. One TASEP seed feeds independent rate queries through the shared
Async preparation/grid APIs. Full fields, roots, failure budgets and independent
analytic endpoints agree in direct, deferred and threaded checks. Native/wide/
sanitizer checks pass 29/39/29; 60 matched-main cases preserve frontend outputs
and physical/control metadata. All 38,070 complete benchmark state rows match.
Seed reuse removes 24–32% of serial work; prepared eight-worker curves achieve
8.9–11.0x combined speedup over independent direct solves. The example and
`docs/parallel-grids.md` in that PR record separate reuse/concurrency costs,
output retention and shared seed-counter semantics. The physical branch,
single-request analytic paths and precision support are unchanged.

Checkpoint #32 validates owning TASEP and finite XYZ requests, including successful
and failed states, deferred/one/multiple workers, and independent error handling
when both models share a scheduler. Substantial filling/eta curves scale
7.3–7.6x with eight workers across native fixed precisions; all complete output
states match direct calls. Scalar kernels and single-request CLIs are unchanged.
Their MPFR paths remain unimplemented; the complex-function and Newton audit
is separate from this measured concurrency checkpoint.

Checkpoint #33 validates owning TB, SU(N) fermion and t–J requests with a compiled
example and complete-state comparisons under deferred and one/four-worker
scheduling. Native/wide/sanitizer checks pass 57/79/57. All 90,888 timed state
rows match direct calls; independent-request sweeps scale 6.4–7.5x at eight
workers, with documented CPU/memory costs. The binary64 t–J benchmark includes
two numerical failures in every variant. A larger 512-site curve has 12 failures;
all states and diagnostics match direct/eight-worker execution, and a bounded
higher-precision probe does not resolve the representative dilute failure.
This convergence limitation is separate from Async correctness. Model kernels,
frontends and precision support are unchanged; no new state family is claimed.

Checkpoint #34 adds `ladder::PreparedSectors<Real>`: one full-sector scan owns its
ordered solved/failing branch prefix and later field/rung queries replay the
same selection rules without new Newton solves. Whole-scan budgets, partial
candidates, analytic product bypasses and exact tie/large-shift handling are
preserved. Copies share immutable storage; dependent Async grids outlive the
preparation handle. Native/wide/sanitizer checks pass 51/67/52, including the
independent ladder oracles and shared-grid held-query lifetime checks.
Matched-main comparisons retain 60 requests / 1,491 identical rows; all 86,088
benchmark state rows converge and match direct calls. Reuse is the primary gain
(roughly 90–270x on these many-parameter curves); subsequent query parallelism
saves only hundredths of a second. Fixed-sector scans and native precision
support remain unchanged. Catalogue/output retention costs are documented.

Checkpoint #35 extends the biquadratic frontend's `--threads` to six ferromagnetic
bound-cluster families and AF module minima. A shared CLI adapter owns its
sequential label generator and per-candidate requests, uses existing continuous
runners, and waits on the final Async result before ranking or output. Numerical
families, budgets, ordering and first-failure rows are preserved; all candidates
are still retained before ranking. Native/wide/sanitizer checks pass 18/18/18,
including held-peer failure completion and complete CLI comparisons. Fresh-main
comparison covers 135 requests / 3,192 rows; standalone long-double bound-pair
fields have documented rounding-scale differences, while worker-count outputs
match exactly. All 45 substantial benchmark runs converge and their 6,309,585
exported rows match main and serial output. Whole-frontend eight-worker scaling
is 5.9–7.4x, with unchanged serial cost. Native precision support is unchanged;
no MPFR claim is made. Single-state and dependent Q-system/singlet modes remain
serial; cheap scans do not justify changing the one-worker default.
