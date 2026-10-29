# Precision and execution validation records

Historical feature checkpoints from October 2026. These preserve the original
configurations, test counts and comparisons; they are not a statement of the
current suite size. See [current capabilities](app-capabilities.md) and
[testing and CI](testing.md) for maintained instructions.

## Validation records

The records below summarize validation at each original feature checkpoint.
At those checkpoints, integration additionally ran the routine native, wide-precision,
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
Gamma-ratio functions used by thermodynamic Hubbard.

The composition audit in #17 adds finite Lieb–Liniger owning scan APIs
and a compiled Richardson/central-spin sweep. Deferred input destruction,
downstream consumers, one/multiple workers, failed solves, cancellation and
independent MPFR precisions are covered. Existing XXX scan and structure-factor
composition tests were rerun; their APIs need no new wrapper. No new precision
support or parallel speedup is claimed by that audit.

Gaudin–Yang's #18 verifies the unreduced complex Bethe equations and
two-body accuracy growth at 128/256/333 bits. Free states, zero-update diagnostics and Newton
updates work with MPFR scalars, using Uni20 CPU LU when MPLAPACK is unavailable. Mixed interacting populations
remain odd, and continuation within each independent solve remains sequential.

Bose–Fermi's #19 validates mixed-species equations, Lieb–Liniger
reductions, free states and missing-energy failures. The shared continuum
Newton driver uses normal Uni20 dispatch, including MPFR CPU LU. Tests resolve a weak-coupling
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
shared-band consumers. That checkpoint did not enable XXZ/XYZ's elliptic models
in MPFR; the [current inventory](app-capabilities.md) records subsequent model-level support.

Sine-Gordon particle dispersion's #23 also uses #20's scalar dispatcher.
Independent algebraic breather masses, precision-sensitive stable-species
selection, unequal-mass threshold minimization and full owning Async queries
pass at 128/256/333 bits. Full MPFR physics runs in scalar-only sanitizers;
1,095 native/binary128 output rows match the baseline. The scattering phase,
Bethe–Yang and finite-volume integral equations required separate precision
checkpoints; their current coverage is listed in the [inventory](app-capabilities.md).

XYZ's #24 prepares the bound coefficient once per branch and reuses it
for both parity copies. Independent elliptic/limiting references, owning Async
readers and failed preparations/queries pass in native/wide/sanitizer builds.
End-to-end binary128 curve time falls from 5.42 to 0.27 seconds with identical
output; an independent 2,097,664-energy parameter sweep falls from 4.480 to
0.585 seconds with eight workers. These are distinct reuse and scheduling
measurements; the CLI stays serial. Commands and limits are in the PR's
`docs/parallel-grids.md`. That scheduling checkpoint did not extend the theta
helper to MPFR.

Merged [Uni20 PR #103](https://github.com/Uni20-dev/uni20/pull/103) supplies
scaled zero-argument rectangular theta constants at native and explicit MPFR
precision. Independent product/Gamma references demonstrate accuracy growth
through 1,024 bits; native/wide/sanitizer checks and upstream CI passed at that
checkpoint. Bethe's [XXZ dispersion](xxz-dispersion.md) uses these constants.

XYZ excitation bands now use the shared complex theta values and derivatives in
[Uni20 #115](https://github.com/Uni20-dev/uni20/pull/115), replacing Bethe's duplicate
evaluator. Explicit MPFR reaches couplings, spinon/bound bands, continuum edges
and owning Async readers. The finite XYZ nonlinear solve also carries explicit
precision through its controls, complex roots and Newton arithmetic, with
independent spin-Hamiltonian and XY references beyond binary128.

Merged [Uni20 PR #104](https://github.com/Uni20-dev/uni20/pull/104) supplies
array-valued tanh-sinh quadrature with explicit precision in endpoint distances,
accumulation and recoverable error diagnostics. Its validation includes
reflected/translated singularities, cutoff sensitivity and accuracy growth
through 1,024 bits. The error estimates are not certified bounds. Each Bethe
consumer still needs model-specific kernel, tolerance and native-output checks;
the [current inventory](app-capabilities.md) and individual model pages record implemented support.

Update the [current inventory](app-capabilities.md) with each integrated app checkpoint. A proposed capability
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
positive energies at eta=1e-7000. Analytic limits, zero-update diagnostics and actual Newton
updates work with MPFR scalars, including Uni20 CPU LU. No new scheduling mechanism is needed.

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

### Bethe–Yang precision

The same-charge sine-Gordon pair solver and its real phase/complex Fourier kernel
carry explicit MPFR precision through constants, quadrature, tail bounds and
failed-state diagnostics. These paths use scalar arithmetic without MPLAPACK.
The existing bounded level grid also accepts owning requests at different
precisions. The vacuum NLIE also supports explicit MPFR/MPC, including independently checked
Fourier tables, free-Dirac energies and interacting D3 reference comparisons.
Both vacuum and selected excited NLIE frontends retain explicit precision in
all controls and success/failure diagnostics. Source-table and free-Dirac level
checks reach 333 bits. Arithmetic precision and Fourier, hole, nonlinear, mesh
and cutoff budgets remain independent; more bits can require more work.
