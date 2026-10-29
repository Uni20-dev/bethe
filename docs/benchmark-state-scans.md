# State and parameter-scan measurements

[Benchmark index and methodology](benchmarks.md) · [Current capabilities](app-capabilities.md)

These are dated October 2026 measurements. Hardware, precision, providers and
timing scope belong to each record; later capabilities do not change what was
measured. Reproduction commands describe those workloads, not recommended defaults.

## Kondo field-response scans

`bethe-kondo-response --field START --field-end END --points N --scale T_B`
evaluates a signed, inclusive field curve. `--threads` overlaps complete,
independent impurity-response solves; each Euler series and Fourier quadrature
remains sequential. The shared grid admits at most 128 active or uncollected
points and emits rows on the calling thread in field order. One worker,
single-point requests and all-zero curves use direct serial evaluation.

The following are whole-process medians of three runs on 28 October 2026,
on an Intel Core i7-12700KF with GCC 13.3 Release (`-O3`, no fast-math),
Uni20 `04163b9b` and oneTBB 2022.3. The host's `long double` is binary80.
Both BLAS/OpenMP limits were one. Runs used default numerical tolerances,
CSV output to a local file and `--no-retain --no-preamble`. Timing includes
startup, scheduling, formatting and file output; memory is peak process RSS
over all three repeats. Worker order was reversed on the second repeat,
and builds did not overlap the timings. This was a shared desktop.

| Field interval, scale 1 | Arithmetic | Points | Direct serial (s) | 4 workers (s) | 8 workers (s) | CPU serial / 8 workers (s) | Peak RSS (MiB) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| [-1,1] | binary64 | 32769 | 0.30 | 0.11 | 0.12 | 0.30 / 0.95 | 5.30 |
| [1,10] | binary64 | 4097 | 6.70 | 1.70 | 0.87 | 6.69 / 6.99 | 5.16 |
| [1,10] | binary80 | 513 | 10.16 | 2.61 | 1.33 | 10.16 / 10.50 | 5.16 |
| [1,10] | binary128 | 65 | 74.02 | 18.93 | 9.64 | 74.01 / 76.50 | 5.32 |

The substantial high-field curves scale about 7.6–7.7x at eight workers.
The cheap low-field curve is dominated by output and scheduling sooner:
four workers are already as fast as eight, and eight use considerably more
CPU time. Serial therefore remains the default. These measurements do not
establish scaling beyond eight workers.

For example, reproduce the high-field binary128 case with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /usr/bin/time -v \
  bethe-kondo-response --field 1 --field-end 10 --points 65 --scale 1 \
  --precision fp128 --threads 8 --no-retain --format csv --no-preamble > response.csv
```

All 336,996 rows from the 36 timed processes agree byte for byte with the
corresponding serial curve, including error estimates and work counters.
Thirty single-point comparisons with merged main `46b636e`, built with the
same compiler and Uni20 revision, preserve every existing numerical field
at one/four requested workers; the new `field` column identifies each row.
Focused native/wide/sanitizer tests additionally cover owning deferred and
threaded queries, mixed numerical failures, both scan directions, a 129-point
window crossing, extreme finite endpoints and non-retained JSON/CSV/TSV.
At that checkpoint, the library's physics, quadrature and precision support
were unchanged; arbitrary-precision Kondo kernels remained separate work.
See the [Kondo guide](kondo.md) for the subsequent explicit-MPFR support.

## Potts momentum scans

The periodic three-state Potts frontend prepares its vacuum reference once,
then evaluates independent charged momenta through the bounded grid:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  bethe-potts-pbc 192 --precision fp64 --threads 8 \
  --no-retain --format csv --no-preamble > potts.csv
```

Each momentum owns a complete auxiliary XXZ solve. Its charge-conjugate rows
reuse that solve, so a full L-site curve performs L charged solves and emits
2L+1 rows including the vacuum. Rows retain increasing momentum order and the
original +1/-1 charge order. Each root iteration remains sequential. The
one-worker default, a selected single momentum and vacuum-only calculations
use direct serial evaluation.

At most 128 charged states are active or awaiting collection; requested
concurrency is also capped by the selected momentum count. Root vectors and
iteration scratch scale with system size, so this is a bound on outstanding
states, not a fixed byte budget. With retention disabled, delivered rows are
discarded. Numerical failures and unavailable reference gaps remain rows;
unexpected exceptions use the shared grid and output completion contracts.
CPU metadata includes the vacuum solve and complete ordered grid with row output.

### Whole-process measurements

Measurements on 28 October 2026 use an Intel Core i7-12700KF, GCC 13.3 Release,
Uni20 `04163b9b`, oneTBB 2022.3 and one BLAS/OpenMP thread. The matched serial
baseline is main `46b636e`; native and wide builds share their dependency and
compiler configuration. Binary80 means this machine's 80-bit extended
`long double`; binary128 is available through the MPLAPACK-enabled build.
These Potts solves use scalar root iterations, not matrix factorizations.

The command above was repeated with 192 binary64 sites, 96 binary80 sites
(`--precision long-double`) and 48 binary128 sites (`--precision fp128`).
Numerical controls were unchanged: at most 10000 sweeps per solve, residual
tolerance 32 epsilon. Whole-process wall/CPU measurements include preparation,
all root solves, formatting and CSV output. Each configuration has three
repeats, reversing execution order on the middle repeat. The baseline omits
`--threads`; the candidate uses 1, 4 or 8.

| Arithmetic | Sites | Baseline serial (s) | New serial (s) | 4 workers (s) | 8 workers (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| binary64 | 192 | 12.76 | 12.74 | 3.30 | 1.72 |
| binary80 | 96 | 11.26 | 11.28 | 2.93 | 1.53 |
| binary128 | 48 | 25.20 | 25.19 | 6.69 | 3.78 |

Eight-worker speedups relative to the new serial path are 7.4x, 7.4x and 6.7x.
CPU work rises modestly; memory includes the scheduler and concurrent root states:

| Arithmetic | Serial CPU (s) | 4-worker CPU (s) | 8-worker CPU (s) | Peak RSS serial / 4 / 8 (MiB) |
| --- | ---: | ---: | ---: | ---: |
| binary64 | 12.74 | 12.99 | 13.15 | 4.50 / 5.00 / 5.00 |
| binary80 | 11.27 | 11.45 | 11.60 | 4.69 / 5.00 / 5.14 |
| binary128 | 25.18 | 25.70 | 26.34 | 4.82 / 5.15 / 5.16 |

All 8,100 rows from the 36 timed processes match the baseline CSV byte for byte
and report convergence. Separate baseline comparisons cover 33 cases across
the three formats, including failed solves, missing vacuum references and
selected charges/momenta. Library tests compare complete returned states,
including roots and labels, under deferred and one/four-worker scheduling.

The serial vacuum solve, ordered collection and unequal momentum costs limit
scaling. These are shared-desktop measurements without CPU affinity, not an
isolated scaling experiment. This scheduling checkpoint preserved the selected
physical state family and used fixed-precision arithmetic. Explicit MPFR support
was added separately; see the [current Potts precision controls](potts.md).
These timing results do not measure MPFR performance.

## TASEP and finite XYZ parameter sweeps

The [complex-root sweep example](async-examples.md#independent-tasep-and-finite-xyz-solves)
compares direct scalar calls with the same calls submitted through the owning
grid. TASEP varies the filling from one particle to half filling, without
particle-hole duplicates; XYZ varies eta across [1/5,4/5] at t=3/4. These are
independent requests. Their coupled root iterations remain sequential, and
the single-request CLIs remain serial.

From a build directory configured with examples enabled:

```sh
cmake --build . --target bethe-async-complex-sweeps-example -j2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./examples/bethe-async-complex-sweeps-example tasep 256 128 0 fp64 > direct.csv
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./examples/bethe-async-complex-sweeps-example tasep 256 128 8 fp64 > async.csv
```

The example's leading comment lines contain timings and worker counts. Compare
the remaining CSV header and **all** rows, including complex roots, quantum
numbers, residuals, statuses and work counters, before interpreting speedup.
The internal timing includes the scheduler and ordered collection, but excludes
printing; external process timing includes both. The grid admits at most 32
active or uncollected states. The example retains collected states for output,
so total memory still grows with the number and size of results.

### Measurement on 28 October 2026

These are whole-process wall medians over three repeats on a shared Intel
Core i7-12700KF workstation, GCC 13.3 Release (`-O3`, no fast-math), oneTBB
2022.3 and Uni20 `04163b9b`. Binary128 uses MPLAPACK 3.0.0 built in Release
mode. BLAS/OpenMP each use one thread; no compilation overlaps the timings.
Worker order is reversed in the middle repeat, with OS-default affinity.
Here `long-double` is x86 binary80 (64 significand bits); these are native
fixed-precision measurements, not MPFR benchmarks.

| Model / precision | Sites | Points | Direct (s) | Async 1 (s) | Async 4 (s) | Async 8 (s) | CPU direct / 8 (s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| TASEP / binary64 | 256 | 128 | 3.05 | 3.07 | 0.81 | 0.42 | 3.05 / 3.23 |
| TASEP / binary80 | 256 | 128 | 12.53 | 12.54 | 3.24 | 1.65 | 12.52 / 12.87 |
| TASEP / binary128 | 128 | 64 | 16.21 | 16.18 | 4.20 | 2.19 | 16.20 / 16.80 |
| XYZ / binary64 | 128 | 1025 | 8.98 | 9.03 | 2.33 | 1.18 | 8.97 / 9.30 |
| XYZ / binary80 | 64 | 513 | 9.71 | 9.72 | 2.50 | 1.30 | 9.70 / 10.14 |
| XYZ / binary128 | 32 | 129 | 8.83 | 8.86 | 2.24 | 1.19 | 8.82 / 9.22 |

All **23,844 complete state rows in 72 processes** match byte for byte, and
every solve converges. Eight-worker speedups range from 7.3 to 7.6; total CPU
work rises by about 3–6%. Peak RSS across these workloads is 14.42 MiB.
One-worker Async overhead is small here; small or analytic-only requests need
not benefit. These shared-workstation measurements establish useful scaling
for these inputs, not a universal performance guarantee. Numerical kernels,
physical state families and precision support are unchanged.

## Nested-model parameter sweeps

The [nested sweep example](async-examples.md#owning-nested-model-requests)
compares direct scalar calls with owning Async grids for TB sizes, SU(N)
repulsive couplings and t–J balanced odd fillings. Each complete state owns its
roots and diagnostics; the coupled scalar iterations remain sequential.

Build `bethe-async-nested-sweeps-example` with `BETHE_BUILD_EXAMPLES=ON` and
run the documented `MODEL SIZE POINTS THREADS PRECISION` commands with zero
threads for the direct control, then one/four/eight threads for Async. Set
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` to prevent nested oversubscription.
Compare all CSV rows after the leading timing comment lines, including nested
roots, labels, residuals, convergence, requested/reached coupling and counters.
Internal timing includes the scheduler and ordered collection, excluding
printing. Whole-process timing includes printing; the retained output means
memory grows with the curve even though admission is bounded at 32 states.

Three-repeat whole-process medians follow, including printing all state fields.
The machine was a shared Intel Core i7-12700KF with GCC 13.3 Release `-O3`,
oneTBB 2022.3 and Uni20 `04163b9b`; no compiler jobs ran during measurement.
Worker order was reversed on the second repeat, with OS-default affinity.
The native configuration disables MPLAPACK and uses normal Uni20 fallback
dispatch; the binary128 configuration enables the Release MPLAPACK 3.0 provider.
No backend is forced by the model or example. All times are seconds; CPU is
user plus system, and peak RSS is the maximum over all repeats/worker counts.

| Model | Arithmetic | SIZE | Points | Direct | Async 1 | Async 4 | Async 8 | CPU direct / 8 | Peak RSS, MiB |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| TB | binary64 | 128 | 449 | 13.61 | 13.63 | 3.56 | 1.94 | 13.60 / 14.57 | 46.37 |
| TB | 80-bit extended | 64 | 193 | 8.73 | 8.72 | 2.31 | 1.26 | 8.73 / 8.88 | 17.83 |
| TB | binary128 | 32 | 81 | 9.33 | 9.32 | 2.49 | 1.34 | 9.33 / 9.56 | 8.43 |
| SU(3) fermions | binary64 | 15 | 4097 | 8.17 | 8.16 | 2.13 | 1.11 | 8.16 / 8.42 | 14.82 |
| SU(3) fermions | 80-bit extended | 9 | 2049 | 8.64 | 8.66 | 2.26 | 1.21 | 8.63 / 9.00 | 10.44 |
| SU(3) fermions | binary128 | 5 | 513 | 9.21 | 9.22 | 2.35 | 1.23 | 9.20 / 9.72 | 6.85 |
| t–J | binary64 | 384 | 96 | 1.59 | 1.58 | 0.41 | 0.24 | 1.59 / 1.73 | 19.24 |
| t–J | 80-bit extended | 256 | 64 | 4.53 | 4.55 | 1.18 | 0.61 | 4.53 / 4.72 | 13.77 |
| t–J | binary128 | 128 | 32 | 5.07 | 5.03 | 1.31 | 0.79 | 5.07 / 5.66 | 7.86 |

The TB curves cover even sizes 128–1024, 64–448 and 32–192 respectively.
The SU(3) curves have three populations of 15, 9 or 5 at total density one,
with c from 1/2 to 4. All TB and SU(3) points converge, as do the two wider
precision t–J curves. These are library sweep benchmarks: the single-request
frontends do not gain a sweep option or parallel Newton iterations.

All 90,888 complete state rows from the 108 timed processes match their direct
counterparts exactly, including unsuccessful iterates. Eight workers give
6.4–7.5x wall-time speedup here, with roughly 2–12% more CPU time. One-worker
overhead is small on these curves. Shorter curves, different population ranges,
providers and shared-machine load can behave differently; the binary64 t–J
curve is substantially shorter than the others.

The 384-site binary64 t–J sweep returns 94 converged states and two failures,
at balanced populations (1,1) and (31,31). Both fail in direct execution too;
the table measures all 96 requests, without dropping failed rows or changing
controls. A larger 512-site trial returned 116 converged states and 12 failures
at balanced populations n=1,5,29,31,...,47. At n=37 and 39 the solver exhausts
its default 10,000-update budget; the others report a stalled line search.
Those failed iterates are not physical energies. The two long runs dominate
that trial's cost, so it is not used for the repeated throughput table.
All 128 full states, including the 12 failures, agree exactly between direct
and eight-worker execution. That single comparison took 135.44 and 66.27 s
respectively (135.42 and 136.19 CPU seconds). Two slow requests limit the gain
despite eight available workers; this is not an eight-worker scaling estimate
from repeated successful curves.

The 512-site two-electron CLI probe also stalls in 80-bit extended precision;
binary128 exhausts a 50-update probe near the same residual, about 0.159.
This leaves a seed/Newton convergence investigation, not evidence that more
workers or arithmetic bits cure the numerical problem. See the
[t–J solver limitations and reproducer](tj.md#library-and-numerical-contract).


## Biquadratic sector and bound-cluster scans

`bethe-biquadratic-obc --threads N` also supports the independent AF
`--sectors` minima and ferromagnetic bound-cluster scans. The latter include
single droplets (`--bound-pairs`, `--bound-triples`, `--bound-quartets`),
`--pair-defects` with any supported real-root count, `--triple-defects`, and
`--two-pair-states`. Real-root excitation scans already use continuous workers.
Selected single states and the dependent Q-system/singlet algorithms remain
on their existing paths; these modes reject an explicit `--threads` option.
One worker remains the default, and scans containing one candidate run directly.

The shared CLI adapter in [`generated-state-scan.hpp`](../apps/generated-state-scan.hpp)
owns its generator and solver until its final Async result completes. A short
serialized step copies the next candidate's labels and assigns its enumeration
ID. Continuous scheduler runners then solve those requests independently,
without a per-candidate coroutine or intermediate result queue. Each runner
writes its own result slot. After completion, the caller performs the existing
ordering, gap subtraction and table output. Numerical failures retain their
original diagnostics; execution exceptions wait for every admitted runner.
The regression suite checks that completion stays pending while a peer solve
is deliberately held open after another throws.

Transient requests are bounded by worker count. Bound-cluster candidates are
still retained until ranking, so result memory scales with candidate and root
counts; `--max-candidates` remains the admission guard. Disabling table retention
does not eliminate this numerical storage. The implementation changes no root
equations, continuation steps, physical family or precision support.

### Validation and measurement

On the Ryzen workstation (Intel Core i7-12700KF), native GCC 13, wide GCC 13 and
Clang 19 ASan/UBSan each pass 18 focused checks. These include the existing
independent small-chain references, complete serial/concurrent output comparisons,
partial scans, export preflight and operation lifetime. A fresh-main comparison
covers 135 requests and 3,192 rows. Labels, ordering, statuses and work counters
are unchanged. Binary64 and binary128 outputs match exactly. In the native
long-double build, 96 standalone bound-pair fields differ at rounding scale
(maximum absolute difference 3.4e-18); all new-build worker counts agree exactly.
A separate GCC 13 probe finds a one-ULP difference between folded and runtime
`acosh(1.5L)`, consistent with this path's sensitivity to compiler inlining.
This cross-build rounding difference is recorded, not hidden by the concurrency
comparison.

Whole-frontend timings on 28 October 2026 use GCC 13.3 Release (`-O3`, no
fast-math), oneTBB 2022.3, Uni20 `04163b9`, and optimized MPLAPACK 3.0.0 for
binary128. BLAS/OpenMP each use one thread. The table gives medians of three
runs, reversing the order on the middle repetition. Compilation had finished;
the workstation was shared, without CPU pinning. Every invocation includes
root tables and JSON output redirected to a host-local file.

| Frontend workload | States | Main serial (s) | New serial (s) | 2 workers (s) | 4 workers (s) | 8 workers (s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 768-site AF sectors, binary64 | 385 | 6.25 | 6.26 | 3.20 | 1.67 | 0.90 |
| 64 sites, pair + 12 real roots, width 15, 80-bit extended | 6,825 | 7.22 | 7.19 | 3.83 | 2.09 | 1.22 |
| 64 sites, pair + 6 real roots, width 11, binary128 | 5,082 | 27.29 | 27.35 | 13.93 | 7.13 | 3.69 |

From one to eight workers, median process CPU times are 6.26 to 6.53 s,
7.18 to 7.56 s, and 27.34 to 28.40 s respectively. Peak RSS across all variants
and repetitions is 18.0–37.7 MiB, 21.6–22.9 MiB and 14.0–16.3 MiB. The sector
scan's larger concurrent Newton matrices explain its higher memory demand.
All 45 measured invocations converge, and all 6,309,585 exported table rows
match their main and serial counterparts exactly after excluding invocation,
revision, timing and concurrency metadata. The standalone bound-pair rounding
differences above do not occur in these three benchmark workloads.

Reproduce with one, two, four or eight workers, enabling binary128 for the last
command. Here native `long double` has a 64-bit significand in the 80-bit
extended format. Run from an installed executable directory or a build directory
with the executables on `PATH`:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
bethe-biquadratic-obc 768 --sectors --precision fp64 \
  --threads 8 --roots --format json > sectors.json
bethe-biquadratic-obc 64 --ferromagnetic --pair-defects all \
  --real-defects 12 --mixed-window 15 --precision long-double \
  --threads 8 --roots --format json > mixed-extended.json
bethe-biquadratic-obc 64 --ferromagnetic --pair-defects all \
  --real-defects 6 --mixed-window 11 --precision fp128 \
  --threads 8 --roots --format json > mixed-binary128.json
```

These substantial scans improve by 5.9–7.4x at eight workers. Small scans
already finish in milliseconds, so this does not justify changing the default.
