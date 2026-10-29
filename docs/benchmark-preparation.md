# Shared-preparation measurements

[Benchmark index and methodology](benchmarks.md) · [Current capabilities](app-capabilities.md)

These are dated October 2026 measurements. Hardware, precision, providers and
timing scope belong to each record; later capabilities do not change what was
measured. Reproduction commands describe those workloads, not recommended defaults.

## XYZ prepared bound branches

`xyz::SpinonDispersion::prepare_bound(s)` evaluates the elliptic coefficient
once and returns an owning immutable branch. Its `energy(Q)` has no theta-function
calls. The CLI uses one preparation for both parity copies and remains serial;
independent parameter curves can compose through the existing Async grid.

Matched whole-process measurements on 28 October 2026 used an i7-12700KF,
GCC 13.3 Release, oneTBB 2022.3 and Uni20 `04163b9b`, with
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`. These are three-repeat medians on
a shared workstation, comparing main `46b636e` with prepared branches.
The numerical output is byte-identical. Timing includes preparation, formatting,
streaming and file output with table retention disabled.

| Precision | Rows | Previous wall / CPU (s) | Prepared wall / CPU (s) | Previous / prepared peak RSS (MiB) |
| --- | ---: | ---: | ---: | ---: |
| binary64 | 786,444 | 2.46 / 2.45 | 1.50 / 1.49 | 4.06 / 4.38 |
| binary128 | 49,164 | 5.42 / 5.41 | 0.27 / 0.27 | 4.22 / 4.53 |

The gain comes from removing redundant special-function work, not parallelism.
Reproduce the binary128 curve from a configured build directory:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./bethe-xyz-dispersion --eta 0.875 --t 0.75 --branch bound \
  --points 4097 --precision fp128 --format csv --no-preamble --no-retain > xyz.csv
```

Use `--points 65537 --precision fp64` for the binary64 case. Both include six
bound branches and their two parity copies. Separate regression comparisons
cover all native formats, selected branches, both theta regimes, tiny gaps,
unavailable gaps and per-momentum overflow.

An opt-in library benchmark evaluates both first and second bound branches
for 256 parameter choices, each at 4,097 reduced momenta: 2,097,664 energies.
Every task prepares its own band/branches and owns its complete curve; the
admission window is 16 parameter curves. It compares every result with direct
serial evaluation before reporting. Scheduler construction and result comparison
are outside the reported interval; preparation, queries and ordered collection
are inside. Reference and actual result arrays are both retained for validation.

```sh
cmake --build . --target bethe_xyz_prepared_benchmark -j2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./tests/bethe_xyz_prepared_benchmark fp128 256 4097
```

The benchmark reports direct serial execution (`workers=0`) and Async execution
with 1/4/8 workers. It waits on the operation's final result; no scheduler-wide
drain is used. It does not add `--threads` to the cheap single-parameter CLI.

On the same workstation, three-repeat median library timings were:

| Execution | Wall (s) | Process CPU (s) |
| --- | ---: | ---: |
| Direct serial | 4.480 | 4.480 |
| Async, 1 worker | 4.470 | 4.470 |
| Async, 4 workers | 1.152 | 4.577 |
| Async, 8 workers | 0.585 | 4.643 |

All 2,097,664 energies match the direct calculation exactly at every worker
count. The whole benchmark's median peak RSS was 73.7 MiB, including both
validation arrays and all four executions; this is not a per-worker memory
measurement. Parameter-curve tasks provide useful granularity even though
individual prepared momentum queries are cheap.


## ASEP shared seeds

At fixed size and filling, ASEP rate queries all begin from the same unit-rate
TASEP solution. `asep::PreparedGap` owns that seed and supports concurrent const
queries with independent continuation scratch. The existing Async preparation
and reader-dependent grid admit rate queries once preparation is ready.
The single-request CLI keeps its
direct solver and analytic fast paths.

Build the [example](../examples/async_asep.cpp) with `BETHE_BUILD_EXAMPLES=ON`,
then run from the build directory:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./examples/bethe-async-asep-example prepared 128 32 1025 8 fp64 > rates.csv
```

The right rate is one; left rates cover [0,0.9] inclusively. `independent`
recomputes the seed per rate; `prepared` constructs it once. Thread count 0
uses direct serial calls and positive counts use Async, capped by the point
count and window 32. All modes use the same numerical controls. The example
retains complete states for subsequent printing, including complex roots and
all diagnostics. The admission bound does not bound that retained output.

### Whole-process measurements

Measurements on 28 October 2026 use an Intel Core i7-12700KF, GCC 13.3 Release,
Uni20 `04163b9b`, oneTBB 2022.3 and one BLAS/OpenMP thread. Binary80 denotes
this machine's 80-bit extended `long double`. Binary128 uses the configured
MPLAPACK 3.0.0 Release provider; normal Uni20 dispatch is unchanged.

Workloads have quarter filling: 128 sites / 1,025 binary64 rates, 64 sites /
1,025 binary80 rates and 64 sites / 65 binary128 rates. Use `long-double` and
`fp128` for the latter precisions. Default budgets are 1,000 seed sweeps,
100 seed Newton updates, 256 continuation attempts and 100 Newton updates per
attempt; both residual tolerances are 128 native epsilons. The kernel,
step-selection rules and precision support are unchanged.

Whole-process wall/CPU times include preparation, scheduler construction,
continuation, ordered collection and CSV printing. Tables give medians over
three repeats, reversing configuration order in the middle repeat. These are
shared-workstation measurements with OS-selected affinity and no overlapping
compiler jobs. They are not an isolated scaling experiment.

| Arithmetic | Independent direct (s) | Independent 8 (s) | Prepared direct (s) | Prepared Async 1 (s) | Prepared 4 (s) | Prepared 8 (s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| binary64 | 12.24 | 1.60 | 8.31 | 8.35 | 2.16 | 1.11 |
| binary80 | 9.87 | 1.31 | 6.98 | 6.98 | 1.80 | 0.94 |
| binary128 | 11.65 | 1.66 | 8.88 | 8.95 | 2.44 | 1.31 |

Seed reuse removes 24–32% of direct serial work in these curves. Scheduling
the prepared queries at eight workers gives another 6.8–7.5x wall-time speedup
relative to prepared direct execution. The combined speedup over independent
direct solves is 8.9–11.0x. Process CPU times distinguish saved preparation
from scheduling overhead:

| Arithmetic | Independent direct CPU (s) | Independent 8 CPU (s) | Prepared direct CPU (s) | Prepared Async 1 CPU (s) | Prepared 4 CPU (s) | Prepared 8 CPU (s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| binary64 | 12.23 | 12.70 | 8.31 | 8.35 | 8.57 | 8.73 |
| binary80 | 9.87 | 10.23 | 6.98 | 6.98 | 7.07 | 7.26 |
| binary128 | 11.64 | 12.10 | 8.87 | 8.95 | 9.10 | 9.21 |

Maximum RSS across these runs is 7.19 MiB. For prepared direct/eight-worker
execution it is 6.41/7.19 MiB in binary64, 6.09/6.41 MiB in binary80 and
6.09/6.25 MiB in binary128. Retained output grows with point count and roots
per state; these figures are not a memory bound for larger sweeps.

All 38,070 complete state rows in 54 timed processes match independent direct
output byte for byte, including every complex root and diagnostic. Separate
matched-main checks cover 60 single-request cases / 291 result and root rows,
including failures, analytic endpoints, tiny bias and rate scaling. All
non-timing/non-provenance metadata and statuses agree exactly. Native, wide and
sanitizer tests also compare full library states and two readers of one
preparation, including numerical seed failure.

Seed work counters in each row describe the shared preparation and must not
be summed as additional CPU work. Timing distinguishes that one-time cost from
the rate queries. An entirely analytic curve should use direct analytic calls;
preparing a TASEP state for it would be unnecessary. The generic grid can also
compose independent geometries or fillings, where each request owns its seed.

## Prepared ladder sweeps

The [ladder example](async-examples.md#prepared-ladder-field-and-rung-sweeps)
separates repeated direct scans from a shared root catalogue and dependent
queries. Reuse and query concurrency are distinct effects; the sequential
preparation remains inside all end-to-end timings. Physical inputs, controls,
full-state diagnostics and output order must match before comparing speed.

Three-repeat whole-process medians on a shared Intel Core i7-12700KF, GCC13.3
Release `-O3`, oneTBB2022.3 and Uni20 `04163b9b` follow. BLAS/OpenMP use one
thread per solve; no compiler jobs overlap measurement. Variant order is
reversed on the second repeat and affinity is left to the OS. Native builds
use normal Uni20 fallback dispatch without MPLAPACK; the binary128 build uses
the configured Release MPLAPACK3 provider. No backend is forced in the code.
All times are seconds and include printing complete states. RSS is the maximum
across worker counts/repeats, including retained results and prepared roots.

| Arithmetic | Rungs / points | Mode | Serial | Async 1 | Async 4 | Async 8 | CPU at 8 | Peak RSS, MiB |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| binary64 | 24 / 1025 | direct | 10.42 | 10.58 | 2.76 | 1.58 | 12.44 | 6.88 |
| binary64 | 24 / 1025 | prepared | 0.12 | 0.13 | 0.05 | 0.04 | 0.19 | 7.19 |
| 80-bit extended | 16 / 2049 | direct | 15.60 | 15.61 | 4.02 | 2.04 | 16.12 | 7.64 |
| 80-bit extended | 16 / 2049 | prepared | 0.09 | 0.09 | 0.05 | 0.05 | 0.15 | 7.81 |
| binary128 | 12 / 513 | direct | 16.07 | 16.04 | 4.13 | 2.13 | 16.80 | 6.56 |
| binary128 | 12 / 513 | prepared | 0.06 | 0.06 | 0.05 | 0.04 | 0.07 | 6.70 |

All 86,088 complete state rows from 72 timed processes converge and match the
direct serial result exactly, including roots, labels, selected populations,
branch diagnostics and scan counters. Independent matched-main CLI checks
retain 1,491 identical rows across 60 requests, including partial budgets,
fixed sectors, both field signs, large shifts and analytic products. Only
revision/invocation/time metadata is excluded from that comparison.

The major benefit is root reuse: prepared serial curves are roughly 90–270x
faster here than repeated serial scans. Those ratios depend strongly on the
number of repeated parameter values and do not describe a single solve.
Repeated direct scans also scale across workers, but retain almost all their
CPU cost. Prepared queries benefit less from parallelism: their whole curves
already take 0.06–0.12 s serially, including the one preparation. Eight-worker
curves take about 0.04–0.05 s, with shared-host variation visible at that scale.
Do not treat this small absolute saving as a reason to parallelize every cheap
query. Four workers can be preferable when CPU efficiency matters.

The native root preparations themselves take about 0.014–0.017 s and binary128
about 0.041 s on these cases. They remain sequential and are included above. Changing
rung count, controls or precision requires another preparation. For a single
request, especially a proven analytic product state, the ordinary direct API
avoids catalogue construction and storage entirely.
