# Parallel observable grids

An observable grid evaluates a scalar quantity at independent logical points.
The model owns the physics and prepares its background; the shared
`bethe::evaluate_observable_grid` adapter owns bounded admission and ordered collection.
Uni20 supplies the scheduler. No output format or model formula belongs in the
adapter.

The first frontends using this path are Hubbard and Lieb–Liniger dispersion.
See the [command-line examples](command-line.md#parallel-observable-grids).
Thermodynamic XXX spectral grids also use this path, with ordered spectrum rows
and final failure/threshold counts; see [the frontend](xxx-structure-factor-thermo.md).
[Thermal Lieb–Liniger](lieb-liniger-thermal.md) also uses this adapter for
independent temperature solves in either ensemble, preserving each point's
iteration and convergence diagnostics.
[Hubbard continuum edges](hubbard-continuum.md) also use this adapter for
independent total momenta, retaining each point's inversion/search diagnostics.
Sine-Gordon Bethe–Yang levels also use this adapter: each level owns an
independent root search and repeated phase integrals.
Other observables can adopt the same path after checking their independence and scratch ownership.

## Using the library

```cpp
#include <bethe/lieb_liniger_thermo.hpp>
#include <bethe/observable_grid.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace model = bethe::lieb_liniger::thermo;
model::Solver<double> const background(4.0, 1.0);
// Validate the physical request and decide how to report background failure
// before opening output. This example samples the Type-I branch on [0, 6].
uni20::async::TbbScheduler scheduler(4);
std::vector<model::Point<double>> points;
bethe::evaluate_observable_grid(
    129,
    [&](std::size_t i) {
      return background.at_momentum(model::Branch::type_i, 6.0 * double(i) / 128.0);
    },
    [&](std::size_t i, auto point) {
      // Calling thread, increasing i: append a row or accumulate a moment here.
      points.push_back(std::move(point));
    },
    {.scheduler = &scheduler, .max_in_flight = 128});
```

The borrowed scheduler must outlive the synchronous call. An omitted scheduler
uses Uni20's active scheduler (serial by default); the adapter never changes it.
`ExecutionOptions::max_in_flight` bounds active and completed-but-uncollected
points and defaults to 128. `batch_size` remains for other batched clients; both
bounds must be positive. Logical
indices are exactly `[0, count)`, including empty and single-point grids.
Frontends keep branches in their declared output order; a row's identity is
its branch and logical point index. Hubbard dispersion flattens its selected
branches into one logical grid, allowing evaluation across branch boundaries
without changing row order. Lieb–Liniger currently evaluates each branch's grid
in turn. In both cases the admission window bounds retained results.

## Ownership, errors and timing

- Evaluation may run concurrently. Share only immutable prepared data, use
  worker-local scratch, and return an owning value. Move-only results work.
  Do not change global precision settings or run dependent continuation steps
  as independent points.
- Collection runs on the calling thread in increasing index order. Tables,
  diagnostics and floating-point reductions belong here, never in workers.
  Ordinary numerical failure/status values are passed through, not dropped.
- As each result is collected, a new point may enter the bounded window. A slow
  early point or collector applies backpressure; there is no whole-batch join.
  Collection may retain more data; disable table retention for streaming output.
- Execution errors and cancellation stop admission and drain every accepted
  task before propagating. Already collected prefixes remain, but a partial
  grid is never reported as successful completion. Numerical failure/status
  values are ordinary results.
- The optional fifth argument, `observe_wait(index, 1)`, scopes waiting for one
  point, excluding collection. It does **not** measure total computation: workers
  may compute while collection runs. Time the entire grid call for total costs.
  The frontends' `CPU time` includes background preparation and complete grid
  calls, including row output; `CPU timing scope` records this in metadata.

`evaluate_observable_grid_async` offers the same window policy with an owning
Async result. Its collector is serialized **on the scheduler**, so use the
synchronous adapter for caller-thread-only CLI output. See
[Async execution](async-execution.md) for ownership and cancellation.

When the background itself is not ready yet,
[`prepared_observable_grid.hpp`](../include/bethe/prepared_observable_grid.hpp)
adds preparation as an Async result and a reader-dependent overload of the same
grid. Several curves can share one prepared solver without copying its meshes;
see [shared backgrounds](async-execution.md#shared-backgrounds) and the complete
Lieb–Liniger example there.
Threading leaves each scalar solver's arithmetic and stopping rules unchanged.
Native-precision tests compare every field against independent scalar calls,
using both a reverse-order scheduler and oneTBB, for binary64, native
`long double` and binary128 when configured. CLI tests also compare row order,
partial failures and non-retained file output across concurrency limits.

## Measuring scaling

The opt-in benchmark is not part of the default build or CTest. From a configured
build directory with tests enabled:

```sh
cmake --build . --target bethe_observable_grid_benchmark -j2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./tests/bethe_observable_grid_benchmark hubbard-half fp128 257 3
```

Model choices are `hubbard-half`, `hubbard-doped` and `lieb-liniger`; the other
arguments are precision, number of points and optional repeat count (default 3).
The fixed workloads are the half-filled holon branch on $`[-\pi,\pi]`$ at U=4,
the doped holon branch at U=4 and n=0.5, and the Lieb–Liniger Type-I branch
on [0,8] at c=4 and n=1. Backgrounds are prepared once per invocation.

Each invocation times 1/2/4/8-worker sweeps, reverses sweep order on alternate
repeats, and checks an exactly matching ordered energy checksum and failure
count. The checksum sums only converged energies; every failed point is counted
and reported, not treated as a zero-energy success. It reports
preparation and total grid wall/process-CPU times separately, including ordered
checksum collection. There is no table output inside those intervals. Use process resource accounting for
peak memory, and time the actual frontend separately when measuring formatting,
streaming and file I/O. Small workloads and serial background preparation can
hide otherwise useful evaluation speedups.

### Initial measurement

On an Intel Core i7-12700KF, GCC 13.3 Release (`-O3`, no fast-math), oneTBB
2022.3, and Uni20 `a25159cc`, the binary128 workloads above gave these median
grid wall times over three repeats on 16 October 2026, using the previous batched adapter.
These are historical baselines, not measurements of continuous execution. Both
`OPENBLAS_NUM_THREADS` and `OMP_NUM_THREADS` were 1; the MPLAPACK 3.0.0
binary128 provider was an optimized Release build. This was a shared workstation,
not an isolated scaling experiment.

| Workload | Points | Preparation (s) | 1 worker (s) | 2 workers (s) | 4 workers (s) | 8 workers (s) |
|---|---:|---:|---:|---:|---:|---:|
| Half-filled Hubbard holon | 257 | 0.000006 | 52.53 | 26.81 | 13.89 | 7.358 |
| Doped Hubbard holon | 65 | 0.7532 | 4.796 | 2.472 | 1.304 | 0.7114 |
| Lieb–Liniger Type I | 4097 | 0.01124 | 2.822 | 1.413 | 0.7690 | 0.4369 |

All these points converged with identical ordered energy checksums. Grid times
include collection but exclude the one-time preparation and table I/O: the
8-worker speedups are 7.14, 6.74 and 6.46 respectively. Including preparation
reduces the small doped-grid speedup to 3.79. Median process CPU times increase
from 52.53 to 55.69 s, 4.796 to 5.133 s, and 2.822 to 3.291 s. Median collection
plus batch-storage overhead is below 0.3 ms for every configuration. Whole-process
peak RSS across each complete sweep was 5.15, 4.66 and 4.68 MiB respectively.

A separate single-pass frontend check included preparation and JSON streaming
with `--no-retain`, comparing complete rows at 1/2/4/8 workers in all three
precisions. In binary128, 1-to-8-worker elapsed times were 7.126 to 1.019 s
(33-point half-filled holon), 19.83 to 3.287 s (257-point doped holon), and
2.789 to 0.4161 s (4097-point Type I). Peak RSS over all frontend checks ranged
from 5.2 to 6.1 MiB. Small binary64 workloads were already only milliseconds
long; their output and scheduling overhead make speedup much less informative.

On the dense Type-I grid, native 80-bit `long double` reported one `mesh_limit`
at p=5.412109375 with the default numerical controls. That row, its missing
energy, and the partial outcome were identical at every worker count; binary64
and binary128 completed this grid. Parallelism preserves these numerical
diagnostics rather than changing the tolerance or silently dropping points.

### Hubbard branch scheduling

A matched frontend comparison on 28 October 2026 measures the additional effect
of removing the wait between Hubbard branches. The baseline is Bethe `46b636e`,
which already uses a continuous window within each branch. The new version puts
all selected branches in one window. Both prepare a doped background once;
this comparison does not measure a change in preparation cost.

The Intel Core i7-12700KF workstation used GCC 13.3 Release, oneTBB 2022.3,
Uni20 `04163b9b` and optimized MPLAPACK 3.0.0. Both builds used the same
dependencies and one BLAS/OpenMP thread. These are medians of three complete
process runs, reversing measurement order on alternate repeats, including
preparation and CSV output with retention disabled.

| Workload / version | Rows | 1 worker (s) | 4 workers (s) | 8 workers (s) | CPU at 8 workers (s) | Peak RSS at 8 workers (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| Half-filled, separate branches | 195 | 34.13 | 8.85 | 4.72 | 35.33 | 7.00 |
| Half-filled, shared window | 195 | 34.22 | 8.82 | 4.48 | 35.25 | 7.50 |
| Doped, separate branches | 387 | 21.59 | 6.13 | 3.53 | 22.46 | 7.19 |
| Doped, shared window | 387 | 21.56 | 6.12 | 3.48 | 22.43 | 7.64 |

The workloads use all three branches at U=4 in binary128: 65 points per branch
at half filling and 129 at density 0.5. All full CSV outputs, including
diagnostics and failed-value fields, match byte for byte between versions at
1/4/8 workers; every point converged. Reproduce each workload with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 bethe-hubbard-dispersion \
  --u 4 --points 65 --precision fp128 --threads 8 \
  --format csv --no-preamble --no-retain > half.csv
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 bethe-hubbard-dispersion \
  --u 4 --density 0.5 --points 129 --precision fp128 --threads 8 \
  --format csv --no-preamble --no-retain > doped.csv
```

Vary `--threads` and compare equally configured builds. The incremental wall-time
reduction here is modest (about 5% and 1% at eight workers), and the smaller
difference should not be generalized from a shared workstation. The substantial
serial-to-parallel benefit was already present in the baseline. The shared
window removes unnecessary branch waits and permits single-momentum queries
across branches to overlap; the separate
[Async example](async-execution.md#shared-backgrounds) demonstrates composable
preparation and dependent curves without introducing another numerical solver.

### Thermodynamic XXX streaming measurement

The thermodynamic XXX frontend now collects directly into its output tables.
With `--no-retain`, it no longer stores a second full spectrum in memory.
The numerical kernel is unchanged. Both implementations produce identical
CSV rows, including density, estimated error, evaluation count and status.

The 28 October 2026 comparison uses a fresh build of Bethe `46b636e` and the streaming frontend,
both with Uni20 `04163b9b`, GCC 13.3 Release and oneTBB 2022.3 on an Intel Core
i7-12700KF (12 physical cores, 20 logical CPUs). It is a shared workstation;
short timings are indicative, not isolated scaling measurements. BLAS and
OpenMP each use one thread. Three repeats alternate the old/new execution
order, reporting medians. The entire process is timed, including CSV formatting.
CSV output goes to `/dev/null`, except the first one-worker pair, which is saved
and compared byte for byte. Native formats group 16 points per task; binary128
uses individual points to balance variable quadrature costs.

Reproduce each row with either executable, selecting `N` as 1, 4 or 8:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /usr/bin/time -v \
  bethe-xxx-structure-factor-thermo --momentum-points 257 --points 1025 \
  --precision fp64 --threads N --format csv --no-preamble --no-retain > /dev/null
```

For binary128 use `--momentum-points 129 --points 257 --precision fp128`.
An exterior-only control uses `--momentum-points 1000 --points 1000`
with `--omega-min -2 --omega-max -1 --precision fp64`: all one million
values are zero, so output and scheduling dominate. This control tests memory scaling,
not useful numerical speedup.

| Grid / arithmetic | Workers | Old / new wall (s) | Old / new CPU (s) | Old / new peak RSS (MiB) |
| --- | ---: | ---: | ---: | ---: |
| 257 × 1025, binary64 | 1 | 1.13 / 1.13 | 1.12 / 1.47 | 16.71 / 4.82 |
| 257 × 1025, binary64 | 4 | 0.56 / 0.49 | 1.28 / 1.87 | 16.72 / 4.66 |
| 257 × 1025, binary64 | 8 | 0.45 / 0.37 | 1.37 / 2.60 | 16.70 / 4.84 |
| 129 × 257, binary128 | 1 | 68.79 / 68.60 | 68.78 / 68.70 | 7.33 / 5.30 |
| 129 × 257, binary128 | 4 | 18.20 / 17.77 | 70.22 / 70.13 | 7.34 / 5.28 |
| 129 × 257, binary128 | 8 | 9.86 / 9.29 | 71.89 / 71.90 | 7.34 / 5.31 |
| 1000 × 1000, exterior binary64 | 1 | 0.79 / 1.13 | 0.79 / 2.26 | 50.31 / 4.69 |
| 1000 × 1000, exterior binary64 | 4 | 1.11 / 1.17 | 1.21 / 4.49 | 50.31 / 4.81 |
| 1000 × 1000, exterior binary64 | 8 | 1.11 / 1.15 | 1.34 / 7.43 | 50.27 / 4.82 |

The physical binary128 grid takes about 7.4 times less wall time at eight
workers than at one, with similar total CPU work. Binary64 physical-grid wall
time also improves over the batched frontend, but scheduling and overlapping
output consume more CPU time. The exterior-only control gains no useful speedup;
one worker is preferable for such cheap grids. Its peak RSS falls from about
50 MiB to about 5 MiB with table retention disabled. The momentum and frequency
axes still occupy storage proportional to their individual lengths.

### Thermal Lieb–Liniger temperature scans

On 28 October 2026, the thermal frontend was compared with a fresh serial
baseline at Bethe `46b636e`, using identical GCC 13.3 Release builds, Uni20
`04163b9b` and oneTBB 2022.3 on an Intel Core i7-12700KF. BLAS/OpenMP each
used one thread. Each row below is the median of three complete-process timings,
including CSV formatting with `--no-retain`; execution order alternated between
repeats. All exported rows, including error estimates and iteration counts,
matched the baseline exactly at every concurrency. The workstation was shared,
so these are workload examples rather than isolated scaling measurements.

Both ensembles use c=4 and a uniform temperature scan from 1 to 2. Reproduce
the binary128 grand-canonical case with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /usr/bin/time -v \
  bethe-lieb-liniger-thermal --c 4 --temperature 1 --temperature-end 2 \
  --points 33 --mu -1 --precision fp128 --threads 8 \
  --format csv --no-preamble --no-retain > /dev/null
```

The binary64 cases use 257 points at `--mu -1`, or 129 points at
`--density 0.1`, with `--precision fp64`. Numerical controls otherwise keep
their precision-dependent defaults. The old baseline has no `--threads` option.

| Ensemble / arithmetic | Points | Old serial (s) | 1 worker (s) | 4 workers (s) | 8 workers (s) | Old / 8-worker CPU (s) | Old / 8-worker peak RSS (MiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Fixed mu, binary64 | 257 | 0.53 | 0.54 | 0.14 | 0.08 | 0.52 / 0.59 | 6.61 / 11.18 |
| Fixed density, binary64 | 129 | 1.88 | 1.91 | 0.48 | 0.26 | 1.87 / 2.01 | 6.69 / 11.68 |
| Fixed mu, binary128 | 33 | 47.75 | 48.03 | 13.37 | 7.45 | 47.74 / 49.70 | 9.99 / 33.80 |

The binary128 scan gains about 6.4 times in wall time at eight workers, with
about 4% more total CPU work than the old serial frontend. Independent solver
workspaces increase peak memory with concurrency; bounded result admission
does not bound total solver scratch independently of worker count. Single-point
requests and very short scans have no corresponding parallel benefit.

### Hubbard continuum frontend

Independent total-momentum queries also scale for continuum edges: the outer
grid can overlap each point's quadrature, inversion and charge-channel extrema
search. The inner algorithms and stopping rules are unchanged. The frontend
uses the shared 128-point admission window and writes ordered rows on the caller.

The following measurements on 28 October 2026 compare the serial frontend at
Bethe `46b636e` with this continuous-grid implementation. Both builds used
Uni20 `04163b9b`, GCC 13.3 Release (`-O3`, no fast-math), oneTBB 2022.3 and
the same optimized MPLAPACK 3.0.0 configuration on an Intel Core i7-12700KF.
These are whole-command wall/process-CPU times and peak RSS, including startup
and CSV formatting. Three repeats alternate command order on a shared
workstation, with no simultaneous compilation. BLAS and OpenMP each use one
thread. All complete CSV rows match the serial baseline exactly at 1/4/8
workers, including the error estimates, work counts and status columns.

From the corresponding build directory, the two workloads are:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./bethe-hubbard-continuum --u 4 --channel two-spinon --points 129 \
  --precision fp128 --threads 8 --format csv --no-preamble --no-retain
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./bethe-hubbard-continuum --u 4 --channel spinon-holon --points 65 \
  --search-tolerance 1e-7 --position-tolerance 1e-5 --tolerance 1e-10 \
  --precision fp64 --threads 8 --format csv --no-preamble --no-retain
```

Use `--threads 1` or `--threads 4` for the other measurements; omit that option
for the old serial frontend. The looser charge controls are part of the
benchmark definition, not a change to the defaults or a certified error bound.

Median wall times:

| Family | Precision | Points | Old serial (s) | New 1 worker (s) | New 4 workers (s) | New 8 workers (s) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Two spinons | binary128 | 129 | 28.55 | 29.64 | 7.65 | 3.95 |
| Spinon–holon | binary64 | 65 | 18.79 | 18.79 | 4.89 | 2.64 |

Relative to the old serial loop, eight workers give 7.23x and 7.12x speedups.
Process CPU time increases from 28.55 to 30.64 s and from 18.78 to 19.65 s;
peak RSS increases from 4.22 to 5.45 MiB and from 4.22 to 6.41 MiB.
The new single-worker two-spinon path costs about 4% more than the old loop,
while the charge workload is unchanged at this timing resolution. Thus the
benefit is concurrent expensive queries; small grids need not become faster.
Per-query numerical scratch grows with active workers even though admission
and retained results are bounded.

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


## Sine-Gordon Bethe–Yang level scans

`bethe-sine-gordon-bethe-yang --threads N` schedules independent interacting
levels with the shared bounded grid. Each level keeps its bisection and phase
quadratures sequential. The window holds at most 128 active or uncollected
results; collection preserves Bethe-number order, complete failure diagnostics
and calling-thread table output. One worker, one level and the analytic free
point use direct serial evaluation. Execution mode, requested concurrency and
whole-grid CPU timing are recorded in metadata.

The CLI regression compares every field under direct and threaded execution
for each configured native precision. It includes complete scans, mixed
successes/failures, exhausted iteration/quadrature/cutoff budgets, both charge
signs, non-retained CSV/TSV/JSON, a 129-level interacting scan crossing the
admission window, and preservation of existing files on invalid requests.
Existing independent energy/phase references remain unchanged.

Benchmarks compare main `46b636e` with this implementation, both using GCC 13.3
Release, Uni20 `04163b9b`, oneTBB 2022.3 and one BLAS/OpenMP thread. They ran on
a shared Intel Core i7-12700KF desktop on 28 October 2026, without concurrent
local compilation. The end-to-end timings include startup, solving and CSV
output; each result is the median of three repeats, with execution order
alternated. Peak RSS is the maximum observed for that variant. All numerical
CSV rows match the baseline exactly at one/four/eight workers, including 2,066
benchmark rows and 495 additional selected/failure/free rows. This comparison
does not establish the physical accuracy of the asymptotic Bethe–Yang levels.

The binary64 workload uses 2,049 levels at tolerance 1e-9; the binary128 workload
uses 17 levels at its default tolerance. Both use mass 1, p=2 and L=100. For
example, compare the following with `--threads 1`, `4` and `8`; omit `--threads`
for the baseline executable:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  bethe-sine-gordon-bethe-yang --p 2 --length 100 --levels 2049 \
  --tolerance 1e-9 --precision fp64 --threads 8 \
  --format csv --no-preamble --no-retain > levels.csv

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  bethe-sine-gordon-bethe-yang --p 2 --length 100 --levels 17 \
  --precision fp128 --threads 8 \
  --format csv --no-preamble --no-retain > levels-wide.csv
```

| Workload | Execution | Wall (s) | CPU (s) | Peak RSS (MiB) |
| --- | --- | ---: | ---: | ---: |
| binary64, 2,049 levels | baseline serial | 6.33 | 6.33 | 3.91 |
| binary64, 2,049 levels | new, 1 worker | 6.32 | 6.32 | 4.50 |
| binary64, 2,049 levels | new, 4 workers | 1.61 | 6.43 | 5.00 |
| binary64, 2,049 levels | new, 8 workers | 0.85 | 6.82 | 5.00 |
| binary128, 17 levels | baseline serial | 30.82 | 30.81 | 4.06 |
| binary128, 17 levels | new, 1 worker | 31.94 | 31.94 | 4.84 |
| binary128, 17 levels | new, 4 workers | 9.43 | 32.33 | 5.16 |
| binary128, 17 levels | new, 8 workers | 5.71 | 32.71 | 4.93 |

The eight-worker wall-time improvement against the baseline is about 7.4x for
binary64 and 5.4x for binary128, with modest additional total CPU work. The new
one-worker binary128 process is about 4% slower in this measurement; the table
includes that cost. The smaller 17-level workload has limited parallelism and
uneven finishing times. These shared-desktop measurements do not establish
scaling beyond eight workers or a benefit for short/free scans. No worker
changes the arithmetic, tolerance or stopping criteria of its scalar solve.

## Finite-volume integral-equation sweeps

The [finite-volume example](async-execution.md#finite-volume-integral-equations)
uses the existing Async grid for independent Lee–Yang and sine-Gordon solves.
Three-repeat medians below were measured on a shared Intel Core i7-12700KF,
GCC 13.3, oneTBB 2022.3 and Uni20 `04163b9b`, with one BLAS/OpenMP thread.
Sweep order was reversed on the second repeat. Times are seconds for solving
and ordered collection; scheduler construction, initial result preallocation
and printing are outside that interval. Every benchmark output row matches
its direct serial counterpart exactly; separate tests compare complete state
diagnostics, including failures, under deferred and threaded execution.

| Model/state | Arithmetic | Points | Direct | Async 1 | Async 4 | Async 8 | CPU at 8 | Peak RSS, MiB |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Lee–Yang vacuum | binary64 | 10,000 | 2.065 | 2.102 | 0.537 | 0.279 | 2.225 | 5.289 |
| Lee–Yang excited | binary64 | 2,048 | 5.571 | 5.571 | 1.425 | 0.738 | 5.897 | 4.219 |
| Sine-Gordon vacuum | binary64 | 256 | 4.891 | 4.896 | 1.262 | 0.650 | 5.107 | 4.062 |
| Sine-Gordon excited | binary64 | 64 | 4.928 | 4.946 | 1.300 | 0.739 | 5.467 | 4.062 |
| Lee–Yang excited | binary128 | 128 | 21.949 | 21.953 | 5.630 | 2.902 | 22.664 | 4.219 |
| Sine-Gordon excited | binary128 | 4 | 33.872 | 33.915 | 9.145 | 9.139 | 34.688 | 4.211 |

The four-point binary128 sine-Gordon example has only four independent tasks;
eight threads cannot improve on four. The longer sweeps show useful scaling,
while one-worker scheduling has a small cost in the cheap Lee–Yang vacuum
case. These are desktop measurements, not guarantees for cluster hardware or
other tolerance/length ranges. Peak RSS covers the complete process and all
retained states, not just active tasks.

To reproduce a row, build `bethe-async-finite-volume-example` with examples
enabled, then run it with the model names and point counts above, precisions
`fp64`/`fp128`, and thread counts `0`, `1`, `4`, `8`:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./examples/bethe-async-finite-volume-example lee-yang-excited 8 128 fp128
```

The model inputs and fixed discretization tolerances are specified in the
[example guide](async-execution.md#finite-volume-integral-equations). Binary128
here improves arithmetic but does not tighten the chosen mesh tolerance.

### Overlapping one vacuum/excited pair

The excited-state CLI option `--threads 2` has a different limit: there are
only two tasks and their costs can be unequal. Three-repeat end-to-end CLI
medians, including JSON output and scheduler startup, are:

| Binary128 request | Main serial | New serial | Two solves | CPU, two solves | Peak RSS, MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| Lee–Yang, L=5 | 10.96 | 10.99 | 10.85 | 10.95 | 5.31 |
| Sine-Gordon, p=2, L=1, tolerance 1e-9 | 19.34 | 19.39 | 15.08 | 19.41 | 5.46 |

Mass is one and the sine-Gordon positive Bethe number is 1/2. The baseline is
main `46b636e`, built with the same compiler and Uni20 revision. All pair
benchmark table rows match exactly. Lee–Yang's roughly 1% gain is negligible
for practical planning; sine-Gordon gains about 22% in elapsed time. The
one-worker default avoids scheduling overhead for short requests. A multi-state
parameter sweep is the better use of a larger CPU allocation.

```sh
bethe-lee-yang-excited --length 5 --precision fp128 --threads 2 --format json
bethe-sine-gordon-excited --p 2 --length 1 --tolerance 1e-9 \
  --precision fp128 --threads 2 --format json
```


Matched-main comparisons cover 30 requests / 120 rows across binary64,
long double and binary128, including partial failures. Physical observables,
statuses and work counters are unchanged. Two long-double Lee–Yang direct-tail
bounds differ between the compiled frontends by less than `5e-19` relative;
these very small diagnostics are far below the requested error tolerance.
Within the candidate, every table field is identical at one/two/eight requested
workers.
