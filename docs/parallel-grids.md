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
The library's physics, quadrature and precision support are unchanged;
arbitrary-precision Kondo kernels remain separate work.

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
isolated scaling experiment. The change neither enlarges the selected physical
state family nor adds arbitrary precision. Explicit MPFR requires a separate
audit of the auxiliary XXZ solver, model constants, twists and diagnostics.

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

## TASEP and finite XYZ parameter sweeps

The [complex-root sweep example](async-execution.md#independent-tasep-and-finite-xyz-solves)
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

The [nested sweep example](async-execution.md#owning-nested-model-requests)
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


## Prepared ladder sweeps

The [ladder example](async-execution.md#prepared-ladder-field-and-rung-sweeps)
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
