# Parallel observable grids

An observable grid evaluates a scalar quantity at independent logical points.
The model owns the physics and prepares its background; the shared
`bethe::evaluate_observable_grid` adapter owns bounded admission and ordered collection.
Uni20 supplies the scheduler. No output format or model formula belongs in the
adapter.

The first frontends using this path are Hubbard and Lieb–Liniger dispersion.
See the [command-line examples](command-line.md#parallel-observable-grids).
Sine-Gordon Bethe–Yang levels also use this adapter: each level owns an
independent root search and repeated phase integrals. Continuum optimizations,
thermal grids and other observables can adopt the same path after checking
their independence and scratch ownership.

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
Frontends keep branches in their declared order and evaluate each branch's grid
in turn; a row's identity is its branch and logical point index.

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
