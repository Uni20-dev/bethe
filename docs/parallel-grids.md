# Parallel observable grids

An observable grid evaluates a scalar quantity at independent logical points.
The model owns the physics and prepares its background; the shared
`bethe::evaluate_observable_grid` adapter owns bounded admission and ordered collection.
Uni20 supplies the scheduler. No output format or model formula belongs in the
adapter.

The first frontends using this path are Hubbard and Lieb–Liniger dispersion.
See the [command-line examples](command-line.md#parallel-observable-grids).
Continuum optimizations, thermal grids and other observables can adopt the
same path after checking their independence and scratch ownership.

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
