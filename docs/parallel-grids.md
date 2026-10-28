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
