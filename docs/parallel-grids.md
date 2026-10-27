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
