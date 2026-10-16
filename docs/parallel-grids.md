# Parallel observable grids

An observable grid evaluates a scalar quantity at independent logical points.
The model owns the physics and prepares its background; the shared
`bethe::evaluate_observable_grid` adapter owns batching and ordered collection.
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
    {&scheduler, 128});
```

The borrowed scheduler must outlive the synchronous call. An omitted scheduler
uses Uni20's active scheduler (serial by default); the adapter never changes it.
`ExecutionOptions::batch_size` must be positive and defaults to 128. Logical
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
- The adapter retains at most one batch of result slots. Collection may retain
  more; for streaming output, disable table retention separately.
- A thrown evaluation error joins the batch's workers and collects none of that
  batch. Earlier batches remain collected. A collection error leaves the
  already-collected prefix; all workers have already finished. These exceptions
  are distinct from a normal point result reporting failed convergence.
- An optional fifth argument, `observe_batch(offset, size)`, returns a scope guard
  on the calling thread. Its lifetime covers only synchronous batch execution,
  ending before collection. The frontends return `run_context::computation()`
  here to measure worker CPU time without concurrent writes to timing metadata.

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
preparation, batch evaluation and total grid wall/process-CPU times separately;
collection plus batch-storage overhead is the difference of wall times. There
is no table output inside those intervals. Use process resource accounting for
peak memory, and time the actual frontend separately when measuring formatting,
streaming and file I/O. Small workloads and serial background preparation can
hide otherwise useful evaluation speedups.

### Initial measurement

On an Intel Core i7-12700KF, GCC 13.3 Release (`-O3`, no fast-math), oneTBB
2022.3, and Uni20 `a25159cc`, the binary128 workloads above gave these median
grid wall times over three repeats on 16 October 2026. Both
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
