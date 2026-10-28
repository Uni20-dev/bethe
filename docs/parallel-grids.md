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
