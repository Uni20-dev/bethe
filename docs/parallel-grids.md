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
Continuum optimizations and other observables can adopt the
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
