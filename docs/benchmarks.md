# Performance measurements

[Async API](async-execution.md) · [CLI parallel calculations](parallel-calculations.md) ·
[Validation history](validation-records.md)

Use measurements for the workload and precision actually tested. A single
vacuum/excited pair has at most two concurrent solves; an eight-worker parameter
sweep is a different calculation. Reusing an expensive background can save more
work than parallelizing cheap prepared queries.

| Records | Workloads |
| --- | --- |
| [Observable grids](benchmark-grids.md) | Hubbard/Lieb–Liniger curves, thermodynamic XXX streaming, thermal grids and continuum edges |
| [Shared preparation](benchmark-preparation.md) | XYZ bound branches, ASEP seeds, ladder sector catalogues |
| [Field theory](benchmark-field-theory.md) | Bethe–Yang levels, finite-volume sweeps and vacuum/excited overlap |
| [State and parameter scans](benchmark-state-scans.md) | Kondo, Potts, TASEP/XYZ, nested models and biquadratic clusters |

## Compare like with like

Check complete results and failure diagnostics before timing. Record the source
revision, CPU, compiler, providers, precision, physical parameters and work
budgets. Separate one-time preparation, solve/collection and serialization time.
Report wall time, process CPU time and peak memory; retained validation arrays
or output tables are not per-worker scratch. Repeat on a quiet machine and
report variation or the aggregation used. Keep BLAS/OpenMP threads at one when
parallelizing independent solves.

The records retain their original commands and measurements. Historical native
precision runs remain native precision runs even where the current model also
supports MPFR. No new timings are claimed by this documentation reorganization.

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

## Matched pipeline benchmark

With tests enabled, build the opt-in comparison from a configured build directory:

```sh
cmake --build . --target bethe_async_structure_factor_benchmark -j2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./tests/bethe_async_structure_factor_benchmark two-stage 128 fp64 4 128 10000 old.tsv
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./tests/bethe_async_structure_factor_benchmark async 128 fp64 4 128 10000 new.tsv
cmp old.tsv new.tsv
```

Arguments are mode (`two-stage`, `sync`, `async`), even site count, precision
(`fp64`, `long-double`, configured `fp128`), scheduler limit, admission/batch
bound, candidate budget and output file. The two-stage control retains the
previous public scanner and prepared form-factor loop under `tests/` only.
Each run prints wall/process-CPU timing and writes complete ordered roots,
diagnostics, spectral lines and reductions using native-precision formatting.
Timing excludes file serialization. Compare the files before interpreting
speedup, repeat measurements, and use a quiet machine for performance claims.

## Monte Carlo task granularity

The opt-in `bethe_monte_carlo_benchmark` isolates task granularity with a
synthetic CPU integrand. It is not a timing of the six-spinon form factor:

```sh
cmake --build . --target bethe_monte_carlo_benchmark -j2
# workers, maximum observations per task, replicates, arithmetic iterations
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./tests/bethe_monte_carlo_benchmark 3 16 4 20000
```

Every replicate contains 768 observations. Setting the final argument to zero
exposes overhead for cheap observations; setting the block size to 768 gives
one task per replicate. Repeat cases and compare their means and covariance
before interpreting timing. The benchmark records wall/process CPU time and
peak simultaneous observations; it does not measure peak memory.

On an i7-12700KF with GCC 13.3, four replicates and three workers, median wall
time over three trials fell from 0.506 s with the previous replicate-task engine
to 0.347 s with 16-sample tasks. Process CPU divided by three times wall time
rose from 68% to 98%. A single replicate reached three simultaneous observations
instead of one. Cheap observations showed an overhead penalty: roughly 0.78 ms
for the previous engine, 1.22 ms with blocks of 16, and 4.85 ms with blocks of one
for all 3,072 draws. These short timings are illustrative, not precise throughput
predictions. A compiler process also ran elsewhere on the host; no CPU affinity
was pinned. The [complete records](validation/monte-carlo-sample-tasks-20261101.json)
retain all repeats, controls and matching numerical results. Full production
speedup depends on observation costs and ordered-window backpressure.

A separate small physical query used k=pi/2, omega=2, R=8, MPFR384, four
replicates of 32 draws, and three workers. With the same candidate binary and
identical spectrum/replicate rows, mean query wall time over two reversed-order
runs was 11.16 s for one task per replicate, 9.42 s for 16-sample tasks, 13.64 s
for four-sample tasks, and 21.20 s for one-sample tasks. Query CPU work remained
26.4–26.9 s. Very small tasks with a six-block ordered window can leave workers
idle behind expensive observations; tuning must consider admission depth as
well as task size. The 16-sample default improved this query by about 16%, with
roughly 95% process-CPU utilization of the three-worker limit. The records
include complete controls and export hashes; they do not claim a speedup for
the much larger R=1040 production calculation.
