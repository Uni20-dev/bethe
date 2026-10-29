# Observable-grid measurements

[Benchmark index and methodology](benchmarks.md) · [Current capabilities](app-capabilities.md)

These are dated October 2026 measurements. Hardware, precision, providers and
timing scope belong to each record; later capabilities do not change what was
measured. Reproduction commands describe those workloads, not recommended defaults.

## Initial measurement

On an Intel Core i7-12700KF, GCC 13.3 Release (`-O3`, no fast-math), oneTBB
2022.3, and Uni20 `a25159cc`, the [binary128 grid workloads](benchmarks.md#measuring-scaling) gave these median
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

## Hubbard branch scheduling

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

## Thermodynamic XXX streaming measurement

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

## Thermal Lieb–Liniger temperature scans

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

## Hubbard continuum frontend

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
