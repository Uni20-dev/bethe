# Finite-volume and Bethe–Yang measurements

[Benchmark index and methodology](benchmarks.md) · [Current capabilities](app-capabilities.md)

These are dated October 2026 measurements. Hardware, precision, providers and
timing scope belong to each record; later capabilities do not change what was
measured. Reproduction commands describe those workloads, not recommended defaults.

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

The [finite-volume example](async-examples.md#finite-volume-integral-equations)
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
[example guide](async-examples.md#finite-volume-integral-equations). Binary128
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
