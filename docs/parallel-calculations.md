# Parallel command-line calculations

[Command-line reference](command-line.md) · [App capabilities](app-capabilities.md) ·
[C++ Async API](async-execution.md) · [Measurements](benchmarks.md)

Concurrency is useful when a tool has enough independent work. Start with one
worker and increase `--threads` for a representative workload; small analytic
queries and sequential preparation may dominate. An admission window bounds
unfinished work, not the size of a retained final spectrum or output table.

## Parallel state scans

`--threads N` is a positive integer, defaulting to 1. It sets a Uni20 oneTBB
scheduler's concurrency limit, not one operating-system thread per state.
The following tools support it:

- `bethe-xxx-structure-factor`: independent states and form factors.
- `bethe-xxx-structure-factor-thermo`: independent spectral points.
- `--excitations` in XXX and XXZ PBC/OBC, Lieb-Liniger PBC/OBC, q-boson PBC,
  quantum-group XXZ OBC, and biquadratic OBC.

The real-root energy tools require an excitation scan for an explicit
`--threads` option, including `--threads 1`. Quantum-group XXZ supports only its
regular-root branch, not the Delta=0 endpoint.

Biquadratic threading additionally supports AF `--sectors` minima and the
ferromagnetic scans `--bound-pairs`, `--bound-triples`, `--bound-quartets`,
`--pair-defects`, `--triple-defects` and `--two-pair-states`. Real-root scans for
either sign include `--ferromagnetic --real-window`. Selected single states,
Q-system, analytic and singlet modes reject explicit threading. Bound-cluster
ranking retains the candidates: `--no-retain` removes output-table storage,
not this numerical result storage. See the [family guide](biquadratic-family.md).

For example, from an installed environment or with the executables on `PATH`:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  bethe-lieb-liniger-pbc 8 --length 8 --c 2 \
  --excitations 10 --padding 2 --threads 4
```

Each state is solved independently; ordered result collection, failure
diagnostics and moment sums are deterministic. Worker threads do not write
tables. Threading does not alter the candidate budget, the family being scanned
or how many converged states are retained. A larger limit can be slower for
small families; benchmark representative sizes and precisions. Metadata records
the requested concurrency limit, not a claim that every worker stayed busy.
Library users can pass `RealExcitationOptions::execution` to the shared scanner
without changing Uni20's active global scheduler.
See the [XXX structure-factor guide](xxx-structure-factor.md#parallel-calculations)
for examples. Set BLAS/OpenMP thread limits to one to avoid nested threading.

## Parallel observable grids

`bethe-xxx-structure-factor-thermo` uses a bounded grid for independent spectral
points, collecting rows in momentum-major order. See its
[streaming and summary conventions](xxx-structure-factor-thermo.md).

`bethe-hubbard-dispersion` (half-filled and doped) and
`bethe-lieb-liniger-dispersion` also accept `--threads N`, defaulting to 1:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  bethe-hubbard-dispersion --u 4 --density 0.5 --points 257 --threads 4
```

The background is prepared once, then independent momentum points are evaluated
in a bounded window that admits a new point after each collected result.
Rows retain their serial branch/momentum order, including
failed points and their diagnostics. Only the calling thread writes tables;
`--no-retain` still permits bounded-memory streaming. `--momentum` accepts the
same option but has only one point per branch. Hubbard places all selected
branches in one work window, so those queries can still overlap; a single
selected branch has no single-point parallelism.

`bethe-hubbard-continuum --threads N` evaluates independent total momenta
through the same bounded grid, for two-spinon and charge-containing edges.
Each point's quadrature, inversion and extrema search remain sequential;
rows preserve momentum order, including failures.

`bethe-lieb-liniger-thermal --threads N` uses the same bounded grid for
temperature scans, at fixed chemical potential or density. Each equilibrium
solve remains sequential; rows preserve the requested temperature order.

The benefit depends on the cost and number of points: background preparation
remains serial and short grids may not benefit. CPU timing covers preparation
and complete grid calls, including row output, since evaluation and collection
can overlap. Export metadata records this as `CPU timing scope`. Set BLAS/OpenMP
thread limits to one to avoid nested threading. The library adapter and its
ownership rules are described in [Parallel observable grids](parallel-grids.md).

`bethe-sine-gordon-bethe-yang` also accepts `--threads N` for its independent
`--levels` calculations. Root searches and phase quadratures stay sequential
inside each level. The same bounded window preserves Bethe-number order and
all failure diagnostics; only the calling thread writes tables. CPU timing
covers the entire level grid including output. One worker, a single level and
the analytic free point use direct serial evaluation; see
[sine-Gordon excitations](sine-gordon-excitations.md).

## Other frontend scans

`bethe-kondo-response` can evaluate a field scan concurrently; `bethe-potts-pbc`
can evaluate charged momenta after one vacuum solve. Their single-point paths
remain direct. See the [Kondo](kondo.md) and [Potts](potts.md) guides for controls.

`bethe-lee-yang-excited` and `bethe-sine-gordon-excited` accept `--threads 2`
for the independent vacuum/excited pair. Higher requests are capped at two
active solves. Their single-scale vacuum tools remain serial. Larger
[parameter sweeps](async-examples.md#finite-volume-integral-equations) are a
library example and can use more workers; they are not the same CLI workload.

## Timing and reproducibility

[CPU/elapsed timing](command-line.md#timing) measures different things. Good
parallel speedup reduces wall time while process CPU consumption can stay the
same or rise. Always compare the complete rows, status values and numerical
controls, not only a checksum of successful energies. [Benchmark records](benchmarks.md)
give their hardware, precision, providers, memory and measurement scope.
