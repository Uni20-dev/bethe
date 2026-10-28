# Async execution

Bethe uses Uni20's scheduler for independent calculations. It creates no thread
pool of its own and never changes the global scheduler. Existing synchronous
scans wait for the same operation exposed by the Async API.
See [app coverage](app-capabilities.md) for the execution and arbitrary-precision
paths currently exposed by each frontend, and the remaining integration work.

## State scans

```cpp
#include <bethe/heisenberg_excitations.hpp>
#include <uni20/async/tbb_scheduler.hpp>

uni20::async::TbbScheduler scheduler(4);
bethe::RealExcitationOptions options;
options.count = 10;
options.execution.scheduler = &scheduler;
options.execution.max_in_flight = 8;
auto pending = bethe::heisenberg::real_excitations_async<double>(
    32, uni20::half_int{1}, options);
auto const& result = pending.get_wait(scheduler);
```

The periodic and open XXX and finite Lieb–Liniger APIs provide
`real_excitations_async`. Other models
using the shared real-root scanner also use continuous execution through their
existing synchronous APIs. After preparing the ground reference, a bounded
number of runners repeatedly solve candidates. There is no barrier between
successive groups of candidates and no task graph proportional to family size.
The top-K collector retains only the requested levels; ties and the first
failed-state diagnostic use enumeration identity, independent of completion
order. Final energy ordering and numerical convergence rules are unchanged.

`max_in_flight` limits active candidates, not scheduler threads. It defaults to
128. `batch_size` remains for clients that still use batches; both must be
positive. Set BLAS/OpenMP thread limits to one when parallelizing many solves.

For finite Lieb–Liniger, include `<bethe/lieb_liniger.hpp>` or
`<bethe/lieb_liniger_open.hpp>` and submit
`lieb_liniger::real_excitations_async(particles, length, c, padding, scan, solver)`
(or `lieb_liniger::open::real_excitations_async` for hard walls). These entry
points own the physical parameters, label-window controls and solver options.
The periodic result also retains its window and physical metadata. The same
finite-family cutoff and top-K selection apply as in the blocking APIs.
Independent MPFR scans can carry different precisions on the same scheduler;
they require the configured MPLAPACK MPFR real-solve provider.

## XXX two-spinon structure factor

`two_spinon_structure_factor_async<Real>(sites, solver, budget, execution)`
owns the complete ground preparation, root solving, form-factor evaluation and
ordered reduction. Each runner retains its candidate through form-factor
evaluation before admitting another. A single immutable
`PreparedRaisingFormFactors<Real>` snapshot is constructed before the runners
start; it is not rebuilt per state. The synchronous entry point waits for this
same operation.

All roots and factor diagnostics remain attached through sorting. Failed roots
never enter the factor kernel, failed factors remain explicit diagnostics, and
final spectral sums follow canonical energy order. The API retains the entire
requested two-spinon family; the admission bound limits transient work, not the
final result.

MPFR callers supply an explicit `SolverOptions<uni20::mpreal>` with finite
`uni20::Precision`. The operation retains that precision in its prepared data,
constants, empty momentum bins and reductions. Independent precisions may run
on the same scheduler without changing ambient defaults.

## Observable grids

`evaluate_observable_grid_async(count, evaluate, collect, execution)` returns
`uni20::async::Async<std::size_t>` containing the completed point count. A bounded
window retains active and completed-but-uncollected points until ordered
collection releases their slots. The collector runs serially on the scheduler.
The synchronous adapter instead collects on its calling thread; cooperative
waits may execute a kernel there. See [parallel grids](parallel-grids.md).

Independent [Kondo field responses](kondo.md) also compose through this API.
Capture the field/scale requests and numerical controls by value so they
survive submission; return the complete `kondo::State<Real>`, including failed
observables and work/error diagnostics. Every call has its own Euler series
and quadrature scratch. No background preparation or global scheduler change
is required. The frontend exposes this as `--field-end`, `--points` and
`--threads`, collecting rows on its calling thread.

The [Potts frontend](potts.md) uses this grid for its charged momenta after
one vacuum solve. Charge conjugation reuses each result for both charge rows.
Independent library calls to `potts::ground_state` and `charged_one_hole` can
also be captured as owning requests in an Async grid: retain the full returned
state, including auxiliary roots, labels, residual and convergence diagnostics.
The root iterations themselves stay sequential.

## Ownership, dependencies and cancellation

[ASEP rate curves](asep.md#reusing-a-seed-across-rate-queries) can share an owning
`asep::PreparedGap` through the reader-dependent grid. Preparation solves one
unit-rate TASEP state; each rate query owns its continuation scratch and returns
its full diagnostics. The compiled [example](../examples/async_asep.cpp) submits
the consumer before preparation is ready, then waits on the final grid result.
Numerical seed failure does not erase an independently valid analytic query.

For curves that depend on a prepared model, see [shared backgrounds](#shared-backgrounds)
below. They use the same bounded grid and add a reader dependency on preparation.

Async entry points own their options and decay-copied callables. References
captured inside callbacks remain borrowed. The scheduler and borrowed data must
outlive completion, even if the result handle is discarded. Evaluators must
permit concurrent calls; prepared data should be immutable and scratch local.
No operation changes ambient arithmetic precision.

Pass `pending.read()` to a dependent Uni20 coroutine to consume a result without
a host wait. Errors propagate through the result dependency. Destroying the
handle does not cancel work; use `ExecutionOptions::stop` with a
`std::stop_source`. A stop request prevents further admission and follow-up work
at checkpoints, then waits for active numerical kernels to finish. Kernels are
not interrupted mid-solve.

Scans throw `ScanCancelled` (derived from `ExecutionCancelled`); grids throw
`ExecutionCancelled`. They do not return a partly completed family as success.
Execution errors take precedence over cancellation and propagate only after
owned work has drained. Numerical nonconvergence remains a normal diagnostic
value. A grid may already have emitted an ordered prefix when it fails.

Argument and budget errors are checked before scientific work. A scheduler
admission failure is an execution error, not a failed root or zero spectral
weight. No scheduler-wide drain is used as the completion boundary for an
individual operation.

## Shared backgrounds

Include `<bethe/prepared_observable_grid.hpp>` to submit preparation and its
dependent grids together. `prepare_background_async(prepare, execution)` runs
the synchronous preparation callable once and returns its owning result as
`Async<Background>`. The prepared-grid overload accepts a reader and calls
`evaluate(background const&, index)` for each point:

```cpp
auto background = bethe::prepare_background_async(
    [] { return model::Solver<double>(4.0, 1.0); }, execution);
auto particle = bethe::evaluate_observable_grid_async(
    background.read(), count,
    [pmax, count](model::Solver<double> const& solver, std::size_t i) {
      return solver.at_momentum(model::Branch::type_i,
                                pmax * double(i) / double(count - 1));
    }, collect_particle, execution);
auto hole = bethe::evaluate_observable_grid_async(
    background.read(), count,
    [pmax, count](model::Solver<double> const& solver, std::size_t i) {
      return solver.at_momentum(model::Branch::type_ii,
                                pmax * double(i) / double(count - 1));
    }, collect_hole, execution);
```

Here `model` is `bethe::lieb_liniger::thermo`, `count` is at least two, and
`pmax` is within the type-II momentum range. Both curves are submitted without
observing the background on the calling thread. They share one scheduler pool;
the waiting coordinators do not reserve a fixed allocation of CPUs.

Each grid holds its reader until all admitted queries have finished, including
failure cleanup. Destroying the original background handle is safe after readers
are enrolled. No prepared meshes are copied per grid or per point. The prepared
type must support concurrent const queries; mutable scratch belongs to each query.
The maximum-in-flight limit applies independently to each grid, while the
scheduler limits shared CPU concurrency. This adds composition; it does not
establish a speedup over an already-prepared synchronous background.

Readiness is separate from numerical convergence. Lieb–Liniger's `at_momentum`
already checks preparation status and returns explicit failed points when the
background did not converge; a zero-energy endpoint remains a successful point.
The generic adapter does not invent a convergence rule for other models. Their
evaluators must preserve preparation diagnostics or reject the failed background.
Preparation exceptions propagate to every dependent grid, including empty grids.

Cancelling one grid stops its queries without cancelling preparation or sibling
grids. A grid waiting for preparation observes cancellation after that dependency
finishes. Preparation itself checks its own stop token before and after its
synchronous kernel; neither preparation nor point kernels are interrupted midway.
Errors in an input dependency take precedence over a dependent grid's cancellation.

The [complete example](../examples/async_background.cpp) submits both branches,
observes all accepted operations even if one fails, and then prints the curves.
Collectors run on the scheduler, so the example collects separate vectors and
prints after completion. It also drains preparation if grid submission fails.
From a separate build directory:

```sh
cmake -DBETHE_BUILD_EXAMPLES=ON /path/to/bethe
cmake --build . --target bethe-async-background-example -j2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./examples/bethe-async-background-example
```

The example is built in CI and is not installed as a frontend executable.

The [doped Hubbard example](../examples/async_hubbard_background.cpp) uses the
same APIs with `bethe::hubbard::thermo::DopedSolver<double>(4.0, 0.5)`. It
submits spinon, holon and charge-particle curves before waiting for preparation,
then collects them at a shared Fermi energy reference. Momentum intervals come
from the prepared solver, so each evaluator owns only its branch selection;
all three readers share the same coarse and fine meshes. The model's const
queries use local inversion and accumulation state, with no mutable mesh cache.

```sh
cmake --build . --target bethe-async-hubbard-background-example -j2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./examples/bethe-async-hubbard-background-example
```

For this model, a failed background produces failed points with the original
background status, even at nominal gapless endpoints. A failed momentum
inversion preserves its own status and missing energies. Tests compare every
point field with serial calls in all configured fixed precisions, using queued
and oneTBB schedulers and exhausting background/query budgets separately.
The ordinary Hubbard frontend also shares one prepared background, but waits
for it before validating branch intervals and opening caller-thread output.
Its queries across selected branches share a single bounded window. The
example demonstrates a client whose dependent computations are submitted
before preparation completes; it does not claim faster preparation itself.

## Ground-state parameter sweeps

The same grid can evaluate independent ground states at different couplings,
densities or system sizes. Each evaluator owns its numerical controls and
returns the model's complete result, including convergence diagnostics. Ordered
collection keeps parameter identities stable even when solves finish out of order.
The operation's final result covers every admitted solve and collected row;
wait on it before reading caller-owned output storage.

The [ground-state sweep example](../examples/async_ground_states.cpp) evaluates
four thermodynamic Lieb–Liniger couplings with a maximum of two active solves:

```sh
cmake -DBETHE_BUILD_EXAMPLES=ON /path/to/bethe
cmake --build . --target bethe-async-ground-states-example -j2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./examples/bethe-async-ground-states-example
```

It uses explicit 160-bit `mpreal` when MPLAPACK MPFR is configured, otherwise
binary64. The printed header identifies the working precision. It waits once
on the grid result, then checks scientific convergence before reporting energies.

Thermodynamic Lieb–Liniger also supports `mpreal` in the shared-background
composition above. Construct `Options<uni20::mpreal>` with a finite
`uni20::Precision`, and construct all model parameters, tolerances and momenta
at that precision. Independent backgrounds can use different precisions on one
scheduler. More bits do not replace mesh refinement or increase iteration
budgets; see [thermodynamic controls](lieb-liniger-thermo.md#c-api).

## Independent TASEP and finite XYZ solves

The [complex-root sweep example](../examples/async_complex_sweeps.cpp) uses
the ordinary owning grid for two independent families:

- TASEP: fixed ring size, particles from 1 through `POINTS`, at unit jump rate.
  The count cannot exceed half the sites, avoiding particle-hole duplicates.
- XYZ: fixed even ring size and elliptic parameter `t=3/4`, with eta equally
  spaced from 1/5 through 4/5. One point selects eta=1/2.

From a build directory configured with `BETHE_BUILD_EXAMPLES=ON`:

```sh
cmake --build . --target bethe-async-complex-sweeps-example -j2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./examples/bethe-async-complex-sweeps-example tasep 128 64 4 fp64
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./examples/bethe-async-complex-sweeps-example xyz 32 129 4 long-double
```

Arguments are model, sites, points, threads and precision (`fp64`,
`long-double`, configured `fp128`). Zero threads calls each solver directly;
positive counts use Async, capped at the smaller of 32 and the point count.
The callable owns its inputs and options; every solve owns its complex roots,
Newton scratch and diagnostics. Completion is the final grid result. Collection
is ordered on the scheduler, and printing happens afterward. The admission
window is 32; this demonstration retains every output state, so output memory
grows with the number of points. It prints all state fields, including roots,
and returns status 2 if any numerical solve fails.

Tests submit both models onto the same scheduler, including deferred and
single-worker execution, after the submission-local request vectors disappear.
They compare every state field with direct calls, retain numerical failures,
and check that an invalid request in either model does not cancel the other
grid. These are fixed-precision native solvers; this example does not enable
MPFR. The existing single-request frontends remain serial. TASEP root updates
and XYZ Newton steps are coupled and remain sequential. For a TASEP rate-only
curve, reuse a unit-rate result and rescale its eigenvalue/gap; the roots do not
depend on the positive jump rate.

## Prepared ladder field and rung sweeps

The [ladder example](../examples/async_ladder.cpp) prepares an owning
`ladder::PreparedSectors<Real>` once and submits a bounded query grid through
`prepare_background_async` and `evaluate_observable_grid_async`. Preparation
remains sequential to preserve the whole-scan budgets and early failure order.
Queries reuse roots, retaining exact partial-state diagnostics and the original
analytic product paths. See the [numerical contract](ladder.md#reusing-roots-across-field-and-rung-queries).

From a build directory configured with `BETHE_BUILD_EXAMPLES=ON`:

```sh
cmake --build . --target bethe-async-ladder-example -j2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./examples/bethe-async-ladder-example prepared 12 257 4 fp64
```

Arguments are `MODE RUNGS POINTS THREADS PRECISION`. `MODE=direct` repeats the
ordinary ground-state scan; `MODE=prepared` reuses one catalogue. Zero threads
runs serially; positive counts use Async with at most 32 active/uncollected
requests. For more than one point, h spans [-2,2] inclusively and J_r=2-h/2;
a single point uses h=0, J_r=2. Precisions are `fp64`, `long-double` and configured
`fp128`. Failed states are printed and cause exit status 2, without loosening
controls. These parameter sweeps belong to the library example, not a new
frontend option.

The example waits on the grid's final result. If consumer admission fails after
preparation was accepted, it also observes that preparation before unwinding.
All output is retained and printed afterward in enumeration order, including
nested roots, labels and every state diagnostic. Thus output memory grows with
POINTS even though admission is bounded. CSV sequence cells use brackets and
semicolons; half-integers are decimal.

Internal timing includes scheduler construction, preparation and ordered
collection but excludes printing. `preparation_seconds` reports the one-time
root work in prepared mode. Compare direct and prepared modes at zero/one/four/
eight workers to separate reuse from concurrency. Once roots are cached, queries
may be too cheap to benefit from parallel execution; measurements belong in
[parallel grids](parallel-grids.md#prepared-ladder-sweeps).

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

## Independent Richardson and central-spin solves

Their synchronous solvers already support explicit MPFR precision, including
scalar-only MPFR configurations without MPLAPACK. Submit independent parameter
or sector points through `evaluate_observable_grid_async`; each point's Newton
iterations and continuation stages stay sequential. A new scheduler or a task
for each continuation step is unnecessary.

The [complete parameter-sweep example](../examples/async_gaudin_sweep.cpp)
captures the level and bath-coupling vectors and both sets of solver options by
value. It constructs the solvers' borrowed spans only inside each evaluation.
Capturing a span by value alone would **not** extend the input array's lifetime.
The bounded collector retains complete states in input order, and the host
waits once on the final point-count result before inspecting or printing them.
The example uses 160-bit MPFR when configured and binary64 otherwise, and reports
the working precision. It is built in CI and is not an installed frontend.

```sh
cmake -DBETHE_BUILD_EXAMPLES=ON /path/to/bethe
cmake --build . --target bethe-async-gaudin-sweep-example -j2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./examples/bethe-async-gaudin-sweep-example
```

Numerical failure remains a state with its original status, reached parameter,
residuals and iteration counts; inspect `converged` before using an energy.
A zero-budget central-spin solve can have no physical energy at all. A grid's
successful completion means every requested evaluation finished, not that every
state converged. Execution exceptions and cancellation follow the shared grid
contract above. This example demonstrates composition, not a measured speedup.

Regression checks compare full states with serial solves, including partial
continuation results, blocked pairing levels and signed bath couplings. Deferred
scheduling exercises destruction of submission-local inputs, while one-worker
and multi-worker runs check composition. MPFR tests interleave 128-, 256- and
333-bit requests and check stored precision without changing ambient defaults.
Existing XXX scan/form-factor tests already cover independent-precision
composition; those APIs require no additional wrapper.

## Finite-volume integral equations

Lee–Yang `ground_state`/`one_particle` and sine-Gordon
`vacuum_energy`/`two_soliton_level` have no shared mutable solver state.
Independent scales, circumferences and allowed levels can use the existing
observable grid. Capture the inputs and numerical options by value, return the
complete state, and wait on the grid's final result before reading collected
states. Coupled fixed-point, hole/root, mesh and contour iterations remain
inside each sequential solve.

The [complete example](../examples/async_finite_volume.cpp) runs any of these
four families. Configure with `BETHE_BUILD_EXAMPLES=ON`, then from the build
directory:

```sh
cmake --build . --target bethe-async-finite-volume-example -j2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./examples/bethe-async-finite-volume-example sine-gordon-excited 8 64 fp64
```

Arguments are model/state, threads, number of points and precision. Model/state
is `lee-yang-vacuum`, `lee-yang-excited`, `sine-gordon-vacuum` or
`sine-gordon-excited`; threads zero is a direct serial comparison, while one or
more uses Async. Native `fp64`, `long-double` and configured `fp128` are
supported. These solvers do not yet support MPFR.

The example uses unit mass, Lee–Yang circumferences `5+5*i/points` at absolute
Y tolerance `1e-10`, or sine-Gordon `1+3*i/points` at p=2 and tolerance `1e-7`.
The excited sine-Gordon state has positive Bethe number 1/2. The grid admits at
most 32 active or uncollected states; this illustrative collector retains all
results. Numerical failures remain in order and cause exit 2 after the complete
sweep. Timing includes solving and ordered collection, excluding scheduler
construction and printing. More arithmetic bits do not improve the deliberately
fixed discretization tolerance automatically.

The two excited-state frontends also offer `--threads 2` for their independent
vacuum/excited pair. Their single-scale vacuum frontends stay serial. Useful
concurrency is across complete solves, not across mutually dependent iterations;
a parameter sweep can scale beyond the pair's two-task ceiling.
