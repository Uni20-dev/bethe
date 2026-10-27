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

The periodic and open XXX APIs provide `real_excitations_async`. Other models
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

## Ownership, dependencies and cancellation

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
