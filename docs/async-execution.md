# Async execution

Bethe uses Uni20's scheduler for independent calculations. It creates no thread
pool of its own and never changes the global scheduler. Existing synchronous
scans wait for the same operation exposed by the Async API.

For C++ build setup, start with [embedding Bethe](library.md#embed-the-library).

See [app coverage](app-capabilities.md) for the execution and arbitrary-precision
paths exposed by each frontend and their numerical limits.

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
their solves use normal Uni20 dispatch, including CPU MPFR LU.

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

## Buildable examples

See [complete composition examples](async-examples.md) for shared Lieb–Liniger,
Hubbard and ladder backgrounds, independent ground-state/nested/complex-root
sweeps, pairing/central-spin queries and finite-volume equations. Their build
commands, supported precisions and output-memory behavior are recorded there.
[Measured benchmarks](benchmarks.md) keep timing evidence separate from contracts.

<a id="ground-state-parameter-sweeps"></a>
[Ground-state parameter sweeps](async-examples.md#ground-state-parameter-sweeps).

<a id="independent-tasep-and-finite-xyz-solves"></a>
[Independent TASEP and finite XYZ solves](async-examples.md#independent-tasep-and-finite-xyz-solves).

<a id="owning-nested-model-requests"></a>
[Owning nested-model requests](async-examples.md#owning-nested-model-requests).

<a id="prepared-ladder-field-and-rung-sweeps"></a>
[Prepared ladder field and rung sweeps](async-examples.md#prepared-ladder-field-and-rung-sweeps).

<a id="matched-pipeline-benchmark"></a>
[Matched pipeline benchmark](benchmarks.md#matched-pipeline-benchmark).

<a id="independent-richardson-and-central-spin-solves"></a>
[Independent Richardson and central-spin solves](async-examples.md#independent-richardson-and-central-spin-solves).

<a id="finite-volume-integral-equations"></a>
[Finite-volume integral equations](async-examples.md#finite-volume-integral-equations).
