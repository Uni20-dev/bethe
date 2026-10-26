# Async execution

Bethe uses Uni20's scheduler for independent calculations. It creates no thread
pool of its own and never changes the global scheduler. Existing synchronous
scans wait for the same operation exposed by the Async API.

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
