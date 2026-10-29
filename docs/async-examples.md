# C++ composition examples

[Library setup](library.md) · [Async contracts](async-execution.md) · [Benchmarks](benchmarks.md)

These are complete buildable examples under `examples/`, enabled with
`BETHE_BUILD_EXAMPLES=ON`. They are checked in CI and are not installed frontend
commands. Each section states its precision choices; several demonstrations
intentionally retain native-only selectors even where the model API supports MPFR.
All wait on final operation results before inspecting caller-owned output.

## Shared-background examples

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

It uses explicit 160-bit `mpreal` when MPFR scalars are configured, otherwise
binary64. The printed header identifies the working precision. It waits once
on the grid result, then checks scientific convergence before reporting energies.

Thermodynamic Lieb–Liniger also supports `mpreal` in the shared-background
[composition](async-execution.md#shared-backgrounds). Construct `Options<uni20::mpreal>` with a finite
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
grid. This example selects only fixed native precisions; the model APIs also support
explicit MPFR. The existing single-request frontends remain serial. TASEP root updates
and XYZ Newton steps are coupled and remain sequential. For a TASEP rate-only
curve, reuse a unit-rate result and rescale its eigenvalue/gap; the roots do not
depend on the positive jump rate.

## Owning nested-model requests

The [nested sweep example](../examples/async_nested_sweeps.cpp) composes
independent TB, SU(N) fermion and t–J solves through the ordinary owning grid.
Each request retains its controls and inputs, including component-population
vectors. Its Newton iterations and coupling continuation stay sequential.
The example changes no state family. Its selector uses native precisions; the
model APIs also support explicit MPFR.

From a build directory configured with `BETHE_BUILD_EXAMPLES=ON`:

```sh
cmake --build . --target bethe-async-nested-sweeps-example -j2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./examples/bethe-async-nested-sweeps-example tb 32 33 4 fp64
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./examples/bethe-async-nested-sweeps-example sun 5 129 4 long-double
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./examples/bethe-async-nested-sweeps-example tj 64 16 4 fp64
```

Arguments are `MODEL SIZE POINTS THREADS PRECISION`. Their physical meaning is:

| Model | SIZE | Sweep |
| --- | --- | --- |
| `tb` | First even site count, at least 4 | Sizes `SIZE, SIZE+2, ...`, last at most 10,000 |
| `sun` | Odd positive population of each of three components | Coupling c from 1/2 to 4 at circumference `3*SIZE`, hence total density one; one point uses c=1 |
| `tj` | Site count, at least 3 | Balanced odd spin populations 1,3,...,`2*POINTS-1`, requiring their sum not to exceed SIZE |

Zero threads uses direct scalar calls; a positive count uses Async, capped at
the smaller of 32 and the number of points. Precisions are `fp64`, `long-double`
and configured `fp128`. These are fixed-precision native paths. The single-request
frontends remain serial.

The example waits on the final grid result and prints complete states afterward.
Collection is ordered on the scheduler; at most 32 states are active or awaiting
collection. Collected output is retained, so its memory grows with the sweep.
CSV sequence cells use brackets/semicolons and complex entries use `real:imag`;
half-integer quantum labels are decimal. The example returns 2 if any numerical
solve fails. Failed TB/t–J rows retain their final-iterate estimates. A partial
SU(N) result retains the **last reached coupling** and its energy; it is not an
energy at the requested coupling. An absent energy remains an empty cell.

Tests compose all three models on one scheduler, including deferred execution
after submission-local vectors disappear and one/four-worker scheduling. They
compare every field against direct calls, including free/XXX limits, population
permutations, spin reversal, numerical budgets and partial continuation. An
unsupported SU(N) shell fails its grid without cancelling the other two.
See [parallel-grid measurements](benchmark-state-scans.md#nested-model-parameter-sweeps).

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
[parallel grids](benchmark-preparation.md#prepared-ladder-sweeps).

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
[contract](async-execution.md#ownership-dependencies-and-cancellation). This example demonstrates composition, not a measured speedup.

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
supported by this example. The Lee–Yang and sine-Gordon vacuum and excited library APIs
also support explicit MPFR controls and owning mixed-precision queries through
`evaluate_observable_grid_async`; the example's precision selector remains
native-only.

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
