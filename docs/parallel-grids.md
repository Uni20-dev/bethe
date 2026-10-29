# Parallel observable grids

An observable grid evaluates a scalar quantity at independent logical points.
The model owns the physics and prepares its background; the shared
`bethe::evaluate_observable_grid` adapter owns bounded admission and ordered collection.
Uni20 supplies the scheduler. No output format or model formula belongs in the
adapter.

Hubbard and Lieb–Liniger dispersion use this path.
See the [command-line examples](parallel-calculations.md#parallel-observable-grids).
Thermodynamic XXX spectral grids also use this path, with ordered spectrum rows
and final failure/threshold counts; see [the frontend](xxx-structure-factor-thermo.md).
[Thermal Lieb–Liniger](lieb-liniger-thermal.md) also uses this adapter for
independent temperature solves in either ensemble, preserving each point's
iteration and convergence diagnostics.
[Hubbard continuum edges](hubbard-continuum.md) also use this adapter for
independent total momenta, retaining each point's inversion/search diagnostics.
Sine-Gordon Bethe–Yang levels also use this adapter: each level owns an
independent root search and repeated phase integrals.
Other observables can adopt the same path after checking their independence and scratch ownership.

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
Frontends keep branches in their declared output order; a row's identity is
its branch and logical point index. Hubbard dispersion flattens its selected
branches into one logical grid, allowing evaluation across branch boundaries
without changing row order. Lieb–Liniger currently evaluates each branch's grid
in turn. In both cases the admission window bounds active and uncollected results;
collected output can still grow with the grid.

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

## Performance measurements

[Benchmarks by workload](benchmarks.md) preserve the dated timings, commands,
precision/provider configurations and whole-result comparisons.

<a id="measuring-scaling"></a>
[Measuring scaling](benchmarks.md#measuring-scaling).

<a id="initial-measurement"></a>
[Initial measurement](benchmark-grids.md#initial-measurement).

<a id="hubbard-branch-scheduling"></a>
[Hubbard branch scheduling](benchmark-grids.md#hubbard-branch-scheduling).

<a id="thermodynamic-xxx-streaming-measurement"></a>
[Thermodynamic XXX streaming measurement](benchmark-grids.md#thermodynamic-xxx-streaming-measurement).

<a id="thermal-liebliniger-temperature-scans"></a>
[Thermal Lieb–Liniger temperature scans](benchmark-grids.md#thermal-liebliniger-temperature-scans).

<a id="hubbard-continuum-frontend"></a>
[Hubbard continuum frontend](benchmark-grids.md#hubbard-continuum-frontend).

<a id="xyz-prepared-bound-branches"></a>
[XYZ prepared bound branches](benchmark-preparation.md#xyz-prepared-bound-branches).

<a id="sine-gordon-betheyang-level-scans"></a>
[Sine-Gordon Bethe–Yang level scans](benchmark-field-theory.md#sine-gordon-betheyang-level-scans).

<a id="finite-volume-integral-equation-sweeps"></a>
[Finite-volume integral-equation sweeps](benchmark-field-theory.md#finite-volume-integral-equation-sweeps).

<a id="overlapping-one-vacuumexcited-pair"></a>
[Overlapping one vacuum/excited pair](benchmark-field-theory.md#overlapping-one-vacuumexcited-pair).

<a id="kondo-field-response-scans"></a>
[Kondo field-response scans](benchmark-state-scans.md#kondo-field-response-scans).

<a id="potts-momentum-scans"></a>
[Potts momentum scans](benchmark-state-scans.md#potts-momentum-scans).

<a id="whole-process-measurements"></a>
[Whole-process measurements](benchmark-state-scans.md#whole-process-measurements).

<a id="asep-shared-seeds"></a>
[ASEP shared seeds](benchmark-preparation.md#asep-shared-seeds).

<a id="whole-process-measurements_1"></a>
<a id="whole-process-measurements-1"></a>
[ASEP whole-process measurements](benchmark-preparation.md#whole-process-measurements).

<a id="tasep-and-finite-xyz-parameter-sweeps"></a>
[TASEP and finite XYZ parameter sweeps](benchmark-state-scans.md#tasep-and-finite-xyz-parameter-sweeps).

<a id="measurement-on-28-october-2026"></a>
[Measurement on 28 October 2026](benchmark-state-scans.md#measurement-on-28-october-2026).

<a id="nested-model-parameter-sweeps"></a>
[Nested-model parameter sweeps](benchmark-state-scans.md#nested-model-parameter-sweeps).

<a id="prepared-ladder-sweeps"></a>
[Prepared ladder sweeps](benchmark-preparation.md#prepared-ladder-sweeps).

<a id="biquadratic-sector-and-bound-cluster-scans"></a>
[Biquadratic sector and bound-cluster scans](benchmark-state-scans.md#biquadratic-sector-and-bound-cluster-scans).

<a id="validation-and-measurement"></a>
[Validation and measurement](benchmark-state-scans.md#validation-and-measurement).
