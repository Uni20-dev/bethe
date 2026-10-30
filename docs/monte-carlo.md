# Monte Carlo execution and physical sampling strategies

`bethe/monte_carlo.hpp` runs independent Monte Carlo replicates over a callable
integrand. The [six-spinon spectral calculation](xxx-six-spinon-structure-factor.md)
is one consumer. The numerical engine manages sample identities, bounded Async
execution, accumulation, covariance and failure records. The physical integrand
owns its proposal distribution, coordinate maps, Jacobians, domain and kernel.

Importance sampling is part of the calculation's design. The engine does not
choose a distribution, optimize it silently, or infer a normalization. A proposal
must expose enough information to evaluate its complete probability density,
including overlapping channels. The spinon integrands use this information for
soft-pair sampling and balanced coarea weights.

## Weighted observations

An integrand is a const-callable object accepting `UniformSample<D> const&` and
a `std::stop_token`. It returns `MonteCarloObservation<Real, Metadata>`:

- `values` contains a vector of **already weighted** contributions, or is absent
  on failure. The engine averages these values without another measure factor.
- `metadata` carries optional model-specific work counters and failure details.
  Its `operator+=` defines their aggregation and must respect the declared work
  budgets. The default metadata type is empty.
- `cancelled` and `reason` identify a terminal cancellation or numerical failure.
  Cancellation remains terminal even if values were supplied.

A value vector can include several observables and their numerical diagnostics.
The engine computes the full covariance of its replicate means. The adapter
still labels which entries are physical observables, arithmetic indicators or
other diagnostics; covariance does not turn them into a combined error bound.

For example, a 50:50 mixture of uniform sampling and density `2*x` on (0,1)
has complete density `q(x)=1/2+x`. This integrand estimates integrals of x squared
and x cubed with the same draws:

```cpp
#include <bethe/monte_carlo.hpp>
#include <uni20/async/tbb_scheduler.hpp>

int main()
{
    auto observe = [](bethe::UniformSample<2> const& identity, std::stop_token) {
        auto u = identity.at_precision(1.0);
        double x = u[0] < 0.5 ? u[1] : std::sqrt(u[1]);
        double q = 0.5 + x; // both proposal channels contribute
        return bethe::MonteCarloObservation<double>{
            std::vector<double>{x*x/q, x*x*x/q}};
    };
    uni20::async::TbbScheduler scheduler(2);
    auto pending = bethe::monte_carlo_async<2>(
        {.replicates=8, .samples_per_replicate=256, .seed=108},
        2, 1.0, observe, {.scheduler=&scheduler, .max_in_flight=2});
    auto result = pending.get_wait(scheduler);
    // Complete results estimate (1/3, 1/4), with a 2x2 covariance of the mean.
    return result.status == bethe::MonteCarloStatus::complete ? 0 : 1;
}
```

The two uniform coordinates select the proposal channel and its position.
Changing to the density of only the selected channel would define a different,
incorrect weighting for this mixture. The supplied engine is ordinary Monte
Carlo, not a quasi-Monte Carlo or adaptive quadrature algorithm.

## Identities, precision and execution

`bethe/sampling.hpp` contains `UniformSample`, `MonteCarloStream` and
`ReplicateStatistics`. A stream is keyed by seed, replicate id and uniform-bit
count. Its integer words represent exact open binary midpoints at every
sufficient precision. A numerical retry retains these identities. It never
redraws a difficult observation. Increasing the uniform-bit count changes the
input grid and is a separate control from increasing arithmetic precision.

`monte_carlo_replicate<D>` evaluates one replicate. `monte_carlo_async<D>` owns
the integrand and schedules a bounded number of replicates using Bethe's shared
observable-grid scheduler. Callers may use `.get_wait(scheduler)` or await the
result in another coroutine. The borrowed scheduler must outlive completion,
including work whose handle was dropped. Cancellation drains admitted work;
dropping a handle does not request cancellation. Callables must support
concurrent const calls without shared mutable numerical scratch.

All arithmetic remains in the exemplar's scalar type and explicit precision.
The engine uses the existing Neumaier accumulator for observations and ordered
replicate reduction for reproducibility. There is no conversion to binary64 for
an integral, variance or diagnostic. Wall time uses native seconds.

## Assessing a proposal

Results retain every available replicate record, its completed sample count,
failed sample index, metadata, mean and largest absolute contribution in each
component. These support inspection of support-hit rates, between-replicate
scatter and unusually influential samples. They do not prove that an unseen
rare region is negligible. Compare independent samples and alternative physical
proposals, and vary cutoffs separately from arithmetic controls.

Complete equal-size replicates supply a mean and, with at least two replicates,
its covariance estimate. Failed or missing replicates cannot form a mean over
survivors. Nonfinite contributions or unrepresentable means/variances fail
explicitly. An integrand may return a genuine support zero; the physical adapter
must distinguish a known zero from a pilot that happened to miss support.

The engine does not certify omitted tails, proposal-grid resolution, contour
convergence or geometry precision. Nor does it supply a quadrature error estimate:
deterministic integration uses its own refinement comparisons even when sharing
the same scheduling and physical kernels.
