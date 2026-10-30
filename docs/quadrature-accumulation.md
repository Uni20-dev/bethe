# Bounded quadrature accumulation

`bethe/quadrature_accumulation.hpp` reduces an indexed set of weighted numerical
observations. It reuses Bethe's bounded [observable grid](async-execution.md)
and preserves deterministic index order. Domain maps, quadrature weights,
importance transformations, physical kernels and refinement choices belong to
the caller. The engine neither rescales a physical result nor estimates its
integration error.

An observation contains an optional fixed-size array of scalar contributions,
optional adapter-specific work counts and failure metadata, a cancellation flag,
and a reason. A successful array can contain physical values and separately
weighted diagnostics. A missing array is a numerical failure. Cancellation is
terminal even if an observation carries a final value.

```cpp
#include <bethe/detail/gauss_legendre.hpp>
#include <bethe/quadrature_accumulation.hpp>
#include <uni20/async/tbb_scheduler.hpp>

auto rule = bethe::detail::gauss_legendre<double>(16);
uni20::async::TbbScheduler scheduler(2);
auto pending = bethe::accumulate_quadrature_async<double, 2>(
    rule.x.size(),
    [rule](std::size_t i, std::stop_token stop) {
        using Observation = bethe::QuadratureObservation<double, 2>;
        if (stop.stop_requested()) return Observation{.cancelled = true};
        double x = (1 + rule.x[i])/2, w = rule.w[i]/2;
        return Observation{std::array{x*x*w, w}};
    }, 1.0, {.scheduler = &scheduler, .max_in_flight = 2});
auto result = pending.get_wait(scheduler);
// A complete result sums to {1/3, 1}, up to numerical roundoff.
```

The exemplar fixes scalar precision, including an explicit finite MPFR context.
Each returned component must have that precision. Callables and numerical state
are owned through completion; references captured by a callable remain borrowed.
The supplied scheduler must outlive all admitted work.

`max_in_flight` bounds active and completed-but-uncollected observations. Results
are accumulated with compensated summation in index order. A failed observation
stops admission, requests cooperative cancellation, and drains accepted work.
Unexpected exceptions also propagate after draining. There is no extra thread
pool and no `run_all()` call.

A reduction contains sums only after complete success. On failure or cancellation
it retains the collected-prefix count and failure information, without presenting
that prefix as an integral. Work counts describe collected observations, including
a collected failure; work discarded during draining is not included. A zero-node
reduction is the empty sum, distinct from a caller's physical support decision.

The [four-spinon observable adapter](xxx-four-spinon-thermodynamics.md) uses this
engine for values, contour and arithmetic diagnostics, and common-shift errors.
Its order and cutoff comparisons are separate from the accumulation mechanism.
