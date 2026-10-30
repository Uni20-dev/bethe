# Bounded spinon-density recovery

`bethe/xxx_spinon_density_recovery.hpp` checks a physical spinon density before
an integrator accepts it. It reuses [the common W2/W4/W6 evaluator](xxx-spinon-density.md);
it introduces no normalization or integration measure.

Increasing arithmetic precision alone does not establish contour convergence.
Each recovery stage therefore evaluates **two increasing contour budgets**, each
at a lower and a higher arithmetic precision. Acceptance requires the contour
and arithmetic screens to pass and the two results to agree. Earlier individually
passing controls also remain part of that comparison: a later successful call
cannot erase a contradictory result.

## Prepared policy

The policy declares every allowed stage in advance. Contexts are prepared once
and support concurrent const calls; each call owns its history and scratch.
For example, in a build with MPFR:

```cpp
#include <bethe/xxx_spinon_density_recovery.hpp>

using namespace bethe::heisenberg;
SpinonRecoveryPolicy policy{
    "1e-8", "0", // relative and absolute density tolerances
    {{120, 190, {{{32, 12}, {64, 20}}}},
     {170, 270, {{{64, 20}, {128, 24}}}}},
    60           // wall seconds per point, checked between kernels
};
XXXSpinonDensityRecovery<uni20::mpreal> recovery(policy);
std::vector<std::string> beta{"-2.1", "-1.3", "-.4", ".3", "1.2", "2"};
auto result = recovery(beta);
if (!result.accepted()) {
    // Retain result.status, result.reason and the attempt history.
    // The enclosing integral is incomplete unless this SAME point is recovered.
} else {
    auto const& density = *result.estimate;
    // density.value is physical W6. Its contour_tail_error and the checked
    // *result.arithmetic_screen remain separate from sampling/cutoff errors.
}
```

A stage is `{low_bits, guard_bits, {{terms, order}, {terms, order}}}`.
Both contour budgets and successive stages must strengthen their controls.
There are at most 64 stages. `XXXSpinonDensityRecovery<double, long double>`
also works where long double has more significand bits than double; its policy
must specify the actual native precisions. Repeating one native type twice is
rejected as a precision check. The same template supports other distinct native
pairs and runtime MPFR stages.

The returned value has the accepted stage's guard precision. An integrator must
explicitly import it into its chosen accumulation context; values from different
MPFR stages must not acquire a hidden common binary64 representation.

## Point identity and terminal results

Decimal strings are parsed afresh at every precision. Exact decimal coincidences
return `unsupported_limit` before evaluation; a collision caused by rounding
distinct inputs can proceed to a higher declared precision. Invalid decimal
syntax returns `invalid_input`.

For generated coordinates, `recover(source, stop_token)` takes a callable
`source(context)` returning the rapidities at that context's precision. It must
reconstruct the same point from retained input identities. For example, an
energy-surface sampler retains its uniform integer words and reconstructs the
geometry; it must never draw another sample because one was difficult. The
callback may run more than once and should have no observable mutable state.

Possible terminal statuses are `accepted`, `exhausted`, `cancelled`, `time_limit`,
`invalid_input`, `unsupported_limit` and `inconsistent_controls`. Failed results
have no accepted estimate. A contour-domain, rounded-collision, arithmetic,
exponent-range or special-function-budget failure is recorded and may advance
through the declared stages. Unexpected exceptions propagate to the caller.

Cancellation and elapsed-time checks occur between numerical kernels. They do
not forcibly interrupt a running special function. Scheduling belongs to the
outer bounded Async calculation; the recovery layer starts no threads, changes
no scheduler and introduces no batch synchronization.

## Meaning of the diagnostics

Each attempt retains the guard estimate, its absolute low/guard difference, the
change in the scalar pair factor, the combined arithmetic screen and the requested
absolute tolerance. The arithmetic screen includes both contexts' diagnostics
and the observed precision differences. This is a conservative numerical screen,
**not an interval certificate**. Contour majorants are also evaluated in ordinary
floating point.

These checks do not establish convergence of a surrounding quadrature or Monte
Carlo estimate, nor of its rapidity cutoff. In particular, a complete density
recovery policy cannot repair geometry whose endpoint distances were already
lost before the rapidities reached the evaluator.

The focused tests exercise native precision pairs, MPFR recovery of nearly
coincident rapidities, exhausted contour and arithmetic budgets, exact input
identities, cancellation and conflicting controls. The independent four-spinon
hypergeometric representation checks the physical value and precision growth.
