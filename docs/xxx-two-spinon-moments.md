# Direct two-spinon moments

[Two-spinon spectrum](xxx-structure-factor-thermo.md) · [Selected sectors](xxx-spinon-sectors.md)

`bethe/xxx_two_spinon_moments.hpp` integrates the two-spinon static structure
factor and first frequency moment directly. It supplies a reference alongside
the four- and six-spinon moment calculations. No displayed frequency grid or
sum-rule rescaling enters the calculation.

## Resolve the spectral edges with the rapidity gap

Let `g` be the positive difference between the two rapidities, `rho=g/(2*pi)`,
and `F(g)=exp(-I(rho))/cosh(g/2)^2`, using the same transition-rate function
as the spectral evaluator. Fold momentum into `0<q<=pi`, and write
`s=sin(q/2)`, `c=cos(q/2)`. Energy on the two-spinon surface is

```math
\omega(g)=\pi J s\sqrt{c^2+s^2/\cosh^2(g/2)}.
```

The upper edge is `g=0`, and the lower edge is approached as `g` tends to
infinity. Changing variables from frequency to this gap absorbs the spectral
edge Jacobian. The longitudinal moments are

```math
m_{0,2}(q)=\int_0^G dg\,F(g)\frac{J s^2}{8\omega(g)},
\qquad
m_{1,2}(q)=\int_0^G dg\,F(g)\frac{J s^2}{8}.
```

The frequency measure `domega/(2*pi)` is already included. The chosen finite
gap cutoff `G` omits the lower-edge tail; it is neither a frequency cutoff
shared between momenta nor an absolute-rapidity box. The common shift remains
unrestricted. At `q=pi` the zeroth moment diverges as `G` grows, so a finite
answer there is explicitly flagged as a truncated estimate.

Integrating the momentum dependence analytically gives the global moments:

```math
I_{0,2}=\int_0^G dg\,F(g)\frac{g}{8\pi^2\tanh(g/2)},
\qquad
I_{1,2}=\int_0^G dg\,F(g)\frac{J}{16}.
```

These coefficients follow from the conservation-law measure, rather than
inserting a known spinon fraction. In particular the first moment's momentum
dependence separates, while its coefficient is still calculated numerically.
The tests independently integrate the spectrum in frequency and the fixed
moments in momentum to check these reductions and their normalization.

## C++ interface and diagnostics

```cpp
#include <bethe/xxx_two_spinon_moments.hpp>
#include <uni20/async/tbb_scheduler.hpp>

using R = double;
bethe::heisenberg::StructureFactorOptions<R> options;
bethe::heisenberg::ThermodynamicTwoSpinonMoments<R> model(1.0, options);
uni20::async::TbbScheduler scheduler(2);
auto pending = model.fixed_moments(1.5, 16, 48.0).evaluate_async(
    {.scheduler = &scheduler, .max_in_flight = 4});
auto result = pending.get_wait(scheduler);
// global_moments(16, 48.0) uses the same owned query interface.
```

The arguments are momentum (for fixed moments), Gauss order per gap panel,
gap cutoff, and an optional total node budget. Panels have width at most two.
Arithmetic precision comes from the supplied scalar/options context; MPFR
requires an explicit finite precision throughout. The owned query survives
destruction of its model, while the scheduler must survive completion.

On success, `result.integral.sums` contains m0, m1 and their two integrated
transition-kernel diagnostics. Global results contain I0 and I1 instead.
The diagnostics exclude gap quadrature and the omitted large-gap region:
vary order and cutoff independently before interpreting the result as
converged. Failed kernels, exponent-range loss and cancellation produce no
partial-prefix estimate.

`sum_rule_denominators` and `sum_rule_fractions` contain the global references
`I0=1/4` and `I1=(2J/3)(log(2)-1/4)`, or the known fixed-q first moment.
There is no invented fixed-q static denominator. At zero momentum both
inelastic moments vanish without kernel evaluations, and the first-moment
fraction is undefined. Channel conversion belongs to the common frontend;
these library values are longitudinal throughout.
