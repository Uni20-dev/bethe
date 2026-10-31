# Four-spinon spectra and direct moments

`bethe/xxx_four_spinon_thermodynamics.hpp` combines physical W4,
[support-adapted geometry](xxx-four-spinon-integration.md), bounded density
recovery and [ordered quadrature accumulation](quadrature-accumulation.md).
It returns the longitudinal response and moments in their physical units.
There is no caller-side normalization correction.

## Choose the observable and domain

`ThermodynamicFourSpinon<Geometry, Low, High>` prepares immutable density
contexts. `Geometry` also supplies the accumulation precision and must retain
the largest density guard precision. For MPFR these three template arguments
can be the same type with separate explicit precision contexts.

| Query | Returned components | Finite domain |
| --- | --- | --- |
| `spectrum(k, omega, order, T)` | S4zz(k,omega) | Log-radius cutoff T |
| `fixed_moments(k, order, T)` | Zeroth and first fixed-momentum moments | Log-radius cutoff T |
| `global_moments(order, gaps)` | Zeroth and first global moments | Three individual gap cutoffs; unrestricted common shift |

Every query accepts an additional maximum node count, defaulting to 100000.
Spectrum and fixed-momentum queries prepare a bounded set of geometry nodes.
Global moments construct nodes by index without storing the three-dimensional
grid. Geometry preparation precedes evaluation; the returned evaluation wall
time excludes that preparation.

Moment measures and sum rules are included internally:

```math
m_r(k)=\int_0^\infty\frac{d\omega}{2\pi}\,\omega^rS_4^{zz}(k,\omega),\qquad
I_r=\int_0^{2\pi}\frac{dk}{2\pi}\,m_r(k).
```

```math
I_0^{\rm exact}=\frac14,\qquad
I_1^{\rm exact}=\frac{2J}{3}\left(\log 2-\frac14\right),\qquad
m_1^{\rm exact}(k)=I_1^{\rm exact}(1-\cos k).
```

Moment results include sum-rule denominators and fractions where defined.
No exact fixed-momentum zeroth-moment denominator is assumed. The first-moment
ratio at k=0 is undefined and is left empty. The k=pi static response diverges
without a cutoff: fixed-momentum results there carry
`unrestricted_static_divergence`, even when their finite-cutoff estimates exist.

For exchange J, spectrum values and diagnostics scale by 1/J at frequency
omega/J. Moment order r scales by J to the power r. The sum-rule denominators
scale consistently, so fractions are dimensionless. These are partial
four-spinon responses, not the complete structure factor.

## A small fresh calculation

With MPFR enabled, the following calculates two global moments using two Gauss
nodes in each of three gap intervals. It is deliberately a tiny finite-domain
example, not a thermodynamic estimate:

```cpp
#include <bethe/xxx_four_spinon_thermodynamics.hpp>
#include <uni20/async/tbb_scheduler.hpp>

using namespace bethe::heisenberg;
using R = uni20::mpreal;
R one(1, uni20::Precision::bits(192));
SpinonRecoveryPolicy recovery{
    "1e-6", "1e-30", {{128, 192, {{{64, 16}, {80, 20}}}}}, 60};
ThermodynamicFourSpinon<R, R> model(recovery, one);
uni20::async::TbbScheduler scheduler(2);
auto pending = model.global_moments(2, {2*one, 2*one, 2*one})
    .evaluate_async({.scheduler = &scheduler, .max_in_flight = 2});
auto result = pending.get_wait(scheduler);
// Check result.estimates before using either component.
```

This particular finite cubature gives approximately 0.00014958616511807964
and 0.00038003614670975796. They are reproduced independently by the 8F7
matrix-element formula and direct positive common-shift integration at 30 and
45 decimal digits. The omitted gap region is large. These numbers check the
integration convention and implementation, not unrestricted spectral fractions.

The complete [composition example](../examples/xxx_four_spinon_moments.cpp)
adds order and cutoff comparisons and can be built with `BETHE_BUILD_EXAMPLES=ON`.
It waits on the final Async result and requires no archive or Python process.

## Keep numerical checks separate

Each component reports:

- Its physical value.
- The weighted contour-tail diagnostic.
- The weighted low/guard arithmetic screen, including applicable common-shift
  arithmetic indicators.
- Common-shift quadrature differences and omitted shift-tail bounds, when the
  zeroth-moment reduction uses its numerical fallback.

These fields exclude the outer quadrature error and omitted domain. Kernel
recovery evaluates the same geometry coordinates at each arithmetic precision;
it does not certify the geometry itself. Repeat the complete geometry at higher
precision to test that separately.

`refined_spectrum`, `refined_fixed_moments` and `refined_global_moments` evaluate
four complete controls: coarse/fine order on a base domain, then coarse/fine on
an extended domain. For example:

```cpp
auto comparison = model.refined_global_moments(
    2, 3, {one, one, one}, {2*one, 2*one, 2*one}).evaluate();
```

The result retains all four controls, their diagnostics, order differences on
both domains, and the signed change between the two fine-order domain estimates.
Only one geometry is retained at a time, and each control uses the same bounded
Async admission. A failed control prevents a refined estimate. Differences are
numerical evidence; they are not certified errors, an extrapolation, or a proof
that the remaining tails are negligible.

The log-radius and gap cutoffs describe different domains. In particular, the
analytical common-shift reduction already integrates the whole shift direction.
Comparing spectra with global moments requires accounting for this distinction.

## Failures and validation

Query construction validates coordinates, precision and refinement controls
before any work is submitted. Refined global moments also check the finest
order on the widest gap domain against the node and span budgets. This
validation does not build the full quadrature grid; geometry remains lazy.
Command-line callers can therefore reject invalid refined requests before
opening their output files.

Resolved exclusions outside the continuum, and the zero-momentum inelastic
response, return `physical_zero`. Thresholds are identified separately and have
no finite point estimate. Unresolved geometry, collisions, exhausted density
recovery and lost exponent range are numerical failures. No failed node is
discarded or replaced with zero. Cancellation drains admitted work and returns
no integral over the surviving prefix.

Results retain precision, work budgets, cutoff definitions, recovery settings,
completed-prefix counts, failure locations and evaluation wall time. Density
attempt counts on failure describe collected observations, not discarded work
finished during draining. The query and originating model can be destroyed
after submission; the caller still owns the scheduler through completion.

Tests cover the independent finite spectral slice, direct-global references,
whole-geometry precision growth, exchange scaling, sum-rule fractions, serial
versus Async agreement, a single-worker refinement chain and explicit failures.
The generic accumulator additionally checks completion while another admitted
observation is deliberately held open. Production command-line integration and
a convergence tutorial are the next release checkpoint.
