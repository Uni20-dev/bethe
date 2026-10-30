# Four-spinon integration geometry

[Spinon densities](xxx-spinon-density.md) · [C++ library](library.md)

A four-spinon spectrum requires two-dimensional integration after imposing
momentum and energy conservation. Direct fixed-momentum moments require three
dimensions. Both use the ordered momentum patches in
`bethe/detail/xxx_four_spinon_geometry.hpp`.

This is an internal integration building block, not a new command-line spectrum
or a convergence controller. It supplies coordinates and measures to a callback;
the physical weight comes from `XXXSpinonDensity`. Kernel failures propagate.

## Integrate the allowed region

Order the four momenta inside `(0,pi)`. Let `a=r*v` and `b=r*(1-v)` be the two
distances closest to an endpoint. They are either both near zero or at opposite
ends. The remaining ordered pair is `h,H-h`:

```math
(a,b,h,H-h),\quad H=Q-r;
\qquad (a,h,H-h,\pi-b),\quad H=Q-\pi+r(1-2v).
```

Ordering gives an interval for `h` and an upper bound on `r`. Support must be
located in **both** remaining coordinates before quadrature. Close to an upper
edge, the allowed fraction `v` can be so narrow that even exact radial bounds
would leave every full-interval fraction node outside support.

First evaluate the maximum energy at the radial endpoint as a function of `v`.
For the same-end patch it increases towards `v=1/2`. The opposite-end patch has
a single maximum. Writing `A=Q-pi`, its interior stationary point obeys
`cos(A*v)=3*cos(A*(1-v))`; when this point lies outside the interval, the maximum
is at `v=0`. Bisection on each side of the maximum gives the allowed fraction
interval, which is then split at changes of the radial bound's formula.

At each fraction node, both radial endpoint energies increase with `r`.
Inverting these bounds gives the occupied radial interval directly. The
quadrature therefore never relies on finding support by sampling zero padding.
Unresolved support boundaries throw instead of becoming outside-support zeros;
points separated from the upper edge by more than the roundoff margin can be
classified as outside support.

At `J=1`, let `e_s` be the energy of the endpoint momenta. Eliminating the hard
pair gives the ordered-root spectral measure

```math
d\mu=\frac{r\,dr\,dv}
{\sqrt{[\pi\sin(H/2)]^2-(\omega-e_s)^2}}.
```

Each emitted `log_measure` includes `4! = 24`, accounting for the labelled
permutations, including the exchanged hard pair. W4 already contains the
reciprocal factorial. No further normalization belongs in the callback.

Both total-momentum sections are integrated. For `0<k<=pi`, the second section
`sum(p)=k+2*pi` is reflected into `Q=2*pi-k` during patch construction, then its
rapidities are reflected back. Momenta above pi use an overall reflection.
The callback therefore sees the actual section momenta, without assuming its
weight is reflection invariant.

## Preserve soft distances

For a radial interval starting at zero, use `r=r_max*exp(-t)` and several panels
in `t`. A finite configured cutoff `T` omits `t>T`. This is a **log-radius
cutoff**, not a uniform rapidity cutoff; geometry supplies no omitted-tail
estimate.

A positive lower radius instead uses a sine-squared map of `log(r)` over the
complete allowed interval. This resolves decades near a moving soft boundary
without adding a cutoff there. Direct moments also split where the hard-pair
interval changes formula.

The code retains `log(a)` and `log(b)` when a soft distance cannot be
exponentiated. The correction to `log(tan(a/2)) = log(a)-log(2)` is omitted
only below the working arithmetic epsilon. Rapidities remain finite instead
of rounding momenta onto zero or pi. The chosen density backend must still
support the resulting coordinates and exponent range.

## Connect a density

This finite-order example calculates a truncated longitudinal spectral slice:

```cpp
#include <bethe/detail/xxx_four_spinon_geometry.hpp>
#include <bethe/xxx_spinon_density.hpp>

using R = long double;
R pi = 4 * uni20::math::atan(R{1});
bethe::heisenberg::detail::FourSpinonGeometry<R> geometry(8, R{16});
bethe::heisenberg::XXXSpinonDensity<R> density;
bethe::detail::CompensatedSum<R> sum;
geometry.spectrum(pi, R{2}, [&](auto const& node) {
    auto point = density.evaluate(node.beta, {64, 16});
    sum.add(uni20::math::exp(uni20::math::log(point.value) + node.log_measure));
});
// Approximately 0.09810221593 at this particular order and cutoff.
```

This checks normalization and geometry, not unrestricted-spectrum accuracy.
A complete integrator must propagate contour/arithmetic diagnostics, compare
quadrature orders, assess the omitted soft region and recover or report failures.

`geometry.moments(k, callback)` supplies direct fixed-momentum nodes. Integrate
`W4*exp(log_measure)/(2*pi)` for the zeroth moment; multiply by `node.energy`
for the first. The frequency measure `1/(2*pi)` belongs to the observable.
Neither that measure nor exchange/channel conversion is part of geometry.
For exchange J, spectra scale by `1/J` with frequency `omega/J`; order-r moments
scale by `J^r`.

MPFR controls require explicit finite precision. Gauss nodes, support inversion,
rapidities and measures stay at that precision; mixed, unset or exact-only
runtime controls are rejected. The prepared rule is immutable and shareable.

## Direct global moments

`bethe/detail/xxx_four_spinon_common_shift.hpp` provides a second geometry for
global moments. Order the rapidities and use their three successive gaps plus
a common shift. The longitudinal density obeys

```math
W_4(\boldsymbol\beta)\prod_j\mathrm{sech}\,\beta_j
=F(a,b,c)\sin^2(P/2),\qquad P=\sum_j2\arctan e^{\beta_j}.
```

`FourSpinonCommonShift` constructs a representative with `P=pi`, preserving
the three gaps. At that point the sine factor is one, so the physical density
and the returned rapidity Jacobian give `F` without division by a small sine.
The common shift can then be integrated independently of the expensive kernel:

```math
A_0(a,b,c)=\int_{-\infty}^{\infty}d\eta\,\sin^2(P/2),\qquad
A_1=\int_{-\infty}^{\infty}d\eta\,E\sin^2(P/2)=\pi^2.
```

Here `J=1`. The first-moment identity follows from `dP/deta=2E/pi` and the increase of
`P` from zero to `4*pi`. The zeroth integral uses a regrouped partial-fraction
expression that avoids subtracting terms growing with the middle gap. When
its arithmetic indicator is too large, positive quadrature supplies the value.
That fallback reports its order difference, arithmetic indicator and an
explicit bound on its omitted exponential shift tails separately. None of
these estimates bounds the omitted **gap** region.

`FourSpinonGapGeometry(order, cutoffs, max_nodes)` emits three-dimensional gap
nodes. It uses geometrically increasing panels with open Gauss nodes, a checked
node budget, and explicit precision. Indexed `node(i)` access permits bounded
parallel evaluation without storing the full tensor grid. Each node contains
the representative
rapidities, `log_measure` with the ordered factor 24 and rapidity Jacobian,
`zeroth_shift`, and `first_shift`. The observable supplies both `1/(2*pi)`
measures. Equivalently, the physical global moments are

```math
I_{0,4}=\frac{6}{\pi^2}\int da\,db\,dc\,F A_0,\qquad
I_{1,4}=6\int da\,db\,dc\,F.
```

Each finite cutoff restricts one gap, while the common shift has already been
integrated over the whole real line. This domain differs from both an absolute
rapidity box and the fixed-momentum geometry's log-radius cutoff. A result
must retain its actual domain and assess its gap tails separately.

The helper verifies physical common-shift invariance pointwise. Independent
positive shift integration checks the analytical zeroth factor and the first
integral `pi^2`; a separable test density checks the ordered and observable
measures together. Coincident-gap limits are supported by the shift integral;
they do not remove the density evaluator's separate coincidence restrictions.

## Independent checks

- Forward momentum and energy reconstruction, ordering and hard-pair velocities
  for both patch types, including MPFR beyond binary128.
- Integrating the energy delta recovers the original hard-pair interval.
- Constant and energy-weighted kernels reproduce independent box-spline section
  volumes, checking the factor of 24 and both momentum sections together.
- The 512-point example agrees with an independently parameterized radial
  quadrature using the separate hypergeometric four-spinon formula.
- Near-edge support tests cover the same-end fraction endpoint and the
  opposite-end interior maximum. MPFR resolves a frequency only `1e-30` below
  the upper edge; an edge unresolved at the working precision fails explicitly.
- Soft-distance tests retain rapidities at `log(r)=-20000`; invalid controls
  and unresolved boundaries fail explicitly.
