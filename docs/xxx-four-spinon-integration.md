# Four-spinon integration geometry

[Spinon densities](xxx-spinon-density.md) · [C++ library](library.md)

A four-spinon spectrum requires two-dimensional integration after imposing
momentum and energy conservation. Direct fixed-momentum moments require three
dimensions. Both use the ordered momentum patches in
`bethe/detail/xxx_four_spinon_geometry.hpp`.

This is an internal integration building block, not a new command-line spectrum
or a convergence controller. It supplies coordinates and measures to a callback;
the physical weight comes from `XXXSpinonDensity`. Kernel failures propagate.
Callbacks may already have run when a failure occurs; discard that partial
accumulation rather than treating it as a completed integral.

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
Inverting these bounds locates the occupied region directly. The quadrature
therefore never relies on finding support by sampling zero padding.
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

A positive lower radius needs a different coordinate. Close to the boundary of
the allowed fraction interval, the two radial endpoints approach each other
quadratically and can round to the same number. Clustering radial nodes at
those endpoints makes this worse as the quadrature order increases. Increasing
the order alone is therefore not a reliable refinement strategy for that map.

Instead use the half-separation of the remaining pair,
`theta=(p_4-p_3)/2` for the same-end patch, and its analogous pair for the
opposite-end patch. The energy and ordered measure become

```math
E(r,\theta)=\frac{\pi}{2}[\sin(rv)+\sin(r(1-v))]
             +\pi\sin(H/2)\cos\theta,
\qquad
d\mu=\frac{r}{E_r}\,d\theta\,dv.
```

Here `E_r` is the radial derivative at fixed angle. It is positive in the
interior of the ordered patch. Invert the lower energy bound to find the
largest allowed radius and its hard-pair angle. For each interior angle,
solve for the distance from that radius using safeguarded Newton iteration.
Trigonometric product identities evaluate the energy difference directly,
retaining distances smaller than an ulp of the full energy or radius. This
avoids subtracting either nearly equal radial endpoints or nearly equal
energies. The angular measure also removes the hard-pair square-root
singularity without dropping any support.

Use ordinary Gauss nodes in this angle and in fraction intervals clipped by
energy support. Unclipped fractions retain the endpoint-resolving map. The
complete positive-radius region is included without an additional cutoff.
Direct moments also split where the hard-pair interval changes formula.

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

## Independent checks

- Forward momentum and energy reconstruction, ordering and hard-pair velocities
  for both patch types, including MPFR beyond binary128.
- Integrating the energy delta recovers the original hard-pair interval.
- Angular and radial coarea integrals agree independently for constant and
  momentum-dependent weights in both patch types.
- Constant and energy-weighted kernels reproduce independent box-spline section
  volumes, checking the factor of 24 and both momentum sections together.
- The 512-point example agrees with an independently parameterized radial
  quadrature using the separate hypergeometric four-spinon formula.
- Near-edge support tests cover the same-end fraction endpoint and the
  opposite-end interior maximum. MPFR resolves a frequency only `1e-30` below
  the upper edge; an edge unresolved at the working precision fails explicitly.
- Previously failing positive-support points remain resolved through order
  256. The native angular measure agrees with an independent 160-bit evaluation
  of the same finite rules; this checks arithmetic, not quadrature convergence.
- Soft-distance tests retain rapidities at `log(r)=-20000`; invalid controls
  and unresolved boundaries fail explicitly.
