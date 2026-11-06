# Momentum input near an endpoint

A tiny momentum can be accurately represented without a large mantissa. The
problem near pi is subtraction: a scalar `pi - 1e-100` rounds to pi in ordinary
floating point. Bethe's shared momentum input keeps the logarithmic coordinate
until a consumer explicitly requests a scalar angle.

The XXX thermodynamic spectrum and moment tools accept three mutually exclusive
forms:

| Option | Meaning |
|---|---|
| `--momentum 1e-100` | A momentum in radians/site; range 0 to 2 pi |
| `--momentum-from-pi 1e-100` | q = pi minus the specified nonnegative distance |
| `--momentum-beta 1000` | b = log(tan(q/2)), with 0 < q < pi |

`--momentum pi` and `--momentum '2*pi'` specify exact endpoints. Zero distance
in `--momentum-from-pi` specifies exact pi; distance `pi` specifies exact zero.
A finite beta always represents an interior point, however large its magnitude.
The ordinary momentum form also preserves reflection from (pi, 2 pi).

For large |b| the distance to the nearest endpoint is approximately
2 exp(-|b|). Thus b = 1000 can be stored in binary64 even though the endpoint
distance cannot. This b describes **external total momentum**, not the
rapidities of the individual spinons being integrated.

## Input and working precision

When MPFR is enabled, decimal input is parsed at a separate guarded precision,
chosen from the input length and the requested working precision. Conversion to
beta precedes conversion to the working scalar type. For example, a long decimal
representation of pi minus 1e-100 can produce a useful binary64 beta even though
binary64 cannot retain that momentum as a single angle. MPFR global precision
settings are not read or changed.

Without MPFR, input conversion uses binary128 when enabled, otherwise native
long double. These have finite exponent range and precision. An ordinary decimal
input indistinguishable from pi is rejected; use an explicit distance or beta.
Nonzero decimal input that underflows to zero is also rejected.

`--momentum-beta-tolerance` sets the warning threshold for estimated **absolute**
beta conversion error; the default is 1e-12. Near an endpoint, an absolute change
in beta is approximately the relative change in endpoint distance. This is
separate from kernel, quadrature, sampling and rapidity-cutoff tolerances. The
reported estimates include a guard-precision roundoff allowance and its
amplification when subtracting an angle from pi or 2 pi; they are not
rigorous interval bounds or an uncertainty assigned to the supplied decimal.

Two conversions are reported separately:

1. High-precision input to the working beta coordinate.
2. The coordinate to a scalar angle, when the selected backend requires one.

Warnings go to stderr. Output metadata retains the original input form and text,
working beta, reflection, input precision, estimated conversion errors and warning
status. Exact endpoints have no finite beta; metadata uses an absent value.

## Current geometry boundary

The public XXX spectrum and moment kernels still consume scalar angles. Their
adapter rejects an interior coordinate that would collapse onto an endpoint.
For instance, `--momentum-from-pi 1e-100 --precision fp64` fails explicitly;
sufficient MPFR geometry precision can retain that scalar displacement.
A representable scalar angle can still lose accuracy in its endpoint distance;
that produces the separate scalar-geometry warning.

The new interface does not itself migrate the integration geometry or eliminate
conditioning inside a kernel. Beta-coordinate consumers can use the preserved
coordinate directly, including endpoint distances beyond the scalar exponent
range. Frequency offsets from continuum boundaries require their own treatment.

## Shared C++ interface

```cpp
#include <bethe/momentum_input.hpp>

auto input = bethe::MomentumInput{
    bethe::MomentumInputKind::distance_from_pi, "1e-100"};
auto q = input.resolve(1.0); // binary64 coordinate, independently parsed input

double beta = *q.coordinate.beta();
double log_distance = *q.coordinate.log_distance();
// Pass beta/log_distance directly to a consumer that supports this coordinate.
// q.require_radians() throws here: binary64 scalar q would equal pi.
```

`MomentumCoordinate<Real>` distinguishes exact endpoints, finite beta and
reflection. `ResolvedMomentum<Real>` carries the conversion diagnostics and an
optional scalar angle. `MomentumInput` owns its original text and is independent
of CLI11, model code and the number of spinons. The shared application helper
`apps/momentum-input.hpp` registers the options, reports metadata and emits
warnings. New frontends should use those facilities rather than parsing through
double or reconstructing pi minus a tiny distance themselves.
