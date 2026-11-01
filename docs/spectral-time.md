# Finite spectra and time correlations

[C++ library](library.md) · [XXX structure factor](xxx-structure-factor.md)

Finite spectral lines give the time correlation directly:

```math
C_{AB}(t)=\sum_n c_n e^{-\mathrm{i}(E_n-E_0)t},\qquad
c_n=\langle0|A|n\rangle\langle n|B|0\rangle.
```

There is no time-stepping error. Accuracy still depends on the retained states,
their energies and coefficients, and the arithmetic used to sum their phases.
The helpers preserve the supplied physical normalization: no broadening, fitted
amplitude or completion of missing states is applied.

## Local and cross spectra

Include `<bethe/local_spectrum.hpp>`. An autocorrelation uses nonnegative
`LocalSpectralLine<Real>::weight`. A cross correlation uses
`CrossSpectralLine<Real>::coefficient`, of type `uni20::complex<Real>`; even a
real coefficient can be negative. Each line also carries a `state_id` and a
finite nonnegative `gap`. An elastic line with gap zero is allowed.

```cpp
#include <bethe/local_spectrum.hpp>
#include <array>

std::array<bethe::LocalSpectralLine<double>, 2> lines{{
    {0, 0.0, 2.0}, {1, 1.0, 3.0}
}};
auto moments = bethe::local_spectral_moments<double>(lines);
auto signal = bethe::local_time_correlation<double>(lines, 0.5);
// moments.weight == 5, moments.first_moment == 3
// signal == 2 + 3*exp(-0.5i), within floating-point accuracy
```

| Function | Result |
| --- | --- |
| `local_spectral_moments(lines, exemplar)` | Real `weight` and `first_moment` |
| `cross_spectral_moments(lines, exemplar)` | Complex `value` and `first_moment` |
| `local_time_correlation(lines, time)` | Complex autocorrelation |
| `cross_time_correlation(lines, time)` | Complex cross correlation |

The moments are sums of coefficients and of gap times coefficient. Consequently
`C(0)` is the zeroth moment and `C'(0)` is `-i` times the first moment. All sums
use Neumaier compensated accumulation, separately for the real and imaginary
components. For an autocorrelation, `C(-t)=conj(C(t))`.
For complex cross coefficients that identity does not hold in general;
conjugating the coefficients and reversing time conjugates the result.

Operator/channel factors belong to the caller. In particular, these local
weights have no implicit system-size or momentum Fourier factor. If a spectral
density is written as `S(omega)=2*pi*sum c_n*delta(omega-gap_n)`, its integral
against `exp(-i*omega*t)/(2*pi)` is exactly the finite sum above.

## Momentum-labelled spectra

Existing `LehmannLine<Real>` data in `<bethe/dynamical_structure_factor.hpp>`
can be used directly:

```cpp
auto signal_by_momentum = bethe::spectral_time_correlation<double>(
    number_of_momenta, lines, time);
```

Here `lines` is a span or container of `LehmannLine<double>`. The result has one
complex value for each momentum index, including zeros for indices with no
supplied lines. The numerical phase sum is shared with the local helpers.
The stored momentum value is metadata; the integer `momentum_index` selects
the output entry, just as in `spectral_moments`.

## Precision and failures

The helpers support fp64, native long double, optional binary128 and MPFR.
For runtime scalars, the `time` argument supplies the finite working precision.
Moment sums instead require an explicit finite-precision `exemplar`:

```cpp
using Real = uni20::mpreal;
auto const precision = uni20::Precision::bits(160);
Real const zero(0, precision), time(1, precision);
std::array<bethe::LocalSpectralLine<Real>, 2> runtime_lines{{
    {0, Real(0, precision), Real(2, precision)},
    {1, Real(1, precision), Real(3, precision)}
}};
auto moments = bethe::local_spectral_moments<Real>(runtime_lines, zero);
auto signal = bethe::local_time_correlation<Real>(runtime_lines, time);
```

The exemplar's precision is used even for an empty spectrum; its numeric value
does not alter the result. Native callers can omit it. Runtime callers must
explicitly convert exact constants, unset values or mixed precisions before
passing them as inputs. No ambient MPFR precision is read or changed.

Invalid gaps, coefficients, times, precision combinations or momentum indices
throw `std::invalid_argument`. Nonfinite arithmetic products, accumulated sums
or phases throw `std::overflow_error`. These helpers do not silently drop bad
lines or replace failed matrix elements with zero. A finite phase can still
lose accuracy at very long times if the supplied energies are insufficiently
accurate; higher summation precision cannot recover missing energy digits.

Independent analytic tests cover zero and nonzero times, complex cross phases,
moment derivatives, signed cancellation and explicit runtime precision. The
tests are small; no calculation archive or model solve is needed to run them.
