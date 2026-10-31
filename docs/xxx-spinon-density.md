# Thermodynamic XXX spinon densities

[C++ library](library.md) · [Two-spinon spectrum](xxx-structure-factor-thermo.md)

`bethe::heisenberg::XXXSpinonDensity<Real>` evaluates the physical longitudinal
weights of two, four or six spinons in the zero-field, zero-temperature XXX
chain. These are the ingredients for thermodynamic spectra and direct moments.
They are not spectra by themselves: the conservation laws still need integration.
The [four-spinon geometry](xxx-four-spinon-integration.md) provides ordered
spectral and direct-moment nodes with explicit measures and cutoffs.

## Physical convention

For `H = J sum S_j . S_(j+1)`, `J > 0`, and `S = sigma/2`, define

```math
S_n^{zz}(k,\omega)=\int_{(0,\pi)^n}d^np\,
W_n(\mathbf p)\delta_{2\pi}\!\left(k-\sum_jp_j\right)
\delta\!\left(\omega-\sum_j\frac{\pi J}{2}\sin p_j\right),
\qquad \delta_{2\pi}(x)=\sum_{s\in\mathbb Z}\delta(x-2\pi s).
```

The evaluator returns W_n, including the longitudinal spin convention,
polarization sum and exactly one `1/n!`. The momenta in this integral are
labelled; do not apply another factorial to the weight. A library integrator
using only ordered momenta must supply its domain multiplicity internally.
Every allowed total-momentum section and both roots of an eliminated pair
must be accounted for by the geometry.

The input is beta, related to a positive spinon momentum by

```math
p=2\arctan e^\beta,\qquad dp=\mathrm{sech}\,\beta\,d\beta.
```

Despite the rapidity argument, the returned weight is against `product(dp)`.
Changing integration variables to beta requires the indicated Jacobian.
W_n is independent of J; the energy delta gives the spectrum its inverse-energy
units. Moments use `domega/(2pi)`, and the global average uses `dk/(2pi)`.
These measures belong to the observable integrator, not a plotting correction.

## Evaluate a density

```cpp
#include <array>
#include <bethe/xxx_spinon_density.hpp>

bethe::heisenberg::XXXSpinonDensity<double> density;
std::array<double, 4> beta{-2.1, -1.3, -0.4, 0.3};
bethe::heisenberg::SpinonDensityOptions options;
options.residue_terms = 64;
options.tail_order = 16;
auto result = density.evaluate(beta, options);
// result.value is physical W4; inspect the diagnostics below before using it.
```

The span length selects 2, 4 or 6 spinons. Other lengths and nonfinite inputs
are rejected. Contexts are immutable and shareable across concurrent calls;
scratch storage belongs to each evaluation. Preparing a context is useful
when many rapidity sets share one precision. There is no internal thread pool.

For MPFR, select finite precision explicitly and parse inputs at that precision:

```cpp
auto p = uni20::Precision::bits(256);
bethe::heisenberg::XXXSpinonDensity<uni20::mpreal> density(p);
std::array beta{density.parse("-2.1"), density.parse("-1.3"),
                density.parse("-.4"), density.parse(".3")};
auto result = density.evaluate(beta, {64, 16});
```

Mixed or unset input precisions and exact-only MPFR arguments are rejected;
use the context parser or an explicit precision conversion. No ambient MPFR
precision is changed. `precision()` and `bits()` identify the context.

## Accuracy and failures

Evaluation produces an estimate, not a convergence verdict. The default
32 residue terms and order-12 tail expansion are starting budgets. Increase
them and compare arithmetic precisions when assessing a point; extra bits
do not control an omitted contour tail or an outer integration error.

| Result field | Meaning |
| --- | --- |
| `value` | Physical W_n |
| `contour_tail_error` | Propagated contour truncation majorant, evaluated in ordinary floating point |
| `arithmetic_indicator` | Propagated roundoff/cancellation indicator, including special-function evaluation errors |
| `log_scale` | Logarithm of the physical coefficient times all pair factors |
| `contour_norm` | Unnormalized polarization sum of squared contour amplitudes |
| `contour_norm_error`, `contour_arithmetic_indicator` | Corresponding contour diagnostics before physical scaling |
| `amplitudes` | Four polarization amplitudes for W4, fifteen for W6, or one for W2 |
| `terms`, `order` | Actual residue and tail budgets used |

The contour majorant is not an interval enclosure. The arithmetic quantity is
a diagnostic, not a rigorous bound. Agreement of two precisions tests arithmetic,
not convergence of the fixed contour approximation. Neither diagnostic covers
rapidity cutoff, sampling, or momentum/frequency integration.

`SpinonDensityError` carries a `SpinonDensityFailure` reason. Unresolved
collisions, a contour tail outside its supported majorant domain, nonfinite
arithmetic and exponent-range failures are not reported as physical zeros.
Native subnormal scale factors, contour norms and weights are rejected because
the relative-roundoff model no longer applies there; use a wider exponent range
or MPFR. Exactly coincident rapidities need a separate limiting calculation.
The initial implementation evaluates separated finite points; it does not automatically
recover them by changing precision or discard failed samples. A higher-level
integrator must retain explicit failure or perform bounded recovery.

## Formula and validation

The pair factor is evaluated through its Barnes-G representation. Gamma
evaluations use Uni20's reusable complex log-Gamma context. Residues share
rational recurrences across polarization components; the six-spinon double
contour separates into four single-contour moments per label. A controlled
rational-ratio tail supplements the finite residue sum. Neumaier summation
retains compensation when a later term is larger than the accumulated sum,
with independent accumulators for the real and imaginary components.

The absolute coefficients include the state-counting factorial at the density
layer. Independent validation uses:

- The existing two-spinon spectrum, reconstructed from W2 and the two labelled
  roots, including exchange scaling.
- A separately implemented four-spinon 8F7 formula with oscillatory pair
  integrals and physical coefficient `G(1/2)^16/3072`.
- High-precision Python fixtures for the same finite contour approximation,
  including individual complex amplitudes, plus independent Hurwitz-zeta tails.
- Permutation checks and precision-growth checks beyond binary128.

The independent Python scripts are validation tools, not runtime dependencies.
`scripts/reference_native_vertex_fixtures.py` regenerates the finite-contour
fixtures with mpmath. `scripts/reference_xxx_four_spinon.py` supplies the
independent physical four-spinon density, including its coefficient internally.
Its `raw_kernel` helper exposes the unnormalized series for algebra checks only.
No sum-rule fraction is used to fit a normalization.
