# Experimental finite-q XXZ vertex spectra

[Experimental APIs](experimental.md) · [Massive two-spinon density](xxz-structure-factor-thermo.md)

This opt-in package evaluates the original **two- and four-spinon complex
vertex amplitudes** for the infinite, zero-field massive XXZ chain. It also
provides fixed-order spectral-angle quadrature for moments and Gaussian
frequency cuts. It complements the supported, unbroadened two-spinon evaluator;
the four-spinon quadrature is an experimental calculation with independent
product, residue and outer-grid controls.

Use the spin convention and exchange of the other XXZ tools:

```math
H=J\sum_j(S_j^xS_{j+1}^x+S_j^yS_{j+1}^y+\Delta S_j^zS_{j+1}^z),
\qquad S^\alpha=\sigma^\alpha/2,\quad J\gt 0,\quad\Delta\gt 1.
```

The observable is the transverse raising response, with physical momentum
k=π−Q in the vertex convention. There is no fitted amplitude, normalization
switch or rescaling to a sum rule. Two and four spinons remain separate sectors.

## Build and run a small example

Configure from a separate build directory:

```sh
cmake /path/to/bethe -DBETHE_BUILD_EXPERIMENTAL=ON -DBETHE_BUILD_EXAMPLES=ON
cmake --build . --target bethe_example_xxz_vertex
./experimental/bethe_example_xxz_vertex
```

The example evaluates both sectors at Δ=1.5 with order 4, displays three
broadened frequencies and the corresponding fixed-momentum and full-zone
moments. This deliberately small grid demonstrates the API, **not convergence**.
Its two-worker executor uses the shared observable-grid final result; it also
checks equality with serial evaluation. There are no raw worker threads or
scheduler-wide drains in the example.

Link consumers to `bethe::experimental`. This target and its examples are absent
by default and are not part of the normal application installation.

## Amplitudes and arithmetic

```cpp
#include <bethe/experimental/xxz_vertex.hpp>

namespace model = bethe::experimental::xxz;
model::VertexEvaluator<double> vertex(1.5);
std::array<double, 4> angles{-1.2, -.2, .6, 1.3};
auto amplitude = vertex(angles);
auto controls = vertex.options();
```

The input angles φ lie in [−π/2,π/2], with ξ=i exp(−iφ), u=−ξ² and
q=−exp(−acosh Δ). `value[vacuum][polarization]` contains complex amplitudes
for both Néel vacua. There is one polarization for two spinons and four for
four spinons. Retaining phases matters: a resolved momentum uses coherent
combinations of the two vacua, rather than a sum of their individual intensities.

`VertexOptions{products, residues}` sets independent truncation counts. Zero
selects a count from Δ and arithmetic precision; `options()` reports the actual
counts. Products must be 4–20000 and residue levels 1–20000. An automatic count
outside those limits is rejected. Refining arithmetic alone does not establish
convergence in either count.

Products retain a separate binary exponent until prefactors cancel. The first
residue levels use literal q-products; later levels use a cancelled recurrence.
Current Uni20 scalar math and scaled theta functions supply kinematics. The
complementary elliptic modulus is retained directly so that a small gap is not
lost by subtracting two numbers near one.

Binary64, configured 80-bit extended precision and binary128 are tested,
alongside explicit MPFR/MPC arithmetic. Construct Δ and
all query coordinates at the same finite precision; exact default exchange 1
and default grid offset 0 are accepted. Other mixed, unset or exact continuous
inputs are rejected. For example, use `uni20::mpreal("1.5", precision)` and
`VertexEvaluator<uni20::mpreal>` for a selected working precision. No ambient
precision is read or changed.

Coincident or unresolved spectral nodes require confluent formulas and are
rejected. Unresolved soft momenta, a physical momentum shift lost when forming π−k,
exhausted scalar range and nonfinite sums
also fail explicitly. These exceptions are not successful zero contributions.

## Moments and broadened cuts

```cpp
#include <bethe/experimental/xxz_vertex_spectrum.hpp>

namespace model = bethe::experimental::xxz;
model::ThermodynamicVertexSpectrum<double> spectrum(1.5);
auto full = spectrum.integrated(4, 12);
std::array<double, 3> omega{1.8, 2.2, 2.6};
auto cut = spectrum.broadened_cut(4, 1.5, omega, .06, 12);
```

`integrated(n, order, offset)` uses staggered periodic spectral-angle grids.
The measure includes the two-vacuum average and the completeness factor:

```math
W_n=\frac{1}{2n!}\int\prod_{j=1}^{n}\frac{d\phi_j}{\pi}
\sum_{i,\ell}|X_{i,\ell}(\phi)|^2.
```

This is the full-zone integral with dk/(2π) and dω/(2π). The returned
`first_moment` additionally weights each node by its total excitation energy.
The complete transverse zeroth-moment sum is 1/2; a selected sector does not
exhaust it. An order-N calculation uses Nⁿ amplitude nodes.

`broadened_cut(n, k, omega, width, order, offset)` resolves momentum first,
retains both coherent vacuum combinations, then applies a Gaussian with
standard deviation `width`. Its density is S itself: the Gaussian integrates
to 2π over the entire frequency axis. The cut's `moments` and `sector_weights`
are computed directly from the momentum measure, independently of the requested
frequency grid. Gaussian tails extend to negative frequencies; clipping them
or integrating a finite plotted interval will not reproduce those moments.

Orders are 2–512, even for cuts. Four-spinon cuts use at most 4N³ amplitude
nodes. The offset shifts staggered free-angle grids and is a useful quadrature
diagnostic; it is not randomized automatically. A cut accepts an optional
executor that invokes every row exactly once and joins all accepted work before
returning, including on failure. Reduction remains in fixed row order. The
example supplies that contract through Bethe's existing Async grid.

## What the diagnostics establish

Every amplitude and moment carries an `arithmetic_indicator`. This is a
conditioning screen based on roundoff and absolute residue sums. It is **not a
certified bound**, a product/residue tail estimate, an outer-quadrature error,
or an estimate of omitted sectors. No single fixed order receives a convergence
status. Increase products, residues, precision and outer order independently;
also vary offsets and inspect near-collision/soft regions. This becomes
especially important approaching Δ=1.

Tests compare eight complex amplitudes with an independent 100-digit residue
implementation, check product/residue refinement and precision ownership, and
compare fixed-order integrated weights and cuts with independent Python
calculations. Two-spinon cuts also agree with a Gaussian convolution of the
supported unbroadened formula. Ising and near-isotropic checks are limited
controls, not a proof of general quadrature convergence.

The reference derivation is
[Abada–Bougourzi–Si-Lakhal, Eqs. (2.6)–(2.7)](../CITATIONS.md#abada-bougourzi-si-lakhal-1997).
The independent scripts, point fixture and fixed-order quadrature controls are in
[`experimental/reference/xxz`](../experimental/reference/xxz/README.md).
They do not use the printed isotropic four-spinon expression to fix a prefactor.
