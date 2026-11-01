# Spectral quadrature and time signals

[C++ library](library.md) · [Thermodynamic XXX spectra](xxx-structure-factor-thermo.md)

`scripts/spectral_time_transform.py` converts sampled densities to time signals
and propagates their sampling covariance. It uses NumPy binary64 arithmetic;
it does not inherit the C++ calculation's arbitrary precision. No private data,
plotting package or compiled Bethe program is required for the analytic example.

From a source checkout:

```sh
python3 -m venv .venv-analysis
. .venv-analysis/bin/activate
python -m pip install -r scripts/requirements-analysis.txt
python scripts/example_spectral_time_transform.py
```

For your own program, put the checkout's `scripts` directory on `PYTHONPATH`,
or place your program beside the module. There is no installed Python package
yet. NumPy 1.24 or later, below 3, is required; CI checks the 1.x minimum and
a current 2.x environment.

## Integration and physical measure

For the density convention used by the thermodynamic tools,

```math
C(t)=\int_{\omega_a}^{\omega_b}\frac{\mathrm{d}\omega}{2\pi}
S(\omega)e^{-\mathrm{i}\omega t}.
```

`weights` returns a complex matrix with one row per time and one column per
supplied frequency:

```python
import numpy as np
from spectral_time_transform import weights

omega = np.array([1.0, 1.1, 1.6, 2.0])
density = np.ones_like(omega)  # Analytic example, not a spin-chain spectrum.
times = np.linspace(0, 20, 81)
W = weights(omega, times, order=32)
signal = W @ density
```

The default representation interpolates `S(omega)` linearly. An `order`-point
Gauss–Legendre rule per cell integrates the exponential against that
interpolant; the default is 24. Negative times are supported. Only the sampled
interval is integrated. `window=(lo, hi)` restricts it further and must lie
within the selected nodes. There is no endpoint filling, tail extrapolation,
broadening or intensity adjustment.

These weights apply to a **density**, not discrete Lehmann line weights. For
finite lines, compute `sum(weight * exp(-1j * gap * time))` directly: their
delta-function measure already cancels the `2*pi`.

## Threshold coordinates and moments

For a square-root lower threshold, supply positive offsets evaluated separately
from the rounded absolute frequencies:

```python
delta = np.geomspace(1e-6, 1.0, 65)
omega = 2.0 + delta
density = 1 / (2*np.sqrt(delta))  # Another analytic example.
W = weights(omega, times, offsets=delta)
signal = W @ density
m0 = (weights(omega, [0.0], offsets=delta) @ density)[0].real
m1 = (weights(omega, [0.0], offsets=delta, energy_power=1) @ density)[0].real
```

This interpolates `G(x)=2*x*S(omega_lower+x*x)` linearly in `x=sqrt(delta)`.
For this example `G` is constant. Supplied offsets are retained at exact grid
nodes rather than reconstructed by subtracting rounded frequencies. Frequencies
and transformed coordinates must still be distinguishable in binary64. A
singular node at offset zero is rejected; its missing interval is not filled.

`energy_power=r` multiplies the **same original interpolant** by `omega**r`
inside the quadrature. Applying ordinary weights to `omega*density` constructs
a different interpolant and generally gives a different first moment.
With consistent quadrature, `C(0)=m0` and `C'(0)=-i*m1`.

## Shared-sample covariance

Let each row of `samples` contain one independent group estimate of the entire
frequency grid. With at least two groups:

```python
from spectral_time_transform import statistics

mean = samples.mean(axis=0)
covariance = np.cov(samples, rowvar=False, ddof=1) / len(samples)
result = statistics(W, mean, covariance)
```

The covariance is that **of the estimated mean**, including off-diagonal
frequency correlations. Returned lists are `real`, `imag`,
`real_standard_error`, `imag_standard_error` and `real_imag_covariance`.
The last is the covariance of the real and imaginary parts of each output row.
Input means/covariances are real; complex cross-spectral statistics would need
an explicit real/imaginary block representation first.

For two transforms of the same data, propagate `W1-W2` to obtain the error of
their difference. A contrast between times similarly uses the difference of
rows. Separate error bars cannot be combined as though they were independent.

For two sectors sampled with shared groups, concatenate their sample columns
and retain the full joint covariance:

```python
joint = np.hstack([sector_a_samples, sector_b_samples])
mean = joint.mean(axis=0)
covariance = np.cov(joint, rowvar=False, ddof=1) / len(joint)
combined = statistics(np.hstack([Wa, Wb]), mean, covariance)
difference = statistics(np.hstack([Wa, -Wb]), mean, covariance)
```

`Wa` and `Wb` can use different frequency grids but must have matching output
times. Cross blocks vanish only for independent estimates. The function checks
dimensions, finiteness, symmetry and positive semidefiniteness, tolerating
roundoff. The covariance eigenvalue check costs `O(nodes^3)`. Returned marginal
statistics avoid storing a dense `times × times` covariance; obtain any required
time contrast through its rows instead. Inputs are not modified.

Ordinary rows use NumPy matrix products. For extreme exponents, covariance
propagation falls back to exact rational products of the supplied binary64
entries, then rounds the standard errors and cross covariance to binary64.
This slower path preserves representable standard errors even when their
squared variances underflow or overflow, and preserves tiny covariances that
large coefficients amplify. It uses Python's standard library, with no extra
dependency. It does not improve the precision of the input measurements or the
symmetry/positive-semidefinite check. An unrepresentable returned value still
raises `OverflowError`; finite means use the usual NumPy transform.

## Separate accuracy controls

1. **Phase quadrature:** increase `order` at fixed grid. Long times may need
   many more points per cell. There is no automatic error estimate or refinement.
2. **Spectral discretization:** refine the frequency grid or compare subsets.
   `coarser_indices(count, stride)` retains both endpoints; pass its output as
   `indices=`. Unselected columns receive zero weights. Increasing phase order
   does not fix an inaccurate interpolant.
3. **Sampling and omitted support:** propagate full covariance and assess
   unsampled frequency/rapidity regions separately. Sampling error excludes
   interpolation bias, a missing threshold interval and omitted sectors.

Invalid inputs raise `ValueError`; nonfinite phase, quadrature or transformed
arithmetic raises `OverflowError`. Failed or absent source values must be
handled before constructing a density array, not zero-filled.

Run the independent analytic and covariance checks with:

```sh
python -m unittest discover -s scripts -p 'test_spectral_time_transform.py' -v
```
