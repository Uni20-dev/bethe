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

## A fresh finite-chain versus thermodynamic example

With the installed Bethe programs on `PATH`, run:

```sh
python scripts/example_xxx_two_spinon_time.py --output xxx-time-example
```

This small calculation compares the longitudinal two-spinon response at
`q=pi/2`, using a 16-site periodic chain and the infinite chain. It calculates
fresh spectra with `bethe-xxx-structure-factor` and
`bethe-xxx-structure-factor-thermo`; no downloaded data or plotting dependency
is needed. `--finite-solver` and `--thermo-solver` accept explicit executable
paths when the programs are not on `PATH`. Both calculations use fp64 and one
solver thread. The example's default grids require only seconds on an ordinary
CPU, not a production run.

`time.csv` contains real and imaginary signals on `0 <= Jt <= 20`:

- `finite_all_*`: every calculated finite-chain two-spinon line at this momentum.
- `finite_window_*`: only those lines within the sampled thermodynamic interval.
- `thermo_window_*`: the integrated thermodynamic two-spinon density on that interval.

The finite lines are summed directly without broadening. The thermodynamic
frequency grid is concentrated near the lower threshold using geometric
segment boundaries, with 129 uniformly spaced points in each segment
(769 unique frequencies by default). The integral uses the threshold-coordinate
interpolant described below. The first sampled point is `omega_lower + 1e-6 J`:
the divergent endpoint and its missing interval are **not filled or estimated**.
The upper boundary is included with its zero density. Both signals omit
higher-spinon sectors, and the finite spectrum is not a thermodynamic
extrapolation; its weights and line positions need not match the infinite chain.

`report.json` records the support, zeroth/first moments, full first-moment sum
rules and calculation provenance. The original solver CSV exports are retained,
including their conventions, build revisions, controls and convergence status.
An existing output directory is refused. A failed source calculation stops the
example, retaining any earlier exports without declaring a completed report.

The time table separately records changes from doubling phase-quadrature order
and coarsening the frequency grid, plus propagation of the solver's pointwise
kernel error estimates. These are diagnostics, not a combined error bar or a
bound on the missing threshold interval. Change these controls separately:

```sh
python scripts/example_xxx_two_spinon_time.py --points 257 --output xxx-time-grid-refined
python scripts/example_xxx_two_spinon_time.py --points 257 \
  --lower-offset 1e-7 --output xxx-time-refined
```

`--sites` accepts multiples of four, keeping the momentum exactly `pi/2`.
`--time-max`, `--time-points` and `--phase-order` control the transform. Long
times require separate phase and frequency-grid checks. No signal is normalized
to another signal or to a sum rule.

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

Ordinary rows use NumPy matrix products. A tiny covariance entry alone does
not select slower arithmetic: an exponent check identifies potentially lost
underflow terms, then compares a conservative allowance for that loss with
each row's variances and real/imaginary cross covariance. Negligible tails
retain the ordinary path. Overflow risk, or underflow loss that cannot be
shown negligible (including active zero or cancelling results), selects exact
rational products of the supplied binary64 entries. The standard errors and
cross covariance are then rounded to binary64.
This slower path preserves representable standard errors even when their
squared variances underflow or overflow, and preserves tiny covariances that
large coefficients amplify. It uses Python's standard library, with no extra
dependency. It does not improve the precision of the input measurements or the
symmetry/positive-semidefinite check. Overflowing returned values still raise
`OverflowError`; values below the subnormal range can round to zero. Finite
means use the usual NumPy transform.

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
