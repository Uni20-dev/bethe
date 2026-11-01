# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Finite-window spectral quadrature and propagation of sampling covariance.

NumPy binary64 analysis, independent of any model, archive or plotting code.
No endpoint extrapolation, broadening or amplitude matching is performed.
"""

from numbers import Integral

import numpy as np
from numpy.polynomial.legendre import leggauss


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _real_array(value, name):
    _require(not np.iscomplexobj(value), f"{name} must be real")
    result = np.asarray(value, dtype=float)
    _require(np.isfinite(result).all(), f"{name} must be finite")
    return result


def _integer(value, minimum, name):
    _require(isinstance(value, Integral) and not isinstance(value, (bool, np.bool_))
             and value >= minimum, f"invalid {name}")
    return int(value)


def weights(frequencies, times, *, offsets=None, indices=None, order=24,
            window=None, energy_power=0):
    """Return W for C(t) = W @ S, including the measure dω/(2π).

    By default S(ω) is interpolated linearly between selected frequency nodes.
    With positive offsets δ=ω-ω_lower, interpolate G(x)=2*x*S(ω_lower+x*x)
    linearly in x=sqrt(δ), preserving supplied offsets at exact grid nodes.
    The phase is integrated by an `order`-point Gauss-Legendre rule per cell.
    Increase order independently of refining the spectral grid at long times.

    `window` must be covered by the selected nodes; omitted columns get zero
    weight. `energy_power` multiplies the SAME interpolant by ω**energy_power,
    rather than interpolating energy-weighted nodal values. A zero-time call
    with energy_power=1 therefore measures the first frequency moment.
    """
    w = _real_array(frequencies, "frequencies")
    t = _real_array(times, "times")
    _require(w.ndim == t.ndim == 1 and len(w) >= 2 and len(t) > 0
             and np.all(w[1:] > w[:-1]), "invalid frequencies or times")
    if indices is None:
        ids = np.arange(len(w))
    else:
        # Preserve original element types: np.asarray([0, True]) would
        # otherwise erase the Boolean by converting it to integer 1.
        selected = np.asarray(indices, dtype=object)
        _require(selected.ndim == 1 and len(selected) >= 2,
                 "invalid selected nodes")
        values = [_integer(value, 0, "selected node") for value in selected]
        _require(all(left < right for left, right in zip(values, values[1:]))
                 and values[-1] < len(w), "invalid selected nodes")
        ids = np.asarray(values, dtype=np.intp)
    order = _integer(order, 2, "phase quadrature order")
    energy_power = _integer(energy_power, 0, "energy power")
    if offsets is None:
        x, jac, threshold = w, np.ones(len(w)), None
    else:
        delta = _real_array(offsets, "offsets")
        _require(delta.shape == w.shape and (delta > 0).all()
                 and np.all(delta[1:] > delta[:-1]), "invalid offsets")
        threshold = w[0] - delta[0]
        scale = np.maximum(1, np.maximum(np.abs(w), np.abs(delta)))
        _require(np.isfinite(threshold) and
                 np.all(np.abs((w - delta) - threshold) <= 8*np.finfo(float).eps*scale),
                 "inconsistent threshold")
        x = np.sqrt(delta)
        jac = 2*x
        _require(np.all(x[1:] > x[:-1]), "unresolved threshold coordinates")
    bounds = w[ids[[0, -1]]] if window is None else _real_array(window, "window")
    _require(bounds.shape == (2,), "window must contain two bounds")
    lo, hi = bounds
    _require(w[ids[0]] <= lo < hi <= w[ids[-1]], "window requires sampled coverage")
    if threshold is None:
        xl, xh = lo, hi
    else:
        def coordinate(value):
            match = np.flatnonzero(w == value)
            return x[match[0]] if len(match) else np.sqrt(value - threshold)
        xl, xh = coordinate(lo), coordinate(hi)
    g, h = leggauss(order)
    u, h = (g+1)/2, h/2
    result = np.zeros((len(t), len(w)), complex)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            for i, j in zip(ids[:-1], ids[1:]):
                a, b = max(x[i], xl), min(x[j], xh)
                if b <= a:
                    continue
                nodes = a+(b-a)*u
                hat = (nodes-x[i])/(x[j]-x[i])
                energy = nodes if threshold is None else threshold+nodes*nodes
                phase = np.exp(-1j*t[:, None]*energy[None, :])
                if energy_power:
                    phase *= energy[None, :]**energy_power
                factor = (b-a)/(2*np.pi)
                result[:, i] += factor*jac[i]*(phase@(h*(1-hat)))
                result[:, j] += factor*jac[j]*(phase@(h*hat))
    except FloatingPointError as error:
        raise OverflowError("nonfinite spectral quadrature arithmetic") from error
    if not np.isfinite(result).all():
        raise OverflowError("nonfinite spectral quadrature weights")
    return result


def statistics(operator, mean, covariance):
    """Transform real sample means with their full covariance-of-the-mean.

    Return real/imaginary means, standard errors and their same-row covariance
    as lists. Off-diagonal input covariance is retained. For shared-data
    differences pass W1-W2; for contrasts between times form the corresponding
    rows before calling. No assumption of independent frequencies/sectors is
    made. Sampling uncertainty does not include grid, tail or phase errors.

    Covariance is checked for symmetry and positive semidefiniteness (an O(n^3)
    eigenvalue check); roundoff-size violations are tolerated, not a physical
    negative variance. Inputs are not mutated.
    """
    W = np.asarray(operator, dtype=complex)
    m = _real_array(mean, "mean")
    C = _real_array(covariance, "covariance")
    _require(W.ndim == 2 and min(W.shape) > 0 and np.isfinite(W).all()
             and m.shape == (W.shape[1],) and C.shape == (len(m), len(m)),
             "invalid transform statistics")
    scale = float(np.max(np.abs(C)))
    tol = 64*np.finfo(float).eps*len(m)
    if scale:
        normalized = C/scale
        _require(np.max(np.abs(normalized-normalized.T)) <= tol, "asymmetric covariance")
        normalized = (normalized + normalized.T)/2
        _require(np.linalg.eigvalsh(normalized)[0] >= -tol,
                 "covariance must be positive semidefinite")
        # Average the original entries: rescaling back from normalized can
        # erase tiny covariances, while halving both entries first erases
        # equal subnormals. Split only sums that could overflow.
        small = np.abs(C) <= np.finfo(float).max/2
        safe_sum = small & small.T
        symmetric = np.empty_like(C)
        symmetric[safe_sum] = (C[safe_sum] + C.T[safe_sum])/2
        symmetric[~safe_sum] = C[~safe_sum]/2 + C.T[~safe_sum]/2
        C = symmetric
    try:
        with np.errstate(over="raise", invalid="raise"):
            y = W@m
            vr = np.einsum("ij,ij->i", W.real@C, W.real)
            vi = np.einsum("ij,ij->i", W.imag@C, W.imag)
            ri = np.einsum("ij,ij->i", W.real@C, W.imag)
            # Scale the negative-variance allowance by this observable's
            # absolute terms, so small physical variances are not zeroed by
            # an unrelated, fixed absolute tolerance.
            allowances = [tol*np.einsum("ij,ij->i", np.abs(part)@np.abs(C), np.abs(part))
                          for part in (W.real, W.imag)]
    except FloatingPointError as error:
        raise OverflowError("nonfinite transformed statistics") from error
    if not all(np.isfinite(x).all() for x in (y, vr, vi, ri, *allowances)):
        raise OverflowError("nonfinite transformed statistics")
    _require(np.all(vr >= -allowances[0]) and np.all(vi >= -allowances[1]),
             "negative transformed variance")
    return dict(real=y.real.tolist(), imag=y.imag.tolist(),
                real_standard_error=np.sqrt(np.maximum(0, vr)).tolist(),
                imag_standard_error=np.sqrt(np.maximum(0, vi)).tolist(),
                real_imag_covariance=ri.tolist())


def coarser_indices(count, stride):
    """Select every stride-th node and always retain both sampled endpoints."""
    count = _integer(count, 2, "node count")
    stride = _integer(stride, 1, "node stride")
    ids = list(range(0, count, stride))
    if ids[-1] != count-1:
        ids.append(count-1)
    return ids
