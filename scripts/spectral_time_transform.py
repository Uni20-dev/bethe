# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Finite-window spectral quadrature and propagation of sampling covariance.

NumPy binary64 analysis, independent of any model, archive or plotting code.
No endpoint extrapolation, broadening or amplitude matching is performed.
"""

from fractions import Fraction
import math
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


def _covariance_range_risks(W, C):
    """Separate possible underflow from possible overflow in quadratic forms."""
    nonzero = C[C != 0]
    if not len(nonzero):
        return np.zeros(len(W), dtype=bool), np.zeros(len(W), dtype=bool)
    ce = np.frexp(nonzero)[1]
    parts = np.stack((W.real, W.imag))
    exponents = np.frexp(parts)[1]
    smallest = np.min(np.where(parts != 0, exponents, 10000), axis=(0, 2))
    largest = np.max(np.where(parts != 0, exponents, -10000), axis=(0, 2))
    # Leave headroom for matrix products, n^2-term allowances and sqrt.
    # This is an arithmetic dispatch check, not a statistical tolerance.
    headroom = 2*(W.shape[1]-1).bit_length()
    return (2*smallest + int(ce.min()) < -900,
            2*largest + int(ce.max()) + headroom > 900)


def _underflow_negligible(W, vr, vi, ri):
    """Check lost underflow terms against each computed observable.

    In (a @ C) @ b, gradual-underflow errors contribute at most about
    2*n*eta*(1+sum(abs(b))), where eta is the minimum binary64 subnormal.
    Use 4*n*(n+1)*eta*max(1,max(abs(b))) as a conservative allowance for
    those errors and their subsequent rounding. Compare exponent bounds,
    avoiding underflow in the allowance itself, with extra mantissa headroom.
    This bounds underflow loss, not the usual relative rounding error or
    covariance-conditioning error. Active zero/cancelling results stay exact.
    """
    n = W.shape[1]
    common = -1074 + 2 + (n-1).bit_length() + n.bit_length()

    def safe(value, right, active):
        largest = np.frexp(np.max(np.abs(right), axis=1))[1]
        lost_exponent = common + np.maximum(0, largest)
        value_exponent = np.frexp(value)[1]
        return ~active | ((value != 0) & (value_exponent > lost_exponent + 56))

    real_active = np.any(W.real != 0, axis=1)
    imag_active = np.any(W.imag != 0, axis=1)
    return (safe(vr, W.real, real_active) & safe(vi, W.imag, imag_active)
            & safe(ri, W.imag, real_active & imag_active))


def _fraction_sqrt(value):
    """Round sqrt of a nonnegative rational without first rounding its square."""
    if not value:
        return 0.
    n, d = value.numerator, value.denominator
    exponent = (n.bit_length()-d.bit_length())//2
    scaled = Fraction(n, d << (2*exponent)) if exponent >= 0 else Fraction(n << (-2*exponent), d)
    return math.ldexp(math.sqrt(float(scaled)), exponent)


def _wide_covariance_row(row, covariance, tolerance):
    """Exact products of binary64 inputs; used only outside the fast range.

    Scaling the whole covariance matrix can erase tiny entries before a large
    operator amplifies them. Rational products also preserve exact correlated
    cancellations without requiring the variance/allowance to fit in binary64.
    """
    real = list(map(Fraction, row.real.tolist()))
    imag = list(map(Fraction, row.imag.tolist()))

    def bilinear(left, right):
        total, absolute = Fraction(), Fraction()
        for a, c_row in zip(left, covariance):
            if a:
                for c, b in zip(c_row, right):
                    if c and b:
                        term = a*c*b
                        total += term
                        absolute += abs(term)
        return total, absolute

    vr, ar = bilinear(real, real)
    vi, ai = bilinear(imag, imag)
    ri, _ = bilinear(real, imag)
    _require(vr >= -tolerance*ar and vi >= -tolerance*ai,
             "negative transformed variance")
    return (_fraction_sqrt(max(0, vr)), _fraction_sqrt(max(0, vi)), float(ri))


def statistics(operator, mean, covariance):
    """Transform real sample means with their full covariance-of-the-mean.

    Return real/imaginary means, standard errors and their same-row covariance
    as lists. Off-diagonal input covariance is retained. For shared-data
    differences pass W1-W2; for contrasts between times form the corresponding
    rows before calling. No assumption of independent frequencies/sectors is
    made. Sampling uncertainty does not include grid, tail or phase errors.

    Covariance is checked for symmetry and positive semidefiniteness (an O(n^3)
    eigenvalue check); roundoff-size violations are tolerated, not a physical
    negative variance. Rows at risk of overflow, or where underflow loss is
    not negligible for a computed quadratic form, use exact products of the
    input binary64 entries before rounding outputs. Inputs are not mutated.
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
    underflow_risk, wide = _covariance_range_risks(W, C)
    ordinary_indices = np.flatnonzero(~wide)
    try:
        with np.errstate(over="raise", invalid="raise"):
            y = W@m
            ordinary = W[ordinary_indices] if wide.any() else W
            vr = np.einsum("ij,ij->i", ordinary.real@C, ordinary.real)
            vi = np.einsum("ij,ij->i", ordinary.imag@C, ordinary.imag)
            ri = np.einsum("ij,ij->i", ordinary.real@C, ordinary.imag)
            # The allowance follows this observable's absolute terms, not an
            # unrelated fixed scale. Extreme rows use exact products below.
            allowances = [tol*np.einsum("ij,ij->i", np.abs(part)@np.abs(C), np.abs(part))
                          for part in (ordinary.real, ordinary.imag)]
        if not all(np.isfinite(x).all() for x in (y, vr, vi, ri, *allowances)):
            raise OverflowError("nonfinite transformed statistics")
        sensitive = underflow_risk[ordinary_indices]
        if sensitive.any():
            keep = np.ones(len(ordinary), dtype=bool)
            keep[sensitive] = _underflow_negligible(
                ordinary[sensitive], vr[sensitive], vi[sensitive], ri[sensitive])
            wide[ordinary_indices[~keep]] = True
            vr, vi, ri = vr[keep], vi[keep], ri[keep]
            allowances = [a[keep] for a in allowances]
        _require(np.all(vr >= -allowances[0]) and np.all(vi >= -allowances[1]),
                 "negative transformed variance")
        real_error, imag_error, covariance_ri = (np.zeros(len(W)) for _ in range(3))
        real_error[~wide] = np.sqrt(np.maximum(0, vr))
        imag_error[~wide] = np.sqrt(np.maximum(0, vi))
        covariance_ri[~wide] = ri
        if wide.any():
            exact = [list(map(Fraction, row)) for row in C.tolist()]
            tolerance = Fraction(tol)
            for i in np.flatnonzero(wide):
                real_error[i], imag_error[i], covariance_ri[i] = _wide_covariance_row(W[i], exact, tolerance)
    except (FloatingPointError, OverflowError) as error:
        raise OverflowError("nonfinite transformed statistics") from error
    return dict(real=y.real.tolist(), imag=y.imag.tolist(),
                real_standard_error=real_error.tolist(),
                imag_standard_error=imag_error.tolist(),
                real_imag_covariance=covariance_ri.tolist())


def coarser_indices(count, stride):
    """Select every stride-th node and always retain both sampled endpoints."""
    count = _integer(count, 2, "node count")
    stride = _integer(stride, 1, "node stride")
    ids = list(range(0, count, stride))
    if ids[-1] != count-1:
        ids.append(count-1)
    return ids
