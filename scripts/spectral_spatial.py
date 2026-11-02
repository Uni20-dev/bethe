# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Explicit linear Fourier/cutoff operators for spectral measures (NumPy fp64).

Compose these with spectral_time_transform.statistics to keep full input
covariance, including sector/time/cutoff contrasts. Measures are always supplied
by the caller; no assumed endpoint values or fitted normalization is inserted.
"""
from numbers import Integral
import numpy as np


def _vector(value, name):
    if np.iscomplexobj(value):
        raise ValueError(f'{name} must be real')
    result = np.asarray(value, dtype=float)
    if result.ndim != 1 or not result.size or not np.isfinite(result).all():
        raise ValueError(f'invalid {name}')
    return result


def spatial_weights(momenta, separations, measure, *, even=False):
    """Weights for C(r)=sum_j measure[j]*exp(+i*q[j]*r)*m0(q[j]).

    q is in radians. For a periodic N-site chain supply measure=1/N at ALL
    N momenta. For full-zone continuum quadrature it includes dq/(2*pi).
    even=True instead uses cos(q*r); half-zone weights must include dq/pi,
    and the caller asserts that the physical spectrum is reflection-even.
    Input values are zeroth frequency moments for equal-time correlations,
    not first moments or unintegrated spectral densities.
    """
    q, mass = _vector(momenta, 'momenta'), _vector(measure, 'momentum measure')
    raw = np.asarray(separations, dtype=object)
    if (raw.ndim != 1 or not raw.size or any(not isinstance(r, Integral)
            or isinstance(r, (bool, np.bool_)) for r in raw)):
        raise ValueError('separations must be lattice integers')
    r = _vector(raw, 'separations')
    if any(int(rounded) != original for rounded, original in zip(r, raw)):
        raise ValueError('separation is not representable at binary64 precision')
    if mass.shape != q.shape or type(even) is not bool:
        raise ValueError('invalid momentum measure or reflection choice')
    with np.errstate(over='raise', invalid='raise'):
        try:
            phase = r[:, None]*q[None, :]
            result = (np.cos(phase) if even else np.exp(1j*phase))*mass
        except FloatingPointError as error:
            raise OverflowError('nonfinite spatial operator') from error
    if not np.isfinite(result).all():
        raise OverflowError('nonfinite spatial operator')
    return result


def cumulative_weights(coordinates, cutoffs, measure):
    """Return measure[j]*1(coordinates[j] <= cutoff[i]), inclusive at equality.

    A coordinate can be rapidity span max(beta)-min(beta), depth max|beta|,
    energy, or a precomputed bin upper edge. These are different selections.
    For raw Monte Carlo observations use the ORIGINAL sample denominator in
    measure (e.g. 1/N), including physical-zero observations. Selecting a
    smaller domain must not divide by the number of retained observations.
    For already integrated bin means supply ones. No overflow-bin filling or
    extrapolation beyond the supplied observations is performed.
    """
    x = _vector(coordinates, 'cutoff coordinates')
    limits = _vector(cutoffs, 'cutoffs')
    mass = _vector(measure, 'observation measure')
    if mass.shape != x.shape or np.any(limits[1:] <= limits[:-1]):
        raise ValueError('invalid observation measure or cutoff order')
    return (x[None, :] <= limits[:, None])*mass
