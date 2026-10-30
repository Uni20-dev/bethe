# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Shared arbitrary-precision XXX pair factor from the Barnes-G identity.

The independent oscillatory-integral oracle remains in reference_xxx_four_spinon.
These special functions are ordinary mpmath evaluations, not interval bounds.
"""
from functools import lru_cache

import mpmath as mp


@lru_cache(maxsize=16)
def _log_constant(precision):
    # Cache by binary precision, never reuse a lower-precision constant.
    with mp.workprec(precision):
        return -2*mp.log(mp.barnesg(mp.mpf('.5')))


def log_barnes_constant():
    """log B, where B=Gamma(1/4)*|A_-(i*pi/2)|=G(1/2)^(-2)."""
    return _log_constant(mp.mp.prec)


def barnes_constant():
    return mp.exp(log_barnes_constant())


def log_transition_rate(rho):
    """log(exp(-I(rho))); rho is a rapidity difference divided by 2*pi."""
    rho = abs(mp.mpf(rho))
    if not mp.isfinite(rho):
        raise ValueError('require a finite real rapidity difference')
    if rho == 0:
        return mp.ninf
    log_ratio = (mp.log(abs(mp.barnesg(1+1j*rho)))
                 - mp.log(abs(mp.barnesg(mp.mpf('.5')+1j*rho))))
    return (2*log_barnes_constant()+mp.log(rho)
            + mp.log(mp.sinh(mp.pi*rho))+mp.log(mp.cosh(mp.pi*rho))
            - 2*mp.log(mp.pi)+4*log_ratio)


def transition_rate(rho):
    return mp.exp(log_transition_rate(rho))
