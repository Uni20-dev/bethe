# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Integration-ready RESEARCH vertex density with separate error diagnostics.

This is the coefficient W_n against labelled product(dp) and the two
unintegrated conservation deltas, in Szz, J=1 conventions. Not a DSF value.
No fitted normalization, finite-chain extrapolation or interval certification.
"""
from dataclasses import dataclass
from itertools import combinations

import mpmath as mp

from reference_xxx_six_spinon_tail import accelerated_contours
from reference_xxx_six_spinon_vertex import validate
from reference_xxx_transition_rate import log_barnes_constant, log_transition_rate


@dataclass(frozen=True)
class DensityEstimate:
    value: object
    contour_tail_error: object
    arithmetic_indicator: object
    log_scale: object
    terms: int
    order: int
    dps: int


@dataclass(frozen=True)
class CheckedDensity:
    estimate: DensityEstimate
    precision_difference: object
    scalar_precision_difference: object
    arithmetic_screen: object
    requested_error: object
    accepted: bool
    reason: str
    base_dps: int

    @property
    def value(self):
        return self.estimate.value

    @property
    def error_screen(self):
        return self.estimate.contour_tail_error+self.arithmetic_screen


def log_density_scale(beta):
    """log of the pair product times the coefficient multiplying sum|G|^2.

    Full labelled momentum cube: one n! and the Szz/frequency conversion.
    Coefficients agree with the independent finite-q derivation in the note.
    """
    beta = validate(beta, 1)
    n = len(beta)
    if n == 2:
        coefficient = -mp.log(4)
    elif n == 4:
        coefficient = -15*mp.log(2)-mp.log(24)-14*mp.log(mp.pi)-8*log_barnes_constant()
    else:
        coefficient = -41*mp.log(2)-mp.log(720)-42*mp.log(mp.pi)-24*log_barnes_constant()
    return coefficient+mp.fsum(log_transition_rate((a-b)/(2*mp.pi)) for a, b in combinations(beta, 2))


def density_at_precision(beta, terms=32, order=12):
    """Evaluate at the caller's precision, retaining tail/roundoff separately.

    Raises for unsupported inputs or a tail outside its majorant domain.
    The arithmetic quantity is a heuristic indicator, NOT a bound on Barnes G.
    Use checked_density and independent precision repeats to assess arithmetic.
    """
    beta = validate(beta, terms)
    contour = accelerated_contours(beta, terms, order)
    log_scale = log_density_scale(beta)
    scale = mp.exp(log_scale)
    value = scale*contour.norm
    pair_count = len(beta)*(len(beta)-1)//2
    arithmetic = (scale*contour.norm_arithmetic_indicator
                  + 64*mp.eps*(1+pair_count+abs(log_scale))*abs(value))
    if not all(mp.isfinite(x) and x >= 0 for x in (value, contour.norm_error, arithmetic)):
        raise ArithmeticError('nonfinite or negative pointwise density/diagnostic')
    return DensityEstimate(value, scale*contour.norm_error, arithmetic,
                           log_scale, terms, order, mp.mp.dps)


def checked_density(beta, terms=32, order=12, *, dps=35, guard_digits=20, rtol='1e-8', atol='0'):
    """Repeat the entire density at two precisions and return an explicit screen.

    Decimal-string inputs are reconstructed independently at each precision.
    A rejected result remains available for diagnostics, never becomes zero.
    Tails at fixed terms/order are NOT tested by the precision difference.
    """
    if (isinstance(dps, bool) or not isinstance(dps, int) or dps < 25
            or isinstance(guard_digits, bool) or not isinstance(guard_digits, int) or guard_digits < 10):
        raise ValueError('require integer dps>=25 and guard_digits>=10')
    beta = tuple(beta)
    with mp.workdps(dps+guard_digits):
        rtol, atol = mp.mpf(rtol), mp.mpf(atol)
        if not mp.isfinite(rtol) or not mp.isfinite(atol) or min(rtol, atol) < 0 or rtol+atol <= 0:
            raise ValueError('require finite nonnegative tolerances, at least one positive')
    with mp.workdps(dps):
        low = density_at_precision(beta, terms, order)
    with mp.workdps(dps+guard_digits):
        high = density_at_precision(beta, terms, order)
        difference = abs(high.value-low.value)
        scalar_difference = abs(mp.expm1(low.log_scale-high.log_scale))*abs(high.value)
        arithmetic = low.arithmetic_indicator+high.arithmetic_indicator+difference+scalar_difference
        tolerance = atol+rtol*abs(high.value)
        accepted = high.contour_tail_error+arithmetic <= tolerance
        reason = ('accepted numerical screens (not an interval certificate)' if accepted
                  else 'contour-tail screen failed' if high.contour_tail_error > tolerance
                  else 'arithmetic/precision screen failed')
        return CheckedDensity(high, difference, scalar_difference, arithmetic, tolerance,
                              accepted, reason, dps)
