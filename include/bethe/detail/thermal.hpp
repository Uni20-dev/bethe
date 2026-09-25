// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <algorithm>
#include <uni20/core/math.hpp>

namespace bethe::detail
{
// Fermionic TBA statistics without exp(+large) or 0*log(0). The particles
// themselves need not be fermions. Callers require finite energy and T>0.
template <uni20::Real Real> struct ThermalFactors
{
    Real log_weight{}, filling{}, entropy{};
};
template <uni20::Real Real> ThermalFactors<Real> thermal_factors(Real energy, Real temperature)
{
  Real const r = std::abs(energy) / temperature, b = std::exp(-r);
  Real const log = std::log1p(b), small = b / (Real{1} + b);
  return {std::max(-energy, Real{0}) + temperature * log, energy >= Real{0} ? small : Real{1} / (Real{1} + b),
          log + (b == Real{0} ? Real{0} : r * small)};
}
} // namespace bethe::detail
