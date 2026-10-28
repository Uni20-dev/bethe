// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <algorithm>
#include <uni20/core/math.hpp>
#include <uni20/core/runtime_precision.hpp>

namespace bethe::detail
{
// Fermionic TBA statistics without exp(+large) or 0*log(0). The particles
// themselves need not be fermions. Callers require finite energy and T>0.
template <uni20::Real Real> struct ThermalFactors
{
    Real log_weight{}, filling{}, entropy{};
};
template <uni20::Real Real> ThermalFactors<Real> thermal_factors(Real const& energy, Real const& temperature)
{
  if constexpr (uni20::has_runtime_precision_v<Real>)
  {
    auto const p = uni20::precision_of(temperature);
    if (p.is_exact() || uni20::precision_of(energy) != p)
      throw std::invalid_argument("thermal factors require matching finite precision");
  }
  Real const zero = uni20::scalar_like(temperature, 0);
  Real const r = uni20::math::abs(energy) / temperature, b = uni20::math::exp(-r);
  Real const log = uni20::math::log1p(b), small = b / (Real{1} + b);
  return {std::max(-energy, zero) + temperature * log, energy >= zero ? small : Real{1} / (Real{1} + b),
          log + (b == zero ? zero : r * small)};
}
} // namespace bethe::detail
