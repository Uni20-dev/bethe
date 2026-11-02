// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/solver.hpp>

namespace bethe::detail
{
/// Internal dimensionless RF(x,y,1) used by the XXZ inverse-dn maps.
/// Duplication followed by the degree-five expansion (DLMF 19.36).
/// Requires 0<=x,y<=1 and at least one positive argument.
/// Keeping one argument exactly unity avoids general-range rescaling/overflow.
template <uni20::Real Real> Real carlson_rf_unit(Real x, Real y)
{
  Real const zero = precision_zero(x), one = uni20::scalar_like(zero, 1);
  require_same_precision(y, zero);
  Real z = one;
  if (x > one || y > one)
    throw std::invalid_argument("this internal RF evaluator requires dimensionless arguments <=1");
  if (!uni20::isfinite(x) || !uni20::isfinite(y) || !uni20::isfinite(z) || x < Real{0} || y < Real{0} || z < Real{0} ||
      x + y <= Real{0} || x + z <= Real{0} || y + z <= Real{0})
    throw std::invalid_argument("RF requires finite nonnegative arguments, at most one zero");
  Real const tolerance = uni20::math::pow(uni20::numeric_limits<Real>::epsilon(zero), one / Real{6}) / Real{8};
  for (std::size_t i = 0; i < 128 + uni20::precision_of(zero).bit_count() / 2; ++i)
  {
    Real const mean = (x + y + z) / Real{3};
    Real const dx = (mean - x) / mean, dy = (mean - y) / mean, dz = (mean - z) / mean;
    if (std::max({uni20::math::abs(dx), uni20::math::abs(dy), uni20::math::abs(dz)}) < tolerance)
    {
      Real const e2 = dx * dy - dz * dz, e3 = dx * dy * dz;
      return (Real{1} + (e2 / Real{24} - Real{1} / Real{10} - Real{3} * e3 / Real{44}) * e2 + e3 / Real{14}) /
             uni20::math::sqrt(mean);
    }
    Real const sx = uni20::math::sqrt(x), sy = uni20::math::sqrt(y), sz = uni20::math::sqrt(z);
    Real const lambda = sx * sy + sx * sz + sy * sz;
    x = (x + lambda) / Real{4};
    y = (y + lambda) / Real{4};
    z = (z + lambda) / Real{4};
  }
  throw std::runtime_error("RF duplication did not converge");
}
} // namespace bethe::detail
