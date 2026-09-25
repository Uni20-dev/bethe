// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <complex>
#include <uni20/core/math.hpp>

namespace bethe::detail
{
/// Principal log(1+z), retaining small complex increments without forming 1+z.
template <uni20::Real Real> std::complex<Real> complex_log1p(std::complex<Real> z)
{
  if (std::abs(z) > Real{0.25}) return std::log(Real{1} + z);
  Real const a = z.real(), b = z.imag();
  return {std::log1p(Real{2} * a + a * a + b * b) / Real{2}, std::atan2(b, Real{1} + a)};
}
} // namespace bethe::detail
