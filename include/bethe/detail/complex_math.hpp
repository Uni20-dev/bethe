// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <complex>
#include <uni20/core/math.hpp>

namespace bethe::detail
{
/// Principal log(1+z), retaining small complex increments without forming 1+z.
template <uni20::Complex C> C complex_log1p(C z)
{
  using Real = decltype(z.real());
  if (uni20::math::abs(z) > Real{1} / Real{4}) return uni20::math::log(Real{1} + z);
  Real const a = z.real(), b = z.imag();
  return {uni20::math::log1p(Real{2} * a + a * a + b * b) / Real{2}, uni20::math::atan2(b, Real{1} + a)};
}
} // namespace bethe::detail
