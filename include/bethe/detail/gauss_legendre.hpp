// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <uni20/core/math.hpp>
#include <vector>

namespace bethe::detail
{
template <uni20::Real Real> struct GaussLegendreRule
{
    std::vector<Real> x, w;
};
// Native-real nodes and weights on [-1,1], shared by Fredholm solvers.
template <uni20::Real Real> GaussLegendreRule<Real> gauss_legendre(std::size_t n)
{
  GaussLegendreRule<Real> result{std::vector<Real>(n), std::vector<Real>(n)};
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t i = 0; i < (n + 1) / 2; ++i)
  {
    Real z = std::cos(pi * (Real(i) + Real{.75}) / (Real(n) + Real{.5}));
    Real derivative{};
    for (int iteration = 0; iteration < 128; ++iteration)
    {
      Real previous{1}, current = z;
      for (std::size_t j = 2; j <= n; ++j)
      {
        Real const next = (Real(2 * j - 1) * z * current - Real(j - 1) * previous) / Real(j);
        previous = current;
        current = next;
      }
      derivative = Real(n) * (z * current - previous) / (z * z - Real{1});
      Real const step = current / derivative;
      if (std::abs(step) <= Real{4} * eps) break;
      z -= step;
    }
    Real const weight = Real{2} / ((Real{1} - z * z) * derivative * derivative);
    result.x[i] = -z;
    result.x[n - 1 - i] = z;
    result.w[i] = result.w[n - 1 - i] = weight;
  }
  return result;
}
} // namespace bethe::detail
