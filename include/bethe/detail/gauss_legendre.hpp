// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <stdexcept>
#include <uni20/core/math.hpp>
#include <uni20/core/runtime_precision.hpp>
#include <vector>

namespace bethe::detail
{
template <uni20::Real Real> struct GaussLegendreRule
{
    std::vector<Real> x, w;
};
// Nodes and weights on [-1,1], shared by Fredholm solvers. Runtime scalars
// require a finite precision; neither constants nor stopping rules use defaults.
template <uni20::Real Real> GaussLegendreRule<Real> gauss_legendre(std::size_t n, uni20::Precision precision)
{
  if (precision.is_exact()) throw std::invalid_argument("quadrature requires finite working precision");
  auto const zero = uni20::scalar_at_precision<Real>(0, precision);
  auto const one = uni20::scalar_at_precision<Real>(1, precision);
  GaussLegendreRule<Real> result{std::vector<Real>(n), std::vector<Real>(n)};
  Real const pi = Real{4} * uni20::math::atan(one), eps = uni20::numeric_limits<Real>::epsilon(precision);
  for (std::size_t i = 0; i < (n + 1) / 2; ++i)
  {
    Real z = uni20::math::cos(pi * (Real(i) + Real{3} * one / Real{4}) / (Real(n) + one / Real{2}));
    Real derivative = zero;
    for (int iteration = 0; iteration < 128; ++iteration)
    {
      Real previous = one, current = z;
      for (std::size_t j = 2; j <= n; ++j)
      {
        Real const next = (Real(2 * j - 1) * z * current - Real(j - 1) * previous) / Real(j);
        previous = current;
        current = next;
      }
      derivative = Real(n) * (z * current - previous) / (z * z - Real{1});
      Real const step = current / derivative;
      if (uni20::math::abs(step) <= Real{4} * eps) break;
      z -= step;
    }
    Real const weight = Real{2} / ((Real{1} - z * z) * derivative * derivative);
    result.x[i] = -z;
    result.x[n - 1 - i] = z;
    result.w[i] = result.w[n - 1 - i] = weight;
  }
  return result;
}

template <uni20::Real Real>
  requires(!uni20::has_runtime_precision_v<Real>)
GaussLegendreRule<Real> gauss_legendre(std::size_t n)
{
  return gauss_legendre<Real>(n, uni20::native_precision<Real>());
}
} // namespace bethe::detail
