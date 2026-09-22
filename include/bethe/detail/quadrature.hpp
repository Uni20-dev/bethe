// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <algorithm>
#include <array>
#include <bethe/solver.hpp>
#include <cmath>
#include <vector>

namespace bethe::detail
{
// Native-precision double-exponential quadrature on a finite interval.
// Independent successive meshes supply an error ESTIMATE, not an interval
// certificate. The caller supplies any infinite-domain tail bound separately.
template <uni20::Real Real, std::size_t N> struct QuadratureResult
{
    std::array<Real, N> value{}, error{};
    bool converged = false;
};
template <uni20::Real Real, std::size_t N, typename Function>
QuadratureResult<Real, N> tanh_sinh(Function const& f, Real lo, Real hi, Real relative_tolerance,
                                    std::size_t& evaluations, std::size_t max_evaluations, std::size_t max_levels,
                                    std::array<Real, N> const& absolute_scales = {})
{
  QuadratureResult<Real, N> out;
  if (lo == hi)
  {
    out.converged = true;
    return out;
  }
  Real const pi = Real{4} * std::atan(Real{1}), width = hi - lo;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  std::array<Real, N> previous{};
  unsigned stable = 0;
  for (std::size_t level = 0; level < max_levels; ++level)
  {
    Real const h = std::ldexp(Real{1}, -int(level));
    std::array<CompensatedSum<Real>, N> sum, absolute;
    auto sample = [&](Real x, Real weight) {
      if (evaluations == max_evaluations) return false;
      ++evaluations;
      auto const y = f(x);
      for (std::size_t j = 0; j < N; ++j)
      {
        if (!uni20::isfinite(y[j])) return false;
        sum[j].add(weight * y[j]);
        absolute[j].add(std::abs(weight * y[j]));
      }
      return true;
    };
    if (!sample(lo + width / Real{2}, width * pi / Real{4})) return out;
    // Stop at an endpoint distance O(epsilon^2). The right endpoint can
    // round to hi earlier; it is never evaluated at a singular endpoint.
    Real const limit = -Real{2} * std::log(eps) + Real{4};
    for (std::size_t k = 1;; ++k)
    {
      Real const t = h * Real(k), argument = pi * std::sinh(t);
      if (argument > limit) break;
      Real const z = std::exp(-argument), distance = width * z / (Real{1} + z);
      Real const weight = width * pi * std::cosh(t) * z / ((Real{1} + z) * (Real{1} + z));
      Real const left = lo + distance, right = hi - distance;
      if (left > lo && !sample(left, weight)) return out;
      if (right < hi && !sample(right, weight)) return out;
    }
    bool ok = level >= 2;
    for (std::size_t j = 0; j < N; ++j)
    {
      out.value[j] = h * sum[j].value();
      Real const roundoff = Real{8} * eps * h * absolute[j].value();
      out.error[j] = std::max(std::abs(out.value[j] - previous[j]), roundoff);
      ok = ok && out.error[j] <= relative_tolerance * std::max(std::abs(out.value[j]), absolute_scales[j]);
    }
    stable = ok ? stable + 1 : 0;
    if (stable == 2)
    {
      out.converged = true;
      return out;
    }
    previous = out.value;
  }
  return out;
}
} // namespace bethe::detail
