// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <uni20/core/math.hpp>
#include <vector>

namespace bethe::detail
{
// Native-precision, row-major, partial-pivot solve for a Newton correction.
// Failure is recoverable and never changes a process-global error policy.
// The RHS may be modified on failure; callers retain their current iterate.
template <uni20::Real Real> bool newton_step(std::vector<Real> a, std::vector<Real>& b)
{
  auto const n = b.size();
  if (n && n > std::numeric_limits<std::size_t>::max() / n) return false;
  if (a.size() != n * n) return false;
  Real scale{};
  for (Real v : a)
  {
    if (!uni20::isfinite(v)) return false;
    scale = std::max(scale, std::abs(v));
  }
  for (Real v : b)
    if (!uni20::isfinite(v)) return false;
  Real const floor = Real{64} * uni20::numeric_limits<Real>::epsilon() * scale;
  for (std::size_t k = 0; k < n; ++k)
  {
    auto p = k;
    for (std::size_t i = k + 1; i < n; ++i)
      if (std::abs(a[i * n + k]) > std::abs(a[p * n + k])) p = i;
    if (!uni20::isfinite(a[p * n + k]) || !(std::abs(a[p * n + k]) > floor)) return false;
    for (std::size_t j = k; j < n; ++j)
      std::swap(a[k * n + j], a[p * n + j]);
    std::swap(b[k], b[p]);
    for (std::size_t i = k + 1; i < n; ++i)
    {
      Real const f = a[i * n + k] / a[k * n + k];
      for (std::size_t j = k + 1; j < n; ++j)
        a[i * n + j] -= f * a[k * n + j];
      b[i] -= f * b[k];
    }
  }
  for (std::size_t i = n; i-- > 0;)
  {
    for (std::size_t j = i + 1; j < n; ++j)
      b[i] -= a[i * n + j] * b[j];
    b[i] /= a[i * n + i];
    if (!uni20::isfinite(b[i])) return false;
  }
  return true;
}
} // namespace bethe::detail
