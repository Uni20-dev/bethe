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

// Overdetermined Newton correction without forming normal equations. Givens
// rotations retain native long-double support and return rank loss instead
// of invoking a terminal LAPACK/error-policy path. Row-major m x n, m>=n.
template <uni20::Real Real> bool least_squares_step(std::vector<Real> a, std::vector<Real>& b, std::size_t n)
{
  auto const m = b.size();
  if (m < n || (m && n > std::numeric_limits<std::size_t>::max() / m) || a.size() != m * n) return false;
  Real scale{};
  for (Real v : a)
  {
    if (!uni20::isfinite(v)) return false;
    scale = std::max(scale, std::abs(v));
  }
  for (Real v : b)
    if (!uni20::isfinite(v)) return false;
  for (std::size_t k = 0; k < n; ++k)
  {
    for (std::size_t i = k + 1; i < m; ++i)
    {
      Real const r = std::hypot(a[k * n + k], a[i * n + k]);
      if (r == Real{0}) continue;
      Real const c = a[k * n + k] / r, s = a[i * n + k] / r;
      a[k * n + k] = r;
      a[i * n + k] = Real{0};
      for (std::size_t j = k + 1; j < n; ++j)
      {
        Real const x = a[k * n + j], y = a[i * n + j];
        a[k * n + j] = c * x + s * y;
        a[i * n + j] = -s * x + c * y;
      }
      Real const x = b[k], y = b[i];
      b[k] = c * x + s * y;
      b[i] = -s * x + c * y;
    }
    if (!uni20::isfinite(a[k * n + k]) ||
        !(std::abs(a[k * n + k]) > Real{64} * uni20::numeric_limits<Real>::epsilon() * scale))
      return false;
  }
  b.resize(n);
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
