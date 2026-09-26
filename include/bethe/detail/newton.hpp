// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <uni20/core/math.hpp>
#include <uni20/linalg/ops/linear_solve.hpp>
#include <uni20/tensor/mdspec_tensor_view.hpp>
#include <vector>

namespace bethe::detail
{
// Native-precision square solve for a row-major Newton correction.
// Failure is recoverable and never changes a process-global error policy.
// The RHS may be modified on failure; callers retain their current iterate.
template <uni20::Real Real> bool newton_step(std::vector<Real> a, std::vector<Real>& b)
{
  auto const n = b.size();
  if (n && n > std::numeric_limits<std::size_t>::max() / n) return false;
  if (a.size() != n * n) return false;
  // Reuse the owned coefficient copy as column-major workspace, allowing
  // LAPACK/MPLAPACK dispatch without another allocation. Long double uses
  // Uni20's native CPU fallback. Neither view outlives its vector storage.
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i + 1; j < n; ++j)
      std::swap(a[i * n + j], a[j * n + i]);
  using Span = stdex::mdspan<Real, stdex::dextents<std::size_t, 2>, stdex::layout_left>;
  using ConstSpan = stdex::mdspan<Real const, stdex::dextents<std::size_t, 2>, stdex::layout_left>;
  using View = uni20::MdspecTensorView<Span, ConstSpan, uni20::HostStorage>;
  View matrix(Span(a.data(), n, n), ConstSpan(a.data(), n, n));
  View rhs(Span(b.data(), n, 1), ConstSpan(b.data(), n, 1));
  return uni20::linalg::solve_inplace_with_info(
             matrix, rhs, {.relative_pivot_tolerance = Real{64} * uni20::numeric_limits<Real>::epsilon()})
      .succeeded();
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
