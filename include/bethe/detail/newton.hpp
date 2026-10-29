// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <uni20/core/math.hpp>
#include <uni20/linalg/ops/least_squares.hpp>
#include <uni20/linalg/ops/linear_solve.hpp>
#include <uni20/tensor/mdspec_tensor_view.hpp>
#include <utility>
#include <vector>

namespace bethe::detail
{
// Square solve for a row-major Newton correction at explicit precision.
// Failure is recoverable and never changes a process-global error policy.
// The RHS may be modified on failure; callers retain their current iterate.
template <uni20::Real Real> bool newton_step(std::vector<Real> a, std::vector<Real>& b, uni20::Precision precision)
{
  auto const n = b.size();
  if (n && n > std::numeric_limits<std::size_t>::max() / n) return false;
  if (a.size() != n * n) return false;
  // Reuse the owned coefficient copy as column-major workspace, allowing
  // LAPACK/MPLAPACK dispatch without another allocation, with Uni20's
  // CPU fallback where needed. Neither view outlives its vector storage.
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i + 1; j < n; ++j)
      std::swap(a[i * n + j], a[j * n + i]);
  using Span = stdex::mdspan<Real, stdex::dextents<std::size_t, 2>, stdex::layout_left>;
  using ConstSpan = stdex::mdspan<Real const, stdex::dextents<std::size_t, 2>, stdex::layout_left>;
  using View = uni20::MdspecTensorView<Span, ConstSpan, uni20::HostStorage>;
  View matrix(Span(a.data(), n, n), ConstSpan(a.data(), n, n));
  View rhs(Span(b.data(), n, 1), ConstSpan(b.data(), n, 1));
  return uni20::linalg::solve_inplace_with_info(
             matrix, rhs, precision,
             {.relative_pivot_tolerance = Real{64} * uni20::numeric_limits<Real>::epsilon(precision)})
      .succeeded();
}

template <uni20::Real Real>
  requires(!uni20::has_runtime_precision_v<Real>)
bool newton_step(std::vector<Real> a, std::vector<Real>& b)
{
  return newton_step(std::move(a), b, uni20::native_precision<Real>());
}

/// Overdetermined Newton correction for row-major m x n coefficients, m>=n.
/// Uni20 owns QR and backend selection; Bethe retains the 64*epsilon diagonal
/// rejection policy (not a rank or condition estimate). No normal equations.
/// Preserve a always, and replace b by the n-component correction only on
/// success. Invalid shapes and numerical failures return false. Invalid or
/// unavailable precision and allocation failures remain ordinary errors.
/// An explicit finite precision permits runtime inputs, converted by Uni20 at
/// that precision; callers must not infer recovered accuracy from conversion.
template <uni20::Real Real>
bool least_squares_step(std::vector<Real> const& a, std::vector<Real>& b, std::size_t n, uni20::Precision precision)
{
  auto const m = b.size();
  if (m < n || (m && n > std::numeric_limits<std::size_t>::max() / m) || a.size() != m * n) return false;
  using Span = stdex::mdspan<Real const, stdex::dextents<std::size_t, 2>, stdex::layout_right>;
  using View = uni20::MdspecTensorView<Span, Span, uni20::HostStorage>;
  View const matrix(Span(a.data(), m, n), Span(a.data(), m, n));
  View const rhs(Span(b.data(), m, 1), Span(b.data(), m, 1));
  auto const result = uni20::linalg::least_squares_with_info(
      matrix, rhs, precision,
      {.relative_diagonal_tolerance = Real{64} * uni20::numeric_limits<Real>::epsilon(precision)});
  if (!result.info.succeeded()) return false;
  std::vector<Real> correction;
  correction.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
    correction.push_back((*result.solution)[i, 0]);
  b = std::move(correction);
  return true;
}

/// Fixed-precision convenience; runtime arithmetic requires an explicit request.
template <uni20::Real Real>
  requires(!uni20::has_runtime_precision_v<Real>)
bool least_squares_step(std::vector<Real> const& a, std::vector<Real>& b, std::size_t n)
{
  return least_squares_step(a, b, n, uni20::native_precision<Real>());
}
} // namespace bethe::detail
