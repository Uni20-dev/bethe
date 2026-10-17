// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <algorithm>
#include <bethe/solver.hpp>
#include <optional>
#include <uni20/linalg/ops/slogdet.hpp>
#include <uni20/tensor/mdspec_tensor_view.hpp>
#include <vector>

namespace bethe::detail
{
template <uni20::Real Real> struct LogDeterminant
{
    std::optional<Real> log_absolute;
    Real pivot_margin = Real{1};
};

/// Column equilibration followed by Uni20's normally dispatched partial-pivot LU.
/// min|U_ii|/max(1,max|U_ij|) is a rejection diagnostic, NOT a condition number
/// or a certified error bound. An unresolved determinant is absent, not weight 0.
/// Runtime scalars require a finite-precision zero, including for an empty matrix.
template <uni20::Real Real, typename Scalar>
LogDeterminant<Real> log_absolute_determinant(std::vector<Scalar> matrix, std::size_t n, Real const& zero = Real{0})
{
  using uni20::math::abs;
  using uni20::math::log;
  if ((n && n > matrix.size() / n) || matrix.size() != n * n)
    throw std::invalid_argument("invalid determinant matrix size");
  auto const working_zero = precision_zero(zero);
  auto const one = uni20::scalar_like(zero, 1);
  LogDeterminant<Real> result{{}, one};
  CompensatedSum<Real> logarithm(working_zero);
  for (std::size_t j = 0; j < n; ++j)
  {
    Real scale = working_zero;
    for (std::size_t i = 0; i < n; ++i)
    {
      Real const magnitude = abs(matrix[i * n + j]);
      require_same_precision(magnitude, zero);
      if (!uni20::isfinite(magnitude)) return {{}, working_zero};
      scale = std::max(scale, magnitude);
    }
    if (scale == Real{0}) return {{}, working_zero};
    logarithm.add(log(scale));
    for (std::size_t i = 0; i < n; ++i)
      matrix[i * n + j] /= scale;
  }
  using Span = stdex::mdspan<Scalar, stdex::dextents<std::size_t, 2>, stdex::layout_right>;
  using ConstSpan = stdex::mdspan<Scalar const, stdex::dextents<std::size_t, 2>, stdex::layout_right>;
  uni20::MdspecTensorView<Span, ConstSpan, uni20::HostStorage> view(Span(matrix.data(), n, n),
                                                                    ConstSpan(matrix.data(), n, n));
  auto const lu = [&] {
#if UNI20_ENABLE_MPFR
    if constexpr (uni20::has_runtime_precision_v<Real>)
      return uni20::linalg::lu_factor_with_info(view, zero.precision());
    else
#endif
      return uni20::linalg::lu_factor_with_info(view);
  }();
  if (!lu.factor) return {{}, working_zero};
  auto const& u = lu.factor->packed();
  Real growth = one, smallest = one;
  for (std::size_t i = 0; i < n; ++i)
  {
    smallest = std::min(smallest, Real(abs(u[i, i])));
    for (std::size_t j = i; j < n; ++j)
      growth = std::max(growth, Real(abs(u[i, j])));
  }
  result.pivot_margin = smallest / growth;
  if (result.pivot_margin <= Real{128} * Real(n) * uni20::numeric_limits<Real>::epsilon(zero)) return result;
  auto const determinant = uni20::linalg::slogdet_with_info(*lu.factor);
  if (!determinant.value || !determinant.info.succeeded()) return result;
  logarithm.add(determinant.value->log_absolute);
  if (uni20::isfinite(logarithm.value())) result.log_absolute = logarithm.value();
  return result;
}
} // namespace bethe::detail
