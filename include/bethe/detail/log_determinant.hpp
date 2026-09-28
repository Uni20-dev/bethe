// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <algorithm>
#include <bethe/solver.hpp>
#include <complex>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

namespace bethe::detail
{
template <uni20::Real Real> struct LogDeterminant
{
    std::optional<Real> log_absolute;
    Real pivot_margin = Real{1};
};

/// Small dense determinant, without forming the determinant itself. Column
/// equilibration followed by complete-pivot LU; scalar arithmetic stays native.
/// The smallest pivot divided by the largest encountered entry is a rejection
/// diagnostic, NOT a condition number or a certified relative error bound.
/// A zero/unresolved determinant is absent, not silently interpreted as weight 0.
template <uni20::Real Real, typename Scalar>
LogDeterminant<Real> log_absolute_determinant(std::vector<Scalar> matrix, std::size_t n)
{
  using std::abs;
  using std::log;
  if ((n && n > matrix.size() / n) || matrix.size() != n * n)
    throw std::invalid_argument("invalid determinant matrix size");
  LogDeterminant<Real> result;
  CompensatedSum<Real> logarithm;
  for (std::size_t j = 0; j < n; ++j)
  {
    Real scale{0};
    for (std::size_t i = 0; i < n; ++i)
    {
      Real const magnitude = abs(matrix[i * n + j]);
      if (!uni20::isfinite(magnitude)) return {{}, Real{0}};
      scale = std::max(scale, magnitude);
    }
    if (scale == Real{0}) return {{}, Real{0}};
    logarithm.add(log(scale));
    for (std::size_t i = 0; i < n; ++i)
      matrix[i * n + j] /= scale;
  }
  Real growth{1};
  for (std::size_t k = 0; k < n; ++k)
  {
    std::size_t row = k, col = k;
    Real largest{0};
    for (std::size_t i = k; i < n; ++i)
      for (std::size_t j = k; j < n; ++j)
      {
        Real const magnitude = abs(matrix[i * n + j]);
        if (!uni20::isfinite(magnitude)) return {{}, Real{0}};
        if (magnitude > largest)
        {
          largest = magnitude;
          row = i;
          col = j;
        }
      }
    growth = std::max(growth, largest);
    result.pivot_margin = std::min(result.pivot_margin, largest / growth);
    if (result.pivot_margin <= Real{128} * Real(n) * uni20::numeric_limits<Real>::epsilon()) return result;
    if (row != k)
      for (std::size_t j = 0; j < n; ++j)
        std::swap(matrix[k * n + j], matrix[row * n + j]);
    if (col != k)
      for (std::size_t i = 0; i < n; ++i)
        std::swap(matrix[i * n + k], matrix[i * n + col]);
    logarithm.add(log(largest));
    for (std::size_t i = k + 1; i < n; ++i)
    {
      Scalar const factor = matrix[i * n + k] / matrix[k * n + k];
      for (std::size_t j = k + 1; j < n; ++j)
        matrix[i * n + j] -= factor * matrix[k * n + j];
    }
  }
  if (uni20::isfinite(logarithm.value())) result.log_absolute = logarithm.value();
  return result;
}
} // namespace bethe::detail
