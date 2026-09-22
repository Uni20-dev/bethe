// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz_common.hpp>
#include <numeric>
#include <uni20/linalg/ops/linear_solve.hpp>

namespace bethe::xxz::detail
{
enum class WronskianStatus
{
  consistent,
  inconsistent,
  ill_conditioned,
  nonfinite
};

/// Numerical consistency of an odd-ring quantum Wronskian, NOT a physical
/// state or ground-state certificate. The generic-q theorem does not cover
/// infinite rapidities or root-of-unity strings. See docs/xxz-wronskian.md.
template <uni20::Real Real> struct WronskianCheck
{
    WronskianStatus status = WronskianStatus::nonfinite;
    Real residual_norm = uni20::numeric_limits<Real>::infinity();
    Real reciprocal_condition = Real{0};
    Real endpoint_ratio = Real{0};
    // Q(t)=sum q[j]*t^(2*j-M); P(t)=sum p[j]*t^(2*j-(N-M))/(2*i).
    // The RHS is (t-t^-1)^N / binomial(N,floor(N/2)).
    std::vector<Real> q, p;
};

/// Independently test (t-t^-1)^N proportional to P^+ Q^- - P^- Q^+.
/// Input is monic Q(x), z=center+scale*x, as in PolynomialBetheSystem.
/// Uses a rank-revealing row selection followed by Uni20's native-precision
/// square solve. ALL N+1 coefficient equations are checked, not just those
/// used to determine P. No normal equations or root extraction are used.
template <uni20::Real Real>
WronskianCheck<Real> check_odd_wronskian(std::size_t sites, std::span<Real const> coefficients, Real delta,
                                         Real center = Real{0}, Real scale = Real{1},
                                         Real tolerance = Real{256} * uni20::numeric_limits<Real>::epsilon())
{
  checked_sites(sites);
  auto const m = coefficients.size();
  if (sites % 2 == 0 || m > sites / 2)
    throw std::invalid_argument("XXZ polynomial Wronskian requires odd N and M <= N/2");
  if (!uni20::isfinite(delta) || delta <= -Real{1} || delta > Real{0})
    throw std::invalid_argument("XXZ polynomial Wronskian requires -1 < Delta <= 0");
  if (!uni20::isfinite(center) || !uni20::isfinite(scale) || scale <= Real{0} || !uni20::isfinite(tolerance) ||
      tolerance <= Real{0})
    throw std::invalid_argument("invalid XXZ Wronskian coordinates or tolerance");
  for (Real c : coefficients)
    if (!uni20::isfinite(c)) throw std::invalid_argument("nonfinite XXZ polynomial coefficient");
  auto const rows = sites + 1, cols = sites - m + 1;
  auto const elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (rows > elements / (cols + 1)) throw std::length_error("XXZ Wronskian workspace is too large");
  WronskianCheck<Real> out;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const s = std::sqrt((Real{1} + delta) / (Real{1} - delta));
  // X=t^2=(s+z)/(s-z). Homogenized Horner transforms Q without dividing
  // by Q(s) or Q(-s), which can vanish at phantom-root configurations.
  Real const unit = std::max({s, std::abs(center), scale});
  Real const a0 = -s / unit - center / unit, a1 = s / unit - center / unit, b = scale / unit;
  auto linear = [](std::vector<Real> const& v, Real a, Real b) {
    std::vector<Real> result(v.size() + 1, Real{0});
    for (std::size_t j = 0; j < v.size(); ++j)
    {
      result[j] += a * v[j];
      result[j + 1] += b * v[j];
    }
    return result;
  };
  out.q = {Real{1}};
  std::vector<Real> power{Real{1}};
  for (std::size_t j = m; j-- > 0;)
  {
    out.q = linear(out.q, a0, a1);
    power = linear(power, b, b);
    for (std::size_t k = 0; k < out.q.size(); ++k)
      out.q[k] += coefficients[j] * power[k];
  }
  Real qnorm = Real{0};
  for (Real value : out.q)
  {
    if (!uni20::isfinite(value)) return out;
    qnorm = std::max(qnorm, std::abs(value));
  }
  if (!(qnorm > Real{0})) return out;
  for (Real& value : out.q)
    value /= qnorm;
  out.endpoint_ratio = std::min(std::abs(out.q.front()), std::abs(out.q.back()));

  // Normalize the binomial RHS without ever forming a potentially huge
  // central binomial coefficient. Exponents are 2*k-N, k=0,...,N.
  std::vector<Real> right(rows);
  auto const middle = sites / 2;
  right[middle] = (sites - middle) % 2 ? -Real{1} : Real{1};
  for (std::size_t j = middle; j > 0; --j)
    right[j - 1] = -right[j] * Real(j) / Real(sites - j + 1);
  for (std::size_t j = middle; j < sites; ++j)
    right[j + 1] = -right[j] * Real(sites - j) / Real(j + 1);
  Real const gamma = std::atan2(std::sqrt((Real{1} - delta) * (Real{1} + delta)), delta);
  uni20::DenseMatrix<Real> matrix(rows, cols), work(rows, cols);
  for (std::size_t k = 0; k < rows; ++k)
    for (std::size_t j = 0; j < cols; ++j)
    {
      Real value = Real{0};
      if (k >= j && k - j <= m)
      {
        auto const l = k - j;
        auto const exponent = 2 * std::int64_t(j) - 2 * std::int64_t(l) - std::int64_t(sites) + 2 * std::int64_t(m);
        value = out.q[l] * std::sin(Real(exponent) * gamma / Real{2});
      }
      matrix[k, j] = work[k, j] = value;
    }

  // Complete pivoting selects independent equations. It is used ONLY for
  // row selection; solve the original, uneliminated square system below.
  std::vector<std::size_t> selected(rows);
  std::iota(selected.begin(), selected.end(), std::size_t{0});
  for (std::size_t k = 0; k < cols; ++k)
  {
    std::size_t row = k, col = k;
    Real pivot = Real{0};
    for (std::size_t i = k; i < rows; ++i)
      for (std::size_t j = k; j < cols; ++j)
        if (std::abs(work[i, j]) > pivot)
        {
          pivot = std::abs(work[i, j]);
          row = i;
          col = j;
        }
    if (!(pivot > Real{64} * eps))
    {
      out.status = WronskianStatus::ill_conditioned;
      return out;
    }
    std::swap(selected[k], selected[row]);
    for (std::size_t j = 0; j < cols; ++j)
      std::swap(work[k, j], work[row, j]);
    for (std::size_t i = 0; i < rows; ++i)
      std::swap(work[i, k], work[i, col]);
    for (std::size_t i = k + 1; i < rows; ++i)
    {
      Real const factor = work[i, k] / work[k, k];
      for (std::size_t j = k + 1; j < cols; ++j)
        work[i, j] -= factor * work[k, j];
    }
  }
  uni20::DenseMatrix<Real> square(cols, cols), rhs(cols, cols + 1);
  Real norm = Real{0};
  for (std::size_t i = 0; i < cols; ++i)
  {
    Real row = Real{0};
    rhs[i, 0] = right[selected[i]];
    for (std::size_t j = 0; j < cols; ++j)
    {
      square[i, j] = matrix[selected[i], j];
      row += std::abs(square[i, j]);
      rhs[i, j + 1] = i == j ? Real{1} : Real{0};
    }
    norm = std::max(norm, row);
  }
  try
  {
    uni20::linalg::solve_inplace(square, rhs);
  }
  catch (std::runtime_error const&)
  {
    out.status = WronskianStatus::ill_conditioned;
    return out;
  }
  Real inverse_norm = Real{0};
  out.p.resize(cols);
  for (std::size_t i = 0; i < cols; ++i)
  {
    out.p[i] = rhs[i, 0];
    Real row = Real{0};
    for (std::size_t j = 0; j < cols; ++j)
    {
      if (!uni20::isfinite(rhs[i, j + 1])) return out;
      row += std::abs(rhs[i, j + 1]);
    }
    if (!uni20::isfinite(row) || !uni20::isfinite(out.p[i])) return out;
    inverse_norm = std::max(inverse_norm, row);
  }
  out.reciprocal_condition = (Real{1} / norm) / inverse_norm;
  out.residual_norm = Real{0};
  for (std::size_t i = 0; i < rows; ++i)
  {
    bethe::detail::CompensatedSum<Real> sum;
    for (std::size_t j = 0; j < cols; ++j)
      sum.add(matrix[i, j] * out.p[j]);
    Real const residual = std::abs(sum.value() - right[i]);
    if (!uni20::isfinite(residual))
    {
      out.residual_norm = uni20::numeric_limits<Real>::infinity();
      return out;
    }
    out.residual_norm = std::max(out.residual_norm, residual);
  }
  if (!(out.reciprocal_condition > Real{64} * eps))
    out.status = WronskianStatus::ill_conditioned;
  else
    out.status = out.residual_norm <= tolerance ? WronskianStatus::consistent : WronskianStatus::inconsistent;
  return out;
}
} // namespace bethe::xxz::detail
