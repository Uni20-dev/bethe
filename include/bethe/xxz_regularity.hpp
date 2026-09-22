// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz_polynomial.hpp>

namespace bethe::xxz::detail
{
enum class RegularityStatus
{
  regular_on_shell,
  off_shell,
  exceptional_or_unresolved,
  ill_conditioned,
  nonfinite
};

/// Numerical checks of the REGULAR Gaudin-Korepin nonzero-vector criterion.
/// This is not an interval proof, a norm value, or a ground-state test.
/// Exceptional/infinite/repeated roots require a different construction.
template <uni20::Real Real> struct PolynomialRegularity
{
    RegularityStatus status = RegularityStatus::nonfinite;
    Real residual_norm = uni20::numeric_limits<Real>::infinity();
    Real endpoint_margin = Real{0}, driving_margin = Real{0};
    Real distinct_roots_margin = Real{0}, scattering_margin = Real{0}, jacobian_margin = Real{0};
};

// Complete-pivot elimination, deliberately NOT a linear solve: singularity
// is an expected diagnostic result, whereas Uni20's solve has a terminal
// singular-matrix policy. Return min|pivot|/(N*max(1,element growth)) after
// normalizing the original largest entry to one. This is a rank-resolution
// indicator, NOT a reciprocal condition number or a rigorous error bound.
template <uni20::Real Real> Real regularity_pivot_margin(uni20::DenseMatrix<Real> matrix)
{
  auto const n = std::size_t(matrix.extent(0));
  if (!n) return Real{1};
  Real scale = Real{0};
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
    {
      if (!uni20::isfinite(matrix[i, j])) return Real{0};
      scale = std::max(scale, std::abs(matrix[i, j]));
    }
  if (!(scale > Real{0})) return Real{0};
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      matrix[i, j] /= scale;
  Real minimum = Real{1}, growth = Real{1};
  for (std::size_t k = 0; k < n; ++k)
  {
    std::size_t row = k, col = k;
    Real pivot = Real{0};
    for (std::size_t i = k; i < n; ++i)
      for (std::size_t j = k; j < n; ++j)
        if (std::abs(matrix[i, j]) > pivot)
        {
          pivot = std::abs(matrix[i, j]);
          row = i;
          col = j;
        }
    if (!(pivot > Real{0})) return Real{0};
    minimum = std::min(minimum, pivot);
    for (std::size_t j = k; j < n; ++j)
      std::swap(matrix[k, j], matrix[row, j]);
    for (std::size_t i = k; i < n; ++i)
      std::swap(matrix[i, k], matrix[i, col]);
    for (std::size_t i = k + 1; i < n; ++i)
    {
      Real const factor = matrix[i, k] / matrix[k, k];
      for (std::size_t j = k + 1; j < n; ++j)
      {
        matrix[i, j] -= factor * matrix[k, j];
        if (!uni20::isfinite(matrix[i, j])) return Real{0};
        growth = std::max(growth, std::abs(matrix[i, j]));
      }
    }
  }
  return (minimum / growth) / Real(n);
}

// Multiplication by g in C[x]/Q. It is invertible exactly when g and Q
// have no common root. A real 2M-by-2M representation preserves fp80/fp128
// without requiring a complex LAPACK backend.
template <uni20::Real Real>
Real quotient_multiplication_margin(std::span<Real const> c, std::vector<std::complex<Real>> g)
{
  auto const m = c.size();
  uni20::DenseMatrix<Real> matrix(2 * m, 2 * m);
  for (std::size_t col = 0; col < m; ++col)
  {
    for (std::size_t row = 0; row < m; ++row)
    {
      matrix[row, col] = matrix[row + m, col + m] = g[row].real();
      matrix[row + m, col] = g[row].imag();
      matrix[row, col + m] = -g[row].imag();
    }
    auto const highest = g.back();
    for (std::size_t j = m; j-- > 0;)
      g[j] = (j ? g[j - 1] : std::complex<Real>{}) - highest * c[j];
  }
  return regularity_pivot_margin(std::move(matrix));
}

// Homogeneous value b^M Q(a/b), relative to a coefficient magnitude bound
// on |x|<=max(1,|a/b|). The unit-radius floor matters for Q(x)=x^M near
// x=0: a small value need not involve cancellation. Neither the endpoint
// coordinate nor powers of 1/b need be formed.
template <uni20::Real Real> Real homogeneous_value_margin(std::span<Real const> c, std::complex<Real> a, Real b)
{
  Real const scale = std::max({std::abs(a.real()), std::abs(a.imag()), b});
  a /= scale;
  b /= scale;
  Real const radius = std::max(std::abs(a), b);
  std::complex<Real> value{1};
  Real bound = Real{1}, power = Real{1};
  for (std::size_t j = c.size(); j-- > 0;)
  {
    power *= b;
    Real const term = c[j] * power;
    value = value * a + term;
    bound = bound * radius + std::abs(term);
  }
  if (!(bound > Real{0}) || !uni20::isfinite(bound)) return Real{0};
  Real const result = std::abs(value) / bound;
  return uni20::isfinite(result) ? result : Real{0};
}

/// On an EXACT solution, nonzero endpoint/driving/scattering factors, distinct
/// roots, and a nonsingular full coefficient Jacobian imply a nonzero Bethe
/// vector by the Gaudin-Korepin scalar-product formula. This routine tests
/// numerical versions of those hypotheses, without extracting roots. No
/// generic-q Wronskian assumption is used. See docs/xxz-regularity.md.
template <uni20::Real Real>
PolynomialRegularity<Real>
check_regular_polynomial(PolynomialBetheSystem<Real> const& system, std::span<Real const> c, Real delta,
                         Real residual_tolerance = Real{32} * uni20::numeric_limits<Real>::epsilon())
{
  if (!uni20::isfinite(residual_tolerance) || residual_tolerance <= Real{0})
    throw std::invalid_argument("regularity residual tolerance must be finite and positive");
  auto const m = system.order;
  auto const elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (m && m > elements / 4 / m) throw std::length_error("XXZ regularity workspace is too large");
  // This validates the coefficient count/values, coordinates and Delta before
  // the numerical-failure handler; invalid input must still throw.
  PolynomialRegularity<Real> out;
  std::vector<std::complex<Real>> scattering;
  try
  {
    scattering = system.scattering_remainder(c, delta);
  }
  catch (std::runtime_error const&)
  {
    return out;
  }
  Real const floor = Real{64} * uni20::numeric_limits<Real>::epsilon();
  if (!m)
  {
    out.status = RegularityStatus::regular_on_shell;
    out.residual_norm = Real{0};
    out.endpoint_margin = out.driving_margin = Real{1};
    out.distinct_roots_margin = out.scattering_margin = out.jacobian_margin = Real{1};
    return out;
  }
  Real const s = std::sqrt((Real{1} + delta) / (Real{1} - delta));
  out.endpoint_margin =
      std::min(homogeneous_value_margin<Real>(c, {s - system.center, Real{0}}, system.coordinate_scale),
               homogeneous_value_margin<Real>(c, {-s - system.center, Real{0}}, system.coordinate_scale));
  out.driving_margin = homogeneous_value_margin<Real>(c, {-system.center, Real{1}}, system.coordinate_scale);
  if (!(out.endpoint_margin > floor) || !(out.driving_margin > floor))
  {
    out.status = RegularityStatus::exceptional_or_unresolved;
    return out;
  }
  std::vector<std::complex<Real>> derivative(m);
  for (std::size_t j = 0; j < m; ++j)
    derivative[j] = Real(j + 1) * (j + 1 == m ? Real{1} : c[j + 1]);
  out.distinct_roots_margin = quotient_multiplication_margin<Real>(c, std::move(derivative));
  out.scattering_margin = quotient_multiplication_margin<Real>(c, scattering);
  if (!(out.distinct_roots_margin > floor) || !(out.scattering_margin > floor))
  {
    out.status = RegularityStatus::exceptional_or_unresolved;
    return out;
  }
  try
  {
    uni20::DenseMatrix<Real> jac(m, m);
    auto const f = system.evaluate(c, delta, &jac);
    out.residual_norm = f.norm;
    if (f.norm > residual_tolerance)
    {
      out.status = RegularityStatus::off_shell;
      return out;
    }
    out.jacobian_margin = regularity_pivot_margin(std::move(jac));
    out.status = out.jacobian_margin > floor ? RegularityStatus::regular_on_shell : RegularityStatus::ill_conditioned;
  }
  catch (std::runtime_error const&)
  {} // Leave nonfinite status; unavailable residuals are infinite.
  return out;
}
} // namespace bethe::xxz::detail
