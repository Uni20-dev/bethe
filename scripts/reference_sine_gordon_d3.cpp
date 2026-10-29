// SPDX-License-Identifier: GPL-3.0-or-later
// Independent real D3 TBA at p=2: Hegedus, 2510.25344, eqs. (2.6)-(2.8).
// Positive-rapidity Gauss mesh and analytic magnon plateau; no contour NLIE.
// Build target bethe_sine_gordon_d3_reference, then pass POINTS BITS LENGTH.
#include <algorithm>
#include <bethe/detail/gauss_legendre.hpp>
#include <bethe/solver.hpp>
#include <iostream>
#include <uni20/common/presentation.hpp>
template <typename Real> Real d3_reference(Real u, std::size_t n)
{
  Real const zero(0, u.precision()), one(1, u.precision());
  Real const pi = 4 * uni20::math::atan(one), log2 = uni20::math::log(2 * one);
  Real const eps = uni20::numeric_limits<Real>::epsilon(u.precision()), b = -uni20::math::log(eps) / Real{2} + Real{4};
  auto const rule = bethe::detail::gauss_legendre<Real>(n, u.precision());
  std::vector<Real> x(n, zero), w(n, zero), bare(n, zero), v(n, zero), a(n, zero), l(n, zero), next(n, zero),
      k(n * n, zero);
  for (std::size_t i = 0; i < n; ++i)
  {
    x[i] = b * (1 + rule.x[i]) / 2;
    w[i] = b * rule.w[i] / 2;
    bare[i] = u * uni20::math::cosh(x[i]) - log2;
  }
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      k[i * n + j] = w[j] / (2 * pi) * (1 / uni20::math::cosh(x[i] - x[j]) + 1 / uni20::math::cosh(x[i] + x[j]));
  for (std::size_t iter = 0; iter < 256; ++iter)
  {
    for (std::size_t j = 0; j < n; ++j)
      l[j] = uni20::math::log1p(uni20::math::exp(-bare[j] - v[j]));
    for (std::size_t i = 0; i < n; ++i)
    {
      bethe::detail::CompensatedSum<Real> sum(zero);
      for (std::size_t j = 0; j < n; ++j)
        sum.add(-k[i * n + j] * l[j]);
      a[i] = uni20::math::log1p(uni20::math::exp(-sum.value())) - log2;
    }
    Real error = zero;
    for (std::size_t i = 0; i < n; ++i)
    {
      bethe::detail::CompensatedSum<Real> sum(zero);
      for (std::size_t j = 0; j < n; ++j)
        sum.add(-2 * k[i * n + j] * a[j]);
      next[i] = sum.value();
      error = std::max(error, uni20::math::abs(next[i] - v[i]));
    }
    if (error < 32 * eps)
    {
      bethe::detail::CompensatedSum<Real> sum(zero);
      for (std::size_t j = 0; j < n; ++j)
        sum.add(-w[j] * u * uni20::math::cosh(x[j]) * l[j] / pi);
      return sum.value();
    }
    v = next;
  }
  throw std::runtime_error("independent D3 reference failed to converge");
}

int main(int argc, char** argv)
{
  using R = uni20::mpreal;
  if (argc != 4)
  {
    std::cerr << "usage: reference POINTS BITS LENGTH\n";
    return 1;
  }
  auto const p = uni20::Precision::bits(std::stoul(argv[2]));
  std::size_t n = std::stoul(argv[1]);
  R u(argv[3], p);
  if (n < 2 || n > 4096 || !uni20::isfinite(u) || u <= R{0})
    throw std::invalid_argument("reference requires 2..4096 points and finite positive scaled length");
  auto y = d3_reference(u, n);
  std::cout << uni20::format_scalar(y) << std::endl;
}
