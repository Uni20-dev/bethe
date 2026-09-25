// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/sine_gordon.hpp>

namespace
{
template <typename Real> class SineGordon : public ::testing::Test {};
TYPED_TEST_SUITE(SineGordon, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(SineGordon, KernelAgainstIndependentFourierTransforms)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const pi = Real{4} * std::atan(Real{1});
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real p : {Real{0.5}, Real{2}})
    for (C z : {C{}, C(Real{0.7}, Real{0.3}), C(Real{4}, Real{-0.6}), C(0, Real{1.2})})
    {
      auto const result = bethe::sine_gordon::scattering_kernel(z, p);
      ASSERT_TRUE(result.converged) << int(result.status);
      // Closed inverse transforms of sech(pi*k/2) and sech²(pi*k/2).
      // The implementation deliberately does not special-case either p.
      C expected;
      if (p == Real{0.5})
        expected = -Real{1} / (Real{2} * pi * std::cosh(z));
      else
        expected = z == C{} ? C(Real{1} / (Real{2} * pi * pi)) : z / (Real{2} * pi * pi * std::sinh(z));
      EXPECT_REAL_NEAR(std::abs(*result.value - expected), Real{0}, Real{4096} * eps);
      EXPECT_LE(result.quadrature_error + result.tail_bound, Real{1024} * eps);
    }
}
TYPED_TEST(SineGordon, AnalyticSymmetriesAndNearFreeCoupling)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real p : {Real{0.3}, Real{0.8}, Real{1.2}, Real{3.7}, Real{1000}})
  {
    C const z(Real{1.3}, Real{0.2});
    auto const a = bethe::sine_gordon::scattering_kernel(z, p);
    auto const b = bethe::sine_gordon::scattering_kernel(-z, p);
    auto const c = bethe::sine_gordon::scattering_kernel(std::conj(z), p);
    ASSERT_TRUE(a.converged);
    ASSERT_TRUE(b.converged);
    ASSERT_TRUE(c.converged);
    EXPECT_REAL_NEAR(std::abs(*a.value - *b.value), Real{0}, Real{32} * eps);
    EXPECT_REAL_NEAR(std::abs(*a.value - std::conj(*c.value)), Real{0}, Real{32} * eps);
  }
  // At z=0 the derivative dG/dp at p=1 is 1/8:
  // integral_0^infinity k/(2*sinh(pi*k)) dk = 1/8.
  for (Real sign : {Real{-1}, Real{1}})
  {
    Real const p = Real{1} + sign * Real{8} * eps;
    bethe::sine_gordon::KernelOptions<Real> options;
    options.tolerance *= std::abs(p - Real{1});
    auto const result = bethe::sine_gordon::scattering_kernel(C{}, p, options);
    ASSERT_TRUE(result.converged);
    EXPECT_REAL_NEAR(result.value->real() / (p - Real{1}), Real{1} / Real{8}, Real{4096} * eps);
  }
}
TYPED_TEST(SineGordon, BudgetsAndValidation)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  using namespace bethe::sine_gordon;
  auto const free = scattering_kernel(C{}, Real{1}, {.max_evaluations = 0, .max_cutoffs = 0});
  ASSERT_TRUE(free.converged);
  EXPECT_EQ(free.value->real(), Real{0});
  EXPECT_EQ(free.value->imag(), Real{0});
  EXPECT_EQ(free.evaluations, 0u);
  auto const cutoff = scattering_kernel(C{}, Real{2}, {.max_cutoffs = 0});
  EXPECT_EQ(cutoff.status, KernelStatus::cutoff_limit);
  EXPECT_FALSE(cutoff.value);
  auto const mesh = scattering_kernel(C{}, Real{2}, {.max_evaluations = 1});
  EXPECT_EQ(mesh.status, KernelStatus::quadrature_limit);
  EXPECT_FALSE(mesh.value);
  auto const levels = scattering_kernel(C{}, Real{2}, {.max_levels = 0});
  EXPECT_EQ(levels.status, KernelStatus::quadrature_limit);
  EXPECT_FALSE(levels.value);
  EXPECT_THROW(scattering_kernel(C{}, Real{0}), std::invalid_argument);
  EXPECT_THROW(scattering_kernel(C(0, 2), Real{0.5}), std::invalid_argument);
  EXPECT_THROW(scattering_kernel(C{}, Real{2}, {.tolerance = Real{0}}), std::invalid_argument);
  EXPECT_THROW(scattering_kernel(C{}, uni20::numeric_limits<Real>::infinity()), std::invalid_argument);
}
TYPED_TEST(SineGordon, NativeReferenceAndOverflowSafeContour)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  C const expected(
      uni20::parse_real<Real>(
          "0.0547765717918015314046546239425095204723847959115512148661684638805025105506472402803056409"),
      uni20::parse_real<Real>(
          "-0.00515484562906219103918414319730092148237776521629576139173090753048764604314854689936674796"));
  auto const result = bethe::sine_gordon::scattering_kernel(
      C(uni20::parse_real<Real>("0.8"), uni20::parse_real<Real>("0.4")), uni20::parse_real<Real>("2.7"));
  ASSERT_TRUE(result.converged);
  EXPECT_REAL_NEAR(std::abs(*result.value - expected), Real{0}, Real{4096} * eps);
  // Separate cosh(k*Im(z)) would overflow in fp64, although their product
  // with the exponentially small Fourier multiplier is well conditioned.
  Real const y = pi - Real{1} / Real{1000};
  auto const integrand = bethe::sine_gordon::detail::fourier_integrand(Real{1000}, C(0, y), Real{2});
  EXPECT_REAL_NEAR(integrand.real(), std::exp(-(pi - y) * Real{1000}) / (Real{2} * pi), Real{128} * eps);
  EXPECT_EQ(integrand.imag(), Real{0});
}
} // namespace
