// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "exclusion_support.hpp"
#include "test_support.hpp"
#include <bethe/tasep.hpp>

namespace
{
template <typename Real> class TASEP : public ::testing::Test {};
TYPED_TEST_SUITE(TASEP, test_support::RealTypes, test_support::PrecisionNames);
TEST(TASEPExact, MarkovGapAtEverySmallFilling)
{
  for (unsigned l = 2; l <= 10; ++l)
    for (unsigned n = 1; n < l; ++n)
    {
      SCOPED_TRACE(::testing::Message() << l << " " << n);
      auto const values = test_support::exclusion_spectrum(l, n);
      auto const state = bethe::tasep::relaxation_gap(l, n, 1.0);
      ASSERT_TRUE(state.converged) << int(state.status);
      EXPECT_NEAR(*state.gap, -values[1].real(), 2e-11);
      EXPECT_NEAR(std::abs(state.eigenvalue->imag()), std::abs(values[1].imag()), 2e-11);
    }
}
TYPED_TEST(TASEP, NativeEquationsScalingAndLargerRings)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t l : {4u, 5u, 8u, 16u, 32u, 64u, 128u, 256u})
    for (std::size_t n :
         {std::size_t{1}, std::size_t{2}, std::size_t{3}, std::max(std::size_t{2}, l / 4), l / 2, l - 1})
    {
      SCOPED_TRACE(::testing::Message() << l << " " << n);
      auto const s = bethe::tasep::relaxation_gap(l, n, Real{1});
      ASSERT_TRUE(s.converged) << int(s.status) << " seed=" << s.seed_iterations << " updates=" << s.iterations
                               << " roots=" << s.roots.size() << " residual=" << uni20::format_scalar(s.residual_norm);
      auto const scaled = bethe::tasep::relaxation_gap(l, n, Real{3});
      ASSERT_TRUE(scaled.converged);
      EXPECT_REAL_NEAR(*scaled.gap, Real{3} * *s.gap, Real{32} * eps);
      C product(-std::pow(Real{2}, Real(l)), 0);
      for (C z : s.roots)
        product *= (z - Real{1}) / (z + Real{1});
      for (C z : s.roots)
      {
        C const lhs =
            std::pow(Real{1} - z, int(s.effective_particles)) * std::pow(Real{1} + z, int(l - s.effective_particles));
        EXPECT_REAL_NEAR(std::abs(lhs / product - C(1)), Real{0}, Real{16384} * Real(l) * eps);
      }
      EXPECT_GT(*s.gap, Real{0});
    }
}
TYPED_TEST(TASEP, BudgetsAndValidation)
{
  using Real = TypeParam;
  using namespace bethe::tasep;
  for (std::size_t n : {0u, 8u})
  {
    auto const state = relaxation_gap(8, n, Real{1});
    EXPECT_TRUE(state.converged);
    EXPECT_FALSE(state.gap);
    EXPECT_FALSE(state.eigenvalue);
    EXPECT_EQ(state.status, Status::stationary_only);
  }
  auto const no_seed = relaxation_gap(8, 3, Real{1}, {.max_seed_iterations = 0});
  EXPECT_EQ(no_seed.status, Status::seed_limit);
  EXPECT_FALSE(no_seed.eigenvalue);
  auto const no_newton = relaxation_gap(8, 3, Real{1}, {.max_iterations = 0});
  EXPECT_EQ(no_newton.status, Status::iteration_limit);
  EXPECT_FALSE(no_newton.eigenvalue);
  EXPECT_THROW(relaxation_gap(1, 0, Real{1}), std::invalid_argument);
  EXPECT_THROW(relaxation_gap(4, 5, Real{1}), std::invalid_argument);
  EXPECT_THROW(relaxation_gap(4, 2, Real{0}), std::invalid_argument);
  EXPECT_THROW(relaxation_gap(4, 2, Real{1}, {.tolerance = Real{0}}), std::invalid_argument);
  EXPECT_THROW(relaxation_gap(8, 3, Real{1}, {.max_sites = 4}), std::length_error);
  auto const overflow = relaxation_gap(2, 1, uni20::numeric_limits<Real>::max());
  EXPECT_FALSE(overflow.converged);
  EXPECT_FALSE(overflow.gap);
  EXPECT_EQ(overflow.status, Status::precision_limit);
  EXPECT_THROW(relaxation_gap(4, 2, uni20::numeric_limits<Real>::quiet_NaN()), std::invalid_argument);
}

TYPED_TEST(TASEP, IndependentNativeReferenceAndJacobian)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  // Independent 90-digit eigensystem of the ten-configuration Markov matrix.
  C const expected(uni20::parse_real<Real>(
                       "-0.710290069123028672493762532365387684613821535715147285599991234743117333842877636225644555"),
                   uni20::parse_real<Real>(
                       "0.326895513054695736451443863914497675039691838878110784103129825870014387978645791366909038"));
  auto const state = bethe::tasep::relaxation_gap(5, 2, Real{1});
  ASSERT_TRUE(state.converged);
  EXPECT_REAL_NEAR(std::abs(*state.eigenvalue - expected), Real{0}, Real{8192} * eps);
  bethe::tasep::detail::System<Real> system{5, 2};
  auto const evaluation = system.evaluate(state.roots, state.log_y, true);
  Real const h = std::cbrt(eps);
  for (std::size_t col = 0; col < 6; ++col)
  {
    auto plus = state.roots, minus = state.roots;
    auto bp = state.log_y, bm = state.log_y;
    C const step = col % 2 ? C(0, h) : C(h, 0);
    if (col < 4)
    {
      plus[col / 2] += step;
      minus[col / 2] -= step;
    }
    else
    {
      bp += step;
      bm -= step;
    }
    auto const fp = system.evaluate(plus, bp, false), fm = system.evaluate(minus, bm, false);
    for (std::size_t row = 0; row < 6; ++row)
      EXPECT_REAL_NEAR(evaluation.jacobian[row * 6 + col], (fp.residual[row] - fm.residual[row]) / (Real{2} * h),
                       Real{10000} * h * h);
  }
}
} // namespace
