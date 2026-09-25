// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "exclusion_support.hpp"
#include "test_support.hpp"
#include <bethe/asep.hpp>

namespace
{
template <typename Real> class ASEP : public ::testing::Test {};
TYPED_TEST_SUITE(ASEP, test_support::RealTypes, test_support::PrecisionNames);
TEST(ASEPExact, MarkovGapEverySmallFilling)
{
  for (unsigned l = 2; l <= 9; ++l)
    for (unsigned n = 1; n < l; ++n)
      for (double q : {0., .1, .5, .9, .99, 1.})
      {
        SCOPED_TRACE(::testing::Message() << l << " " << n << " q=" << q);
        auto const values = test_support::exclusion_spectrum(l, n, 1, q);
        auto const state = bethe::asep::relaxation_gap(l, n, 1., q);
        ASSERT_TRUE(state.converged) << int(state.status) << " reached=" << state.reached_ratio;
        EXPECT_NEAR(*state.gap, -values[1].real(), 2e-11);
        EXPECT_NEAR(state.eigenvalue->imag(), std::abs(values[1].imag()), 2e-11);
      }
}
TYPED_TEST(ASEP, OriginalEquationsScalingAndSymmetries)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t l : {5u, 8u, 16u, 32u})
    for (std::size_t n : {std::size_t{2}, l / 2})
      for (Real q : {Real{0.3}, Real{0.9}, Real{0.99}})
      {
        SCOPED_TRACE(::testing::Message() << l << " " << n << " " << uni20::format_scalar(q));
        auto const state = bethe::asep::relaxation_gap(l, n, Real{1}, q);
        ASSERT_TRUE(state.converged) << int(state.status) << " reached=" << uni20::format_scalar(state.reached_ratio);
        auto const scaled = bethe::asep::relaxation_gap(l, l - n, Real{3} * q, Real{3});
        ASSERT_TRUE(scaled.converged);
        EXPECT_REAL_NEAR(*scaled.gap, Real{3} * *state.gap, Real{32768} * eps);
        std::vector<C> z;
        for (std::size_t j = 0; j < state.scaled_roots.size(); ++j)
          z.push_back(Real{1} + (Real{1} - q) * state.scaled_roots[j] +
                      (j == state.wave_index ? state.wave_base : C{}));
        C eigenvalue{};
        for (C root : z)
          eigenvalue += q * root + Real{1} / root - Real{1} - q;
        EXPECT_REAL_NEAR(std::abs(eigenvalue - *state.eigenvalue), Real{0}, Real{32768} * eps);
        for (std::size_t i = 0; i < z.size(); ++i)
        {
          C rhs = z.size() % 2 ? C(1) : C(-1);
          for (std::size_t j = 0; j < z.size(); ++j)
            if (i != j)
              rhs *= (q * z[i] * z[j] - (Real{1} + q) * z[i] + Real{1}) /
                     (q * z[i] * z[j] - (Real{1} + q) * z[j] + Real{1});
          EXPECT_REAL_NEAR(std::abs(std::pow(z[i], int(l)) / rhs - C(1)), Real{0}, Real{1000000} * eps);
        }
      }
}
TYPED_TEST(ASEP, SymmetricEndpointAndTinyBias)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  for (std::size_t n : {1u, 2u, 3u, 4u})
  {
    auto const symmetric = bethe::asep::relaxation_gap(5, n, Real{2}, Real{2});
    ASSERT_TRUE(symmetric.converged);
    EXPECT_TRUE(symmetric.analytic);
    EXPECT_REAL_NEAR(*symmetric.gap, Real{8} * std::sin(pi / Real{5}) * std::sin(pi / Real{5}), Real{128} * eps);
    EXPECT_EQ(symmetric.eigenvalue->imag(), Real{0});
  }
  Real const q = Real{1} - Real{8} * eps;
  auto const near = bethe::asep::relaxation_gap(5, 2, Real{1}, q);
  ASSERT_TRUE(near.converged) << int(near.status) << " reached=" << uni20::format_scalar(near.reached_ratio);
  EXPECT_FALSE(near.analytic);
  // First-order drift of the one-magnon descendant: sin(2*pi/L)*(L-2N)/(L-2).
  EXPECT_REAL_NEAR(near.eigenvalue->imag() / (Real{1} - q), std::sin(Real{2} * pi / Real{5}) / Real{3},
                   Real{32768} * eps);
  C const tiny(eps * eps, eps * eps);
  EXPECT_REAL_NEAR(std::abs(bethe::detail::complex_log1p(tiny) / tiny - C(1)), Real{0}, Real{32} * eps);
}
TYPED_TEST(ASEP, BudgetsAndValidation)
{
  using Real = TypeParam;
  using namespace bethe::asep;
  auto const empty = relaxation_gap(8, 0, Real{1}, Real{1});
  EXPECT_TRUE(empty.converged);
  EXPECT_FALSE(empty.gap);
  EXPECT_EQ(empty.status, Status::stationary_only);
  auto const no_steps = relaxation_gap(8, 3, Real{1}, Real{0.5}, {.max_continuation_steps = 0});
  EXPECT_EQ(no_steps.status, Status::continuation_limit);
  EXPECT_FALSE(no_steps.eigenvalue);
  auto const no_newton = relaxation_gap(8, 3, Real{1}, Real{0.5}, {.max_iterations = 0});
  EXPECT_EQ(no_newton.status, Status::iteration_limit);
  EXPECT_FALSE(no_newton.gap);
  Options<Real> no_seed;
  no_seed.seed_options.max_iterations = 0;
  auto const failed_seed = relaxation_gap(8, 3, Real{1}, Real{0.5}, no_seed);
  EXPECT_EQ(failed_seed.status, Status::seed_limit);
  EXPECT_FALSE(failed_seed.eigenvalue);
  auto const overflow = relaxation_gap(2, 1, uni20::numeric_limits<Real>::max(), Real{0});
  EXPECT_EQ(overflow.status, Status::precision_limit);
  EXPECT_FALSE(overflow.gap);
  EXPECT_THROW(relaxation_gap(8, 9, Real{1}, Real{1}), std::invalid_argument);
  EXPECT_THROW(relaxation_gap(8, 3, Real{0}, Real{0}), std::invalid_argument);
  EXPECT_THROW(relaxation_gap(8, 3, Real{-1}, Real{1}), std::invalid_argument);
  EXPECT_THROW(relaxation_gap(8, 3, Real{1}, Real{1}, {.tolerance = Real{0}}), std::invalid_argument);
}
TYPED_TEST(ASEP, IndependentNativeReferenceAndJacobian)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  // Independent 90-digit diagonalization of the ten-configuration generator.
  C const expected(uni20::parse_real<Real>(
                       "-1.04039558904641505636142345412607413799326810221158347669691888089548565042267332323361234"),
                   uni20::parse_real<Real>(
                       "0.159168821610180253202429215239399086838149271486788653791649970802963326547367324815530809"));
  auto const reference = bethe::asep::relaxation_gap(5, 2, Real{1}, Real{0.5});
  ASSERT_TRUE(reference.converged);
  EXPECT_REAL_NEAR(std::abs(*reference.eigenvalue - expected), Real{0}, Real{8192} * eps);
  for (Real q : {Real{0.3}, Real{0.9}, Real{1} - Real{8} * eps})
  {
    auto const state = bethe::asep::relaxation_gap(8, 3, Real{1}, q);
    ASSERT_TRUE(state.converged);
    bethe::asep::detail::System<Real> system{8, 3, state.wave_index, q, Real{1} - q, state.wave_base};
    auto const e = system.evaluate(state.scaled_roots, true);
    Real const h = std::cbrt(eps);
    for (std::size_t col = 0; col < 6; ++col)
    {
      auto plus = state.scaled_roots, minus = plus;
      C const step = col % 2 ? C(0, h) : C(h, 0);
      plus[col / 2] += step;
      minus[col / 2] -= step;
      auto const fp = system.evaluate(plus), fm = system.evaluate(minus);
      for (std::size_t row = 0; row < 6; ++row)
        EXPECT_REAL_NEAR(e.jacobian[row * 6 + col], (fp.residual[row] - fm.residual[row]) / (Real{2} * h),
                         Real{10000} * h * h);
    }
  }
}
} // namespace
