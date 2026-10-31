// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxx_thermodynamic_structure_factor.hpp>
#include <bethe/xxx_two_spinon_moments.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
namespace m = bethe::heisenberg;
std::array<double, 2> frequency_moments(double q, double exchange, double gap)
{
  m::StructureFactorOptions<double> options;
  options.tolerance = 1e-10;
  m::ThermodynamicTwoSpinonStructureFactor<double> spectrum(exchange, options);
  auto const [lower, upper] = spectrum.boundaries(q);
  double const pi = 4 * std::atan(1.);
  double const bottom = std::sqrt(lower * lower + (upper - lower) * (upper + lower) / std::pow(std::cosh(gap / 2), 2));
  auto const rule = bethe::detail::gauss_legendre<double>(96);
  std::array<double, 2> sum{};
  for (std::size_t i = 0; i < rule.x.size(); ++i)
  {
    double const a = pi * (1 + rule.x[i]) / 4;
    double const omega = bottom + (upper - bottom) * std::pow(std::sin(a), 2);
    double const measure = (upper - bottom) * std::sin(2 * a) * rule.w[i] / 8;
    auto point = spectrum(q, omega);
    EXPECT_TRUE(point.value && point.error);
    if (point.value)
    {
      sum[0] += measure * *point.value;
      sum[1] += measure * omega * *point.value;
    }
  }
  return sum;
}
} // namespace

TEST(TwoSpinonMoments, IndependentFrequencyIntegralsAndExchange)
{
  m::StructureFactorOptions<double> options;
  options.tolerance = 1e-10;
  double const pi = 4 * std::atan(1.);
  for (double j : {1., 3.})
    for (double q : {pi / 4, pi / 2, pi})
    {
      auto const result = m::ThermodynamicTwoSpinonMoments<double>(j, options).fixed_moments(q, 32, 6.).evaluate();
      ASSERT_EQ(result.status, m::TwoSpinonMomentStatus::estimated);
      ASSERT_TRUE(result.integral.sums);
      auto const reference = frequency_moments(q, j, 6.);
      for (int i = 0; i < 2; ++i)
        EXPECT_NEAR((*result.integral.sums)[i], reference[i], 3e-8 * reference[i]);
      EXPECT_FALSE(result.sum_rule_denominators[0]);
      EXPECT_TRUE(result.sum_rule_fractions[1]);
      EXPECT_EQ(result.unrestricted_static_divergence, q == pi);
    }
}

TEST(TwoSpinonMoments, GlobalMeasureMatchesIndependentMomentumIntegration)
{
  m::StructureFactorOptions<double> options;
  options.tolerance = 1e-10;
  m::ThermodynamicTwoSpinonMoments<double> model(1., options);
  auto const global = model.global_moments(8, 4.).evaluate();
  ASSERT_TRUE(global.integral.sums);
  double const pi = 4 * std::atan(1.);
  auto const rule = bethe::detail::gauss_legendre<double>(32);
  std::array<double, 2> sum{};
  for (std::size_t i = 0; i < rule.x.size(); ++i)
  {
    auto const fixed = model.fixed_moments(pi * (1 + rule.x[i]) / 2, 8, 4.).evaluate();
    ASSERT_TRUE(fixed.integral.sums);
    for (int n = 0; n < 2; ++n)
      sum[n] += rule.w[i] / 2 * (*fixed.integral.sums)[n];
  }
  for (int n = 0; n < 2; ++n)
    EXPECT_NEAR(sum[n], (*global.integral.sums)[n], 2e-12);
  EXPECT_DOUBLE_EQ(*global.sum_rule_denominators[0], .25);
  EXPECT_NEAR(*global.sum_rule_fractions[0], 4 * sum[0], 1e-11);
  EXPECT_NEAR(*global.sum_rule_fractions[1], sum[1] / (2. / 3 * (std::log(2.) - .25)), 1e-11);
}

TEST(TwoSpinonMoments, OwnedAsyncCutoffAndEndpointControls)
{
  m::StructureFactorOptions<double> options;
  options.tolerance = 1e-10;
  m::ThermodynamicTwoSpinonMoments<double> model(1., options);
  uni20::async::TbbScheduler scheduler(2);
  auto pending = m::ThermodynamicTwoSpinonMoments<double>(1., options)
                     .global_moments(8, 4.)
                     .evaluate_async({.scheduler = &scheduler, .max_in_flight = 2});
  auto result = pending.get_wait(scheduler);
  auto direct = model.global_moments(8, 4.).evaluate({.max_in_flight = 1});
  EXPECT_EQ(result.integral.sums, direct.integral.sums);
  auto extended = model.global_moments(8, 8.).evaluate();
  ASSERT_TRUE(extended.integral.sums && result.integral.sums);
  for (int i = 0; i < 2; ++i)
    EXPECT_GT((*extended.integral.sums)[i], (*result.integral.sums)[i]);
  double const pi = 4 * std::atan(1.);
  for (double q : {0., 2 * pi})
  {
    auto zero = model.fixed_moments(q, 8, 4.).evaluate();
    ASSERT_EQ(zero.status, m::TwoSpinonMomentStatus::physical_zero);
    ASSERT_TRUE(zero.integral.sums);
    EXPECT_EQ((*zero.integral.sums)[0], 0.);
    EXPECT_EQ((*zero.integral.sums)[1], 0.);
    EXPECT_EQ(zero.integral.work_units, 0u);
    EXPECT_FALSE(zero.sum_rule_fractions[1]);
  }
  std::stop_source stop;
  stop.request_stop();
  auto cancelled = model.fixed_moments(pi / 2, 8, 4.).evaluate({.stop = stop.get_token()});
  EXPECT_EQ(cancelled.status, m::TwoSpinonMomentStatus::cancelled);
  EXPECT_FALSE(cancelled.integral.sums);
  options.max_evaluations = 1;
  auto failed = m::ThermodynamicTwoSpinonMoments<double>(1., options).global_moments(8, 4.).evaluate();
  EXPECT_EQ(failed.status, m::TwoSpinonMomentStatus::numerical_failure);
  EXPECT_FALSE(failed.integral.sums);
  EXPECT_FALSE(failed.sum_rule_fractions[0]);
  EXPECT_THROW(model.global_moments(8, 4., 8), std::length_error);
  EXPECT_THROW(model.fixed_moments(-1., 8, 4.), std::invalid_argument);
}

#if UNI20_HAS_MPREAL
TEST(TwoSpinonMomentsMpfr, RetainsExplicitPrecision)
{
  using R = uni20::mpreal;
  for (unsigned bits : {160u, 256u})
  {
    auto p = uni20::Precision::bits(bits);
    R const one(1, p), pi = 4 * uni20::math::atan(one);
    m::StructureFactorOptions<R> options(p);
    options.tolerance = R("1e-25", p);
    m::ThermodynamicTwoSpinonMoments<R> model(one, options);
    auto result = model.fixed_moments(pi / 2, 4, 4 * one).evaluate();
    ASSERT_TRUE(result.integral.sums);
    for (auto const& v : *result.integral.sums)
      EXPECT_EQ(v.precision(), p);
    auto scaled = m::ThermodynamicTwoSpinonMoments<R>(3 * one, options).fixed_moments(pi / 2, 4, 4 * one).evaluate();
    ASSERT_TRUE(scaled.integral.sums);
    R const tol = 256 * uni20::numeric_limits<R>::epsilon(one);
    EXPECT_REAL_NEAR((*scaled.integral.sums)[0], (*result.integral.sums)[0], tol);
    EXPECT_REAL_NEAR((*scaled.integral.sums)[1], 3 * (*result.integral.sums)[1], tol);
    EXPECT_THROW(model.global_moments(4, R(4, uni20::Precision::bits(bits + 16))), std::invalid_argument);
  }
}
#endif
