// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxx_four_spinon_thermodynamics.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using namespace bethe::heisenberg;
using Native = ThermodynamicFourSpinon<long double, double, long double>;
SpinonRecoveryPolicy native_policy()
{
  return {"1e-4", "1e-24", {{53, uni20::numeric_limits<long double>::digits(1.L), {{{64, 16}, {80, 20}}}}}, 60};
}
constexpr bool wider_long_double = std::numeric_limits<long double>::digits > std::numeric_limits<double>::digits;

TEST(FourSpinonThermodynamics, SpectralSupportAndNodeBudgetAreExplicit)
{
  if (!wider_long_double) GTEST_SKIP();
  Native model(native_policy(), 1.L);
  long double const pi = 4 * std::atan(1.L);
  auto outside = model.spectrum(pi / 2, 20.L, 2, 2.L).evaluate();
  ASSERT_EQ(outside.status, FourSpinonIntegralStatus::physical_zero);
  ASSERT_TRUE(outside.estimates);
  EXPECT_EQ((*outside.estimates)[0].value, 0.L);
  EXPECT_EQ(outside.density_attempts, 0u);
  auto edge = model.spectrum(pi / 2, pi / 2, 2, 2.L).evaluate();
  EXPECT_EQ(edge.status, FourSpinonIntegralStatus::threshold);
  EXPECT_FALSE(edge.estimates);
  auto budget = model.spectrum(pi, 2.L, 8, 16.L, 2).evaluate();
  EXPECT_EQ(budget.status, FourSpinonIntegralStatus::numerical_failure);
  EXPECT_FALSE(budget.estimates);
  EXPECT_EQ(budget.completed_nodes, 0u);
  EXPECT_THROW(model.spectrum(pi, 2.L, 0, 16.L), std::invalid_argument);
}

TEST(FourSpinonThermodynamics, HigherOrderSpectrumGeometryAndRefinementEntry)
{
  if (!wider_long_double) GTEST_SKIP();
  Native model(native_policy(), 1.L);
  long double const pi = 4 * std::atan(1.L);
  // These interior points failed during geometry construction with the old
  // radial map. Check both refinement orders and domains without evaluating
  // the expensive physical density at every node.
  for (auto const& [q, omega] :
       std::array<std::pair<long double, long double>, 3>{{{.5L, 3.117081L}, {pi, 4.4428825L}, {pi, 4.220739L}}})
    for (std::size_t order : {64u, 128u})
      for (long double cutoff : {16.L, 20.L})
      {
        SCOPED_TRACE(::testing::Message() << q << ", " << omega << ", " << order << ", " << cutoff);
        auto query = model.spectrum(q, omega, order, cutoff, 1000000);
        EXPECT_GT(query.size(), 0u);
      }

  // Exercise the refinement controller's first high-order geometry as well.
  // Cancellation deliberately stops before density work or the later controls;
  // their geometries are checked directly above, not claimed as integrated.
  std::stop_source stop;
  stop.request_stop();
  auto result =
      model.refined_spectrum(.5L, 3.117081L, 64, 128, 16.L, 20.L, 1000000).evaluate({nullptr, 8, 2, stop.get_token()});
  ASSERT_TRUE(result.controls[0]);
  EXPECT_GT(result.controls[0]->total_nodes, 0u);
  EXPECT_EQ(result.controls[0]->density_attempts, 0u);
  EXPECT_EQ(result.status, FourSpinonIntegralStatus::cancelled);
  EXPECT_FALSE(result.estimates);
  EXPECT_FALSE(result.controls[1]);
}

TEST(FourSpinonThermodynamics, FixedMomentsAndSpectraApplyExchangeInternally)
{
  if (!wider_long_double) GTEST_SKIP();
  Native one(native_policy(), 1.L), two(native_policy(), 2.L);
  long double const pi = 4 * std::atan(1.L);
  auto a = one.fixed_moments(pi / 2, 2, 2.L).evaluate();
  ASSERT_TRUE(a.estimates) << a.reason;
  auto b = two.fixed_moments(pi / 2, 2, 2.L).evaluate();
  ASSERT_TRUE(b.estimates) << b.reason;
  EXPECT_EQ(a.completed_nodes, a.total_nodes);
  EXPECT_FALSE(a.sum_rule_fractions[0]);
  ASSERT_TRUE(a.sum_rule_fractions[1]);
  ASSERT_TRUE(b.sum_rule_fractions[1]);
  EXPECT_REAL_NEAR(*a.sum_rule_fractions[1], *b.sum_rule_fractions[1], 1e-16L);
  for (int i = 0; i < 2; ++i)
  {
    auto const& x = (*a.estimates)[i];
    auto const& y = (*b.estimates)[i];
    long double const factor = i ? 2.L : 1.L;
    EXPECT_REAL_NEAR(y.value, factor * x.value, 1e-16L);
    EXPECT_REAL_NEAR(y.contour_diagnostic, factor * x.contour_diagnostic, 1e-16L);
    EXPECT_REAL_NEAR(y.arithmetic_diagnostic, factor * x.arithmetic_diagnostic, 1e-16L);
  }
  auto s1 = one.spectrum(pi, 2.L, 2, 2.L).evaluate();
  auto s2 = two.spectrum(pi, 4.L, 2, 2.L).evaluate();
  ASSERT_TRUE(s1.estimates) << s1.reason;
  ASSERT_TRUE(s2.estimates) << s2.reason;
  EXPECT_REAL_NEAR((*s1.estimates)[0].value, 2 * (*s2.estimates)[0].value, 1e-16L);
  EXPECT_TRUE(one.fixed_moments(pi, 2, 2.L).evaluate().unrestricted_static_divergence);
}

TEST(FourSpinonThermodynamics, AsyncOwnsItsQueryAndMatchesOrderedSerialReduction)
{
  if (!wider_long_double) GTEST_SKIP();
  Native model(native_policy(), 1.L);
  auto serial = model.global_moments(2, {2.L, 2.L, 2.L}).evaluate();
  ASSERT_TRUE(serial.estimates) << serial.reason;
  ASSERT_TRUE(serial.sum_rule_fractions[0]);
  ASSERT_TRUE(serial.sum_rule_fractions[1]);
  EXPECT_REAL_NEAR(*serial.sum_rule_fractions[0], 4 * (*serial.estimates)[0].value, 1e-18L);
  EXPECT_REAL_NEAR(*serial.sum_rule_denominators[1], 2.L / 3 * (std::log(2.L) - .25L), 1e-18L);
  auto zero = model.fixed_moments(0.L, 2, 2.L).evaluate();
  ASSERT_TRUE(zero.estimates);
  EXPECT_FALSE(zero.sum_rule_fractions[1]); // 0/0 is undefined.
  // Independent 8F7 density and direct positive common-shift integration;
  // two Gauss points in each gap on (0,2), not a converged global response.
  std::array<long double, 2> const reference{0.0001495861651180796436601650318071853595051L,
                                             0.0003800361467097579563009399447519737662089L};
  for (int i = 0; i < 2; ++i)
  {
    auto const& v = (*serial.estimates)[i];
    long double const diagnostic =
        v.contour_diagnostic + v.arithmetic_diagnostic + v.shift_quadrature_diagnostic + v.shift_tail_bound;
    EXPECT_REAL_NEAR(v.value, reference[i], diagnostic + 1e-20L);
    EXPECT_LT(diagnostic, 1e-6L * reference[i]);
  }
  uni20::async::TbbScheduler scheduler(2);
  auto pending = Native(native_policy(), 1.L).global_moments(2, {2.L, 2.L, 2.L}).evaluate_async({&scheduler, 8, 2});
  auto parallel = pending.get_wait(scheduler);
  ASSERT_TRUE(parallel.estimates) << parallel.reason;
  EXPECT_EQ(serial.total_nodes, 8u);
  EXPECT_EQ(parallel.completed_nodes, 8u);
  for (int i = 0; i < 2; ++i)
  {
    EXPECT_EQ((*serial.estimates)[i].value, (*parallel.estimates)[i].value);
    EXPECT_EQ((*serial.estimates)[i].contour_diagnostic, (*parallel.estimates)[i].contour_diagnostic);
    EXPECT_EQ((*serial.estimates)[i].arithmetic_diagnostic, (*parallel.estimates)[i].arithmetic_diagnostic);
  }
}

TEST(FourSpinonThermodynamics, CancellationAndFailedRecoveryDoNotReturnPrefixEstimates)
{
  if (!wider_long_double) GTEST_SKIP();
  Native model(native_policy(), 1.L);
  std::stop_source stop;
  stop.request_stop();
  auto cancelled = model.global_moments(2, {2.L, 2.L, 2.L}).evaluate({nullptr, 8, 2, stop.get_token()});
  EXPECT_EQ(cancelled.status, FourSpinonIntegralStatus::cancelled);
  EXPECT_FALSE(cancelled.estimates);
  EXPECT_EQ(cancelled.completed_nodes, 0u);
  auto policy = native_policy();
  policy.relative_tolerance = "1e-40";
  policy.absolute_tolerance = "0";
  uni20::async::TbbScheduler scheduler(2);
  auto failed = Native(policy, 1.L).global_moments(2, {2.L, 2.L, 2.L}).evaluate({&scheduler, 8, 2});
  EXPECT_EQ(failed.status, FourSpinonIntegralStatus::numerical_failure);
  EXPECT_FALSE(failed.estimates);
  EXPECT_TRUE(failed.failed_node);
  EXPECT_TRUE(failed.recovery_failure);
  EXPECT_LT(failed.completed_nodes, failed.total_nodes);
}

TEST(FourSpinonWeightedValue, SubnormalsAreNotRelativeAccuracyEstimates)
{
  EXPECT_THROW(detail::spinon_weighted_value(1., -740.), SpinonDensityError);
  EXPECT_EQ(detail::spinon_weighted_value(0., -740.), 0.);
}

TEST(FourSpinonThermodynamics, RefinementKeepsOrderAndCutoffChangesSeparate)
{
  if (!wider_long_double) GTEST_SKIP();
  Native model(native_policy(), 1.L);
  uni20::async::TbbScheduler scheduler(1);
  auto pending =
      model.refined_global_moments(2, 3, {1.L, 1.L, 1.L}, {2.L, 2.L, 2.L}).evaluate_async({&scheduler, 8, 2});
  auto result = pending.get_wait(scheduler);
  ASSERT_TRUE(result.estimates) << result.reason;
  ASSERT_TRUE(result.quadrature_difference);
  ASSERT_TRUE(result.cutoff_change);
  for (auto const& control : result.controls)
  {
    ASSERT_TRUE(control);
    ASSERT_TRUE(control->estimates);
    EXPECT_EQ(control->completed_nodes, control->total_nodes);
  }
  for (int i = 0; i < 2; ++i)
  {
    auto value = [&](std::size_t n) { return (*result.controls[n]->estimates)[i].value; };
    EXPECT_EQ((*result.estimates)[i].value, value(3));
    EXPECT_EQ((*result.quadrature_difference)[i], std::abs(value(3) - value(2)));
    EXPECT_EQ((*result.cutoff_change)[i], value(3) - value(1));
  }
  EXPECT_THROW(model.refined_global_moments(2, 3, {1.L, 1.L, 1.L}, {2.L, 2.L, 2.L}, 1).evaluate(),
               std::invalid_argument);
  auto bad = native_policy();
  bad.relative_tolerance = "1e-40";
  bad.absolute_tolerance = "0";
  auto failed = Native(bad, 1.L).refined_global_moments(2, 3, {1.L, 1.L, 1.L}, {2.L, 2.L, 2.L}).evaluate();
  EXPECT_FALSE(failed.estimates);
  EXPECT_TRUE(failed.controls[0]);
  EXPECT_FALSE(failed.controls[1]);
  EXPECT_THROW(model.refined_spectrum(1.L, 2.L, 4, 2, 2.L, 4.L), std::invalid_argument);
}

#if UNI20_HAS_MPREAL
TEST(FourSpinonThermodynamicsMpfr, GlobalMomentsRetainGeometryPrecisionAndIndependentMeasures)
{
  using R = uni20::mpreal;
  SpinonRecoveryPolicy policy{"1e-6", "1e-30", {{128, 192, {{{64, 16}, {80, 20}}}}}, 60};
  std::optional<std::array<R, 2>> previous;
  for (int bits : {192, 256})
  {
    auto const p = uni20::Precision::bits(bits);
    R const one(1, p);
    ThermodynamicFourSpinon<R, R> model(policy, one);
    auto result = model.global_moments(2, {2 * one, 2 * one, 2 * one}).evaluate();
    ASSERT_TRUE(result.estimates) << result.reason;
    std::array<R, 2> const reference{R("0.0001495861651180796436601650318071853595051", p),
                                     R("0.0003800361467097579563009399447519737662089", p)};
    std::array<R, 2> current;
    for (int i = 0; i < 2; ++i)
    {
      auto const& value = (*result.estimates)[i];
      EXPECT_EQ(value.value.precision(), p);
      EXPECT_REAL_NEAR(value.value, reference[i],
                       R(value.contour_diagnostic + value.arithmetic_diagnostic + value.shift_quadrature_diagnostic +
                         value.shift_tail_bound + R("1e-40", p)));
      current[i] = value.value;
      if (previous) EXPECT_LT(uni20::math::abs(current[i] - (*previous)[i].at(p)), R("1e-45", p) * current[i]);
    }
    previous = current;
  }
}

TEST(FourSpinonThermodynamicsMpfr, PhysicalSliceWithIndependentKernelReference)
{
  using R = uni20::mpreal;
  auto const p = uni20::Precision::bits(192);
  R const one(1, p), pi = 4 * uni20::math::atan(one);
  SpinonRecoveryPolicy policy{"1e-6", "1e-30", {{128, 192, {{{64, 16}, {80, 20}}}}}, 60};
  ThermodynamicFourSpinon<R, R> model(policy, one);
  uni20::async::TbbScheduler scheduler(2);
  auto result = model.spectrum(pi, R(2, p), 8, R(16, p)).evaluate({&scheduler, 8, 2});
  ASSERT_TRUE(result.estimates) << result.reason;
  ASSERT_EQ(result.status, FourSpinonIntegralStatus::estimated);
  EXPECT_EQ(result.completed_nodes, 512u);
  auto const& value = (*result.estimates)[0];
  EXPECT_EQ(value.value.precision(), p);
  // Independent 8F7 density and independently parameterized finite quadrature.
  R const reference("0.0981022159286235520475", p);
  EXPECT_REAL_NEAR(value.value, reference, R(value.contour_diagnostic + value.arithmetic_diagnostic + R("2e-10", p)));
  EXPECT_LT(value.contour_diagnostic + value.arithmetic_diagnostic, R("1e-7", p));
  EXPECT_THROW(model.global_moments(2, {one, one.at(uni20::Precision::bits(160)), one}), std::invalid_argument);
}
#endif
} // namespace
