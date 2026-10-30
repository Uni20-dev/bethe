// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxx_spinon_density_recovery.hpp>

using namespace bethe::heisenberg;

namespace
{
std::vector<std::string> const two{"-2", "1"};
std::vector<std::string> const bulk{"-2.1", "-1.3", "-.4", ".3", "1.2", "2"};
std::vector<std::string> const collision{"-2", "-.6", "-.59999999999999999999", ".4", "1.3", "2.1"};

SpinonRecoveryPolicy native_policy()
{
  return {"1e-6",
          "0",
          {{uni20::numeric_limits<double>::digits(1.),
            uni20::numeric_limits<long double>::digits(1.L),
            {{{32, 12}, {64, 20}}}}},
          60};
}
} // namespace

TEST(SpinonRecovery, DecimalIdentityAndMalformedInputs)
{
  using bethe::heisenberg::detail::spinon_decimal_key;
  EXPECT_EQ(spinon_decimal_key("+0.6000"), spinon_decimal_key("6e-1"));
  EXPECT_EQ(spinon_decimal_key("-0e99"), spinon_decimal_key("0.0"));
  EXPECT_NE(spinon_decimal_key("-.6"), spinon_decimal_key("-.59999999999999999999"));
  for (auto s : {"", "+", "-", "nan", "inf", "0x1p0", "1/2", ".", "1e", "1e+", "1e+-1", "1e--1", "1..2",
                 "1e999999999999999999999999", "0.1e-9223372036854775808", "10e9223372036854775807"})
    EXPECT_THROW(spinon_decimal_key(s), std::invalid_argument) << s;
}

TEST(SpinonRecovery, NativeAcceptanceAndExplicitExhaustion)
{
  if (uni20::numeric_limits<long double>::digits(1.L) <= uni20::numeric_limits<double>::digits(1.))
    GTEST_SKIP() << "no distinct native low/guard pair";
  XXXSpinonDensityRecovery<double, long double> recovery(native_policy());
  auto result = recovery(two);
  ASSERT_TRUE(result.accepted()) << result.reason;
  ASSERT_TRUE(result.estimate);
  ASSERT_TRUE(result.arithmetic_screen);
  EXPECT_EQ(result.attempts.size(), 2u);
  EXPECT_GT(*result.arithmetic_screen, 0.L);
  for (auto const& a : result.attempts)
  {
    ASSERT_TRUE(a.checked);
    EXPECT_TRUE(a.checked->accepted());
    EXPECT_GT(a.checked->precision_difference, 0.L);
  }
  XXXSpinonDensity<long double> direct;
  std::array b{direct.parse("-2"), direct.parse("1")};
  EXPECT_EQ(result.estimate->value, direct.evaluate(b, {64, 20}).value);
  auto strict = native_policy();
  strict.relative_tolerance = "1e-60";
  auto failed = XXXSpinonDensityRecovery<double, long double>(strict)(two);
  EXPECT_EQ(failed.status, SpinonRecoveryStatus::exhausted);
  EXPECT_EQ(failed.attempts.size(), 2u);
  EXPECT_FALSE(failed.estimate);
  EXPECT_FALSE(failed.arithmetic_screen);
  EXPECT_EQ(failed.attempts.back().failure, SpinonDensityFailure::arithmetic_precision);
  std::vector<std::string> tiny{"0", "1e-400"};
  auto unrepresentable = recovery(tiny);
  EXPECT_EQ(unrepresentable.status, SpinonRecoveryStatus::exhausted);
  EXPECT_FALSE(unrepresentable.estimate);
  EXPECT_EQ(unrepresentable.attempts.front().failure, SpinonDensityFailure::exponent_range);
}

TEST(SpinonRecovery, ValidatePoliciesBeforeKernels)
{
  using Driver = XXXSpinonDensityRecovery<double, long double>;
  auto p = native_policy();
  if (p.stages.front().guard_bits <= p.stages.front().bits) GTEST_SKIP();
  auto invalid = p;
  invalid.stages.clear();
  EXPECT_THROW(Driver{invalid}, std::invalid_argument);
  invalid = p;
  invalid.stages[0].bits = 52;
  EXPECT_THROW(Driver{invalid}, std::invalid_argument);
  invalid = p;
  invalid.stages[0].guard_bits = 53;
  EXPECT_THROW(Driver{invalid}, std::invalid_argument);
  invalid = p;
  invalid.stages[0].contours[1] = invalid.stages[0].contours[0];
  EXPECT_THROW(Driver{invalid}, std::invalid_argument);
  invalid = p;
  invalid.stages.push_back(p.stages.front());
  EXPECT_THROW(Driver{invalid}, std::invalid_argument);
  invalid = p;
  invalid.relative_tolerance = "-1";
  EXPECT_THROW(Driver{invalid}, std::invalid_argument);
  invalid.relative_tolerance = "0";
  EXPECT_THROW(Driver{invalid}, std::invalid_argument);
  invalid = p;
  invalid.max_wall_seconds = std::numeric_limits<double>::infinity();
  EXPECT_THROW(Driver{invalid}, std::invalid_argument);
}

TEST(SpinonRecovery, InvalidInputCancellationAndBudgetAreNotEstimates)
{
  auto p = native_policy();
  if (p.stages.front().guard_bits <= p.stages.front().bits) GTEST_SKIP();
  XXXSpinonDensityRecovery<double, long double> recovery(p);
  for (auto const& b : {std::vector<std::string>{"nan", "1"}, std::vector<std::string>{"1"}})
  {
    auto result = recovery(b);
    EXPECT_EQ(result.status, SpinonRecoveryStatus::invalid_input);
    EXPECT_TRUE(result.attempts.empty());
    EXPECT_FALSE(result.estimate);
  }
  std::vector<std::string> same{".6", "6e-1"};
  EXPECT_EQ(recovery(same).status, SpinonRecoveryStatus::unsupported_limit);
  std::stop_source stop;
  stop.request_stop();
  auto cancelled = recovery(two, stop.get_token());
  EXPECT_EQ(cancelled.status, SpinonRecoveryStatus::cancelled);
  EXPECT_TRUE(cancelled.attempts.empty());
  EXPECT_FALSE(cancelled.estimate);
  p.max_wall_seconds = 1e-300;
  auto timed = XXXSpinonDensityRecovery<double, long double>(p)(two);
  EXPECT_EQ(timed.status, SpinonRecoveryStatus::time_limit);
  EXPECT_FALSE(timed.estimate);
  EXPECT_THROW(recovery.recover([](auto const& c) -> std::vector<decltype(c.parse("0"))> {
    throw std::logic_error("a programming error is not a numerical retry");
  }),
               std::logic_error);
}

TEST(SpinonRecovery, ContradictoryPassingControlsAreTerminal)
{
  auto p = native_policy();
  if (p.stages.front().guard_bits <= p.stages.front().bits) GTEST_SKIP();
  XXXSpinonDensityRecovery<double, long double> recovery(p);
  int calls = 0;
  // Inject contradictory controls by changing the input between contour
  // budgets. Each low/guard pair agrees in isolation; the final check must not
  // discard the earlier passing result. A production source never does this.
  auto result = recovery.recover([&](auto const& context) {
    bool const first = calls++ < 2;
    return std::vector{context.parse("-2"), context.parse(first ? "1" : "2")};
  });
  EXPECT_EQ(result.status, SpinonRecoveryStatus::inconsistent_controls);
  ASSERT_EQ(result.conflicting_attempts.size(), 1u);
  EXPECT_EQ(result.conflicting_attempts[0], 0u);
  EXPECT_FALSE(result.estimate);
  EXPECT_FALSE(result.arithmetic_screen);
}

#if UNI20_HAS_FLOAT128
TEST(SpinonRecovery, NativeBinary128Guard)
{
  using R = uni20::float128;
  int bits = uni20::numeric_limits<long double>::digits(1.L), guard = uni20::numeric_limits<R>::digits(R{1});
  if (bits >= guard) GTEST_SKIP() << "no distinct wider native guard";
  SpinonRecoveryPolicy p{"1e-12", "0", {{bits, guard, {{{32, 12}, {64, 20}}}}}, 60};
  auto result = XXXSpinonDensityRecovery<long double, R>(p)(two);
  ASSERT_TRUE(result.accepted()) << result.reason;
  EXPECT_GT(*result.arithmetic_screen, R{0});
  EXPECT_LT(*result.arithmetic_screen / result.estimate->value, uni20::parse_real<R>("1e-12"));
}
#endif

#if UNI20_HAS_MPREAL
namespace
{
using MpDriver = XXXSpinonDensityRecovery<uni20::mpreal>;
SpinonRecoveryPolicy mp_policy()
{
  return {"1e-8", "0", {{120, 190, {{{32, 12}, {64, 20}}}}, {170, 270, {{{64, 20}, {128, 24}}}}}, 60};
}
} // namespace

TEST(SpinonRecoveryMpfr, BulkAndNearCollisionRequireDifferentStages)
{
  MpDriver recovery(mp_policy());
  auto ordinary = recovery(bulk);
  ASSERT_TRUE(ordinary.accepted()) << ordinary.reason;
  EXPECT_EQ(ordinary.attempts.size(), 2u);
  EXPECT_EQ(ordinary.estimate->value.precision(), uni20::Precision::bits(190));
  auto recovered = recovery(collision);
  ASSERT_TRUE(recovered.accepted()) << recovered.reason;
  EXPECT_EQ(recovered.attempts.size(), 4u);
  EXPECT_NE(recovered.attempts.front().failure, SpinonDensityFailure::none);
  EXPECT_EQ(recovered.estimate->value.precision(), uni20::Precision::bits(270));
  auto limited = mp_policy();
  limited.stages.resize(1);
  auto failed = MpDriver(limited)(collision);
  EXPECT_EQ(failed.status, SpinonRecoveryStatus::exhausted);
  EXPECT_EQ(failed.attempts.size(), 2u);
  EXPECT_FALSE(failed.estimate);
}

TEST(SpinonRecoveryMpfr, RoundedCollisionReconstructsOriginalInputs)
{
  auto p = mp_policy();
  p.stages[0].bits = 60;
  p.stages[0].guard_bits = 100;
  auto recovered = MpDriver(p)(collision);
  ASSERT_TRUE(recovered.accepted()) << recovered.reason;
  EXPECT_EQ(recovered.attempts.size(), 4u);
  EXPECT_EQ(recovered.attempts.front().failure, SpinonDensityFailure::rounded_collision);
  EXPECT_FALSE(recovered.attempts.front().checked);
}

TEST(SpinonRecoveryMpfr, TailFailureAdvancesOnlyDeclaredBudgets)
{
  auto p = mp_policy();
  p.stages = {{120, 190, {{{2, 1}, {4, 2}}}}, {120, 190, {{{32, 12}, {64, 20}}}}};
  auto result = MpDriver(p)(bulk);
  ASSERT_TRUE(result.accepted()) << result.reason;
  ASSERT_EQ(result.attempts.size(), 4u);
  EXPECT_TRUE(result.attempts.front().failure == SpinonDensityFailure::tail_domain ||
              result.attempts.front().failure == SpinonDensityFailure::contour_tail);
  p.stages.resize(1);
  auto failed = MpDriver(p)(bulk);
  EXPECT_EQ(failed.status, SpinonRecoveryStatus::exhausted);
  EXPECT_FALSE(failed.estimate);
}

TEST(SpinonRecoveryMpfr, IndependentFourSpinonNormalizationAndPrecisionGrowth)
{
  auto p = mp_policy();
  p.stages = {{160, 240, {{{64, 16}, {128, 24}}}}};
  auto source = [](auto const& c) {
    auto pi = bethe::detail::pi_like(c.parse("1"));
    std::vector<decltype(pi)> beta;
    for (auto s : {"-.6", "-.1", ".3", ".8"})
      beta.push_back(2 * pi * c.parse(s));
    return beta;
  };
  auto first = MpDriver(p).recover(source);
  ASSERT_TRUE(first.accepted()) << first.reason;
  // Independent 8F7 contour sums at 60/85 decimal digits; regenerated by
  // reference_xxx_spinon_recovery.py. The older oscillatory-quadrature fixture
  // was only suitable for its original 1e-32 comparison allowance.
  auto expected = uni20::mpreal("0.044857368060020152745531125322084827034805784163101572458260824621140092618698817",
                                uni20::Precision::bits(240));
  EXPECT_LE(uni20::math::abs(first.estimate->value - expected),
            first.estimate->contour_tail_error + *first.arithmetic_screen);
  EXPECT_LT(*first.arithmetic_screen, expected * uni20::mpreal("1e-30", expected.precision()));
  p.stages[0].bits = 240;
  p.stages[0].guard_bits = 360;
  auto higher = MpDriver(p).recover(source);
  ASSERT_TRUE(higher.accepted()) << higher.reason;
  auto import = [](auto const& x) { return uni20::mpreal(x, uni20::Precision::bits(360)); };
  EXPECT_LT(*higher.arithmetic_screen, import(*first.arithmetic_screen));
  EXPECT_LE(uni20::math::abs(higher.estimate->value - import(first.estimate->value)),
            higher.estimate->contour_tail_error + import(first.estimate->contour_tail_error) +
                *higher.arithmetic_screen + import(*first.arithmetic_screen));
}
#endif
