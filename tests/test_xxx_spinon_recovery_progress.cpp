// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxx_spinon_density_recovery.hpp>

using namespace bethe::heisenberg;
namespace
{
struct SavedControl
{};
using Driver = XXXSpinonDensityRecovery<double, long double>;
std::vector<std::string> const beta{"-2", "1"};

SpinonRecoveryPolicy policy()
{
  return {"1e-6",
          "0",
          {{uni20::numeric_limits<double>::digits(1.),
            uni20::numeric_limits<long double>::digits(1.L),
            {{{32, 12}, {64, 20}}}}},
          60};
}

bool native_pair() { return uni20::numeric_limits<long double>::digits(1.L) > 53; }

TEST(SpinonRecoveryProgress, PartialPrefixAndCompletedReplayDoNotRepeatKernels)
{
  if (!native_pair()) GTEST_SKIP();
  Driver recovery(policy());
  Driver::Result saved;
  int calls = 0, publications = 0;
  auto source = [&](auto const& context) {
    ++calls;
    return std::vector{context.parse("-2"), context.parse("1")};
  };
  EXPECT_THROW(recovery.recover(source, {},
                                [&](auto const& progress) {
                                  saved = progress;
                                  ++publications;
                                  throw SavedControl{};
                                }),
               SavedControl);
  ASSERT_EQ(saved.attempts.size(), 1u);
  EXPECT_EQ(calls, 2);
  EXPECT_EQ(publications, 1);
  EXPECT_GT(saved.attempts[0].completed_wall_seconds, 0);
  auto result = recovery.recover(source, saved, [&](auto const& progress) {
    ++publications;
    EXPECT_EQ(progress.attempts.size(), 2u);
  });
  ASSERT_TRUE(result.accepted()) << result.reason;
  EXPECT_EQ(calls, 4);
  EXPECT_EQ(publications, 2);
  EXPECT_GE(result.wall_seconds, saved.wall_seconds);
  auto fresh = recovery(beta);
  ASSERT_TRUE(fresh.estimate);
  EXPECT_EQ(result.estimate->value, fresh.estimate->value);
  EXPECT_EQ(result.arithmetic_screen, fresh.arithmetic_screen);
  auto replay =
      recovery.recover(source, result, [&](auto const&) { ADD_FAILURE() << "completed replay published work"; });
  EXPECT_EQ(replay.estimate->value, result.estimate->value);
  EXPECT_EQ(calls, 4);
}

TEST(SpinonRecoveryProgress, DecimalIdentityCancellationAndUntrustedTopLevelFlags)
{
  if (!native_pair()) GTEST_SKIP();
  Driver recovery(policy());
  std::stop_source stop;
  auto partial = recovery.resume(beta, {}, [&](auto const&) { stop.request_stop(); }, stop.get_token());
  ASSERT_EQ(partial.attempts.size(), 1u);
  EXPECT_EQ(partial.status, SpinonRecoveryStatus::cancelled);
  EXPECT_FALSE(partial.estimate);
  EXPECT_EQ(partial.decimal_inputs, beta);
  // A forged top-level acceptance flag does not replace a missing control.
  auto fake = partial;
  fake.status = SpinonRecoveryStatus::accepted;
  fake.estimate = partial.attempts[0].checked->estimate;
  auto checked = recovery.restore(fake);
  EXPECT_FALSE(checked.accepted());
  EXPECT_FALSE(checked.estimate);
  int publications = 0;
  auto done = recovery.resume(beta, fake, [&](auto const&) { ++publications; });
  EXPECT_TRUE(done.accepted());
  EXPECT_EQ(publications, 1);
  std::vector<std::string> other{"-2.0", "1"};
  EXPECT_THROW(recovery.resume(other, partial, [](auto const&) {}), std::invalid_argument);
  other = {"-2", "2"};
  EXPECT_THROW(recovery.resume(other, partial, [](auto const&) {}), std::invalid_argument);
  fake.policy.stages.clear();
  EXPECT_THROW(recovery.resume(beta, fake, [](auto const&) {}), std::invalid_argument);
  EXPECT_THROW(recovery.recover([](auto const& c) { return std::vector{c.parse("-2"), c.parse("1")}; }, fake,
                                [](auto const&) {}),
               std::invalid_argument);
}

TEST(SpinonRecoveryProgress, PolicyOrderingDiagnosticsAndBudgetCannotBeReset)
{
  if (!native_pair()) GTEST_SKIP();
  Driver recovery(policy());
  auto saved = recovery(beta);
  ASSERT_TRUE(saved.accepted());
  auto reject = [&](auto mutate) {
    auto bad = saved;
    mutate(bad);
    EXPECT_THROW(recovery.restore(bad), std::invalid_argument);
  };
  reject([](auto& x) { x.policy.max_wall_seconds *= 2; });
  reject([](auto& x) { x.policy.relative_tolerance = "1e-5"; });
  reject([](auto& x) { ++x.policy.stages[0].guard_bits; });
  reject([](auto& x) { ++x.policy.stages[0].contours[0].tail_order; });
  reject([](auto& x) { x.wall_seconds = 0; });
  reject([](auto& x) { x.wall_seconds = std::numeric_limits<double>::infinity(); });
  reject([](auto& x) { x.attempts[1].completed_wall_seconds = -1; });
  reject([](auto& x) { x.attempts[0].stage = 1; });
  reject([](auto& x) { x.attempts[0].control = 1; });
  reject([](auto& x) { x.attempts[0].checked.reset(); });
  reject([](auto& x) { x.attempts[0].failure = SpinonDensityFailure::tail_domain; });
  reject([](auto& x) { x.attempts[0].checked->failure = SpinonDensityFailure::arithmetic_precision; });
  reject([](auto& x) { x.attempts[0].checked->requested_error *= 2; });
  reject([](auto& x) { x.attempts[0].checked->arithmetic_screen *= 2; });
  reject([](auto& x) { x.attempts[0].checked->low_arithmetic_indicator = -1; });
  reject([](auto& x) { x.attempts[0].checked->estimate.value = 0; });
  reject([](auto& x) { ++x.attempts[0].checked->estimate.terms; });
  reject([](auto& x) { x.attempts[0].checked->estimate.amplitudes.clear(); });
  reject([](auto& x) { x.decimal_inputs = {"-2", "1", "2", "3"}; });
  reject([](auto& x) { x.decimal_inputs = {"-2", "-2.0"}; });
  saved.wall_seconds = policy().max_wall_seconds;
  auto expired = recovery.resume(beta, saved, [](auto const&) { ADD_FAILURE() << "expired history ran control"; });
  EXPECT_EQ(expired.status, SpinonRecoveryStatus::time_limit);
  EXPECT_FALSE(expired.estimate);
  EXPECT_FALSE(expired.arithmetic_screen);
}

TEST(SpinonRecoveryProgress, EarlierPassingControlRemainsAContradictionAfterRestart)
{
  if (!native_pair()) GTEST_SKIP();
  auto p = policy();
  auto second = p.stages[0];
  second.contours = {{{96, 22}, {128, 24}}};
  p.stages.push_back(second);
  Driver recovery(p);
  int calls = 0;
  Driver::Result saved;
  EXPECT_THROW(recovery.recover(
                   [&](auto const& c) {
                     if (++calls == 3)
                       throw SpinonDensityError(SpinonDensityFailure::tail_domain, "injected control failure");
                     return std::vector{c.parse("-2"), c.parse("1")};
                   },
                   {},
                   [&](auto const& progress) {
                     if (progress.attempts.size() == 2)
                     {
                       saved = progress;
                       throw SavedControl{};
                     }
                   }),
               SavedControl);
  ASSERT_EQ(saved.attempts.size(), 2u);
  ASSERT_TRUE(saved.attempts[0].checked->accepted());
  EXPECT_EQ(saved.attempts[1].failure, SpinonDensityFailure::tail_domain);
  // Deliberately violate input stability to inject later contradictory values;
  // this is an acceptance test, not a legitimate production source.
  auto inconsistent = recovery.recover([](auto const& c) { return std::vector{c.parse("-2"), c.parse("2")}; }, saved,
                                       [](auto const&) {});
  EXPECT_EQ(inconsistent.status, SpinonRecoveryStatus::inconsistent_controls);
  EXPECT_FALSE(inconsistent.estimate);
  EXPECT_EQ(inconsistent.conflicting_attempts, std::vector<std::size_t>{0});
  auto replay = recovery.recover(
      [&](auto const& c) {
        ADD_FAILURE() << "inconsistent history ran another kernel";
        return std::vector{c.parse("-2"), c.parse("2")};
      },
      inconsistent, [](auto const&) {});
  EXPECT_EQ(replay.status, SpinonRecoveryStatus::inconsistent_controls);
  EXPECT_EQ(replay.conflicting_attempts, inconsistent.conflicting_attempts);
}

TEST(SpinonRecoveryProgress, FourAndSixSpinonEvidenceRetainsItsPolarizations)
{
  if (!native_pair()) GTEST_SKIP();
  Driver recovery(policy());
  for (std::vector<std::string> const point : {std::vector<std::string>{"-2.1", "-.4", ".3", "2"},
                                               std::vector<std::string>{"-2.1", "-1.3", "-.4", ".3", "1.2", "2"}})
  {
    SCOPED_TRACE(point.size());
    Driver::Result saved;
    EXPECT_THROW(recovery.resume(point, {},
                                 [&](auto const& progress) {
                                   saved = progress;
                                   throw SavedControl{};
                                 }),
                 SavedControl);
    ASSERT_EQ(saved.attempts.size(), 1u);
    ASSERT_TRUE(saved.attempts.front().checked);
    EXPECT_EQ(saved.attempts.front().checked->estimate.amplitudes.size(), point.size() == 4 ? 4u : 15u);
    int publications = 0;
    auto resumed = recovery.resume(point, saved, [&](auto const&) { ++publications; });
    auto fresh = recovery(point);
    EXPECT_EQ(publications, 1);
    EXPECT_EQ(resumed.status, fresh.status);
    EXPECT_EQ(resumed.density_failure, fresh.density_failure);
    EXPECT_EQ(resumed.conflicting_attempts, fresh.conflicting_attempts);
    ASSERT_EQ(resumed.attempts.size(), fresh.attempts.size());
    for (std::size_t i = 0; i < resumed.attempts.size(); ++i)
    {
      auto const& a = resumed.attempts[i];
      auto const& b = fresh.attempts[i];
      EXPECT_EQ(a.failure, b.failure);
      ASSERT_TRUE(a.checked && b.checked);
      EXPECT_EQ(a.checked->estimate.value, b.checked->estimate.value);
      EXPECT_EQ(a.checked->arithmetic_screen, b.checked->arithmetic_screen);
      EXPECT_EQ(a.checked->estimate.amplitudes, b.checked->estimate.amplitudes);
    }
    auto restored = recovery.restore(resumed);
    EXPECT_EQ(restored.status, fresh.status);
    EXPECT_EQ(restored.density_failure, fresh.density_failure);
  }
}

#if UNI20_HAS_FLOAT128
TEST(SpinonRecoveryProgress, NativeBinary128EvidenceRoundTrip)
{
  using R = uni20::float128;
  using D = XXXSpinonDensityRecovery<long double, R>;
  int const low = uni20::numeric_limits<long double>::digits(1.L);
  int const high = uni20::numeric_limits<R>::digits(R{1});
  if (low >= high) GTEST_SKIP();
  SpinonRecoveryPolicy p{"1e-12", "0", {{low, high, {{{32, 12}, {64, 20}}}}}, 60};
  D recovery(p);
  D::Result partial;
  EXPECT_THROW(recovery.resume(beta, {},
                               [&](auto const& x) {
                                 partial = x;
                                 throw SavedControl{};
                               }),
               SavedControl);
  ASSERT_EQ(partial.attempts.size(), 1u);
  auto done = recovery.resume(beta, partial, [](auto const&) {});
  ASSERT_TRUE(done.accepted()) << done.reason;
  EXPECT_EQ(recovery.restore(done).estimate->value, done.estimate->value);
}
#endif

#if UNI20_HAS_MPREAL
TEST(SpinonRecoveryProgressMpfr, MixedStagePrecisionsSurviveRestart)
{
  using R = uni20::mpreal;
  using D = XXXSpinonDensityRecovery<R>;
  SpinonRecoveryPolicy p{"1e-8", "0", {{80, 128, {{{32, 12}, {64, 20}}}}, {160, 256, {{{64, 20}, {128, 24}}}}}, 60};
  D recovery(p);
  D::Result saved;
  int calls = 0;
  EXPECT_THROW(recovery.recover(
                   [&](auto const& c) {
                     if (++calls == 1)
                       throw SpinonDensityError(SpinonDensityFailure::tail_domain, "injected control failure");
                     return std::vector{c.parse("-2"), c.parse("1")};
                   },
                   {},
                   [&](auto const& progress) {
                     if (progress.attempts.size() == 2)
                     {
                       saved = progress;
                       throw SavedControl{};
                     }
                   }),
               SavedControl);
  ASSERT_EQ(saved.attempts.size(), 2u);
  auto result = recovery.recover(
      [&](auto const& c) {
        ++calls;
        return std::vector{c.parse("-2"), c.parse("1")};
      },
      saved, [](auto const&) {});
  ASSERT_TRUE(result.accepted()) << result.reason;
  EXPECT_EQ(calls, 7);
  EXPECT_EQ(result.attempts[1].checked->estimate.value.precision(), uni20::Precision::bits(128));
  EXPECT_EQ(result.attempts[1].checked->low_value.precision(), uni20::Precision::bits(128));
  EXPECT_EQ(result.estimate->value.precision(), uni20::Precision::bits(256));
  auto invalid = result;
  invalid.attempts[1].checked->low_value = invalid.attempts[1].checked->low_value.at(uni20::Precision::bits(256));
  EXPECT_THROW(recovery.restore(invalid), std::invalid_argument);
  invalid = result;
  invalid.attempts[1].checked->estimate.log_scale = R(0);
  EXPECT_THROW(recovery.restore(invalid), std::invalid_argument);
  EXPECT_EQ(recovery.restore(result).estimate->value, result.estimate->value);
}
#endif
} // namespace
