// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxx_six_spinon_moments.hpp>
#include <uni20/async/tbb_scheduler.hpp>

#if UNI20_HAS_MPREAL
namespace
{
using R = uni20::mpreal;
using Moments = bethe::heisenberg::ThermodynamicSixSpinonMoments<R, R>;
using Kind = bethe::heisenberg::detail::SpinonProposalKind;
auto proposal(R const& one)
{
  return bethe::heisenberg::detail::SpinonJointProposal<R>(8 * one, one / 2,
                                                           {{Kind::background, one / 2, one, one},
                                                            {Kind::one_pair, one / 4, 4 * one, one / 2},
                                                            {Kind::two_pairs, one / 4, 8 * one, one}});
}
auto policy()
{
  return bethe::heisenberg::SpinonRecoveryPolicy{
      "1e-7", "0", {{128, 192, {{{32, 12}, {64, 20}}}}, {224, 320, {{{64, 20}, {128, 24}}}}}, 60};
}
} // namespace
TEST(SixSpinonMomentsMpfr, IndependentPhysicalGlobalObservation)
{
  R one(1, uni20::Precision::bits(384));
  Moments model(proposal(one), policy(), one);
  bethe::UniformSample<7>::Words words;
  std::array<unsigned, 7> raw{4, 32, 80, 112, 180, 220, 240};
  for (int i = 0; i < 7; ++i)
    words[i] = {raw[i]};
  auto result = model.global_observation({8, words});
  ASSERT_TRUE(result.values) << result.reason;
  // Independent positive-PDF/global measure and 85-dps Python W6, 128 terms/order24.
  std::array<char const*, 2> refs{"0.0000039404529116124723724485546556222387959072271381721554310630449464825055242",
                                  "0.0000182628663917981358830789296704061304908479313649550212512434727095286472474"};
  std::array<char const*, 2> tails{"9.13713663096281799777082620216005582236418835824615052921959851972246693727e-55",
                                   "4.23480014196118164189326281531325148328548307801851813280160987194681328954e-54"};
  for (int order = 0; order < 2; ++order)
  {
    auto const& v = *result.values;
    R expected(refs[order], one.precision()), tail(tails[order], one.precision());
    EXPECT_LE(uni20::math::abs(v[order] - expected), v[2 + order] + v[4 + order] + tail);
    EXPECT_LT(v[2 + order] + v[4 + order], expected * R("1e-7", one.precision()));
  }
}

TEST(SixSpinonMomentsMpfr, JointMomentsExchangeScalingAndOwnedAsync)
{
  R one(1, uni20::Precision::bits(384));
  Moments model(proposal(one), policy(), one);
  bethe::MonteCarloOptions options{2, 4, 42};
  auto a = model.global(options, {.max_in_flight = 1});
  std::string failures;
  for (auto const& row : a.integral.replicates)
    if (row && !row->means) failures += row->reason + "; ";
  ASSERT_EQ(a.integral.status, bethe::MonteCarloStatus::complete) << a.integral.reason << failures;
  ASSERT_TRUE(a.integral.mean && a.integral.covariance && a.fraction);
  EXPECT_EQ(a.integral.completed_samples, 8u);
  EXPECT_EQ(a.integral.metadata.supported_configurations, 8u);
  EXPECT_GT((*a.integral.mean)[0], 0);
  EXPECT_GT((*a.integral.mean)[1], 0);
  EXPECT_GT((*a.integral.covariance)[1], 0);
  uni20::async::TbbScheduler scheduler(2);
  auto pending =
      Moments(proposal(one), policy(), one).global_async(options, {.scheduler = &scheduler, .max_in_flight = 2});
  auto b = pending.get_wait(scheduler);
  EXPECT_EQ(a.integral.mean, b.integral.mean);
  EXPECT_EQ(a.integral.covariance, b.integral.covariance);
  auto scaled =
      Moments(proposal(one), policy(), 3 * one).global(options, {.scheduler = &scheduler, .max_in_flight = 2});
  ASSERT_TRUE(scaled.integral.mean && scaled.integral.covariance && scaled.fraction) << scaled.integral.reason;
  R tol("1e-100", one.precision());
  EXPECT_EQ((*scaled.integral.mean)[0], (*a.integral.mean)[0]);
  EXPECT_REAL_NEAR((*scaled.integral.mean)[1], 3 * (*a.integral.mean)[1], tol * (*a.integral.mean)[1]);
  EXPECT_REAL_NEAR((*scaled.integral.covariance)[7], 9 * (*a.integral.covariance)[7],
                   tol * (*a.integral.covariance)[7]);
  EXPECT_EQ((*scaled.fraction)[0], (*a.fraction)[0]);
  EXPECT_REAL_NEAR((*scaled.fraction)[1], (*a.fraction)[1], tol * (*a.fraction)[1]);
  RecordProperty("I0_R8_smoke", uni20::format_scalar((*a.integral.mean)[0], {.precision = 16}));
  RecordProperty("I1_R8_smoke", uni20::format_scalar((*a.integral.mean)[1], {.precision = 16}));
}

TEST(SixSpinonMomentsMpfr, FailureAndCancellationLeaveNoMomentOrFraction)
{
  R one(1, uni20::Precision::bits(384));
  auto p = policy();
  p.stages = {{128, 192, {{{2, 1}, {4, 2}}}}};
  Moments model(proposal(one), p, one);
  auto result = model.global({2, 2, 42});
  EXPECT_EQ(result.integral.status, bethe::MonteCarloStatus::numerical_failure);
  EXPECT_FALSE(result.integral.mean);
  EXPECT_FALSE(result.fraction);
  EXPECT_GT(result.integral.metadata.density_attempts, 0u);
  std::stop_source stop;
  stop.request_stop();
  auto cancelled = model.global({2, 2, 42}, {.stop = stop.get_token()});
  EXPECT_EQ(cancelled.integral.status, bethe::MonteCarloStatus::cancelled);
  EXPECT_FALSE(cancelled.integral.mean);
  EXPECT_FALSE(cancelled.fraction);
  EXPECT_EQ(cancelled.integral.metadata.density_attempts, 0u);
  EXPECT_THROW((Moments(proposal(one), policy(), -one)), std::invalid_argument);
}
#endif
