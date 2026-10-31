// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxx_six_spinon_moments.hpp>
#include <uni20/async/tbb_scheduler.hpp>

TEST(SixSpinonFixedMoments, NativeZeroMomentumDoesNotStartDensityWork)
{
  using namespace bethe::heisenberg;
  if (std::numeric_limits<long double>::digits <= std::numeric_limits<double>::digits) GTEST_SKIP();
  detail::SpinonJointProposal<long double> proposal(8.L, .5L,
                                                    {{detail::SpinonProposalKind::background, 1.L, 1.L, 1.L}});
  SpinonRecoveryPolicy policy{
      "1e-5", "0", {{53, std::numeric_limits<long double>::digits, {{{32, 12}, {64, 20}}}}}, 60};
  ThermodynamicSixSpinonMoments<long double, double, long double> model(proposal, policy, 1.L);
  auto const pi = 4 * std::atan(1.L);
  for (auto k : {0.L, 2 * pi})
  {
    auto result = model.fixed_moments(k).evaluate({2, 4});
    EXPECT_EQ(result.status, SixSpinonMomentStatus::physical_zero);
    ASSERT_TRUE(result.integral.mean);
    EXPECT_EQ((*result.integral.mean)[0], 0.L);
    EXPECT_EQ(result.integral.metadata.density_attempts, 0u);
    EXPECT_FALSE(result.first_moment_fraction);
  }
  EXPECT_THROW(model.fixed_moments(-1.L), std::invalid_argument);
  EXPECT_THROW(model.fixed_moments(pi / 2).evaluate({0, 4}), std::invalid_argument);
}

#if UNI20_HAS_MPREAL
namespace
{
using namespace bethe::heisenberg;
using R = uni20::mpreal;
using Moments = ThermodynamicSixSpinonMoments<R, R>;
using Kind = detail::SpinonProposalKind;
auto proposal(R const& one)
{
  return detail::SpinonJointProposal<R>(8 * one, one / 2,
                                        {{Kind::background, one / 2, one, one},
                                         {Kind::one_pair, one / 4, 4 * one, one / 2},
                                         {Kind::two_pairs, one / 4, 8 * one, one}});
}
auto policy()
{
  return SpinonRecoveryPolicy{
      "1e-7", "0", {{128, 192, {{{32, 12}, {64, 20}}}}, {224, 320, {{{64, 20}, {128, 24}}}}}, 60};
}
} // namespace

TEST(SixSpinonFixedMomentsMpfr, IndependentPhysicalObservation)
{
  R const one(1, uni20::Precision::bits(384)), pi = 4 * uni20::math::atan(one);
  Moments model(proposal(one), policy(), one);
  bethe::UniformSample<6>::Words words;
  std::array<unsigned, 6> raw{2, 32, 80, 112, 180, 220};
  for (int i = 0; i < 6; ++i)
    words[i] = {raw[i]};
  auto result = model.fixed_moments(pi / 2).observation({8, words});
  ASSERT_TRUE(result.values) << result.reason;
  EXPECT_EQ(result.metadata.supported_configurations, 1u);
  EXPECT_GT(result.metadata.density_attempts, 0u);
  // Independent positive full-chart PDF and 85-dps Python W6, 128 terms/order24.
  std::array<char const*, 2> refs{
      "0.00000000767156467774964093309507201352870341152693541997552765266214319862735702515",
      "0.0000000372313292235237059410754393305114985867864419548104062719697283781970906354"};
  std::array<char const*, 2> tails{"4.29715733570045984090388156480214325150014481245391474082259888501292082374e-57",
                                   "2.08547911946529488074290956250998683989856803008784896995997947822634507801e-56"};
  for (int i = 0; i < 2; ++i)
  {
    auto const& v = *result.values;
    R const expected(refs[i], one.precision()), tail(tails[i], one.precision());
    EXPECT_LE(uni20::math::abs(v[i] - expected), v[i + 2] + v[i + 4] + tail);
    EXPECT_LT(v[i + 2] + v[i + 4], expected * R("1e-7", one.precision()));
  }
}

TEST(SixSpinonFixedMomentsMpfr, DirectMomentsCovarianceScalingAndOwnedAsync)
{
  R const one(1, uni20::Precision::bits(384)), pi = 4 * uni20::math::atan(one);
  Moments model(proposal(one), policy(), one);
  bethe::MonteCarloOptions options{2, 4, 42};
  auto serial = model.fixed_moments(pi / 2).evaluate(options, {.max_in_flight = 1});
  std::string failures;
  for (auto const& row : serial.integral.replicates)
    if (row && !row->means) failures += row->reason + "; ";
  ASSERT_EQ(serial.status, SixSpinonMomentStatus::estimated) << serial.reason << failures;
  ASSERT_TRUE(serial.integral.mean && serial.integral.covariance && serial.first_moment_fraction);
  EXPECT_EQ(serial.integral.completed_samples, 8u);
  EXPECT_GT(serial.integral.metadata.supported_configurations, 0u);
  EXPECT_GT((*serial.integral.covariance)[1], 0);
  R const exact = 2 * one / 3 * (uni20::math::log(2 * one) - one / 4), tol("1e-100", one.precision());
  EXPECT_REAL_NEAR(*serial.first_moment_denominator, exact, tol);
  EXPECT_REAL_NEAR(*serial.first_moment_fraction, (*serial.integral.mean)[1] / exact, tol);
  uni20::async::TbbScheduler scheduler(2);
  auto pending = Moments(proposal(one), policy(), one)
                     .fixed_moments(pi / 2)
                     .evaluate_async(options, {.scheduler = &scheduler, .max_in_flight = 2});
  auto parallel = pending.get_wait(scheduler);
  EXPECT_EQ(serial.integral.mean, parallel.integral.mean);
  EXPECT_EQ(serial.integral.covariance, parallel.integral.covariance);
  auto scaled = Moments(proposal(one), policy(), 3 * one)
                    .fixed_moments(pi / 2)
                    .evaluate(options, {.scheduler = &scheduler, .max_in_flight = 2});
  ASSERT_TRUE(scaled.integral.mean && scaled.integral.covariance && scaled.first_moment_fraction);
  EXPECT_EQ((*scaled.integral.mean)[0], (*serial.integral.mean)[0]);
  EXPECT_REAL_NEAR((*scaled.integral.mean)[1], 3 * (*serial.integral.mean)[1], tol);
  EXPECT_REAL_NEAR((*scaled.integral.covariance)[7], 9 * (*serial.integral.covariance)[7], tol);
  EXPECT_REAL_NEAR(*scaled.first_moment_fraction, *serial.first_moment_fraction, tol);
}

TEST(SixSpinonFixedMomentsMpfr, EndpointsMissingSupportAndFailuresStayDistinct)
{
  R const one(1, uni20::Precision::bits(384)), pi = 4 * uni20::math::atan(one);
  Moments model(proposal(one), policy(), one);
  for (R const k : {0 * one, 2 * pi})
  {
    auto result = model.fixed_moments(k).evaluate({2, 4, 42});
    EXPECT_EQ(result.status, SixSpinonMomentStatus::physical_zero);
    ASSERT_TRUE(result.integral.mean);
    EXPECT_EQ((*result.integral.mean)[0], 0);
    EXPECT_EQ((*result.integral.mean)[1], 0);
    EXPECT_EQ(result.integral.completed_samples, 0u);
    EXPECT_FALSE(result.first_moment_denominator);
    EXPECT_FALSE(result.first_moment_fraction);
  }
  auto endpoint = model.fixed_moments(pi).evaluate({1, 2, 42});
  EXPECT_TRUE(endpoint.unrestricted_static_divergence);
  ASSERT_TRUE(endpoint.first_moment_denominator);
  R const exact = 4 * one / 3 * (uni20::math::log(2 * one) - one / 4);
  EXPECT_REAL_NEAR(*endpoint.first_moment_denominator, exact, R("1e-100", one.precision()));
  auto p = policy();
  p.stages = {{128, 192, {{{2, 1}, {4, 2}}}}};
  auto failed = Moments(proposal(one), p, one).fixed_moments(pi / 2).evaluate({2, 4, 42});
  EXPECT_EQ(failed.status, SixSpinonMomentStatus::numerical_failure);
  EXPECT_FALSE(failed.integral.mean);
  EXPECT_FALSE(failed.first_moment_fraction);
  EXPECT_GT(failed.integral.metadata.density_attempts, 0u);
  std::stop_source stop;
  stop.request_stop();
  for (R const k : {0 * one, pi / 2})
  {
    auto cancelled = model.fixed_moments(k).evaluate({2, 4, 42}, {.stop = stop.get_token()});
    EXPECT_EQ(cancelled.status, SixSpinonMomentStatus::cancelled);
    EXPECT_FALSE(cancelled.integral.mean);
    EXPECT_FALSE(cancelled.first_moment_fraction);
    EXPECT_EQ(cancelled.integral.metadata.density_attempts, 0u);
  }
  detail::SpinonJointProposal<R> narrow(one / 20, one / 2, {{Kind::background, one, one, one}});
  auto empty = Moments(narrow, policy(), one).fixed_moments(pi / 2).evaluate({2, 4, 42});
  EXPECT_EQ(empty.status, SixSpinonMomentStatus::no_support_sampled);
  EXPECT_FALSE(empty.integral.mean);
  EXPECT_FALSE(empty.integral.covariance);
  EXPECT_FALSE(empty.first_moment_fraction);
  EXPECT_EQ(empty.integral.metadata.density_attempts, 0u);
  EXPECT_THROW(model.fixed_moments(-one), std::invalid_argument);
  EXPECT_THROW(model.fixed_moments(pi.at(uni20::Precision::bits(256))), std::invalid_argument);
}
#endif
