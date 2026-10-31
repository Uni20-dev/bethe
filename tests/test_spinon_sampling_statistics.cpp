// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/sampling.hpp>

namespace
{
namespace sampling = bethe;
namespace math = uni20::math;
template <class R> class SpinonSampling : public ::testing::Test {};
TYPED_TEST_SUITE(SpinonSampling, test_support::RealTypes, test_support::PrecisionNames);

template <class R> void endpoints(R const& one)
{
  auto const bits = unsigned(uni20::numeric_limits<R>::digits(one) - 1);
  std::size_t const count = (std::size_t(bits) + 63) / 64;
  unsigned const first = (bits - 1) % 64 + 1;
  sampling::UniformSample<2>::Words words;
  words[0].resize(count, 0);
  words[1].resize(count, std::numeric_limits<std::uint64_t>::max());
  if (first < 64) words[1][0] = (std::uint64_t{1} << first) - 1;
  sampling::UniformSample<2> const inputs(bits, words);
  auto const u = inputs.at_precision(one);
  R const epsilon = uni20::numeric_limits<R>::epsilon(one);
  EXPECT_EQ(u[0], epsilon / R{2});
  EXPECT_EQ(u[1], one - epsilon / R{2});
  EXPECT_GT(u[0], R{0});
  EXPECT_LT(u[1], one);
  EXPECT_EQ(inputs.words(), words);
}

template <class R> void covariance(R const& one)
{
  sampling::ReplicateStatistics<R> forward(4, 100, 2, one), reverse(4, 100, 2, one);
  for (int i = 0; i < 4; ++i)
  {
    R const a = R(1 + 2 * i) * one, b = R(7 - 2 * i) * one;
    forward.accept(i, 100, {a, R{2} * a});
    reverse.accept(3 - i, 100, {b, R{2} * b});
  }
  auto const a = forward.summary(), b = reverse.summary();
  ASSERT_EQ(a.status, sampling::ReplicateStatus::complete);
  ASSERT_TRUE(a.mean);
  ASSERT_TRUE(a.covariance);
  EXPECT_EQ(a.mean, b.mean);
  EXPECT_EQ(a.covariance, b.covariance);
  EXPECT_EQ((*a.mean)[0], R{4});
  EXPECT_EQ((*a.mean)[1], R{8});
  R const eps = uni20::numeric_limits<R>::epsilon(one);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      EXPECT_REAL_NEAR((*a.covariance)[2 * i + j], R(5 * (i + 1) * (j + 1)) * one / R{3}, R{16} * eps);
  // Fully correlated frequencies: Var(2X-Y) must cancel, not add two SEs.
  R const contrast = R{4} * (*a.covariance)[0] - R{4} * (*a.covariance)[1] + (*a.covariance)[3];
  EXPECT_REAL_NEAR(contrast, one - one, R{32} * eps);

  sampling::ReplicateStatistics<R> offset(3, 1, 1, one);
  R const anchor = math::ldexp(one, int(uni20::numeric_limits<R>::digits(one)) - 4);
  offset.accept(0, 1, {anchor});
  offset.accept(1, 1, {anchor + R{2}});
  offset.accept(2, 1, {anchor + R{4}});
  auto const shifted = offset.summary();
  EXPECT_EQ((*shifted.mean)[0], anchor + R{2});
  EXPECT_REAL_NEAR((*shifted.covariance)[0], R{4} * one / R{3}, R{8} * eps);
}

template <class R> void cancelling_mean(R const& one)
{
  R const zero = one - one;
  R const large = math::ldexp(one, int(uni20::numeric_limits<R>::digits(one)) + 2);
  sampling::ReplicateStatistics<R> stats(4, 1, 1, one);
  std::array<R, 4> rows{zero, one, large, -large};
  for (std::size_t i = 0; i < rows.size(); ++i)
    stats.accept(i, 1, {rows[i]});
  auto const result = stats.summary();
  ASSERT_TRUE(result.mean);
  // The small positive contribution precedes the large cancelling pair.
  // Kahan reduction loses it; the exact mean is representable at every precision.
  EXPECT_EQ((*result.mean)[0], one / R{4});
}

template <class R> void cancelling_covariance(R const& one)
{
  R const zero = one - one;
  R const large = math::ldexp(one, int(uni20::numeric_limits<R>::digits(one)) + 2);
  std::array<std::array<R, 2>, 7> rows{
      {{zero, zero}, {large, large}, {-large, -large}, {one, one}, {-one, -one}, {large, -large}, {-large, large}}};
  sampling::ReplicateStatistics<R> stats(rows.size(), 1, 2, one);
  for (std::size_t i = 0; i < rows.size(); ++i)
    stats.accept(i, 1, {rows[i][0], rows[i][1]});
  auto const result = stats.summary();
  ASSERT_TRUE(result.mean && result.covariance);
  EXPECT_EQ((*result.mean)[0], zero);
  EXPECT_EQ((*result.mean)[1], zero);
  // Means cancel exactly even with Kahan. Cross-products sum to 2, despite
  // their large positive/negative terms: Cov(mean_x,mean_y)=2/(7*6).
  EXPECT_EQ((*result.covariance)[1], one / R{21});
  EXPECT_EQ((*result.covariance)[2], one / R{21});
}

template <class R> void nonzero_anchor(R const& one)
{
  R const large = math::ldexp(one, int(uni20::numeric_limits<R>::digits(one)) + 2);
  std::array<R, 3> rows{-large, one, large};
  do
  {
    sampling::ReplicateStatistics<R> stats(3, 1, 1, one);
    for (std::size_t i = 0; i < rows.size(); ++i)
      stats.accept(i, 1, {rows[i]});
    auto const result = stats.summary();
    ASSERT_TRUE(result.mean);
    EXPECT_EQ((*result.mean)[0], one / R{3});
  }
  while (std::next_permutation(rows.begin(), rows.end()));
}

TYPED_TEST(SpinonSampling, NonzeroAnchorDoesNotEraseSignedMean) { nonzero_anchor(TypeParam{1}); }
TYPED_TEST(SpinonSampling, SignedCancellationRetainsSmallMean) { cancelling_mean(TypeParam{1}); }
TYPED_TEST(SpinonSampling, SignedCancellationRetainsCrossCovariance) { cancelling_covariance(TypeParam{1}); }
TYPED_TEST(SpinonSampling, ExactOpenGridEndpoints) { endpoints(TypeParam{1}); }
TYPED_TEST(SpinonSampling, OrderedCovarianceOfTheMean) { covariance(TypeParam{1}); }
TYPED_TEST(SpinonSampling, FailedOrMissingReplicatesCannotProduceSuccessfulMean)
{
  using R = TypeParam;
  sampling::ReplicateStatistics<R> stats(3, 100, 1, R{1});
  stats.accept(0, 100, {R{2}});
  auto partial = stats.summary();
  EXPECT_EQ(partial.status, sampling::ReplicateStatus::incomplete);
  EXPECT_EQ(partial.received, 1u);
  EXPECT_FALSE(partial.mean);
  EXPECT_FALSE(partial.covariance);
  EXPECT_THROW(stats.accept(1, 99, {R{1}}), std::invalid_argument);
  EXPECT_THROW(stats.accept(1, 100, {R{1}, R{2}}), std::invalid_argument);
  stats.reject(1, "unrecovered sample");
  stats.accept(2, 100, {R{3}});
  auto const failed = stats.summary();
  EXPECT_EQ(failed.status, sampling::ReplicateStatus::failed);
  EXPECT_EQ(failed.received, 3u);
  EXPECT_EQ(failed.failed, 1u);
  EXPECT_FALSE(failed.mean);
  EXPECT_FALSE(failed.covariance);
  ASSERT_EQ(failed.failures.size(), 1u);
  EXPECT_NE(failed.failures[0].find("replicate 1"), std::string::npos);
  EXPECT_THROW(stats.accept(0, 100, {R{5}}), std::invalid_argument);
  EXPECT_THROW(stats.reject(1, "again"), std::invalid_argument);
  EXPECT_THROW(stats.reject(3, "outside"), std::out_of_range);
  sampling::ReplicateStatistics<R> single(1, 100, 1, R{1});
  single.accept(0, 100, {R{0}}); // genuine support zeros count
  auto const result = single.summary();
  EXPECT_EQ(result.status, sampling::ReplicateStatus::complete);
  ASSERT_TRUE(result.mean);
  EXPECT_EQ((*result.mean)[0], R{0});
  EXPECT_FALSE(result.covariance); // one replicate cannot estimate an SE
}

TEST(SpinonSamplingStream, ReplayReplicateIdentityAndResolution)
{
  sampling::MonteCarloStream<5> a(1234, 8, 52), replay(1234, 8, 52), other(1234, 9, 52);
  for (int i = 0; i < 20; ++i)
  {
    auto const x = a.next(), y = replay.next(), z = other.next();
    EXPECT_EQ(x.words(), y.words());
    EXPECT_NE(x.words(), z.words());
    for (auto const u : x.at_precision(1.0))
    {
      EXPECT_GT(u, 0);
      EXPECT_LT(u, 1);
    }
  }
  sampling::MonteCarloStream<5> wide(1234, 8, 64);
  auto const x = wide.next();
  EXPECT_THROW(x.at_precision(1.0), std::invalid_argument);
  EXPECT_THROW((sampling::MonteCarloStream<1>(0, 0, 0)), std::invalid_argument);
  EXPECT_THROW((sampling::UniformSample<1>(8, {{{256}}})), std::invalid_argument);
  EXPECT_THROW((sampling::UniformSample<1>(65, {{{1}}})), std::invalid_argument);
}

TEST(SpinonSamplingStatistics, UnderflowDoesNotReportZeroSamplingError)
{
  sampling::ReplicateStatistics<double> stats(2, 1, 1, 1.0);
  stats.accept(0, 1, {0});
  stats.accept(1, 1, {2 * std::numeric_limits<double>::min()});
  EXPECT_THROW(stats.summary(), std::runtime_error);
}

TEST(SpinonSamplingStatistics, CommonOffsetSurvivesRawSumOverflow)
{
  double const common = std::numeric_limits<double>::max() / 2;
  sampling::ReplicateStatistics<double> stats(4, 1, 1, 1.);
  for (std::size_t i = 0; i < 4; ++i)
    stats.accept(i, 1, {common});
  auto const result = stats.summary();
  ASSERT_TRUE(result.mean && result.covariance);
  EXPECT_EQ((*result.mean)[0], common);
  EXPECT_EQ((*result.covariance)[0], 0.);
}

#if UNI20_HAS_MPREAL
TEST(SpinonSamplingMpfr, NonzeroAnchorDoesNotEraseSignedMean)
{
  for (unsigned bits : {160u, 333u})
    nonzero_anchor(uni20::mpreal(1, uni20::Precision::bits(bits)));
}
TEST(SpinonSamplingMpfr, SignedCancellationRetainsSmallMean)
{
  for (unsigned bits : {160u, 333u})
    cancelling_mean(uni20::mpreal(1, uni20::Precision::bits(bits)));
}
TEST(SpinonSamplingMpfr, SignedCancellationRetainsCrossCovariance)
{
  for (unsigned bits : {160u, 333u})
    cancelling_covariance(uni20::mpreal(1, uni20::Precision::bits(bits)));
}
TEST(SpinonSamplingMpfr, IdenticalRationalInputAtBothPrecisions)
{
  using R = uni20::mpreal;
  auto const p = uni20::Precision::bits(256), high = uni20::Precision::bits(2048);
  R const one(1, p), higher(1, high);
  endpoints(one);
  endpoints(higher);
  covariance(one);
  sampling::MonteCarloStream<5> stream(1234, 91, 255);
  auto const source = stream.next();
  auto const low = source.at_precision(one), guard = source.at_precision(higher);
  for (int i = 0; i < 5; ++i)
  {
    EXPECT_EQ(R(low[i], high), guard[i]);
    EXPECT_EQ(low[i].precision(), p);
    EXPECT_EQ(guard[i].precision(), high);
  }
  EXPECT_THROW(source.at_precision(R(1, uni20::Precision::bits(128))), std::invalid_argument);
  EXPECT_THROW(source.at_precision(R{1}), std::invalid_argument);
  sampling::ReplicateStatistics<R> stats(2, 10, 1, one);
  EXPECT_THROW(stats.accept(0, 10, {higher}), std::invalid_argument);
  stats.accept(0, 10, {R(2, p)});
  stats.accept(1, 10, {R(4, p)});
  auto const result = stats.summary();
  EXPECT_EQ((*result.mean)[0].precision(), p);
  EXPECT_EQ((*result.covariance)[0].precision(), p);
  EXPECT_EQ((*result.mean)[0], 3);
  EXPECT_EQ((*result.covariance)[0], 1);
}
#endif
} // namespace
