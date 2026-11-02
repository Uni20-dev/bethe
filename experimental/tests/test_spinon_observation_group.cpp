// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/experimental/spinon_observation_group.hpp>
#include <bethe/xxx_six_spinon_structure_factor.hpp>

namespace x = bethe::experimental::xxx;
namespace h = bethe::heisenberg;
namespace
{
template <uni20::Real R> auto proposal(R const& one)
{
  using Kind = h::detail::SpinonProposalKind;
  return h::detail::SpinonJointProposal<R>{8 * one,
                                           one / 2,
                                           {{Kind::background, one / 2, one, one},
                                            {Kind::one_pair, one / 4, 4 * one, one / 2},
                                            {Kind::two_pairs, one / 4, 8 * one, one}}};
}
template <uni20::Real R> auto requests(R const& one, bool masks = true)
{
  R const zero = uni20::scalar_like(one, 0), pi = 4 * uni20::math::atan(one);
  std::vector<x::SpinonSectionRequest<R>> rows;
  for (int s = 0; s < 3; ++s)
  {
    rows.push_back({"whole-" + std::to_string(s), pi / 2, 2 * one, zero, 8 * one, s});
    if (masks)
    {
      rows.push_back({"compact-" + std::to_string(s), pi / 2, 2 * one, zero, 4 * one, s});
      rows.push_back({"remainder-" + std::to_string(s), pi / 2, 2 * one, 4 * one, 8 * one, s});
    }
  }
  return rows;
}
template <uni20::Real R> auto unit_density(std::size_t& evaluated)
{
  return [&](auto, auto, std::span<R const>, R const& log_weight, std::stop_token) {
    ++evaluated;
    bethe::MonteCarloObservation<R, h::SpinonSampleMetadata> result;
    R const weight = uni20::math::exp(log_weight);
    result.values = std::vector<R>{weight, weight / 16, weight / 32};
    result.metadata.supported_configurations = 1;
    return result;
  };
}

template <uni20::Real R> void check_geometry(R const& one)
{
  R const zero = uni20::scalar_like(one, 0), pi = 4 * uni20::math::atan(one);
  R const tolerance = 32768 * uni20::numeric_limits<R>::epsilon(one);
  x::SixSpinonObservationGroup<R> group(proposal(one), one, requests(one));
  bethe::MonteCarloStream<5> stream(42, 0, 24);
  bethe::ReplicateStatistics<R> statistics(24, 1, 9, one);
  std::size_t supported = 0, zeros = 0;
  for (std::size_t i = 0; i < 24; ++i)
  {
    auto identity = stream.next();
    auto prepared = group.prepare(std::to_string(i), identity);
    ASSERT_TRUE(prepared.draw());
    EXPECT_EQ(prepared.identity().words(), identity.words());
    ASSERT_EQ(prepared.sections().size(), 3u); // Three masks share each solve.
    ASSERT_EQ(prepared.rows().size(), 9u);
    for (std::size_t s = 0; s < prepared.sections().size(); ++s)
    {
      auto const& geometry = prepared.sections()[s].geometry;
      ASSERT_TRUE(geometry) << prepared.sections()[s].reason;
      ASSERT_TRUE(prepared.sections()[s].reason.empty());
      if (!geometry->inside)
      {
        ++zeros;
        continue;
      }
      ++supported;
      ASSERT_EQ(geometry->beta[0].size(), 6u);
      ASSERT_EQ(geometry->beta[1].size(), 6u);
      for (std::size_t j = 0; j < 4; ++j)
      {
        EXPECT_EQ(geometry->beta[0][j], prepared.draw()->beta[j]);
        EXPECT_EQ(geometry->beta[1][j], prepared.draw()->beta[j]);
      }
      EXPECT_EQ(geometry->beta[0][4], geometry->beta[1][5]);
      EXPECT_EQ(geometry->beta[0][5], geometry->beta[1][4]);
      for (auto const& root : geometry->beta)
      {
        R energy = zero, momentum = zero;
        for (auto const& beta : root)
        {
          energy += pi / (2 * uni20::math::cosh(beta));
          momentum += 2 * uni20::math::atan(uni20::math::exp(beta));
        }
        EXPECT_REAL_NEAR(energy, 2 * one, tolerance);
        EXPECT_REAL_NEAR(momentum, pi / 2 + 2 * s * pi, 8 * tolerance);
      }
    }
    std::size_t evaluated = 0;
    auto value = group.observe(prepared, unit_density<R>(evaluated));
    ASSERT_TRUE(value.values) << value.reason;
    EXPECT_EQ(evaluated, value.metadata.supported_configurations);
    EXPECT_LE(evaluated, 6u); // No repeat density call for the three masks.
    for (std::size_t s = 0; s < 3; ++s)
      for (std::size_t field = 0; field < 3; ++field)
      {
        auto base = 9 * field + 3 * s;
        EXPECT_EQ((*value.values)[base], (*value.values)[base + 1] + (*value.values)[base + 2]);
      }
    statistics.accept(i, 1, {value.values->begin(), value.values->begin() + 9});
  }
  EXPECT_GT(supported, 0u);
  EXPECT_GT(zeros, 0u);
  auto result = statistics.summary();
  ASSERT_TRUE(result.mean && result.covariance);
  EXPECT_EQ(result.received, 24u); // Includes resolved zeros in the denominator.
  for (std::size_t s = 0; s < 3; ++s)
  {
    auto a = 3 * s, b = a + 1, c = a + 2;
    EXPECT_REAL_NEAR((*result.mean)[a], (*result.mean)[b] + (*result.mean)[c], tolerance);
    auto const& covariance = *result.covariance;
    R const reconstructed = covariance[9 * b + b] + covariance[9 * c + c] + 2 * covariance[9 * b + c];
    EXPECT_REAL_NEAR(covariance[9 * a + a], reconstructed, tolerance * (one + covariance[9 * a + a]));
  }
}

template <uni20::Real G, uni20::Real Low, uni20::Real High>
void physical_match(G const& one, h::SpinonRecoveryPolicy policy)
{
  h::ThermodynamicSixSpinonStructureFactor<G, Low, High> reference(proposal(one), policy, one);
  h::XXXSpinonDensityRecovery<Low, High> recovery(policy);
  x::SixSpinonObservationGroup<G> group(proposal(one), one, requests(one, false));
  auto const point = reference.at(2 * uni20::math::atan(one), 2 * one);
  bethe::MonteCarloStream<5> stream(42, 0, 24);
  std::size_t supported = 0;
  for (std::size_t i = 0; i < 12; ++i)
  {
    auto const identity = stream.next();
    auto const prepared = group.prepare(std::to_string(i), identity);
    auto const actual = group.observation(prepared, recovery);
    auto const expected = point.observation(identity);
    ASSERT_EQ(actual.values.has_value(), expected.values.has_value()) << actual.reason << " / " << expected.reason;
    EXPECT_EQ(actual.metadata.failure, expected.failure);
    EXPECT_EQ(actual.metadata.supported_configurations, expected.supported_roots);
    EXPECT_EQ(actual.metadata.density_attempts, expected.density_attempts);
    if (!actual.values) continue;
    supported += actual.metadata.supported_configurations;
    for (std::size_t f = 0; f < 3; ++f)
    {
      G sum = uni20::scalar_like(one, 0);
      for (std::size_t s = 0; s < 3; ++s)
        sum += (*actual.values)[3 * f + s];
      G const tol = 256 * uni20::numeric_limits<G>::epsilon(one) * (one + uni20::math::abs(sum));
      EXPECT_REAL_NEAR(sum, (*expected.values)[f], tol);
    }
  }
  EXPECT_GT(supported, 0u);
}
} // namespace

template <typename R> class SpinonObservationGroup : public ::testing::Test {};
TYPED_TEST_SUITE(SpinonObservationGroup, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(SpinonObservationGroup, LabelledRootsMasksAndJointCovariance) { check_geometry(TypeParam{1}); }

TEST(SpinonObservationGroup, FreshPhysicalObservationsMatchPublicPoint)
{
  int const bits = uni20::numeric_limits<long double>::digits();
  if (bits <= 53) GTEST_SKIP();
  physical_match<long double, double, long double>(1.L, {"1e-5", "0", {{53, bits, {{{32, 12}, {64, 20}}}}}, 60});
}

TEST(SpinonObservationGroup, FailedRootsCancellationAndOwnershipNeverBecomeZero)
{
  using Group = x::SixSpinonObservationGroup<double>;
  Group group(proposal(1.), 1., requests(1., false)), other(proposal(1.), 1., requests(1., false));
  bethe::MonteCarloStream<5> stream(42, 0, 24);
  bool checked = false;
  for (int i = 0; i < 24 && !checked; ++i)
  {
    auto prepared = group.prepare(std::to_string(i), stream.next());
    if (std::none_of(prepared.rows().begin(), prepared.rows().end(),
                     [](auto const& row) { return row.support == x::SpinonSectionSupport::selected; }))
      continue;
    checked = true;
    std::size_t evaluated = 0;
    EXPECT_THROW(other.observe(prepared, unit_density<double>(evaluated)), std::invalid_argument);
    auto failed = group.observe(prepared, [&](auto, auto, auto, auto, auto) {
      Group::Observation result;
      result.metadata.failure = h::SpinonObservationFailure::density;
      result.metadata.density_status = h::SpinonRecoveryStatus::exhausted;
      result.reason = "deliberate failed density";
      return result;
    });
    EXPECT_FALSE(failed.values);
    EXPECT_EQ(failed.metadata.failure, h::SpinonObservationFailure::density);
    bethe::ReplicateStatistics<double> statistics(2, 1, 3, 1.);
    statistics.accept(0, 1, {1, 2, 3});
    statistics.reject(1, failed.reason);
    EXPECT_FALSE(statistics.summary().mean);
    std::stop_source stop;
    auto cancelled = group.observe(
        prepared,
        [&](auto a, auto b, auto beta, auto weight, auto token) {
          stop.request_stop();
          return unit_density<double>(evaluated)(a, b, beta, weight, token);
        },
        stop.get_token());
    EXPECT_TRUE(cancelled.cancelled);
    EXPECT_FALSE(cancelled.values);
    auto copy = prepared;
    auto moved = std::move(copy);
    EXPECT_THROW(copy.rows(), std::logic_error);
    EXPECT_EQ(moved.id(), prepared.id());
    auto const& row = *std::find_if(prepared.rows().begin(), prepared.rows().end(),
                                    [](auto const& r) { return r.support == x::SpinonSectionSupport::selected; });
    double depth = 0;
    for (auto b : prepared.sections()[row.section].geometry->beta[0])
      depth = std::max(depth, std::abs(b));
    auto r = requests(1., false)[row.section];
    r.outer_cutoff = depth;
    Group ambiguous(proposal(1.), 1., {r});
    auto unresolved = ambiguous.prepare("boundary", prepared.identity());
    EXPECT_EQ(unresolved.rows()[0].support, x::SpinonSectionSupport::geometry_failure);
    EXPECT_FALSE(ambiguous.observe(unresolved, unit_density<double>(evaluated)).values);
  }
  EXPECT_TRUE(checked);
  auto r = requests(1., false);
  r[1].id = r[0].id;
  EXPECT_THROW(Group(proposal(1.), 1., r), std::invalid_argument);
  r = requests(1., false);
  r[0].outer_cutoff = 9;
  EXPECT_THROW(Group(proposal(1.), 1., r), std::invalid_argument);
  r = requests(1., false);
  r[0].momentum = 0;
  EXPECT_THROW(Group(proposal(1.), 1., r), std::invalid_argument);
}

#if UNI20_HAS_MPREAL
TEST(SpinonObservationGroupMpfr, GeometryAndPhysicalRecoveryKeepExplicitPrecision)
{
  using R = uni20::mpreal;
  auto const ambient = mpfr_get_default_prec();
  R const one(1, uni20::Precision::bits(160));
  check_geometry(one);
  physical_match<R, R, R>(one, {"1e-5", "0", {{80, 128, {{{32, 12}, {64, 20}}}}}, 60});
  auto r = requests(one);
  r[0].outer_cutoff = R(8, uni20::Precision::bits(80));
  EXPECT_THROW((x::SixSpinonObservationGroup<R>(proposal(one), one, r)), std::invalid_argument);
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
#endif
