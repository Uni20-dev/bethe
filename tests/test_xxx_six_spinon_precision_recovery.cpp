// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxx_six_spinon_precision_recovery.hpp>

namespace
{
namespace model = bethe::heisenberg;
using Failure = model::SpinonObservationFailure;
using Status = model::SpinonRecoveryStatus;

struct Calls
{
    std::vector<int> stages;
    bethe::UniformSample<5>::Words expected;
};

template <class R> struct MockPoint
{
    model::SpinonSpectrumObservation<R> result;
    std::shared_ptr<Calls> calls;
    int stage;
    bool throws = false;
    std::stop_source* cancel_after = nullptr;
    R exemplar = R{1};
    std::optional<model::SixSpinonSpectrumStatus> boundary;
    auto boundary_status() const { return boundary; }
    auto observation(bethe::UniformSample<5> const& identity, std::stop_token) const
    {
      if (identity.bits() != 52 || identity.words() != calls->expected)
        throw std::logic_error("precision recovery changed the sample identity");
      // Reconstruct each stage's coordinates from the identity, not another
      // stage's rounded rapidities. The typed points need no common scalar.
      auto coordinates = identity.at_precision(exemplar);
      for (auto const& x : coordinates)
        if (!(x > R{0} && x < R{1})) throw std::logic_error("uniform left open support");
      calls->stages.push_back(stage);
      if (throws) throw std::logic_error("invalid model input");
      if (cancel_after) cancel_after->request_stop();
      return result;
    }
};

template <class R>
auto failure(Failure kind, std::size_t roots, std::size_t attempts, std::optional<Status> density = {})
{
  model::SpinonSpectrumObservation<R> r;
  r.failure = kind;
  r.supported_roots = roots;
  r.density_attempts = attempts;
  r.density_status = density;
  r.reason = "controlled failure";
  return r;
}

template <class R> auto success(R value, std::size_t roots = 2, std::size_t attempts = 4)
{
  model::SpinonSpectrumObservation<R> r;
  auto const zero = uni20::scalar_like(value, 0);
  r.values = std::array<R, 3>{value, zero, zero};
  r.supported_roots = roots;
  r.density_attempts = attempts;
  return r;
}

class SixSpinonPrecisionRecovery : public testing::Test {
  protected:
    bethe::UniformSample<5> identity = bethe::MonteCarloStream<5>(123, 8, 52).next();
    std::shared_ptr<Calls> calls = std::make_shared<Calls>();
    MockPoint<double> a{success(1.25), calls, 0};
    MockPoint<long double> b{success(2.5L), calls, 1}, c{success(3.75L), calls, 2};
    void SetUp() override { calls->expected = identity.words(); }
};

TEST_F(SixSpinonPrecisionRecovery, AcceptsFirstStageAndPreservesZeroSupport)
{
  auto run = model::make_six_spinon_precision_recovery(1.L, a, b, c);
  auto first = run.observation(identity);
  EXPECT_EQ(calls->stages, (std::vector<int>{0}));
  ASSERT_TRUE(first.observation.values);
  EXPECT_EQ((*first.observation.values)[0], 1.25L);
  EXPECT_EQ(first.attempted_stages, 1u);
  calls->stages.clear();
  a.result = success(0., 0, 0);
  auto zero = model::make_six_spinon_precision_recovery(1.L, a, b, c).observation(identity);
  ASSERT_TRUE(zero.observation.values);
  EXPECT_EQ((*zero.observation.values)[0], 0.L);
  EXPECT_EQ(zero.attempted_stages, 1u);
}

TEST_F(SixSpinonPrecisionRecovery, ReconstructsIdenticalDrawAndSeparatesFinalFromInclusiveWork)
{
  a.result = failure<double>(Failure::geometry, 1, 0);
  b.result = failure<long double>(Failure::density, 2, 1, Status::exhausted);
  auto all = model::make_six_spinon_precision_recovery(1.L, a, b, c).observation(identity);
  EXPECT_EQ(calls->stages, (std::vector<int>{0, 1, 2}));
  ASSERT_TRUE(all.observation.values);
  EXPECT_EQ((*all.observation.values)[0], 3.75L);
  EXPECT_EQ(all.attempted_stages, 3u);
  EXPECT_EQ(all.evaluated_roots, 5u);
  EXPECT_EQ(all.density_attempts, 5u);
  EXPECT_EQ(all.observation.supported_roots, 2u);
  EXPECT_EQ(all.observation.density_attempts, 4u);
  EXPECT_EQ(all.attempts[0].failure, Failure::geometry);
  EXPECT_EQ(all.attempts[1].density_status, Status::exhausted);
}

TEST_F(SixSpinonPrecisionRecovery, RetriesBoundaryDisagreementsBeforeCallingTheNativePoint)
{
  for (auto boundary : {model::SixSpinonSpectrumStatus::outside_continuum, model::SixSpinonSpectrumStatus::upper_edge,
                        model::SixSpinonSpectrumStatus::lower_threshold})
  {
    SCOPED_TRACE(int(boundary));
    calls->stages.clear();
    a.boundary = boundary;
    a.result = success(0., 0, 0);
    a.throws = boundary == model::SixSpinonSpectrumStatus::lower_threshold;
    auto result = model::make_six_spinon_precision_recovery(1.L, a, b).observation(identity);
    EXPECT_EQ(calls->stages, (std::vector<int>{1}));
    EXPECT_EQ(result.attempted_stages, 2u);
    EXPECT_EQ(result.attempts[0].failure, Failure::geometry);
    ASSERT_TRUE(result.observation.values);
    EXPECT_EQ((*result.observation.values)[0], 2.5L);
    EXPECT_EQ(result.density_attempts, 4u);
  }
  // Matching physical boundary zeros remain successful zeros.
  calls->stages.clear();
  a.boundary = b.boundary = model::SixSpinonSpectrumStatus::upper_edge;
  a.throws = false;
  a.result = success(0., 0, 0);
  auto zero = model::make_six_spinon_precision_recovery(1.L, a, b).observation(identity);
  ASSERT_TRUE(zero.observation.values);
  EXPECT_EQ((*zero.observation.values)[0], 0.L);
  EXPECT_EQ(calls->stages, (std::vector<int>{0}));
}

TEST_F(SixSpinonPrecisionRecovery, DoesNotRetryTerminalStatuses)
{
  for (auto status : {Status::cancelled, Status::time_limit, Status::invalid_input, Status::unsupported_limit,
                      Status::inconsistent_controls})
  {
    SCOPED_TRACE(int(status));
    calls->stages.clear();
    a.result = failure<double>(Failure::density, 1, 1, status);
    auto result = model::make_six_spinon_precision_recovery(1.L, a, b, c).observation(identity);
    EXPECT_EQ(calls->stages, (std::vector<int>{0}));
    EXPECT_FALSE(result.observation.values);
    EXPECT_EQ(result.observation.density_status, status);
    EXPECT_EQ(result.density_attempts, 1u);
  }
  calls->stages.clear();
  a.result = failure<double>(Failure::cancelled, 0, 0);
  auto result = model::make_six_spinon_precision_recovery(1.L, a, b, c).observation(identity);
  EXPECT_FALSE(result.observation.values);
  EXPECT_EQ(calls->stages, (std::vector<int>{0}));
}

TEST_F(SixSpinonPrecisionRecovery, RetainsAllFailuresWhenExhausted)
{
  a.result = failure<double>(Failure::exponent_range, 1, 2);
  b.result = failure<long double>(Failure::geometry, 1, 0);
  c.result = failure<long double>(Failure::density, 2, 1, Status::exhausted);
  auto failed = model::make_six_spinon_precision_recovery(1.L, a, b, c).observation(identity);
  EXPECT_FALSE(failed.observation.values);
  EXPECT_EQ(failed.attempted_stages, 3u);
  EXPECT_EQ(failed.observation.failure, Failure::density);
  EXPECT_EQ(failed.evaluated_roots, 4u);
  EXPECT_EQ(failed.density_attempts, 3u);
}

TEST_F(SixSpinonPrecisionRecovery, CancelsBeforeAndBetweenStages)
{
  std::stop_source cancelled;
  cancelled.request_stop();
  auto stopped = model::make_six_spinon_precision_recovery(1.L, a, b, c).observation(identity, cancelled.get_token());
  EXPECT_TRUE(calls->stages.empty());
  EXPECT_EQ(stopped.attempted_stages, 0u);
  EXPECT_EQ(stopped.observation.failure, Failure::cancelled);
  EXPECT_FALSE(stopped.observation.values);
  std::stop_source between;
  a.result = failure<double>(Failure::exponent_range, 1, 2);
  a.cancel_after = &between;
  auto result = model::make_six_spinon_precision_recovery(1.L, a, b, c).observation(identity, between.get_token());
  EXPECT_EQ(calls->stages, (std::vector<int>{0}));
  EXPECT_EQ(result.observation.failure, Failure::cancelled);
  EXPECT_EQ(result.attempted_stages, 1u);
  EXPECT_EQ(result.evaluated_roots, 1u);
  EXPECT_EQ(result.density_attempts, 2u);
}

TEST_F(SixSpinonPrecisionRecovery, PropagatesExceptions)
{
  a.throws = true;
  auto run = model::make_six_spinon_precision_recovery(1.L, a, b, c);
  EXPECT_THROW((void)run.observation(identity), std::logic_error);
  EXPECT_EQ(calls->stages, (std::vector<int>{0}));
}

TEST_F(SixSpinonPrecisionRecovery, RejectsInsufficientOutputPrecision)
{
  if constexpr (std::numeric_limits<long double>::digits <= std::numeric_limits<double>::digits)
    GTEST_SKIP() << "long double has no wider significand on this platform";
  auto run = model::make_six_spinon_precision_recovery(1., b);
  EXPECT_THROW((void)run.observation(identity), std::invalid_argument);
}

TEST_F(SixSpinonPrecisionRecovery, RejectsWorkOverflowAndNonfiniteValues)
{
  a.result = failure<double>(Failure::geometry, std::numeric_limits<std::size_t>::max(), 0);
  b.result = success(2.5L, 1, 1);
  auto run = model::make_six_spinon_precision_recovery(1.L, a, b);
  EXPECT_THROW((void)run.observation(identity), std::length_error);
  b.result = success(std::numeric_limits<long double>::infinity());
  auto nonfinite = model::make_six_spinon_precision_recovery(1.L, b);
  EXPECT_THROW((void)nonfinite.observation(identity), std::range_error);
}

#if UNI20_HAS_MPREAL
TEST_F(SixSpinonPrecisionRecovery, DetectsActualRoundedLowerBoundary)
{
  if constexpr (std::numeric_limits<long double>::digits <= std::numeric_limits<double>::digits)
    GTEST_SKIP() << "requires a distinct native comparison format";
  using Mp = uni20::mpreal;
  using Kind = model::detail::SpinonProposalKind;
  auto make = []<class R, class Low>(R const& one, Low const&) {
    model::detail::SpinonJointProposal<R> proposal(8 * one, one / 2,
                                                   {{Kind::background, one / 2, one, one},
                                                    {Kind::one_pair, one / 4, 4 * one, one / 2},
                                                    {Kind::two_pairs, one / 4, 8 * one, one}});
    int const high = uni20::numeric_limits<R>::digits(one);
    int const low = uni20::has_runtime_precision_v<R> ? high - 64 : std::numeric_limits<Low>::digits;
    model::SpinonRecoveryPolicy policy{"1e-5", "0", {{low, high, {{{32, 12}, {64, 20}}}}}, 60};
    return model::ThermodynamicSixSpinonStructureFactor<R, Low, R>(std::move(proposal), policy, one);
  };
  auto native = make(1.L, 1.);
  Mp const one(1, uni20::Precision::bits(256));
  auto reference = make(one, one);
  auto const pi = Mp{4} * uni20::math::atan(one);
  auto const native_pi = 4.L * uni20::math::atan(1.L);
  bool checked = false;
  for (int j = 1; j <= 32 && !checked; ++j)
  {
    long double const q = static_cast<long double>(j) / 32;
    auto const native_lower = native_pi / 2 * uni20::math::sin(q);
    Mp const qref(q, one.precision()), rounded_lower(native_lower, one.precision());
    auto const true_lower = pi / Mp{2} * uni20::math::sin(qref);
    if (rounded_lower <= true_lower) continue;
    // Interior in the reference format, but the native conversion lands
    // exactly on its lower threshold. No density evaluation is required.
    Mp const omega = rounded_lower + uni20::math::ldexp(one, -240);
    auto native_point = native.at(q, static_cast<long double>(omega));
    auto reference_point = reference.at(qref, omega);
    ASSERT_EQ(native_point.boundary_status(), model::SixSpinonSpectrumStatus::lower_threshold);
    ASSERT_FALSE(reference_point.boundary_status());
    auto final = MockPoint<Mp>{success(Mp(7, one.precision())), calls, 1};
    final.exemplar = one;
    final.boundary = reference_point.boundary_status();
    auto result = model::make_six_spinon_precision_recovery(one, native_point, final).observation(identity);
    EXPECT_EQ(result.attempted_stages, 2u);
    EXPECT_EQ(result.attempts[0].failure, Failure::geometry);
    ASSERT_TRUE(result.observation.values);
    EXPECT_EQ((*result.observation.values)[0], Mp(7, one.precision()));
    checked = true;
  }
  EXPECT_TRUE(checked) << "no native upward-rounded lower boundary in the fixture grid";
}

TEST_F(SixSpinonPrecisionRecovery, RejectsInexactButAcceptsExactNativeSubnormal)
{
  using Mp = uni20::mpreal;
  auto const precision = uni20::Precision::bits(53);
  auto tiny = uni20::math::ldexp(Mp(3, precision), -1075);
  auto point = MockPoint<Mp>{success(tiny), calls, 0};
  point.exemplar = Mp(1, precision);
  auto inexact = model::make_six_spinon_precision_recovery(1., point);
  EXPECT_THROW((void)inexact.observation(identity), std::range_error);
  point.result = success(uni20::math::ldexp(Mp(1, precision), -1074));
  auto exact = model::make_six_spinon_precision_recovery(1., point).observation(identity);
  ASSERT_TRUE(exact.observation.values);
  EXPECT_EQ((*exact.observation.values)[0], std::numeric_limits<double>::denorm_min());
}
#endif
} // namespace

TEST_F(SixSpinonPrecisionRecovery, CallbackRetainsFinalAndPerStageAccounting)
{
  a.result = failure<double>(Failure::geometry, 1, 0);
  b.result = failure<long double>(Failure::density, 2, 1, Status::exhausted);
  auto callback = model::make_six_spinon_precision_recovery(1.L, a, b, c);
  auto sample = callback(identity);
  ASSERT_TRUE(sample.values);
  EXPECT_EQ((*sample.values)[0], 3.75L);
  EXPECT_EQ(sample.metadata.supported_configurations, 2u);
  EXPECT_EQ(sample.metadata.density_attempts, 4u);
  EXPECT_EQ(sample.metadata.failure, Failure::none);
  EXPECT_FALSE(sample.cancelled);
  auto const& stages = sample.metadata.stages;
  for (auto const& stage : stages)
    EXPECT_EQ(stage.attempts, 1u);
  EXPECT_EQ(stages[0].accepted, 0u);
  EXPECT_EQ(stages[1].accepted, 0u);
  EXPECT_EQ(stages[2].accepted, 1u);
  EXPECT_EQ(stages[0].retried, 1u);
  EXPECT_EQ(stages[1].retried, 1u);
  EXPECT_EQ(stages[2].retried, 0u);
  EXPECT_EQ(stages[0].evaluated_roots + stages[1].evaluated_roots + stages[2].evaluated_roots, 5u);
  EXPECT_EQ(stages[0].density_attempts + stages[1].density_attempts + stages[2].density_attempts, 5u);
  EXPECT_EQ(stages[0].failures.begin()->first,
            (model::SixSpinonPrecisionStageWork::FailureKey{Failure::geometry, {}, {}, "controlled failure"}));
  EXPECT_EQ(stages[1].failures.begin()->first, (model::SixSpinonPrecisionStageWork::FailureKey{
                                                   Failure::density, Status::exhausted, {}, "controlled failure"}));
  EXPECT_EQ(stages[0].failures.begin()->second.evaluated_roots, 1u);
  EXPECT_EQ(stages[1].failures.begin()->second.density_attempts, 1u);
  EXPECT_TRUE(stages[2].failures.empty());
  auto aggregate = sample.metadata;
  aggregate += sample.metadata;
  EXPECT_EQ(aggregate.supported_configurations, 4u);
  EXPECT_EQ(aggregate.density_attempts, 8u);
  EXPECT_EQ(aggregate.stages[0].failures.begin()->second.attempts, 2u);
  EXPECT_EQ(aggregate.stages[2].accepted, 2u);
}

TEST_F(SixSpinonPrecisionRecovery, CallbackCancellationDoesNotInventAnAttempt)
{
  std::stop_source stop;
  a.result = failure<double>(Failure::geometry, 1, 0);
  a.cancel_after = &stop;
  auto sample = model::make_six_spinon_precision_recovery(1.L, a, b)(identity, stop.get_token());
  EXPECT_TRUE(sample.cancelled);
  EXPECT_FALSE(sample.values);
  EXPECT_EQ(sample.metadata.failure, Failure::cancelled);
  EXPECT_EQ(sample.metadata.stages[0].attempts, 1u);
  EXPECT_EQ(sample.metadata.stages[0].retried, 0u);
  EXPECT_EQ(sample.metadata.stages[1].attempts, 0u);
  EXPECT_TRUE(sample.metadata.stages[1].failures.empty());
}

TEST_F(SixSpinonPrecisionRecovery, CallbackWorksWithModelAssembly)
{
  if (uni20::numeric_limits<long double>::digits(1.L) <= uni20::numeric_limits<double>::digits(1.)) GTEST_SKIP();
  using Kind = model::detail::SpinonProposalKind;
  model::detail::SpinonJointProposal<long double> proposal(8.L, .5L, {{Kind::background, 1.L, 1.L, 1.L}});
  model::SpinonRecoveryPolicy policy{
      "1e-5", "0", {{53, uni20::numeric_limits<long double>::digits(1.L), {{{32, 12}, {64, 20}}}}}, 60};
  model::ThermodynamicSixSpinonStructureFactor<long double, double, long double> spectrum(proposal, policy, 1.L);
  a.result = failure<double>(Failure::geometry, 1, 0);
  auto callback = model::make_six_spinon_precision_recovery(1.L, a, b);
  auto result = spectrum.at(1.5L, 2.L).evaluate_observations_async({1, 1, 123, 8, 52}, callback).get_wait();
  ASSERT_EQ(result.spectrum.status, model::SixSpinonSpectrumStatus::estimated);
  EXPECT_EQ(result.spectrum.value, 2.5L);
  EXPECT_FALSE(result.spectrum.sampling_error); // one replicate has no SE
  EXPECT_EQ(result.metadata.stages[0].retried, 1u);
  EXPECT_EQ(result.metadata.stages[1].accepted, 1u);
  ASSERT_TRUE(result.replicate_metadata[0]);
  EXPECT_EQ(result.replicate_metadata[0]->stages[0].evaluated_roots, 1u);
  EXPECT_EQ(result.spectrum.supported_roots, 2u);
}

TEST(SixSpinonPrecisionMetadata, OverflowCannotWrapReportedWork)
{
  model::SixSpinonPrecisionMetadata<1> all, extra;
  all.supported_configurations = std::numeric_limits<std::size_t>::max();
  extra.supported_configurations = 1;
  EXPECT_THROW(all += extra, std::length_error);
  all = {};
  extra = {};
  all.stages[0].attempts = std::numeric_limits<std::size_t>::max();
  extra.stages[0].attempts = 1;
  EXPECT_THROW(all += extra, std::length_error);
}
