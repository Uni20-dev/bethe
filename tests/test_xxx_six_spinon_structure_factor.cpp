// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "xxx_six_spinon_spectrum_references.hpp"
#include <bethe/xxx_six_spinon_structure_factor.hpp>
#include <uni20/async/tbb_scheduler.hpp>

using namespace bethe::heisenberg;

namespace
{
using Kind = bethe::heisenberg::detail::SpinonProposalKind;
namespace sampling = bethe::heisenberg::detail;

template <class R> sampling::SpinonJointProposal<R> proposal(R const& one)
{
  return {8 * one,
          one / 2,
          {{Kind::background, one / 2, one, one},
           {Kind::one_pair, one / 4, 4 * one, one / 2},
           {Kind::two_pairs, one / 4, 8 * one, one}}};
}

SpinonRecoveryPolicy native_policy()
{
  return {"1e-5",
          "0",
          {{uni20::numeric_limits<double>::digits(1.),
            uni20::numeric_limits<long double>::digits(1.L),
            {{{32, 12}, {64, 20}}}}},
          60};
}
} // namespace

TEST(SixSpinonSpectrum, BoundaryLimitsDoNotStartSampling)
{
  if (uni20::numeric_limits<long double>::digits(1.L) <= uni20::numeric_limits<double>::digits(1.)) GTEST_SKIP();
  ThermodynamicSixSpinonStructureFactor<long double, double, long double> f(proposal(1.L), native_policy(), 1.L);
  auto const pi = 4 * std::atan(1.L);
  for (auto q : {0.L, 2 * pi})
  {
    auto result = f.at(q, 2.L).evaluate({2, 8});
    EXPECT_EQ(result.status, SixSpinonSpectrumStatus::outside_continuum);
    ASSERT_TRUE(result.value);
    EXPECT_EQ(*result.value, 0.L);
    EXPECT_EQ(result.density_attempts, 0u);
    EXPECT_TRUE(result.replicates.empty());
  }
  EXPECT_EQ(f.at(pi / 2, 0.L).evaluate({2, 8}).status, SixSpinonSpectrumStatus::outside_continuum);
  EXPECT_EQ(f.at(pi, 0.L).evaluate({2, 8}).status, SixSpinonSpectrumStatus::lower_threshold);
  auto upper = f.at(pi, 3 * pi).evaluate({2, 8});
  EXPECT_EQ(upper.status, SixSpinonSpectrumStatus::upper_edge);
  EXPECT_EQ(*upper.value, 0.L);
  auto above = f.at(pi / 2, 10.L).evaluate({2, 8});
  EXPECT_EQ(above.status, SixSpinonSpectrumStatus::outside_continuum);
  EXPECT_THROW(f.at(-1.L, 2.L), std::invalid_argument);
  EXPECT_THROW(f.at(pi / 2, 2.L).evaluate({0, 8}), std::invalid_argument);
  EXPECT_THROW(f.at(pi / 2, 2.L).evaluate({2, 8, 0, 0, 100}), std::invalid_argument);
}

namespace
{
struct ObservedMetadata : SpinonSampleMetadata
{
    std::size_t observations = 0;
    ObservedMetadata& operator+=(ObservedMetadata const& other)
    {
      SpinonSampleMetadata::operator+=(other);
      observations += other.observations;
      return *this;
    }
};
using ObservedSample = bethe::MonteCarloObservation<long double, ObservedMetadata>;

ObservedSample known_observation(bethe::UniformSample<5> const& identity, std::stop_token)
{
  auto const u = identity.at_precision(1.L);
  ObservedSample result;
  result.values = std::vector{u[0], u[0] / 8, u[0] / 16};
  result.metadata.supported_configurations = 2;
  result.metadata.density_attempts = 4;
  result.metadata.observations = 1;
  return result;
}
} // namespace

TEST(SixSpinonObservedSpectrum, ReferenceBoundariesBypassObserverAndKeepPhysicalStatuses)
{
  if (uni20::numeric_limits<long double>::digits(1.L) <= uni20::numeric_limits<double>::digits(1.)) GTEST_SKIP();
  ThermodynamicSixSpinonStructureFactor<long double, double, long double> f(proposal(1.L), native_policy(), 1.L);
  auto const pi = 4 * std::atan(1.L);
  auto never = [](auto const&, std::stop_token) -> ObservedSample { throw std::logic_error("boundary sampled"); };
  for (auto const& point : {f.at(0.L, 2.L), f.at(pi, 0.L), f.at(pi, 3 * pi)})
  {
    auto ordinary = point.evaluate({2, 8});
    auto observed = point.evaluate_observations_async({2, 8}, never).get_wait();
    EXPECT_EQ(observed.spectrum.status, ordinary.status);
    EXPECT_EQ(observed.spectrum.value, ordinary.value);
    EXPECT_EQ(observed.spectrum.sampling_error, ordinary.sampling_error);
    EXPECT_EQ(observed.spectrum.reason, ordinary.reason);
    EXPECT_TRUE(observed.spectrum.replicates.empty());
    EXPECT_TRUE(observed.replicate_metadata.empty());
    EXPECT_EQ(observed.metadata.observations, 0u);
    std::stop_source stop;
    stop.request_stop();
    auto cancelled = point.evaluate_observations_async({2, 8}, never, {.stop = stop.get_token()}).get_wait();
    EXPECT_EQ(cancelled.spectrum.status, SixSpinonSpectrumStatus::cancelled);
    EXPECT_FALSE(cancelled.spectrum.value);
  }
}

TEST(SixSpinonObservedSpectrum, AssemblesSamplingErrorsAndReplicateAccounting)
{
  if (uni20::numeric_limits<long double>::digits(1.L) <= uni20::numeric_limits<double>::digits(1.)) GTEST_SKIP();
  ThermodynamicSixSpinonStructureFactor<long double, double, long double> f(proposal(1.L), native_policy(), 1.L);
  auto const pi = 4 * std::atan(1.L);
  SpinonSamplingOptions options{3, 11, 42, 9, 52};
  auto expected = bethe::monte_carlo_async<5>(options, 3, 1.L, known_observation).get_wait();
  uni20::async::TbbScheduler scheduler(2);
  auto pending = f.at(pi / 2, 2.L)
                     .evaluate_observations_async(options, known_observation,
                                                  {.scheduler = &scheduler, .batch_size = 4, .max_in_flight = 3});
  auto observed = pending.get_wait(scheduler);
  auto const& r = observed.spectrum;
  ASSERT_EQ(r.status, SixSpinonSpectrumStatus::estimated);
  ASSERT_TRUE(r.value && r.sampling_error && r.contour_diagnostic && r.arithmetic_diagnostic);
  EXPECT_EQ(*r.value, (*expected.mean)[0]);
  EXPECT_EQ(*r.contour_diagnostic, (*expected.mean)[1]);
  EXPECT_EQ(*r.arithmetic_diagnostic, (*expected.mean)[2]);
  EXPECT_EQ(*r.sampling_error, uni20::math::sqrt((*expected.covariance)[0]));
  EXPECT_EQ(r.completed_samples, 33u);
  EXPECT_EQ(r.supported_roots, 66u);
  EXPECT_EQ(r.density_attempts, 132u);
  EXPECT_EQ(observed.metadata.observations, 33u);
  ASSERT_EQ(observed.replicate_metadata.size(), 3u);
  for (std::size_t i = 0; i < 3; ++i)
  {
    ASSERT_TRUE(r.replicates[i]);
    ASSERT_TRUE(observed.replicate_metadata[i]);
    EXPECT_EQ(r.replicates[i]->id, 9 + i);
    EXPECT_EQ(r.replicates[i]->completed_samples, 11u);
    EXPECT_EQ(observed.replicate_metadata[i]->observations, 11u);
    EXPECT_EQ((*r.replicates[i]->means)[0], (*expected.replicates[i]->means)[0]);
  }
}

TEST(SixSpinonObservedSpectrum, NoSupportFailureCancellationAndExceptionsStayDistinct)
{
  if (uni20::numeric_limits<long double>::digits(1.L) <= uni20::numeric_limits<double>::digits(1.)) GTEST_SKIP();
  ThermodynamicSixSpinonStructureFactor<long double, double, long double> f(proposal(1.L), native_policy(), 1.L);
  auto point = f.at(2 * std::atan(1.L), 2.L);
  auto empty = [](auto const&, std::stop_token) {
    ObservedSample sample;
    sample.values = std::vector{0.L, 0.L, 0.L};
    sample.metadata.observations = 1;
    return sample;
  };
  auto unsupported = point.evaluate_observations_async({2, 3}, empty).get_wait();
  EXPECT_EQ(unsupported.spectrum.status, SixSpinonSpectrumStatus::no_support_sampled);
  EXPECT_FALSE(unsupported.spectrum.value);
  EXPECT_FALSE(unsupported.spectrum.sampling_error);
  EXPECT_EQ(unsupported.metadata.observations, 6u);
  auto failure = [](auto const&, std::stop_token) {
    ObservedSample sample;
    sample.metadata.failure = SpinonObservationFailure::density;
    sample.metadata.density_status = SpinonRecoveryStatus::exhausted;
    sample.metadata.observations = 1;
    sample.reason = "controlled exhausted density";
    return sample;
  };
  auto failed = point.evaluate_observations_async({2, 3}, failure).get_wait();
  EXPECT_EQ(failed.spectrum.status, SixSpinonSpectrumStatus::numerical_failure);
  EXPECT_FALSE(failed.spectrum.value);
  EXPECT_EQ(failed.metadata.observations, 2u);
  for (auto const& row : failed.spectrum.replicates)
  {
    ASSERT_TRUE(row);
    EXPECT_EQ(row->failed_sample, 0u);
    EXPECT_EQ(row->failure, SpinonObservationFailure::density);
    EXPECT_EQ(row->density_status, SpinonRecoveryStatus::exhausted);
  }
  std::stop_source stop;
  stop.request_stop();
  auto never = [](auto const&, std::stop_token) -> ObservedSample { throw std::logic_error("deliberate exception"); };
  auto cancelled = point.evaluate_observations_async({2, 3}, never, {.stop = stop.get_token()}).get_wait();
  EXPECT_EQ(cancelled.spectrum.status, SixSpinonSpectrumStatus::cancelled);
  EXPECT_EQ(cancelled.metadata.observations, 0u);
  EXPECT_FALSE(cancelled.spectrum.value);
  auto invalid = point.evaluate_observations_async({2, 3}, never);
  EXPECT_THROW(invalid.get_wait(), std::logic_error);
}

#if UNI20_HAS_MPREAL
namespace
{
using R = uni20::mpreal;
using Spectrum = ThermodynamicSixSpinonStructureFactor<R, R>;
SpinonRecoveryPolicy mp_policy()
{
  return {"1e-7", "0", {{128, 192, {{{32, 12}, {64, 20}}}}, {224, 320, {{{64, 20}, {128, 24}}}}}, 60};
}
} // namespace

TEST(SixSpinonSpectrumMpfr, FreshReplicatesSerialAsyncAndExchangeScaling)
{
  R const one(1, uni20::Precision::bits(384)), pi = 4 * uni20::math::atan(one);
  Spectrum spectrum(proposal(one), mp_policy(), one);
  auto point = spectrum.at(pi / 2, 2 * one);
  SpinonSamplingOptions options{3, 8, 42, 0, 52};
  auto serial = point.evaluate(options, {.max_in_flight = 1});
  std::string failures;
  for (auto const& row : serial.replicates)
    if (row && !row->means) failures += " replicate " + std::to_string(row->id) + ": " + row->reason;
  ASSERT_EQ(serial.status, SixSpinonSpectrumStatus::estimated) << serial.reason << failures;
  ASSERT_TRUE(serial.value);
  ASSERT_TRUE(serial.sampling_error);
  EXPECT_GT(*serial.value, 0);
  EXPECT_GT(*serial.sampling_error, 0);
  RecordProperty("finite_cutoff_S6zz", uni20::format_scalar(*serial.value, {.precision = 16}));
  RecordProperty("sampling_error", uni20::format_scalar(*serial.sampling_error, {.precision = 16}));
  RecordProperty("contour_diagnostic", uni20::format_scalar(*serial.contour_diagnostic, {.precision = 16}));
  RecordProperty("arithmetic_diagnostic", uni20::format_scalar(*serial.arithmetic_diagnostic, {.precision = 16}));
  RecordProperty("serial_wall_seconds", std::to_string(serial.wall_seconds));
  EXPECT_EQ(serial.completed_samples, 24u);
  EXPECT_GT(serial.supported_roots, 0u);
  uni20::async::TbbScheduler scheduler(2);
  auto pending = point.evaluate_async(options, {.scheduler = &scheduler, .batch_size = 1, .max_in_flight = 2});
  auto concurrent = pending.get_wait(scheduler);
  ASSERT_EQ(concurrent.status, SixSpinonSpectrumStatus::estimated) << concurrent.reason;
  EXPECT_EQ(concurrent.value, serial.value);
  EXPECT_EQ(concurrent.sampling_error, serial.sampling_error);
  EXPECT_EQ(concurrent.contour_diagnostic, serial.contour_diagnostic);
  EXPECT_EQ(concurrent.arithmetic_diagnostic, serial.arithmetic_diagnostic);
  EXPECT_EQ(concurrent.completed_samples, serial.completed_samples);
  for (std::size_t i = 0; i < options.replicates; ++i)
    EXPECT_EQ(concurrent.replicates[i]->means, serial.replicates[i]->means);
  // The custom policy seam uses the same physical assembly and stream order.
  auto observe = [point](auto const& identity, std::stop_token stop) {
    auto value = point.observation(identity, stop);
    bethe::MonteCarloObservation<R, ObservedMetadata> sample;
    if (value.values) sample.values.emplace(value.values->begin(), value.values->end());
    static_cast<SpinonSampleMetadata&>(sample.metadata) = {value.supported_roots, value.density_attempts, value.failure,
                                                           value.density_status, value.density_failure};
    sample.metadata.observations = 1;
    sample.cancelled = value.failure == SpinonObservationFailure::cancelled;
    sample.reason = std::move(value.reason);
    return sample;
  };
  auto observed_pending = point.evaluate_observations_async(
      options, std::move(observe), {.scheduler = &scheduler, .batch_size = 3, .max_in_flight = 2});
  auto observed = observed_pending.get_wait(scheduler);
  EXPECT_EQ(observed.spectrum.status, serial.status);
  EXPECT_EQ(observed.spectrum.value, serial.value);
  EXPECT_EQ(observed.spectrum.sampling_error, serial.sampling_error);
  EXPECT_EQ(observed.spectrum.contour_diagnostic, serial.contour_diagnostic);
  EXPECT_EQ(observed.spectrum.arithmetic_diagnostic, serial.arithmetic_diagnostic);
  EXPECT_EQ(observed.spectrum.supported_roots, serial.supported_roots);
  EXPECT_EQ(observed.spectrum.density_attempts, serial.density_attempts);
  EXPECT_EQ(observed.metadata.observations, 24u);
  for (std::size_t i = 0; i < options.replicates; ++i)
  {
    EXPECT_EQ(observed.spectrum.replicates[i]->means, serial.replicates[i]->means);
    EXPECT_EQ(observed.replicate_metadata[i]->observations, 8u);
  }
  // The temporary point/evaluator can die before the owning Async completes.
  auto scaled_pending = Spectrum(proposal(one), mp_policy(), 3 * one)
                            .at(pi / 2, 6 * one)
                            .evaluate_async(options, {.scheduler = &scheduler, .batch_size = 1, .max_in_flight = 2});
  auto triple = scaled_pending.get_wait(scheduler);
  ASSERT_EQ(triple.status, SixSpinonSpectrumStatus::estimated) << triple.reason;
  R const tolerance("1e-100", one.precision());
  EXPECT_REAL_NEAR(3 * *triple.value, *serial.value, tolerance * *serial.value);
  EXPECT_REAL_NEAR(3 * *triple.sampling_error, *serial.sampling_error, tolerance * *serial.sampling_error);
}

TEST(SixSpinonSpectrumMpfr, FailedSamplesAndCancellationCannotBecomeSuccessfulMeans)
{
  R const one(1, uni20::Precision::bits(384)), pi = 4 * uni20::math::atan(one);
  auto policy = mp_policy();
  policy.stages = {{128, 192, {{{2, 1}, {4, 2}}}}};
  Spectrum spectrum(proposal(one), policy, one);
  auto point = spectrum.at(pi / 2, 2 * one);
  auto result = point.evaluate({3, 8, 42, 0, 52});
  EXPECT_EQ(result.status, SixSpinonSpectrumStatus::numerical_failure);
  EXPECT_FALSE(result.value);
  bool failed = false;
  for (auto const& row : result.replicates)
    if (row && row->failed_sample)
    {
      failed = true;
      EXPECT_EQ(row->failure, SpinonObservationFailure::density);
      EXPECT_EQ(row->density_status, SpinonRecoveryStatus::exhausted);
      EXPECT_FALSE(row->means);
    }
  EXPECT_TRUE(failed);
  std::stop_source stop;
  stop.request_stop();
  auto cancelled = point.evaluate({3, 8, 42, 0, 52}, {.stop = stop.get_token()});
  EXPECT_EQ(cancelled.status, SixSpinonSpectrumStatus::cancelled);
  EXPECT_FALSE(cancelled.value);
  EXPECT_EQ(cancelled.completed_samples, 0u);
  EXPECT_EQ(cancelled.density_attempts, 0u);
}

TEST(SixSpinonSpectrumMpfr, IndependentPhysicalObservationsCoverAllSections)
{
  R const one(1, uni20::Precision::bits(384)), pi = 4 * uni20::math::atan(one);
  Spectrum spectrum(proposal(one), mp_policy(), one);
  for (auto const& row : six_spinon_spectrum_reference::points)
  {
    SCOPED_TRACE(row.section);
    SCOPED_TRACE(row.omega);
    bethe::UniformSample<5>::Words words;
    words[0] = {32}; // background component: selector (32+1/2)/256
    for (std::size_t j = 0; j < 4; ++j)
      words[j + 1] = {row.words[j]};
    auto result = spectrum.at(pi / 2, R(row.omega, one.precision())).observation({8, words});
    ASSERT_TRUE(result.values) << result.reason;
    EXPECT_EQ(result.supported_roots, 2u);
    auto const& value = *result.values;
    R const reference(row.value, one.precision()), tail(row.tail, one.precision());
    EXPECT_LE(uni20::math::abs(value[0] - reference), value[1] + value[2] + tail);
    EXPECT_LT((value[1] + value[2]) / reference, R("1e-7", one.precision()));
  }
}

TEST(SixSpinonSpectrumMpfr, NoSupportedDrawIsNotAZeroSpectrum)
{
  R const one(1, uni20::Precision::bits(384)), pi = 4 * uni20::math::atan(one);
  Spectrum spectrum(proposal(one), mp_policy(), one);
  // Almost all unconditioned draws miss this tiny upper-edge neighbourhood.
  // Retain an explicit sampling status, not an apparently exact physical zero.
  auto result = spectrum.at(pi, 3 * pi - R("1e-30", one.precision())).evaluate({2, 1, 42, 0, 52});
  EXPECT_EQ(result.status, SixSpinonSpectrumStatus::no_support_sampled);
  EXPECT_EQ(result.completed_samples, 2u);
  EXPECT_EQ(result.supported_roots, 0u);
  EXPECT_FALSE(result.value);
  EXPECT_FALSE(result.sampling_error);
}

TEST(SixSpinonSpectrumMpfr, GeometryPrecisionIsDistinctFromDensityRecovery)
{
  auto const& row = six_spinon_spectrum_reference::points[std::size(six_spinon_spectrum_reference::points) - 1];
  bethe::UniformSample<5>::Words words;
  words[0] = {32};
  for (std::size_t j = 0; j < 4; ++j)
    words[j + 1] = {row.words[j]};
  std::optional<R> previous;
  for (unsigned bits : {384u, 512u})
  {
    R const one(1, uni20::Precision::bits(bits)), pi = 4 * uni20::math::atan(one);
    Spectrum spectrum(proposal(one), mp_policy(), one);
    auto result = spectrum.at(pi / 2, R(row.omega, one.precision())).observation({8, words});
    ASSERT_TRUE(result.values) << result.reason;
    auto const& value = (*result.values)[0];
    EXPECT_EQ(value.precision(), one.precision());
    if (previous)
      EXPECT_LT(uni20::math::abs(value - previous->at(one.precision())), R("1e-100", one.precision()) * value);
    previous = value;
  }
  R const too_narrow(1, uni20::Precision::bits(256));
  EXPECT_THROW((Spectrum(proposal(too_narrow), mp_policy(), too_narrow)), std::invalid_argument);
}
#endif
