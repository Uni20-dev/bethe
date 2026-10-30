// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <atomic>
#include <bethe/monte_carlo.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
struct Counters
{
    std::size_t uniform = 0, linear = 0;
    Counters& operator+=(Counters const& b)
    {
      uniform += b.uniform;
      linear += b.linear;
      return *this;
    }
};

template <class R> auto mixture_integrand(R one)
{
  return [one](bethe::UniformSample<2> const& identity, std::stop_token) {
    auto u = identity.at_precision(one);
    bool const linear = u[0] >= one / R{2};
    R const x = linear ? uni20::math::sqrt(u[1]) : u[1];
    R const q = one / R{2} + x; // COMPLETE mixture, including both channels
    R const weight = x * x / q;
    return bethe::MonteCarloObservation<R, Counters>{std::vector<R>{weight, R{3} * weight, one},
                                                     {std::size_t(!linear), std::size_t(linear)}};
  };
}

template <class R> auto run_known(R one)
{
  bethe::MonteCarloOptions options{8, 256, 108, 0, 52};
  uni20::async::TbbScheduler scheduler(2);
  auto serial = bethe::monte_carlo_async<2>(options, 3, one, mixture_integrand(one), {.max_in_flight = 1});
  auto a = serial.get_wait();
  auto concurrent = bethe::monte_carlo_async<2>(options, 3, one, mixture_integrand(one),
                                                {.scheduler = &scheduler, .max_in_flight = 3});
  auto b = concurrent.get_wait(scheduler);
  EXPECT_EQ(a.status, bethe::MonteCarloStatus::complete);
  EXPECT_EQ(a.mean, b.mean);
  EXPECT_EQ(a.covariance, b.covariance);
  EXPECT_EQ(a.completed_samples, 2048u);
  EXPECT_EQ(a.metadata.uniform + a.metadata.linear, 2048u);
  EXPECT_GT(a.metadata.uniform, 0u);
  EXPECT_GT(a.metadata.linear, 0u);
  EXPECT_TRUE(a.mean && a.covariance);
  if (a.mean && a.covariance)
  {
    auto const& v = *a.mean;
    auto const& c = *a.covariance;
    R const eps = uni20::numeric_limits<R>::epsilon(one);
    EXPECT_LT(uni20::math::abs(v[0] - one / R{3}), R{6} * uni20::math::sqrt(c[0]));
    EXPECT_REAL_NEAR(v[1], R{3} * v[0], R{32} * eps);
    EXPECT_EQ(v[2], one);
    EXPECT_REAL_NEAR(c[1], R{3} * c[0], R{128} * eps * c[0]);
    EXPECT_REAL_NEAR(c[4], R{9} * c[0], R{1024} * eps * c[0]);
    EXPECT_EQ(c[8], one - one);
    for (auto const& row : a.replicates)
    {
      EXPECT_GT(row->maximum_absolute_observation[0], (*row->means)[0]);
      EXPECT_EQ(row->maximum_absolute_observation[2], one);
    }
  }
  return a;
}
} // namespace

TEST(MonteCarlo, ImportanceMixtureVectorCovarianceAndScheduling) { run_known(1.); }

TEST(MonteCarlo, FailuresKeepIdentityAndCannotFormSurvivorMean)
{
  auto observe = [](auto const& identity, std::stop_token) {
    auto u = identity.at_precision(1.);
    bethe::MonteCarloObservation<double> result;
    if (u[0] < .5)
      result.values = std::vector{u[0]};
    else
      result.reason = "deliberate failure";
    return result;
  };
  auto pending = bethe::monte_carlo_async<1>({3, 16, 19}, 1, 1., observe);
  auto result = pending.get_wait();
  EXPECT_EQ(result.status, bethe::MonteCarloStatus::numerical_failure);
  EXPECT_FALSE(result.mean);
  EXPECT_FALSE(result.covariance);
  for (auto const& row : result.replicates)
  {
    ASSERT_TRUE(row);
    ASSERT_TRUE(row->failed_sample);
    bethe::MonteCarloStream<1> stream(19, row->id, 52);
    for (std::size_t i = 0; i < *row->failed_sample; ++i)
      EXPECT_LT(stream.next().at_precision(1.)[0], .5);
    EXPECT_GE(stream.next().at_precision(1.)[0], .5);
    EXPECT_EQ(row->completed_samples, *row->failed_sample);
    EXPECT_EQ(row->reason, "deliberate failure");
  }
}

TEST(MonteCarlo, InFlightCancellationDrainsAndRetainsCompletedObservations)
{
  auto stop = std::make_shared<std::stop_source>();
  auto calls = std::make_shared<std::atomic<unsigned>>(0);
  uni20::async::TbbScheduler scheduler(2);
  auto pending = bethe::monte_carlo_async<1>({8, 64, 7}, 1, 1.,
                                             [stop, calls](auto const&, std::stop_token) {
                                               calls->fetch_add(1);
                                               stop->request_stop();
                                               return bethe::MonteCarloObservation<double>{std::vector{1.}};
                                             },
                                             {.scheduler = &scheduler, .max_in_flight = 2, .stop = stop->get_token()});
  auto result = pending.get_wait(scheduler);
  EXPECT_EQ(result.status, bethe::MonteCarloStatus::cancelled);
  EXPECT_FALSE(result.mean);
  EXPECT_GT(result.completed_samples, 0u);
  EXPECT_EQ(result.completed_samples, calls->load());
  EXPECT_LE(result.completed_samples, 2u);
}

TEST(MonteCarlo, InvalidContractsPropagateAndNumericalRangeFailsExplicitly)
{
  auto observe = [](auto const&, std::stop_token) { return bethe::MonteCarloObservation<double>{std::vector{1.}}; };
  EXPECT_THROW((bethe::monte_carlo_async<1>({0, 2}, 1, 1., observe)), std::invalid_argument);
  EXPECT_THROW((bethe::monte_carlo_async<1>({2, 2}, 0, 1., observe)), std::invalid_argument);
  auto wrong = bethe::monte_carlo_async<1>({2, 2}, 2, 1., observe);
  EXPECT_THROW(wrong.get_wait(), std::invalid_argument);
  auto nonfinite = bethe::monte_carlo_async<1>({2, 2}, 1, 1., [](auto const&, std::stop_token) {
    return bethe::MonteCarloObservation<double>{std::vector{std::numeric_limits<double>::infinity()}};
  });
  auto result = nonfinite.get_wait();
  EXPECT_EQ(result.status, bethe::MonteCarloStatus::numerical_failure);
  EXPECT_FALSE(result.mean);
  EXPECT_EQ(result.completed_samples, 0u);
  auto throwing = bethe::monte_carlo_async<1>({2, 2}, 1, 1.,
                                              [](auto const&, std::stop_token) -> bethe::MonteCarloObservation<double> {
                                                throw std::logic_error("programmer error");
                                              });
  EXPECT_THROW(throwing.get_wait(), std::logic_error);
}

#if UNI20_HAS_MPREAL
TEST(MonteCarloMpfr, KeepsPrecisionAndVectorRelations)
{
  for (unsigned bits : {160u, 256u})
  {
    uni20::mpreal one(1, uni20::Precision::bits(bits));
    auto result = run_known(one);
    ASSERT_TRUE(result.mean && result.covariance);
    for (auto const& value : *result.mean)
      EXPECT_EQ(value.precision(), one.precision());
    for (auto const& value : *result.covariance)
      EXPECT_EQ(value.precision(), one.precision());
  }
}
#endif
