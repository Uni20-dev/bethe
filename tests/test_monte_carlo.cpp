// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <atomic>
#include <bethe/monte_carlo.hpp>
#include <condition_variable>
#include <mutex>
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

TEST(MonteCarlo, HeldBlockDrainsBeforeExceptionalCompletion)
{
  struct State
  {
      std::mutex mutex;
      std::condition_variable cv;
      unsigned entered = 0;
      bool held = false, thrown = false, release = false;
  };
  auto state = std::make_shared<State>();
  uni20::async::TbbScheduler scheduler(3);
  auto pending = bethe::monte_carlo_async<1>(
      {1, 8, 7}, 1, 1.,
      [state](auto const&, std::stop_token) -> bethe::MonteCarloObservation<double> {
        std::unique_lock lock(state->mutex);
        if (state->entered++ == 0)
        {
          state->held = true;
          state->cv.notify_all();
          if (!state->cv.wait_for(lock, std::chrono::seconds(10), [&] { return state->release; }))
            throw std::runtime_error("held replicate timeout");
          return {std::vector{1.}};
        }
        state->thrown = true;
        state->cv.notify_all();
        throw std::logic_error("deliberate concurrent failure");
      },
      {.scheduler = &scheduler, .batch_size = 1, .max_in_flight = 2});
  {
    std::unique_lock lock(state->mutex);
    EXPECT_TRUE(state->cv.wait_for(lock, std::chrono::seconds(10), [&] { return state->held && state->thrown; }));
    auto reader = pending.read();
    EXPECT_FALSE(reader.await_ready());
    state->release = true;
    state->cv.notify_all();
  }
  EXPECT_THROW(pending.get_wait(scheduler), std::logic_error);
}
TEST(MonteCarlo, CancellationWithValueIsTerminalOnSingleWorker)
{
  uni20::async::TbbScheduler scheduler(1);
  auto pending = bethe::monte_carlo_async<1>({2, 1, 7}, 1, 1.,
                                             [](auto const&, std::stop_token) {
                                               bethe::MonteCarloObservation<double> r{std::vector{1.}};
                                               r.cancelled = true;
                                               r.reason = "cancelled after evaluation";
                                               return r;
                                             },
                                             {.scheduler = &scheduler, .max_in_flight = 2});
  auto result = pending.get_wait(scheduler);
  EXPECT_EQ(result.status, bethe::MonteCarloStatus::cancelled);
  EXPECT_FALSE(result.mean);
  EXPECT_FALSE(result.covariance);
  EXPECT_EQ(result.completed_samples, 0u);
  for (auto const& row : result.replicates)
  {
    ASSERT_TRUE(row);
    EXPECT_EQ(row->failed_sample, 0u);
  }
}

TEST(MonteCarlo, SingleReplicateUsesMultipleWorkers)
{
  struct State
  {
      std::mutex mutex;
      std::condition_variable cv;
      unsigned calls = 0;
      bool overlap = false;
  };
  auto state = std::make_shared<State>();
  uni20::async::TbbScheduler scheduler(3);
  auto pending = bethe::monte_carlo_async<1>({1, 8, 23}, 1, 1.,
                                             [state](auto const&, std::stop_token) {
                                               std::unique_lock lock(state->mutex);
                                               if (++state->calls == 1)
                                                 state->overlap = state->cv.wait_for(lock, std::chrono::seconds(2),
                                                                                     [&] { return state->calls >= 2; });
                                               else
                                                 state->cv.notify_all();
                                               return bethe::MonteCarloObservation<double>{std::vector{1.}};
                                             },
                                             {.scheduler = &scheduler, .batch_size = 2, .max_in_flight = 4});
  auto result = pending.get_wait(scheduler);
  EXPECT_TRUE(state->overlap);
  EXPECT_EQ(result.status, bethe::MonteCarloStatus::complete);
  EXPECT_EQ(result.completed_samples, 8u);
  ASSERT_TRUE(result.mean);
  EXPECT_EQ(result.mean->front(), 1.);
}

namespace
{
template <class R> void check_sample_blocks(R one)
{
  bethe::MonteCarloOptions const options{3, 37, 9182, 5, 52};
  auto observe = mixture_integrand(one);
  std::vector<bethe::MonteCarloReplicate<R, Counters>> serial;
  for (std::size_t i = 0; i < options.replicates; ++i)
    serial.push_back(bethe::monte_carlo_replicate<2>(options, i, 3, one, observe));
  uni20::async::TbbScheduler scheduler(3);
  for (std::size_t grain : {1, 4, 16, 128})
  {
    SCOPED_TRACE(grain);
    auto result = bethe::monte_carlo_async<2>(options, 3, one, observe,
                                              {.scheduler = &scheduler, .batch_size = grain, .max_in_flight = 5})
                      .get_wait(scheduler);
    ASSERT_EQ(result.status, bethe::MonteCarloStatus::complete);
    EXPECT_EQ(result.completed_samples, options.replicates * options.samples_per_replicate);
    for (std::size_t i = 0; i < serial.size(); ++i)
    {
      ASSERT_TRUE(result.replicates[i]);
      auto const& actual = *result.replicates[i];
      EXPECT_EQ(actual.id, serial[i].id);
      EXPECT_EQ(actual.means, serial[i].means);
      EXPECT_EQ(actual.maximum_absolute_observation, serial[i].maximum_absolute_observation);
      EXPECT_EQ(actual.metadata.uniform, serial[i].metadata.uniform);
      EXPECT_EQ(actual.metadata.linear, serial[i].metadata.linear);
    }
  }
}
} // namespace

TEST(MonteCarlo, SampleBlocksPreserveSerialDrawsAndArithmetic)
{
  check_sample_blocks(1.);
  check_sample_blocks(1.L);
}

#if UNI20_HAS_MPREAL
TEST(MonteCarloMpfr, SampleBlocksPreserveExplicitPrecision)
{
  check_sample_blocks(uni20::mpreal(1, uni20::Precision::bits(160)));
}
#endif

TEST(MonteCarlo, StreamDiscardPreservesWideUniformIdentities)
{
  for (unsigned bits : {1u, 52u, 64u, 65u, 129u})
    for (std::size_t count : {0, 1, 7, 64, 513})
    {
      bethe::MonteCarloStream<3> serial(81, 19, bits), blocked(81, 19, bits);
      for (std::size_t i = 0; i < count; ++i)
        (void)serial.next();
      blocked.discard(count);
      for (int i = 0; i < 4; ++i)
        EXPECT_EQ(serial.next().words(), blocked.next().words());
    }
}

TEST(MonteCarlo, BlockFailuresRetainTheFirstLogicalFailure)
{
  auto observe = [](auto const& identity, std::stop_token) {
    auto const x = identity.at_precision(1.)[0];
    return x < .7 ? bethe::MonteCarloObservation<double>{std::vector{x}}
                  : bethe::MonteCarloObservation<double>{{}, {}, false, "deliberate failure"};
  };
  bethe::MonteCarloOptions options{3, 31, 101};
  uni20::async::TbbScheduler scheduler(3);
  for (std::size_t grain : {1, 3, 8})
  {
    auto result = bethe::monte_carlo_async<1>(options, 1, 1., observe,
                                              {.scheduler = &scheduler, .batch_size = grain, .max_in_flight = 5})
                      .get_wait(scheduler);
    EXPECT_EQ(result.status, bethe::MonteCarloStatus::numerical_failure);
    EXPECT_FALSE(result.mean);
    for (std::size_t i = 0; i < options.replicates; ++i)
    {
      auto serial = bethe::monte_carlo_replicate<1>(options, i, 1, 1., observe);
      ASSERT_TRUE(result.replicates[i]);
      EXPECT_EQ(result.replicates[i]->failed_sample, serial.failed_sample);
      EXPECT_EQ(result.replicates[i]->completed_samples, serial.completed_samples);
      EXPECT_EQ(result.replicates[i]->maximum_absolute_observation, serial.maximum_absolute_observation);
    }
  }
}

TEST(MonteCarlo, BlockCancellationRetainsCompletedObservations)
{
  auto stop = std::make_shared<std::stop_source>();
  auto calls = std::make_shared<std::atomic<unsigned>>(0);
  uni20::async::TbbScheduler scheduler(3);
  auto result = bethe::monte_carlo_async<1>(
                    {2, 32, 19}, 1, 1.,
                    [stop, calls](auto const&, std::stop_token) {
                      ++*calls;
                      stop->request_stop();
                      return bethe::MonteCarloObservation<double>{std::vector{1.}};
                    },
                    {.scheduler = &scheduler, .batch_size = 2, .max_in_flight = 3, .stop = stop->get_token()})
                    .get_wait(scheduler);
  EXPECT_EQ(result.status, bethe::MonteCarloStatus::cancelled);
  EXPECT_FALSE(result.mean);
  EXPECT_GT(calls->load(), 0u);
  EXPECT_EQ(result.completed_samples, calls->load());
  EXPECT_LE(result.completed_samples, 3u);
}

TEST(MonteCarlo, CancellationAtBlockBoundaryClassifiesUnfinishedReplicate)
{
  for (unsigned workers : {1u, 2u})
    for (std::size_t samples : {4u, 8u})
    {
      SCOPED_TRACE(workers);
      SCOPED_TRACE(samples);
      std::stop_source stop;
      unsigned calls = 0;
      uni20::async::TbbScheduler scheduler(workers);
      auto result = bethe::monte_carlo_async<1>(
                        {2, samples, 19}, 1, 1.,
                        [&](auto const&, std::stop_token) {
                          if (++calls == 4) stop.request_stop();
                          return bethe::MonteCarloObservation<double>{std::vector{1.}};
                        },
                        {.scheduler = &scheduler, .batch_size = 4, .max_in_flight = 1,
                         .stop = stop.get_token()})
                        .get_wait(scheduler);
      EXPECT_EQ(result.status, bethe::MonteCarloStatus::cancelled);
      EXPECT_FALSE(result.mean);
      EXPECT_FALSE(result.covariance);
      EXPECT_EQ(calls, 4u);
      ASSERT_TRUE(result.replicates[0]);
      auto const& row = *result.replicates[0];
      EXPECT_EQ(row.completed_samples, 4u);
      EXPECT_EQ(row.cancelled, samples == 8);
      if (samples == 8)
      {
        EXPECT_EQ(row.failed_sample, 4u);
        EXPECT_EQ(row.reason, "replicate cancelled");
        EXPECT_FALSE(row.means);
      }
      else
      {
        EXPECT_FALSE(row.failed_sample);
        ASSERT_TRUE(row.means);
        EXPECT_EQ(row.means->front(), 1.);
      }
      EXPECT_FALSE(result.replicates[1]); // Never admitted; no fabricated row.
    }
}

TEST(MonteCarlo, CancelledCountsIncludeCompletedSpeculativeBlocks)
{
  bethe::MonteCarloOptions options{1, 12, 19};
  std::vector<std::uint64_t> words;
  bethe::MonteCarloStream<1> stream(options.seed, 0, options.uniform_bits);
  for (std::size_t i = 0; i < options.samples_per_replicate; ++i)
    words.push_back(stream.next().words()[0][0]);
  std::stop_source stop;
  std::mutex mutex;
  std::condition_variable cv;
  bool first_entered = false;
  unsigned later_calls = 0;
  uni20::async::TbbScheduler scheduler(3);
  auto result = bethe::monte_carlo_async<1>(
                    options, 1, 1.,
                    [&](auto const& identity, std::stop_token) {
                      auto found = std::ranges::find(words, identity.words()[0][0]);
                      if (found == words.end()) throw std::logic_error("unexpected sample identity");
                      auto index = std::size_t(found - words.begin());
                      std::unique_lock lock(mutex);
                      if (index == 0)
                      {
                        first_entered = true;
                        cv.notify_all();
                        if (!cv.wait_for(lock, std::chrono::seconds(10), [&] { return stop.stop_requested(); }))
                          throw std::runtime_error("later blocks did not complete");
                      }
                      else
                      {
                        if (!cv.wait_for(lock, std::chrono::seconds(10), [&] { return first_entered; }))
                          throw std::runtime_error("first block did not enter");
                        if (++later_calls == 8)
                        {
                          stop.request_stop();
                          cv.notify_all();
                        }
                      }
                      return bethe::MonteCarloObservation<double, Counters>{std::vector{double(index)}, {1, 0}};
                    },
                    {.scheduler = &scheduler, .batch_size = 4, .max_in_flight = 3, .stop = stop.get_token()})
                    .get_wait(scheduler);
  EXPECT_EQ(result.status, bethe::MonteCarloStatus::cancelled);
  EXPECT_FALSE(result.mean);
  EXPECT_FALSE(result.covariance);
  ASSERT_TRUE(result.replicates[0]);
  auto const& row = *result.replicates[0];
  EXPECT_TRUE(row.cancelled);
  EXPECT_EQ(row.failed_sample, 1u); // First missing index, not the completed count.
  EXPECT_EQ(row.completed_samples, 9u); // Sample 0 and both later four-sample blocks.
  EXPECT_EQ(row.metadata.uniform, 9u);
  EXPECT_EQ(row.maximum_absolute_observation, std::vector{11.});
  EXPECT_FALSE(row.means);
  EXPECT_EQ(result.completed_samples, 9u);
  EXPECT_EQ(result.metadata.uniform, 9u);
}

TEST(MonteCarlo, ObservationStorageIsBoundedByAdmittedBlocks)
{
  struct Counts
  {
      std::atomic<int> live{0}, peak{0};
  } counts;
  struct Payload
  {
      Counts& counts;
      explicit Payload(Counts& c) : counts(c)
      {
        int n = ++counts.live, old = counts.peak.load();
        while (old < n && !counts.peak.compare_exchange_weak(old, n))
        {}
      }
      ~Payload() { --counts.live; }
  };
  struct Metadata
  {
      std::shared_ptr<Payload> payload;
      Metadata& operator+=(Metadata const&) { return *this; }
  };
  uni20::async::TbbScheduler scheduler(3);
  auto pending = bethe::monte_carlo_async<1>({4, 257, 19}, 1, 1.,
                                             [&](auto const&, std::stop_token) {
                                               return bethe::MonteCarloObservation<double, Metadata>{
                                                   std::vector{1.}, {std::make_shared<Payload>(counts)}};
                                             },
                                             {.scheduler = &scheduler, .batch_size = 7, .max_in_flight = 5});
  auto result = pending.get_wait(scheduler);
  EXPECT_EQ(result.status, bethe::MonteCarloStatus::complete);
  EXPECT_EQ(result.completed_samples, 4u * 257);
  EXPECT_LE(counts.peak, 35);
  EXPECT_EQ(counts.live, 0);
}

TEST(MonteCarlo, TerminalReplicatesDoNotStartNewObservationBlocks)
{
  bethe::MonteCarloOptions options{3, 65, 71, 9, 52};
  // Make every identity in the first replicate fail, while the other two
  // replicates remain healthy. Verify work performed, not merely the sums.
  std::vector<std::uint64_t> failing;
  bethe::MonteCarloStream<1> stream(options.seed, options.first_replicate, options.uniform_bits);
  for (std::size_t i = 0; i < options.samples_per_replicate; ++i)
    failing.push_back(stream.next().words()[0][0]);
  uni20::async::TbbScheduler scheduler(3);
  for (std::size_t window : {1, 3})
  {
    SCOPED_TRACE(window);
    std::atomic<std::size_t> failed_calls = 0, healthy_calls = 0;
    auto result = bethe::monte_carlo_async<1>(options, 1, 1.,
                                              [&](auto const& identity, std::stop_token) {
                                                if (std::ranges::find(failing, identity.words()[0][0]) != failing.end())
                                                {
                                                  ++failed_calls;
                                                  return bethe::MonteCarloObservation<double>{
                                                      {}, {}, false, "terminal replicate"};
                                                }
                                                ++healthy_calls;
                                                return bethe::MonteCarloObservation<double>{std::vector{1.}};
                                              },
                                              {.scheduler = &scheduler, .batch_size = 4, .max_in_flight = window})
                      .get_wait(scheduler);
    EXPECT_EQ(result.status, bethe::MonteCarloStatus::numerical_failure);
    EXPECT_FALSE(result.mean);
    EXPECT_GE(failed_calls, 1u);
    EXPECT_LE(failed_calls, window); // At most the already-admitted window can run.
    EXPECT_EQ(healthy_calls, 2 * options.samples_per_replicate);
    ASSERT_TRUE(result.replicates[0]);
    EXPECT_EQ(result.replicates[0]->failed_sample, 0u);
    EXPECT_EQ(result.replicates[0]->completed_samples, 0u);
    for (std::size_t i = 1; i < options.replicates; ++i)
    {
      ASSERT_TRUE(result.replicates[i]);
      EXPECT_EQ(result.replicates[i]->completed_samples, options.samples_per_replicate);
      ASSERT_TRUE(result.replicates[i]->means);
      EXPECT_EQ(result.replicates[i]->means->front(), 1.);
    }
  }
}
