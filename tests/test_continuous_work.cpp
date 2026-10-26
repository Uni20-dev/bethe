// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <array>
#include <bethe/detail/continuous_work.hpp>
#include <condition_variable>
#include <oneapi/tbb/info.h>
#include <uni20/async/tbb_scheduler.hpp>

using bethe::detail::continuous_work;

TEST(ContinuousWork, CrossesLogicalBatchBoundaryBeforeSlowPointFinishes)
{
  if (oneapi::tbb::info::default_concurrency() < 2)
    GTEST_SKIP() << "cross-batch concurrency requires at least two available CPUs";
  uni20::async::TbbScheduler scheduler(2);
  std::mutex mutex;
  std::condition_variable ready;
  bool later_batch_finished = false, crossed = false;
  // Logical batches [0,1] and [2,3]. A join after [0,1] necessarily fails
  // this test, even if both points within that batch execute in parallel.
  auto completed = continuous_work(scheduler, 2, 4, [&](std::size_t i, std::stop_token) {
    if (i == 0)
    {
      std::unique_lock lock(mutex);
      crossed = ready.wait_for(lock, std::chrono::seconds(10), [&] { return later_batch_finished; });
    }
    if (i == 2)
    {
      std::lock_guard lock(mutex);
      later_batch_finished = true;
      ready.notify_all();
    }
  });
  EXPECT_EQ(completed, 4u);
  EXPECT_TRUE(crossed);
}

TEST(ContinuousWork, EveryIndexOnceWithBoundedActiveCallbacks)
{
  uni20::async::TbbScheduler scheduler(4);
  std::array<std::atomic<int>, 257> visits{};
  std::atomic<int> active{0}, peak{0};
  auto completed = continuous_work(scheduler, 3, visits.size(), [&](std::size_t i, std::stop_token) {
    auto current = ++active;
    auto old = peak.load();
    while (old < current && !peak.compare_exchange_weak(old, current))
    {}
    ++visits[i];
    --active;
  });
  EXPECT_EQ(completed, visits.size());
  EXPECT_LE(peak.load(), 3);
  for (auto const& visits_to_point : visits)
    EXPECT_EQ(visits_to_point.load(), 1);
}

TEST(ContinuousWork, SerialCancellationAndExceptionsStopAdmission)
{
  uni20::async::TbbScheduler scheduler(1);
  std::stop_source stop;
  int calls = 0;
  EXPECT_EQ(continuous_work(
                scheduler, 1, 20,
                [&](auto, std::stop_token token) {
                  ++calls;
                  stop.request_stop();
                  EXPECT_TRUE(token.stop_requested());
                },
                stop.get_token()),
            1u);
  EXPECT_EQ(calls, 1);
  EXPECT_EQ(continuous_work(scheduler, 1, 20, [&](auto, auto) { ++calls; }, stop.get_token()), 0u);
  EXPECT_THROW(continuous_work(scheduler, 1, 20,
                               [&](auto, auto) {
                                 ++calls;
                                 throw std::runtime_error("publication failed");
                               }),
               std::runtime_error);
  EXPECT_EQ(calls, 2);
  EXPECT_EQ(continuous_work(scheduler, 1, 0, [&](auto, auto) { ++calls; }), 0u);
  EXPECT_THROW(continuous_work(scheduler, 0, 1, [&](auto, auto) {}), std::invalid_argument);
}

TEST(ContinuousWork, ErrorDrainsActiveWorkBeforePropagating)
{
  if (oneapi::tbb::info::default_concurrency() < 2)
    GTEST_SKIP() << "concurrent cancellation requires at least two available CPUs";
  uni20::async::TbbScheduler scheduler(2);
  std::mutex mutex;
  std::condition_variable_any ready;
  bool second_started = false, drained = false;
  EXPECT_THROW(continuous_work(scheduler, 2, 20,
                               [&](std::size_t i, std::stop_token token) {
                                 std::unique_lock lock(mutex);
                                 if (i == 0)
                                 {
                                   EXPECT_TRUE(
                                       ready.wait_for(lock, std::chrono::seconds(10), [&] { return second_started; }));
                                   throw std::runtime_error("checkpoint publication failed");
                                 }
                                 EXPECT_EQ(i, 1u);
                                 second_started = true;
                                 ready.notify_all();
                                 ready.wait_for(lock, token, std::chrono::seconds(10), [] { return false; });
                                 drained = token.stop_requested();
                               }),
               std::runtime_error);
  EXPECT_TRUE(drained);
}
