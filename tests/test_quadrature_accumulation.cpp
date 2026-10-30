// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include <bethe/quadrature_accumulation.hpp>
#include <condition_variable>
#include <uni20/async/tbb_scheduler.hpp>

TEST(QuadratureAccumulation, OwnedVectorReductionAndEmptyPrecision)
{
  uni20::async::TbbScheduler scheduler(2);
  auto evaluate = [](std::size_t i, std::stop_token) {
    return bethe::QuadratureObservation<double, 2>{{{{double(i), double(i * i)}}}, 2};
  };
  auto pending = bethe::accumulate_quadrature_async<double, 2>(100, evaluate, 1., {&scheduler, 8, 3});
  auto result = pending.get_wait(scheduler);
  ASSERT_EQ(result.status, bethe::QuadratureStatus::complete);
  ASSERT_TRUE(result.sums);
  EXPECT_EQ((*result.sums)[0], 4950.);
  EXPECT_EQ((*result.sums)[1], 328350.);
  EXPECT_EQ(result.completed_nodes, 100u);
  EXPECT_EQ(result.work_units, 200u);
  auto empty = bethe::accumulate_quadrature_async<double, 2>(0, evaluate, 1., {&scheduler});
  ASSERT_TRUE(empty.get_wait(scheduler).sums);
  EXPECT_EQ((*empty.get_wait(scheduler).sums)[0], 0.);
}

TEST(QuadratureAccumulation, FailureCannotCompleteWhileAdmittedObservationIsHeld)
{
  struct State
  {
      std::mutex mutex;
      std::condition_variable cv;
      bool held = false, failed = false, release = false;
  };
  // The caller deliberately does not join the scheduler while observing the
  // held state; reserve room for both blocked/background observations.
  uni20::async::TbbScheduler scheduler(3);
  for (bool unexpected : {false, true})
  {
    auto state = std::make_shared<State>();
    auto pending = bethe::accumulate_quadrature_async<double, 1, int>(
        8,
        [state, unexpected](std::size_t i, std::stop_token) {
          std::unique_lock lock(state->mutex);
          if (i == 1)
          {
            state->held = true;
            state->cv.notify_all();
            if (!state->cv.wait_for(lock, std::chrono::seconds(10), [&] { return state->release; }))
              throw std::runtime_error("held observation timeout");
            return bethe::QuadratureObservation<double, 1, int>{{{{1.}}}};
          }
          if (!state->cv.wait_for(lock, std::chrono::seconds(10), [&] { return state->held; }))
            throw std::runtime_error("second observation not admitted");
          state->failed = true;
          state->cv.notify_all();
          if (unexpected) throw std::runtime_error("injected unexpected evaluator failure");
          return bethe::QuadratureObservation<double, 1, int>{{}, 3, 42, false, "injected numerical failure"};
        },
        1., {&scheduler, 8, 2});
    {
      std::unique_lock lock(state->mutex);
      bool const started =
          state->cv.wait_for(lock, std::chrono::seconds(10), [&] { return state->held && state->failed; });
      EXPECT_TRUE(started);
      auto reader = pending.read();
      EXPECT_FALSE(reader.await_ready());
      state->release = true;
      state->cv.notify_all();
    }
    if (unexpected)
    {
      EXPECT_THROW(pending.get_wait(scheduler), std::runtime_error);
    }
    else
    {
      auto result = pending.get_wait(scheduler);
      EXPECT_EQ(result.status, bethe::QuadratureStatus::numerical_failure);
      EXPECT_FALSE(result.sums);
      EXPECT_EQ(result.failed_node, 0u);
      EXPECT_EQ(result.failure_metadata, 42);
      EXPECT_EQ(result.completed_nodes, 0u);
    }
  }
}

TEST(QuadratureAccumulation, CancelledAndNonfiniteObservationsDoNotBecomeSums)
{
  uni20::async::TbbScheduler scheduler(1);
  auto evaluate = [](std::size_t, std::stop_token) {
    return bethe::QuadratureObservation<double, 1>{{{{std::numeric_limits<double>::infinity()}}}};
  };
  std::stop_source stop;
  stop.request_stop();
  auto cancelled = bethe::accumulate_quadrature_async<double, 1>(4, evaluate, 1., {&scheduler, 8, 2, stop.get_token()});
  EXPECT_EQ(cancelled.get_wait(scheduler).status, bethe::QuadratureStatus::cancelled);
  EXPECT_FALSE(cancelled.get_wait(scheduler).sums);
  auto nonfinite = bethe::accumulate_quadrature_async<double, 1>(4, evaluate, 1., {&scheduler, 8, 2});
  EXPECT_EQ(nonfinite.get_wait(scheduler).status, bethe::QuadratureStatus::numerical_failure);
  EXPECT_FALSE(nonfinite.get_wait(scheduler).sums);
  auto contradictory = bethe::accumulate_quadrature_async<double, 1>(
      1,
      [](std::size_t, std::stop_token) {
        return bethe::QuadratureObservation<double, 1>{{{{1.}}}, 0, {}, true, "cancelled with a value"};
      },
      1., {&scheduler});
  EXPECT_EQ(contradictory.get_wait(scheduler).status, bethe::QuadratureStatus::cancelled);
  EXPECT_FALSE(contradictory.get_wait(scheduler).sums);
}

#if UNI20_HAS_MPREAL
TEST(QuadratureAccumulationMpfr, PrecisionComesFromTheExplicitContext)
{
  using R = uni20::mpreal;
  auto const p = uni20::Precision::bits(160);
  R const one(1, p);
  uni20::async::TbbScheduler scheduler(1);
  auto evaluate = [one](std::size_t i, std::stop_token) {
    return bethe::QuadratureObservation<R, 1>{{{{uni20::scalar_like(one, i)}}}};
  };
  auto pending = bethe::accumulate_quadrature_async<R, 1>(4, evaluate, one, {&scheduler, 8, 2});
  auto result = pending.get_wait(scheduler);
  ASSERT_TRUE(result.sums);
  EXPECT_EQ((*result.sums)[0], 6);
  EXPECT_EQ((*result.sums)[0].precision(), p);
  auto empty = bethe::accumulate_quadrature_async<R, 1>(0, evaluate, one, {&scheduler});
  auto zero = empty.get_wait(scheduler);
  ASSERT_TRUE(zero.sums);
  EXPECT_EQ((*zero.sums)[0].precision(), p);
}
#endif
