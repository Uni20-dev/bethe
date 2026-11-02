// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxx_four_spinon_thermodynamics.hpp>
#include <bethe/xxx_six_spinon_moments.hpp>
#include <bethe/xxx_six_spinon_structure_factor.hpp>
#include <bethe/xxx_two_spinon_moments.hpp>
#include <condition_variable>
#include <oneapi/tbb/global_control.h>
#include <uni20/async/debug_scheduler.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using namespace uni20::async;
class RejectSubmission : public DebugScheduler {
  public:
    explicit RejectSubmission(std::size_t n) : reject(n) {}
    std::size_t calls = 0, reject;
    void schedule(AsyncTask&& task) override
    {
      if (++calls == reject) throw std::bad_alloc{};
      DebugScheduler::schedule(std::move(task));
    }
};

template <class Submit> void check_owner(Submit submit)
{
  RejectSubmission first(1);
  EXPECT_THROW((void)submit(first), std::bad_alloc);
  EXPECT_TRUE(first.done());

  // Only one orchestration task may be admitted by the public call. Later
  // observation admission failures belong to its returned result and drain
  // any blocks accepted before the failing admission.
  for (std::size_t reject : {2u, 3u})
  {
    SCOPED_TRACE(reject);
    RejectSubmission scheduler(reject);
    try
    {
      auto pending = submit(scheduler);
      EXPECT_EQ(scheduler.calls, 1u);
      EXPECT_THROW((void)pending.get_wait(scheduler), std::bad_alloc);
    }
    catch (std::bad_alloc const&)
    {
      ADD_FAILURE() << "public call admitted a dependency before failing";
    }
    EXPECT_TRUE(scheduler.done());
    if (!scheduler.done()) scheduler.run_all(); // Drain only a failing old-code test.
  }
}

struct ActiveObserver
{
    std::mutex mutex;
    std::condition_variable cv;
    bool started = false, rejected = false, stopped = false, release = false;
};

class RejectAfterObserverStarts : public TbbScheduler {
    std::shared_ptr<ActiveObserver> state_;
    std::atomic<std::size_t> calls_{0};

  public:
    explicit RejectAfterObserverStarts(std::shared_ptr<ActiveObserver> state)
        : TbbScheduler(3), state_(std::move(state))
    {}
    void schedule(AsyncTask&& task) override
    {
      if (++calls_ == 3)
      {
        std::unique_lock lock(state_->mutex);
        if (!state_->cv.wait_for(lock, std::chrono::seconds(5), [&] { return state_->started; }))
          throw std::runtime_error("first observer did not start");
        state_->rejected = true;
        state_->cv.notify_all();
        throw std::bad_alloc{};
      }
      TbbScheduler::schedule(std::move(task));
    }
};

template <class Submit> void check_admission_stop(Submit submit)
{
  // The coordinator deliberately blocks admission until an observer starts.
  // Permit both background workers even when CPU affinity limits the default
  // pool to one. This scoped test budget does not require three physical CPUs.
  oneapi::tbb::global_control workers(oneapi::tbb::global_control::max_allowed_parallelism, 3);
  auto state = std::make_shared<ActiveObserver>();
  RejectAfterObserverStarts scheduler(state);
  std::stop_source external;
  auto wait_for_stop = [state](std::stop_token stop) {
    std::stop_callback callback(stop, [state] {
      std::lock_guard lock(state->mutex);
      state->stopped = true;
      state->cv.notify_all();
    });
    std::unique_lock lock(state->mutex);
    state->started = true;
    state->cv.notify_all();
    if (!state->cv.wait_for(lock, std::chrono::seconds(10), [&] { return state->stopped || state->release; }))
      throw std::runtime_error("observer did not receive stop or test release");
  };
  auto pending = submit(scheduler, wait_for_stop, external.get_token());
  {
    std::unique_lock lock(state->mutex);
    EXPECT_TRUE(state->cv.wait_for(lock, std::chrono::seconds(5), [&] { return state->rejected; }));
    EXPECT_TRUE(state->cv.wait_for(lock, std::chrono::seconds(5), [&] { return state->stopped; }));
    // Bound an old-code failure without leaking the deliberately held worker.
    state->release = true;
    state->cv.notify_all();
  }
  EXPECT_THROW((void)pending.get_wait(scheduler), std::bad_alloc);
  EXPECT_FALSE(external.stop_requested()); // failure is local to this operation
}
} // namespace

TEST(SamplingOwnership, QuadratureAdmissionFailureStopsActiveObserver)
{
  check_admission_stop([](auto& scheduler, auto wait_for_stop, std::stop_token external) {
    return bethe::accumulate_quadrature_async<double, 1>(
        4,
        [wait_for_stop](std::size_t, std::stop_token stop) {
          wait_for_stop(stop);
          return bethe::QuadratureObservation<double, 1>{{{{1.}}}};
        },
        1., {&scheduler, 1, 2, external});
  });
}

TEST(SamplingOwnership, MonteCarloAdmissionFailureStopsActiveObserver)
{
  check_admission_stop([](auto& scheduler, auto wait_for_stop, std::stop_token external) {
    return bethe::monte_carlo_async<1>(
        {2, 2, 42}, 1, 1.,
        [wait_for_stop](auto const&, std::stop_token stop) {
          wait_for_stop(stop);
          return bethe::MonteCarloObservation<double>{std::vector{1.}};
        },
        {&scheduler, 1, 2, external});
  });
}

TEST(SamplingOwnership, MonteCarloAdmissionFailuresRetainOneOwner)
{
  check_owner([](auto& scheduler) {
    return bethe::monte_carlo_async<1>(
        {2, 2, 42}, 1, 1.,
        [](auto const&, std::stop_token) { return bethe::MonteCarloObservation<double>{std::vector{1.}}; },
        {&scheduler, 1, 2});
  });
}

TEST(SamplingOwnership, ThrowingObserverMoveCannotAdmitWork)
{
  struct ThrowMove
  {
      ThrowMove() = default;
      ThrowMove(ThrowMove&&) { throw std::bad_alloc{}; }
      bethe::MonteCarloObservation<double> operator()(bethe::UniformSample<1> const&, std::stop_token) const
      {
        return {std::vector{1.}};
      }
  };
  RejectSubmission scheduler(99);
  EXPECT_THROW((void)bethe::monte_carlo_async<1>({2, 2}, 1, 1., ThrowMove{}, {&scheduler}), std::bad_alloc);
  EXPECT_EQ(scheduler.calls, 0u);
  EXPECT_TRUE(scheduler.done());
}

TEST(SamplingOwnership, QuadratureAdmissionFailuresRetainOneOwner)
{
  check_owner([](auto& scheduler) {
    return bethe::accumulate_quadrature_async<double, 1>(
        4, [](std::size_t, std::stop_token) { return bethe::QuadratureObservation<double, 1>{{{{1.}}}}; }, 1.,
        {&scheduler, 1, 2});
  });
}

TEST(SamplingOwnership, TwoSpinonMomentsOwnTheirQuadratureDependency)
{
  using namespace bethe::heisenberg;
  StructureFactorOptions<double> options;
  options.tolerance = 1e-10;
  ThermodynamicTwoSpinonMoments<double> model(1., options);
  check_owner([&](auto& scheduler) { return model.global_moments(8, 4.).evaluate_async({&scheduler, 1, 2}); });
}

TEST(SamplingOwnership, FourSpinonMomentsOwnTheirQuadratureDependency)
{
  using namespace bethe::heisenberg;
  if (uni20::numeric_limits<long double>::digits() <= uni20::numeric_limits<double>::digits()) GTEST_SKIP();
  SpinonRecoveryPolicy policy{
      "1e-5", "0", {{53, uni20::numeric_limits<long double>::digits(), {{{32, 12}, {64, 20}}}}}, 60};
  ThermodynamicFourSpinon<long double, double, long double> model(policy, 1.L);
  check_owner([&](auto& scheduler) { return model.global_moments(2, {2.L, 2.L, 2.L}).evaluate_async({&scheduler, 1, 2}); });
}

TEST(SamplingOwnership, SpectrumAndMomentsOwnTheirSamplingDependency)
{
  using namespace bethe::heisenberg;
  if (uni20::numeric_limits<long double>::digits() <= uni20::numeric_limits<double>::digits()) GTEST_SKIP();
  bethe::heisenberg::detail::SpinonJointProposal<long double> proposal(
      8.L, .5L, {{bethe::heisenberg::detail::SpinonProposalKind::background, 1.L, 1.L, 1.L}});
  SpinonRecoveryPolicy policy{
      "1e-5", "0", {{53, uni20::numeric_limits<long double>::digits(), {{{32, 12}, {64, 20}}}}}, 60};
  ThermodynamicSixSpinonStructureFactor<long double, double, long double> spectrum(proposal, policy, 1.L);
  ThermodynamicSixSpinonMoments<long double, double, long double> moments(proposal, policy, 1.L);
  auto point = spectrum.at(4 * std::atan(1.L) / 2, 2.L);
  check_owner([&](auto& scheduler) { return point.evaluate_async({2, 2, 42}, {&scheduler, 1, 2}); });
  check_owner([&](auto& scheduler) {
    return point.evaluate_observations_async({2, 2, 42},
                                             [](auto const&, std::stop_token) {
                                               return bethe::MonteCarloObservation<long double, SpinonSampleMetadata>{
                                                   std::vector{1.L, 0.L, 0.L}};
                                             },
                                             {&scheduler, 1, 2});
  });
  check_owner(
      [&](auto& scheduler) { return moments.fixed_moments(1.L).evaluate_async({2, 2, 42}, {&scheduler, 1, 2}); });
  check_owner([&](auto& scheduler) { return moments.global_async({2, 2, 42}, {&scheduler, 1, 2}); });
}
