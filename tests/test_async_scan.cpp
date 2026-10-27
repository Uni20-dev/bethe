// SPDX-License-Identifier: GPL-3.0-or-later
#include "batched_scan_reference.hpp"
#include "test_support.hpp"
#include <atomic>
#include <bethe/heisenberg_excitations.hpp>
#include <condition_variable>
#include <memory>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using namespace uni20::async;
struct State
{
    std::vector<uni20::half_int> quantum_numbers;
    double energy = 0;
    bool converged = true;
    std::shared_ptr<void> payload;
};
using Scan = bethe::RealExcitationScan<State>;
auto ground() { return State{{}, -1}; }
auto submit(bethe::ExecutionOptions execution, auto solve, std::size_t slots = 20, std::size_t keep = 3)
{
  return bethe::detail::scan_real_combinations_async<Scan>(slots, 1, 0, {keep, 100000, execution}, std::move(solve),
                                                           ground);
}
AsyncTask consume(ReadBuffer<Scan> input, WriteBuffer<std::size_t> output)
{
  auto const& scan = co_await input;
  co_await output = scan.converged_count;
}
} // namespace

TEST(AsyncScan, RefillsPastBlockedFirstCandidateAndRetainsFirstFailureByIdentity)
{
  TbbScheduler scheduler(2);
  std::mutex mutex;
  std::condition_variable ready;
  bool later = false, crossed = false;
  auto result = submit({&scheduler, 128, 2}, [&](auto labels) {
    auto id = labels[0].twice() / 2;
    if (id == 0)
    {
      std::unique_lock lock(mutex);
      crossed = ready.wait_for(lock, std::chrono::seconds(5), [&] { return later; });
    }
    if (id == 3)
    {
      std::lock_guard lock(mutex);
      later = true;
      ready.notify_all();
    }
    return State{labels, double(id), id != 0 && id != 2};
  });
  auto const& value = result.get_wait(scheduler);
  EXPECT_TRUE(crossed);
  EXPECT_EQ(value.converged_count, 18u);
  ASSERT_TRUE(value.first_unconverged);
  EXPECT_EQ(value.first_unconverged->quantum_numbers.front().twice(), 0);
}

TEST(AsyncScan, BoundedLivePayloadsForLargeDiscardingScan)
{
  struct Counts
  {
      std::atomic<int> live{0}, peak{0};
  } counts;
  struct Payload
  {
      Counts& counts;
      Payload(Counts& c) : counts(c)
      {
        auto n = ++counts.live, old = counts.peak.load();
        while (old < n && !counts.peak.compare_exchange_weak(old, n))
        {}
      }
      ~Payload() { --counts.live; }
  };
  TbbScheduler scheduler(4);
  std::vector<std::atomic<int>> visits(10000);
  {
    auto result = submit(
        {&scheduler, 128, 3},
        [&](auto labels) {
          auto id = labels[0].twice() / 2;
          ++visits[id];
          return State{labels, double(10000 - id), true, std::make_shared<Payload>(counts)};
        },
        10000, 5);
    auto const& value = result.get_wait(scheduler);
    ASSERT_EQ(value.levels.size(), 5u);
    EXPECT_EQ(value.converged_count, 10000u);
    EXPECT_LE(counts.peak, 8); // top 5 plus 3 active results, independent of 10000
    EXPECT_EQ(counts.live, 5);
    for (auto const& n : visits)
      EXPECT_EQ(n, 1);
  }
  scheduler.run_all(); // only flush completed frame destruction for lifetime assertion
  EXPECT_EQ(counts.live, 0);
}

TEST(AsyncScan, OwnsSubmissionStateAndComposesWithoutHostWait)
{
  for (auto order : {DebugSchedulerOrder::fifo, DebugSchedulerOrder::reverse, DebugSchedulerOrder::random})
  {
    DebugScheduler scheduler({.order = order, .random_seed = 41});
    auto owned = std::make_shared<int>(7);
    std::weak_ptr<int> weak = owned;
    Async<std::size_t> count;
    {
      auto result = submit({&scheduler, 128, 3}, [owned](auto labels) { return State{labels, double(*owned)}; });
      scheduler.schedule(consume(result.read(), count.write()));
      owned.reset();
    } // both submission-local result handle and original configuration disappear
    EXPECT_FALSE(weak.expired());
    EXPECT_EQ(count.get_wait(scheduler), 20u);
    scheduler.run_all();
    EXPECT_TRUE(weak.expired());
  }
}

TEST(AsyncScan, DiscardedOutputDoesNotCancelOwnedWork)
{
  DebugScheduler scheduler;
  int calls = 0;
  auto owned = std::make_shared<int>(9);
  std::weak_ptr<int> weak = owned;
  {
    auto ignored = submit({&scheduler, 128, 2}, [owned, &calls](auto labels) {
      ++calls;
      return State{labels, double(*owned)};
    });
  }
  owned.reset();
  EXPECT_FALSE(weak.expired());
  scheduler.run_all(); // no output exists; explicitly drive the borrowed scheduler
  EXPECT_EQ(calls, 20);
  EXPECT_TRUE(weak.expired());
}

TEST(AsyncScan, CancellationBeforePreparationAndDuringScan)
{
  DebugScheduler scheduler;
  for (bool before : {false, true})
  {
    std::stop_source stop;
    int calls = 0, preparations = 0;
    if (before) stop.request_stop();
    auto result = bethe::detail::scan_real_combinations_async<Scan>(
        20, 1, 0, {3, 100, {&scheduler, 128, 2, stop.get_token()}},
        [&](auto labels) {
          ++calls;
          stop.request_stop();
          return State{labels};
        },
        [&] {
          ++preparations;
          return ground();
        });
    EXPECT_THROW((void)result.get_wait(scheduler), bethe::ScanCancelled);
    EXPECT_EQ(calls, before ? 0 : 1);
    EXPECT_EQ(preparations, before ? 0 : 1);
  }
}

TEST(AsyncScan, CancellationAndExceptionJoinActiveSolves)
{
  TbbScheduler scheduler(2);
  for (bool cancel : {false, true})
  {
    std::stop_source stop;
    std::mutex mutex;
    std::condition_variable ready;
    bool sibling_started = false, failure = false;
    std::atomic<bool> sibling_done{false};
    std::atomic<int> active{0};
    auto result = submit({&scheduler, 128, 2, stop.get_token()}, [&](auto labels) {
      struct Guard
      {
          std::atomic<int>& a;
          Guard(std::atomic<int>& a) : a(a) { ++a; }
          ~Guard() { --a; }
      } guard(active);
      auto id = labels[0].twice() / 2;
      std::unique_lock lock(mutex);
      if (id == 0)
      {
        EXPECT_TRUE(ready.wait_for(lock, std::chrono::seconds(5), [&] { return sibling_started; }));
        failure = true;
        ready.notify_all();
        if (!cancel) throw std::bad_alloc{};
        stop.request_stop();
      }
      else if (id == 1)
      {
        sibling_started = true;
        ready.notify_all();
        EXPECT_TRUE(ready.wait_for(lock, std::chrono::seconds(5), [&] { return failure; }));
        sibling_done = true;
      }
      return State{labels};
    });
    if (cancel)
      EXPECT_THROW((void)result.get_wait(scheduler), bethe::ScanCancelled);
    else
      EXPECT_THROW((void)result.get_wait(scheduler), std::bad_alloc);
    EXPECT_TRUE(sibling_done);
    EXPECT_EQ(active, 0);
  }
}

TEST(AsyncScan, PreflightAndGroundFailures)
{
  DebugScheduler scheduler;
  int calls = 0;
  auto solve = [&](auto labels) {
    ++calls;
    return State{labels};
  };
  EXPECT_THROW(submit({&scheduler, 128, 0}, solve), std::invalid_argument);
  EXPECT_THROW(bethe::detail::scan_real_combinations_async<Scan>(20, 2, 0, {3, 10, {&scheduler}}, solve, ground),
               std::length_error);
  auto result = bethe::detail::scan_real_combinations_async<Scan>(20, 1, 0, {3, 100, {&scheduler}}, solve,
                                                                  []() -> State { throw std::bad_alloc{}; });
  EXPECT_THROW((void)result.get_wait(scheduler), std::bad_alloc);
  EXPECT_EQ(calls, 0);
}

TEST(AsyncScan, InitialAdmissionFailureDoesNotLeaveWork)
{
  struct Reject : DebugScheduler
  {
      void schedule(AsyncTask&&) override { throw std::bad_alloc{}; }
  } scheduler;
  auto owner = std::make_shared<int>(0);
  std::weak_ptr<int> weak = owner;
  EXPECT_THROW(submit({&scheduler},
                      [owner](auto labels) {
                        ++*owner;
                        return State{labels};
                      }),
               std::bad_alloc);
  EXPECT_EQ(*owner, 0);
  owner.reset();
  EXPECT_TRUE(weak.expired());
}

TEST(AsyncScan, RunnerAdmissionFailureJoinsAcceptedWork)
{
  struct RejectRunner : DebugScheduler
  {
      std::atomic<bool> finished{false};
      TbbScheduler inner{2};
      void execute_batch_impl(LightweightTaskBatch const& batch) override
      {
        // Model a batch implementation that fails after another runner starts.
        // Its contract must join accepted callbacks before reporting failure.
        std::mutex mutex;
        std::condition_variable ready;
        bool started = false;
        inner.execute_batch(2, [&](std::size_t i) {
          if (i == 0)
          {
            {
              std::lock_guard lock(mutex);
              started = true;
              ready.notify_all();
            }
            batch(0);
            finished = true;
          }
          else
          {
            std::unique_lock lock(mutex);
            EXPECT_TRUE(ready.wait_for(lock, std::chrono::seconds(5), [&] { return started; }));
            throw std::bad_alloc{};
          }
        });
      }
  } scheduler;
  auto owner = std::make_shared<int>(1);
  std::weak_ptr<int> weak = owner;
  auto result = submit({&scheduler}, [owner](auto labels) { return State{labels, double(*owner)}; });
  owner.reset();
  EXPECT_THROW((void)result.get_wait(scheduler), std::bad_alloc);
  EXPECT_TRUE(scheduler.finished);
  scheduler.run_all();
  EXPECT_TRUE(weak.expired());
}

TEST(AsyncScan, UnconvergedGroundRemainsNumericalOutcomeWithAbsentGaps)
{
  DebugScheduler scheduler;
  auto result = bethe::detail::scan_real_combinations_async<Scan>(
      5, 1, 0, {5, 100, {&scheduler}}, [](auto labels) { return State{labels}; },
      [] { return State{{uni20::half_int{0}}, -1, false}; });
  auto const& value = result.get_wait(scheduler);
  EXPECT_FALSE(value.converged());
  EXPECT_FALSE(value.family_converged());
  EXPECT_EQ(value.converged_count, 4u); // failed reference reused, not solved a second time
  ASSERT_TRUE(value.first_unconverged);
  for (auto const& level : value.levels)
    EXPECT_FALSE(level.gap.has_value());
}

TEST(AsyncScan, HeapAndProjectedEnergyMatchFrozenScanner)
{
  DebugScheduler serial;
  TbbScheduler parallel(4);
  for (auto order : {bethe::detail::EnergyOrder::ascending, bethe::detail::EnergyOrder::descending})
    for (std::size_t keep : {1, 5, 10000})
    {
      auto solve = [](auto labels) { return State{labels, 1e20, labels[1].twice() % 7 != 0}; };
      auto project = [](State const& s) {
        return s.quantum_numbers.empty() ? -1. : double(s.quantum_numbers[0].twice() % 5);
      };
      auto a = scan_reference::scan_real_combinations<Scan>(40, 3, 0, {keep, 10000, {&serial}}, solve, ground, order,
                                                            project);
      auto b = bethe::detail::scan_real_combinations<Scan>(40, 3, 0, {keep, 10000, {&parallel, 128, 3}}, solve, ground,
                                                           order, project);
      EXPECT_EQ(a.candidate_count, 9880u);
      EXPECT_EQ(a.converged_count, b.converged_count);
      ASSERT_TRUE(a.first_unconverged && b.first_unconverged);
      EXPECT_EQ(a.first_unconverged->quantum_numbers, b.first_unconverged->quantum_numbers);
      ASSERT_EQ(a.levels.size(), b.levels.size());
      for (std::size_t i = 0; i < a.levels.size(); ++i)
      {
        EXPECT_EQ(a.levels[i].state.quantum_numbers, b.levels[i].state.quantum_numbers);
        EXPECT_EQ(a.levels[i].gap, b.levels[i].gap);
      }
    }
}

template <typename T> class AsyncScanPrecision : public ::testing::Test {};
TYPED_TEST_SUITE(AsyncScanPrecision, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(AsyncScanPrecision, PublicXXXMatchesFrozenScanner)
{
  using R = TypeParam;
  DebugScheduler serial;
  TbbScheduler parallel(4);
  constexpr std::size_t n = 32, m = n / 2 - 1, slots = n - m;
  auto check = [&]<bool open>() {
    using S = std::conditional_t<open, bethe::heisenberg::open::RealState<R>, bethe::heisenberg::RealState<R>>;
    auto solve = [n](auto const& labels) {
      if constexpr (open)
        return bethe::heisenberg::open::solve_real<R>(n, labels);
      else
        return bethe::heisenberg::solve_real<R>(n, labels);
    };
    auto reference = [n] {
      if constexpr (open)
        return bethe::heisenberg::open::ground_state<R>(n);
      else
        return bethe::heisenberg::ground_state<R>(n);
    };
    for (std::size_t keep : {5, 1000})
    {
      auto a = scan_reference::scan_real_combinations<bethe::RealExcitationScan<S>>(
          slots, m, open ? 2 : -std::int64_t(slots - 1), {keep, 1000, {&serial}}, solve, reference);
      auto pending = [&] {
        if constexpr (open)
          return bethe::heisenberg::open::real_excitations_async<R>(n, uni20::half_int{1},
                                                                    {keep, 1000, {&parallel, 128, 3}});
        else
          return bethe::heisenberg::real_excitations_async<R>(n, uni20::half_int{1}, {keep, 1000, {&parallel, 128, 3}});
      }();
      auto const& b = pending.get_wait(parallel);
      EXPECT_EQ(b.spin, uni20::half_int{1});
      EXPECT_EQ(a.converged_count, b.converged_count);
      EXPECT_EQ(a.candidate_count, b.candidate_count);
      ASSERT_EQ(a.levels.size(), b.levels.size());
      for (std::size_t i = 0; i < a.levels.size(); ++i)
      {
        auto const& x = a.levels[i].state;
        auto const& y = b.levels[i].state;
        EXPECT_EQ(x.quantum_numbers, y.quantum_numbers);
        EXPECT_EQ(x.rapidities, y.rapidities);
        EXPECT_EQ(x.energy, y.energy);
        EXPECT_EQ(x.residual_norm, y.residual_norm);
        EXPECT_EQ(x.iterations, y.iterations);
        EXPECT_EQ(a.levels[i].gap, b.levels[i].gap);
      }
    }
  };
  check.template operator()<false>();
  check.template operator()<true>();
}
