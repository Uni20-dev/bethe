// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include "two_stage_dsf_reference.hpp"
#include <atomic>
#include <condition_variable>
#include <memory>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
namespace model = bethe::heisenberg;
using namespace uni20::async;
using State = model::RealState<double>;
using Factor = model::FormFactor<double>;
using Result = model::TwoSpinonStructureFactor<double>;

struct Fixture
{
    std::vector<model::QuantumNumbers> labels;
    Fixture()
    {
      model::QuantumNumbers q{uni20::half_int{-2}, uni20::half_int{-1}, uni20::half_int{0}};
      do
      {
        labels.push_back(q);
      }
      while (bethe::detail::advance_combination<uni20::half_int>(q, uni20::half_int{2}));
    }
    std::size_t id(model::QuantumNumbers const& q) const
    {
      return std::lower_bound(labels.begin(), labels.end(), q) - labels.begin();
    }
    State solve(model::QuantumNumbers q) const
    {
      auto const i = id(q);
      State result;
      result.rapidities = {double(i)}; // synthetic identity; never passed to a numerical kernel
      result.quantum_numbers = std::move(q);
      result.energy = double(i % 3); // ties, different from enumeration order
      result.momentum_index = i % 8;
      result.converged = true;
      return result;
    }
};
State ground()
{
  State s;
  s.quantum_numbers = {uni20::half_int{99}};
  s.energy = -100;
  s.converged = true;
  return s;
}
Factor factor(State const&, State const& state)
{
  return {.weight = state.rapidities[0] + 1,
          .log_weight = state.rapidities[0],
          .pivot_margin = state.rapidities[0] + .5,
          .status = model::FormFactorStatus::converged};
}
// Adapt the synthetic two-state kernel to the owning prepare-once contract.
template <typename Solve, typename Ground, typename Weight>
auto launch(std::size_t n, std::size_t budget, bethe::ExecutionOptions execution, Solve solve, Ground ground_state,
            Weight weight)
{
  return model::detail::submit_two_spinon_structure_factor<double>(
      n, budget, execution, std::move(solve), std::move(ground_state),
      [weight = std::move(weight)](State const& reference) {
        return [weight, reference](State const& state) { return weight(reference, state); };
      },
      0.);
}
template <typename Solve, typename Weight = decltype(&factor)>
auto submit(bethe::ExecutionOptions execution, Solve solve, Weight weight = factor)
{
  return launch(8, 100, execution, std::move(solve), ground, std::move(weight));
}
AsyncTask consume(ReadBuffer<Result> input, WriteBuffer<std::size_t> output)
{
  auto const& result = co_await input;
  co_await output = result.lines.size();
}
} // namespace

TEST(AsyncStructureFactor, LaterFactorsPassBlockedRootAndBlockedFactor)
{
  Fixture fixture;
  for (bool block_factor : {false, true})
  {
    TbbScheduler scheduler(2);
    std::mutex mutex;
    std::condition_variable ready;
    bool later = false, crossed = false;
    auto wait = [&] {
      std::unique_lock lock(mutex);
      crossed = ready.wait_for(lock, std::chrono::seconds(5), [&] { return later; });
    };
    auto pending = submit(
        {&scheduler, 128, 2},
        [&](auto q) {
          if (!block_factor && fixture.id(q) == 0) wait();
          return fixture.solve(std::move(q));
        },
        [&](auto const& g, auto const& s) {
          auto const id = fixture.id(s.quantum_numbers);
          if (block_factor && id == 0) wait();
          if (id == 3)
          {
            std::lock_guard lock(mutex);
            later = true;
            ready.notify_all();
          }
          return factor(g, s);
        });
    auto const& result = pending.get_wait(scheduler);
    EXPECT_TRUE(crossed);
    ASSERT_TRUE(result.converged());
    for (std::size_t i = 0; i < result.scan.levels.size(); ++i)
    {
      auto id = fixture.id(result.scan.levels[i].state.quantum_numbers);
      EXPECT_EQ(result.form_factors[i].weight, double(id + 1));
      EXPECT_EQ(result.lines[i].state_id, i + 1);
      EXPECT_EQ(result.lines[i].weight, double(id + 1));
    }
  }
}

TEST(AsyncStructureFactor, CapacityHeldThroughBlockedFormFactor)
{
  Fixture fixture;
  TbbScheduler scheduler(4);
  std::mutex mutex;
  std::condition_variable ready;
  bool release = false;
  int blocked = 0;
  std::atomic<int> solves{0}, live{0}, peak{0};
  std::vector<std::atomic<int>> visits(10);
  auto pending = submit(
      {&scheduler, 128, 2},
      [&](auto q) {
        ++solves;
        ++visits[fixture.id(q)];
        auto n = ++live, old = peak.load();
        while (old < n && !peak.compare_exchange_weak(old, n))
        {}
        return fixture.solve(std::move(q));
      },
      [&](auto const& g, auto const& s) {
        if (fixture.id(s.quantum_numbers) < 2)
        {
          std::unique_lock lock(mutex);
          ++blocked;
          ready.notify_all();
          EXPECT_TRUE(ready.wait_for(lock, std::chrono::seconds(5), [&] { return release; }));
        }
        --live;
        return factor(g, s);
      });
  {
    std::unique_lock lock(mutex);
    EXPECT_TRUE(ready.wait_for(lock, std::chrono::seconds(5), [&] { return blocked == 2; }));
    EXPECT_EQ(solves, 2); // no solved-root queue can accumulate behind the factors
    release = true;
    ready.notify_all();
  }
  EXPECT_TRUE(pending.get_wait(scheduler).converged());
  EXPECT_EQ(peak, 2);
  EXPECT_EQ(live, 0);
  for (auto const& count : visits)
    EXPECT_EQ(count, 1);
}

TEST(AsyncStructureFactor, FailuresStayAlignedAndNeverEnterFactorKernel)
{
  Fixture fixture;
  DebugScheduler scheduler;
  int calls = 0;
  auto pending = submit(
      {&scheduler},
      [&](auto q) {
        auto s = fixture.solve(std::move(q));
        s.converged = fixture.id(s.quantum_numbers) != 2;
        if (fixture.id(s.quantum_numbers) == 1) s.energy = -101; // invalid negative gap
        return s;
      },
      [&](auto const& g, auto const& s) {
        EXPECT_TRUE(g.converged && s.converged);
        ++calls;
        return fixture.id(s.quantum_numbers) == 3 ? Factor{} : factor(g, s);
      });
  auto const& result = pending.get_wait(scheduler);
  EXPECT_FALSE(result.converged());
  EXPECT_EQ(calls, 9);
  ASSERT_EQ(result.scan.levels.size(), 9u);
  ASSERT_EQ(result.form_factors.size(), 9u);
  ASSERT_TRUE(result.scan.first_unconverged);
  EXPECT_EQ(fixture.id(result.scan.first_unconverged->quantum_numbers), 2u);
  EXPECT_EQ(result.lines.size(), 7u);
  for (std::size_t i = 0; i < result.form_factors.size(); ++i)
  {
    auto id = fixture.id(result.scan.levels[i].state.quantum_numbers);
    if (id == 1 || id == 3)
    {
      EXPECT_EQ(result.form_factors[i].status, model::FormFactorStatus::precision_limit);
      EXPECT_FALSE(result.form_factors[i].weight);
    }
    else
      EXPECT_EQ(result.form_factors[i].weight, double(id + 1));
  }
  for (auto const& line : result.lines)
  {
    auto id = fixture.id(result.scan.levels[line.state_id - 1].state.quantum_numbers);
    EXPECT_EQ(line.weight, double(id + 1));
  }
  auto failed_ground = launch(
      8, 100, {&scheduler}, [&](auto q) { return fixture.solve(std::move(q)); },
      [] {
        auto s = ground();
        s.converged = false;
        return s;
      },
      [&](auto const&, auto const&) {
        ++calls;
        return Factor{};
      });
  auto const& failed = failed_ground.get_wait(scheduler);
  EXPECT_EQ(calls, 9);
  EXPECT_EQ(failed.scan.converged_count, 10u);
  EXPECT_TRUE(failed.lines.empty());
  EXPECT_FALSE(failed.first_moment_fraction);
  for (auto const& f : failed.form_factors)
    EXPECT_EQ(f.status, model::FormFactorStatus::roots_unconverged);
}

TEST(AsyncStructureFactor, CancellationAtStageBoundaries)
{
  Fixture fixture;
  DebugScheduler scheduler;
  for (int stage : {0, 1, 2})
  {
    std::stop_source stop;
    int solves = 0, factors = 0, preparations = 0;
    if (stage == 0) stop.request_stop();
    auto pending = launch(
        8, 100, {&scheduler, 128, 2, stop.get_token()},
        [&](auto q) {
          ++solves;
          if (stage == 1) stop.request_stop();
          return fixture.solve(std::move(q));
        },
        [&] {
          ++preparations;
          return ground();
        },
        [&](auto const& g, auto const& s) {
          ++factors;
          stop.request_stop();
          return factor(g, s);
        });
    EXPECT_THROW((void)pending.get_wait(scheduler), bethe::ScanCancelled);
    EXPECT_EQ(preparations, stage == 0 ? 0 : 1);
    EXPECT_EQ(solves, stage == 0 ? 0 : 1);
    EXPECT_EQ(factors, stage == 2 ? 1 : 0);
  }
}

TEST(AsyncStructureFactor, FactorExceptionOrCancellationJoinsActiveSibling)
{
  Fixture fixture;
  for (bool cancel : {false, true})
  {
    TbbScheduler scheduler(2);
    std::stop_source stop;
    std::mutex mutex;
    std::condition_variable ready;
    bool sibling = false, failure = false;
    std::atomic<bool> done = false;
    std::atomic<int> active = 0;
    auto pending = submit(
        {&scheduler, 128, 2, stop.get_token()}, [&](auto q) { return fixture.solve(std::move(q)); },
        [&](auto const& g, auto const& s) {
          struct Guard
          {
              std::atomic<int>& n;
              Guard(std::atomic<int>& n) : n(n) { ++n; }
              ~Guard() { --n; }
          } guard(active);
          auto id = fixture.id(s.quantum_numbers);
          std::unique_lock lock(mutex);
          if (id == 0)
          {
            EXPECT_TRUE(ready.wait_for(lock, std::chrono::seconds(5), [&] { return sibling; }));
            failure = true;
            ready.notify_all();
            if (!cancel) throw std::bad_alloc{};
            stop.request_stop();
          }
          else if (id == 1)
          {
            sibling = true;
            ready.notify_all();
            EXPECT_TRUE(ready.wait_for(lock, std::chrono::seconds(5), [&] { return failure; }));
            done = true;
          }
          return factor(g, s);
        });
    if (cancel)
      EXPECT_THROW((void)pending.get_wait(scheduler), bethe::ScanCancelled);
    else
      EXPECT_THROW((void)pending.get_wait(scheduler), std::bad_alloc);
    EXPECT_TRUE(done);
    EXPECT_EQ(active, 0);
  }
}

TEST(AsyncStructureFactor, OwnsContextAndComposesAfterSubmissionHandleDisappears)
{
  Fixture fixture;
  for (auto order : {DebugSchedulerOrder::fifo, DebugSchedulerOrder::reverse, DebugSchedulerOrder::random})
  {
    DebugScheduler scheduler({.order = order, .random_seed = 31});
    auto owner = std::make_shared<int>(0);
    std::weak_ptr<int> weak = owner;
    Async<std::size_t> count;
    {
      auto pending = submit({&scheduler}, [fixture, owner](auto q) {
        ++*owner;
        return fixture.solve(std::move(q));
      });
      scheduler.schedule(consume(pending.read(), count.write()));
    }
    owner.reset();
    EXPECT_FALSE(weak.expired());
    EXPECT_EQ(count.get_wait(scheduler), 10u);
    scheduler.run_all(); // frame destruction only, not production completion contract
    EXPECT_TRUE(weak.expired());
  }
}

TEST(AsyncStructureFactor, DiscardedHandleDoesNotCancel)
{
  Fixture fixture;
  DebugScheduler scheduler;
  auto owner = std::make_shared<int>(0);
  std::weak_ptr<int> weak = owner;
  int factors = 0;
  {
    auto ignored = submit(
        {&scheduler},
        [fixture, owner](auto q) {
          ++*owner;
          return fixture.solve(std::move(q));
        },
        [&](auto const& g, auto const& s) {
          ++factors;
          return factor(g, s);
        });
  }
  owner.reset();
  EXPECT_FALSE(weak.expired());
  scheduler.run_all(); // explicitly drive work after voluntarily discarding its handle
  EXPECT_EQ(factors, 10);
  EXPECT_TRUE(weak.expired());
}

TEST(AsyncStructureFactor, PreflightPreparationAndInitialAdmissionFailures)
{
  DebugScheduler scheduler;
  int calls = 0;
  auto solve = [&](auto) {
    ++calls;
    return State{};
  };
  auto prepare = [&] {
    ++calls;
    return ground();
  };
  auto request = [&](std::size_t n, std::size_t budget, bethe::ExecutionOptions execution) {
    return launch(n, budget, execution, solve, prepare, factor);
  };
  EXPECT_THROW(request(5, 100, {&scheduler}), std::invalid_argument);
  EXPECT_THROW(request(8, 9, {&scheduler}), std::length_error);
  EXPECT_THROW(request(8, 100, {&scheduler, 128, 0}), std::invalid_argument);
  EXPECT_EQ(calls, 0);
  struct Reject : DebugScheduler
  {
      void schedule(AsyncTask&&) override { throw std::bad_alloc{}; }
  } reject;
  EXPECT_THROW(request(8, 100, {&reject}), std::bad_alloc);
  EXPECT_EQ(calls, 0);
  auto pending = launch(8, 100, {&scheduler}, solve, []() -> State { throw std::bad_alloc{}; }, factor);
  EXPECT_THROW((void)pending.get_wait(scheduler), std::bad_alloc);
  EXPECT_EQ(calls, 0);
}

TEST(AsyncStructureFactor, ActiveSchedulerIsResolvedAtSubmissionAndOneRunnerBatchIsUsed)
{
  struct Recording : DebugScheduler
  {
      std::vector<std::size_t> batches;
      void execute_batch_impl(LightweightTaskBatch const& batch) override
      {
        batches.push_back(batch.size());
        for (std::size_t i = 0; i < batch.size(); ++i)
          batch(i);
      }
  } first, other;
  Fixture fixture;
  auto pending = [&] {
    ScopedScheduler scope(&first);
    return submit({nullptr, 1, 3}, [fixture](auto q) { return fixture.solve(std::move(q)); });
  }();
  ScopedScheduler scope(&other);
  EXPECT_TRUE(pending.get_wait(first).converged());
  EXPECT_EQ(first.batches, (std::vector<std::size_t>{3}));
  EXPECT_TRUE(other.batches.empty());
  EXPECT_EQ(get_global_scheduler(), &other);
}

TEST(AsyncStructureFactor, FinalReductionFailureReachesObserver)
{
  Fixture fixture;
  DebugScheduler scheduler;
  int factors = 0;
  auto pending = submit(
      {&scheduler}, [fixture](auto q) { return fixture.solve(std::move(q)); },
      [&](auto const& g, auto const& s) {
        ++factors;
        auto f = factor(g, s);
        f.weight = std::numeric_limits<double>::quiet_NaN();
        return f;
      });
  EXPECT_THROW((void)pending.get_wait(scheduler), std::invalid_argument);
  EXPECT_EQ(factors, 10); // all owned work has completed before final reduction
}

TEST(AsyncStructureFactor, PreparesOneOwnedSnapshotBeforeAnyCandidate)
{
  Fixture fixture;
  TbbScheduler scheduler(4);
  std::atomic<int> preparations{0}, factors{0};
  std::weak_ptr<int> lifetime;
  auto pending = model::detail::submit_two_spinon_structure_factor<double>(
      8, 100, {&scheduler, 128, 3},
      [&](auto q) {
        EXPECT_EQ(preparations, 1);
        return fixture.solve(std::move(q));
      },
      ground,
      [&](State const& reference) {
        ++preparations;
        auto owned = std::make_shared<int>(42);
        lifetime = owned;
        return [owned, reference, &factors, workspace = std::make_unique<int>(7)](State const& state) {
          EXPECT_EQ(*owned, 42);
          EXPECT_EQ(*workspace, 7);
          ++factors;
          return factor(reference, state);
        };
      },
      0.);
  ASSERT_TRUE(pending.get_wait(scheduler).converged());
  EXPECT_EQ(preparations, 1);
  EXPECT_EQ(factors, 10);
  scheduler.run_all(); // complete frame destruction before checking ownership
  EXPECT_TRUE(lifetime.expired());
}

TEST(AsyncStructureFactor, SnapshotPreparationFailureOrCancellationPreventsSolves)
{
  DebugScheduler scheduler;
  for (bool cancel : {false, true})
  {
    std::stop_source stop;
    int solves = 0;
    auto pending = model::detail::submit_two_spinon_structure_factor<double>(
        8, 100, {&scheduler, 128, 3, stop.get_token()},
        [&](auto) {
          ++solves;
          return State{};
        },
        ground,
        [&](State const&) {
          if (!cancel) throw std::runtime_error("snapshot preparation failed");
          stop.request_stop();
          return [](State const&) { return Factor{}; };
        },
        0.);
    if (cancel)
      EXPECT_THROW((void)pending.get_wait(scheduler), bethe::ScanCancelled);
    else
      EXPECT_THROW((void)pending.get_wait(scheduler), std::runtime_error);
    EXPECT_EQ(solves, 0);
  }
}

namespace
{
template <typename R>
void equal(model::TwoSpinonStructureFactor<R> const& a, model::TwoSpinonStructureFactor<R> const& b)
{
  auto state = [](auto const& x, auto const& y) {
    EXPECT_EQ(x.rapidities, y.rapidities);
    EXPECT_EQ(x.quantum_numbers, y.quantum_numbers);
    EXPECT_EQ(x.energy, y.energy);
    EXPECT_EQ(x.residual_norm, y.residual_norm);
    EXPECT_EQ(x.iterations, y.iterations);
    EXPECT_EQ(x.converged, y.converged);
    EXPECT_EQ(x.sz, y.sz);
    EXPECT_EQ(x.momentum_index, y.momentum_index);
    EXPECT_EQ(x.momentum, y.momentum);
    EXPECT_EQ(x.spin_reversed, y.spin_reversed);
  };
  EXPECT_EQ(a.scan.spin, b.scan.spin);
  EXPECT_EQ(a.scan.candidate_count, b.scan.candidate_count);
  EXPECT_EQ(a.scan.converged_count, b.scan.converged_count);
  EXPECT_EQ(a.converged(), b.converged());
  state(a.scan.ground_state, b.scan.ground_state);
  ASSERT_EQ(a.scan.first_unconverged.has_value(), b.scan.first_unconverged.has_value());
  if (a.scan.first_unconverged) state(*a.scan.first_unconverged, *b.scan.first_unconverged);
  ASSERT_EQ(a.scan.levels.size(), b.scan.levels.size());
  ASSERT_EQ(a.form_factors.size(), b.form_factors.size());
  for (std::size_t i = 0; i < a.scan.levels.size(); ++i)
  {
    state(a.scan.levels[i].state, b.scan.levels[i].state);
    EXPECT_EQ(a.scan.levels[i].gap, b.scan.levels[i].gap);
    EXPECT_EQ(a.form_factors[i].status, b.form_factors[i].status);
    EXPECT_EQ(a.form_factors[i].weight, b.form_factors[i].weight);
    EXPECT_EQ(a.form_factors[i].log_weight, b.form_factors[i].log_weight);
    EXPECT_EQ(a.form_factors[i].pivot_margin, b.form_factors[i].pivot_margin);
  }
  ASSERT_EQ(a.lines.size(), b.lines.size());
  for (std::size_t i = 0; i < a.lines.size(); ++i)
  {
    EXPECT_EQ(a.lines[i].state_id, b.lines[i].state_id);
    EXPECT_EQ(a.lines[i].momentum_index, b.lines[i].momentum_index);
    EXPECT_EQ(a.lines[i].momentum, b.lines[i].momentum);
    EXPECT_EQ(a.lines[i].gap, b.lines[i].gap);
    EXPECT_EQ(a.lines[i].weight, b.lines[i].weight);
  }
  ASSERT_EQ(a.moments.size(), b.moments.size());
  for (std::size_t i = 0; i < a.moments.size(); ++i)
  {
    EXPECT_EQ(a.moments[i].weight, b.moments[i].weight);
    EXPECT_EQ(a.moments[i].first_moment, b.moments[i].first_moment);
  }
  EXPECT_EQ(a.integrated_weight, b.integrated_weight);
  EXPECT_EQ(a.weight_fraction, b.weight_fraction);
  EXPECT_EQ(a.first_moment_fraction, b.first_moment_fraction);
}
} // namespace

template <typename T> class AsyncStructureFactorPrecision : public ::testing::Test {};
TYPED_TEST_SUITE(AsyncStructureFactorPrecision, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(AsyncStructureFactorPrecision, CompleteResultsExactlyMatchTwoStagePath)
{
  using R = TypeParam;
  DebugScheduler serial;
  TbbScheduler parallel(4);
  for (std::size_t n : {2, 6, 16, 32})
  {
    auto old = dsf_reference::two_spinon_structure_factor<R>(n, {}, 1000, {&serial, 7});
    ASSERT_TRUE(old.converged());
    for (std::size_t window : {1, 3, 128})
    {
      auto pending = model::two_spinon_structure_factor_async<R>(n, {}, 1000, {&parallel, 7, window});
      equal(old, pending.get_wait(parallel));
    }
    equal(old, model::two_spinon_structure_factor<R>(n, {}, 1000, {&serial, 1, 2}));
  }
}
TYPED_TEST(AsyncStructureFactorPrecision, NumericalFailuresExactlyMatchTwoStagePath)
{
  using R = TypeParam;
  DebugScheduler serial;
  TbbScheduler parallel(4);
  for (auto options :
       {model::SolverOptions<R>{.max_iterations = 0}, model::SolverOptions<R>{.residual_tolerance = R{.1}}})
  {
    auto old = dsf_reference::two_spinon_structure_factor<R>(8, options, 100, {&serial, 3});
    auto pending = model::two_spinon_structure_factor_async<R>(8, options, 100, {&parallel, 3, 2});
    EXPECT_FALSE(old.converged());
    equal(old, pending.get_wait(parallel));
  }
}
