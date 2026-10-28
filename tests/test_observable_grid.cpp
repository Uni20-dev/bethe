// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include <array>
#include <atomic>
#include <bethe/hubbard_doped.hpp>
#include <bethe/lieb_liniger_thermo.hpp>
#include <bethe/observable_grid.hpp>
#include <condition_variable>
#include <memory>
#include <thread>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using namespace uni20::async;
template <typename E, typename C>
void run(bool async, std::size_t count, E evaluate, C collect, bethe::ExecutionOptions options)
{
  if (async)
  {
    auto pending = bethe::evaluate_observable_grid_async(count, std::move(evaluate), std::move(collect), options);
    EXPECT_EQ(pending.get_wait(*options.scheduler), count);
  }
  else
    bethe::evaluate_observable_grid(count, evaluate, collect, options);
}

struct Counts
{
    std::atomic<int> live{0}, peak{0};
};
struct Payload
{
    Counts& counts;
    std::size_t id;
    Payload(Counts& c, std::size_t i) : counts(c), id(i)
    {
      auto n = ++counts.live, old = counts.peak.load();
      while (old < n && !counts.peak.compare_exchange_weak(old, n))
      {}
    }
    ~Payload() { --counts.live; }
};
AsyncTask downstream(ReadBuffer<std::size_t> input, WriteBuffer<std::size_t> out)
{
  auto n = co_await input;
  co_await out = n + 1;
}
} // namespace

TEST(ObservableGrid, OrderedBoundedMoveOnlyResultsAndCallerThreadCollection)
{
  TbbScheduler scheduler(4);
  auto const caller = std::this_thread::get_id();
  for (bool async : {false, true})
  {
    Counts counts;
    std::size_t next = 0;
    std::array<std::atomic<int>, 1000> visits{};
    run(
        async, visits.size(),
        [&](std::size_t i) {
          ++visits[i];
          return std::make_unique<Payload>(counts, i);
        },
        [&](std::size_t i, auto value) {
          if (!async) EXPECT_EQ(std::this_thread::get_id(), caller);
          EXPECT_EQ(i, next++);
          EXPECT_EQ(value->id, i);
        },
        {&scheduler, 128, 3});
    EXPECT_EQ(next, visits.size());
    EXPECT_LE(counts.peak, 3);
    EXPECT_EQ(counts.live, 0);
    for (auto const& n : visits)
      EXPECT_EQ(n, 1);
  }
}

TEST(ObservableGrid, ContinuousAdmissionPastSlowEndOfWindow)
{
  TbbScheduler scheduler(2);
  // The async coordinator can progress on either scheduler participant. A
  // synchronous caller doing get_wait may itself execute the artificial gated
  // kernel, so demanding that same caller collect concurrently would be invalid.
  for (bool async : {true})
  {
    std::mutex mutex;
    std::condition_variable changed;
    bool later = false, crossed = false;
    std::size_t next = 0;
    run(
        async, 12,
        [&](std::size_t i) {
          if (i == 1)
          {
            std::unique_lock lock(mutex);
            crossed = changed.wait_for(lock, std::chrono::seconds(5), [&] { return later; });
          }
          if (i == 2)
          {
            std::lock_guard lock(mutex);
            later = true;
            changed.notify_all();
          }
          return i;
        },
        [&](std::size_t i, auto value) {
          EXPECT_EQ(i, next++);
          EXPECT_EQ(value, i);
        },
        {&scheduler, 2, 2});
    EXPECT_TRUE(crossed); // A join of [0,1] before admitting 2 cannot pass.
    EXPECT_EQ(next, 12u);
  }
}

TEST(ObservableGrid, SynchronousRefillsBeforePreviousWindowCompletes)
{
  struct Recording : DebugScheduler
  {
      Recording() : DebugScheduler({.order = DebugSchedulerOrder::fifo}) {}
      std::size_t admitted = 0;
      WaitPredicate const* waiting = nullptr;
      void help_while_waiting(WaitPredicate const& ready) override
      {
        // DebugScheduler normally drains its ready wave. Stop between tasks
        // when the awaited value is ready, to test the adapter's own admission
        // policy independently of that scheduler's coarser stepping policy.
        struct Restore
        {
            WaitPredicate const*& slot;
            WaitPredicate const* old;
            ~Restore() { slot = old; }
        } guard{waiting, waiting};
        waiting = &ready;
        DebugScheduler::help_while_waiting(ready);
      }
      void resume_task(BasicTask& task) override
      {
        if (waiting && (*waiting)())
          enqueue_task(std::move(task));
        else
          DebugScheduler::resume_task(task);
      }
      void schedule(AsyncTask&& task) override
      {
        ++admitted;
        DebugScheduler::schedule(std::move(task));
      }
  } scheduler;
  std::size_t collected = 0;
  bethe::evaluate_observable_grid(
      12,
      [&](std::size_t i) {
        // FIFO executes point 0 first. Its result is collected and point 2 is
        // admitted BEFORE point 1 runs, which a whole-window join cannot do.
        if (i == 1)
        {
          EXPECT_GE(scheduler.admitted, 3u);
          EXPECT_EQ(collected, 1u);
        }
        return i;
      },
      [&](auto i, auto value) {
        EXPECT_EQ(i, collected++);
        EXPECT_EQ(i, value);
      },
      {&scheduler, 2, 2});
  EXPECT_EQ(collected, 12u);
}

TEST(ObservableGrid, SlowCollectorRetainsCreditAndBoundsAdmission)
{
  TbbScheduler scheduler(4);
  for (bool async : {false, true})
  {
    std::atomic<int> started{0};
    Counts counts;
    run(
        async, 80,
        [&](std::size_t i) {
          ++started;
          return std::make_unique<Payload>(counts, i);
        },
        [&](std::size_t i, auto value) {
          if (i == 0)
          {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            EXPECT_LE(started, 3);
            EXPECT_LE(counts.peak, 3);
          }
          EXPECT_EQ(value->id, i);
        },
        {&scheduler, 128, 3});
    EXPECT_EQ(started, 80);
    EXPECT_EQ(counts.live, 0);
  }
}

TEST(ObservableGrid, SerialSchedulerOrdersEmptyAndFailedNumericalValues)
{
  for (auto order : {DebugSchedulerOrder::fifo, DebugSchedulerOrder::reverse, DebugSchedulerOrder::random})
    for (bool async : {false, true})
    {
      DebugScheduler scheduler({.order = order, .random_seed = 42});
      for (std::size_t count : {0, 1, 17})
      {
        std::size_t next = 0;
        run(
            async, count, [](std::size_t i) { return i % 2 ? std::optional<int>{} : std::optional<int>(int(i)); },
            [&](std::size_t i, auto value) {
              EXPECT_EQ(i, next++);
              EXPECT_EQ(value.has_value(), i % 2 == 0);
              if (value) EXPECT_EQ(*value, int(i));
            },
            {&scheduler, 128, 3});
        EXPECT_EQ(next, count);
      }
    }
}

TEST(ObservableGrid, SingleWorkerTbbAndWaitObservation)
{
  TbbScheduler scheduler(1);
  auto caller = std::this_thread::get_id();
  bool observing = false;
  std::size_t observed = 0, collected = 0;
  struct Scope
  {
      bool& b;
      ~Scope() { b = false; }
  };
  bethe::evaluate_observable_grid(
      17, [](auto i) { return i; },
      [&](auto i, auto v) {
        EXPECT_FALSE(observing);
        EXPECT_EQ(caller, std::this_thread::get_id());
        EXPECT_EQ(i, collected++);
        EXPECT_EQ(i, v);
      },
      {&scheduler, 128, 3},
      [&](auto i, auto n) {
        EXPECT_FALSE(observing);
        EXPECT_EQ(i, observed++);
        EXPECT_EQ(n, 1u);
        observing = true;
        return Scope{observing};
      });
  EXPECT_EQ(observed, 17u);
  EXPECT_EQ(collected, 17u);
  run(true, 17, [](auto i) { return i; }, [](auto i, auto v) { EXPECT_EQ(i, v); }, {&scheduler, 128, 3});
}

TEST(ObservableGrid, CancellationPreventsSuccessfulPartialCompletion)
{
  TbbScheduler scheduler(2);
  for (bool async : {false, true})
    for (bool pre_cancel : {false, true})
    {
      std::stop_source stop;
      Counts counts;
      std::size_t collected = 0;
      if (pre_cancel) stop.request_stop();
      EXPECT_THROW(run(
                       async, 100, [&](auto i) { return std::make_unique<Payload>(counts, i); },
                       [&](auto i, auto) {
                         EXPECT_EQ(i, collected++);
                         stop.request_stop();
                       },
                       {&scheduler, 128, 3, stop.get_token()}),
                   bethe::ExecutionCancelled);
      EXPECT_EQ(collected, pre_cancel ? 0u : 1u);
      EXPECT_EQ(counts.live, 0);
    }
}

TEST(ObservableGrid, EvaluatorAndCollectorFailuresDrainOwnedWork)
{
  TbbScheduler scheduler(4);
  for (bool async : {false, true})
    for (bool fail_collect : {false, true})
    {
      Counts counts;
      std::mutex mutex;
      std::condition_variable changed;
      bool entered = false, release = false;
      std::atomic<bool> finished = false;
      auto fail = [&] {
        std::unique_lock lock(mutex);
        EXPECT_TRUE(changed.wait_for(lock, std::chrono::seconds(5), [&] { return entered; }));
        release = true;
        changed.notify_all();
        throw std::runtime_error("injected failure");
      };
      EXPECT_THROW(run(
                       async, 100,
                       [&](auto i) {
                         auto value = std::make_unique<Payload>(counts, i);
                         if (i == 1)
                         {
                           std::unique_lock lock(mutex);
                           entered = true;
                           changed.notify_all();
                           // A synchronous get_wait may run this kernel on the
                           // collector's own thread; do not wait for that same
                           // thread to enter collection while it is still here.
                           if (!fail_collect || async)
                             EXPECT_TRUE(changed.wait_for(lock, std::chrono::seconds(5), [&] { return release; }));
                           lock.unlock();
                           std::this_thread::sleep_for(std::chrono::milliseconds(10));
                           finished = true;
                         }
                         if (!fail_collect && i == 0) fail();
                         return value;
                       },
                       [&](auto, auto) {
                         if (fail_collect) fail();
                       },
                       {&scheduler, 128, 4}),
                   std::runtime_error);
      EXPECT_TRUE(finished); // the active sibling finished before error propagation
      EXPECT_EQ(counts.live, 0);
    }
}

TEST(ObservableGrid, DiscardDestroysPayloadBeforeReleasedStorageOwners)
{
  DebugScheduler scheduler;
  Counts counts;
  using Value = std::optional<std::unique_ptr<Payload>>;
  auto result = std::make_unique<Async<Value>>(std::make_unique<Payload>(counts, 0));
  auto retained_storage = result->read();
  (void)retained_storage.get_wait(scheduler);
  retained_storage.release();
  // A released buffer can retain storage after signalling readiness. It must
  // not retain a discarded payload past the enclosing operation's completion.
  bethe::detail::discard_observable_result(*result, scheduler);
  result.reset();
  EXPECT_EQ(counts.live, 0);
}

TEST(ObservableGrid, TakeDestroysMovedFromValueBeforeReleasedStorageOwners)
{
  struct Tracked
  {
      int& live;
      bool fail_move;
      Tracked(int& n, bool fail) : live(n), fail_move(fail) { ++live; }
      Tracked(Tracked&& other) : live(other.live), fail_move(other.fail_move)
      {
        if (fail_move) throw std::runtime_error("injected move failure");
        ++live;
      }
      ~Tracked() { --live; }
  };
  DebugScheduler scheduler;
  int live = 0;
  using Value = std::optional<Tracked>;
  for (bool fail_move : {false, true})
  {
    auto result = std::make_unique<Async<Value>>(std::in_place, live, fail_move);
    auto retained_storage = result->read();
    (void)retained_storage.get_wait(scheduler);
    retained_storage.release();
    if (fail_move)
      EXPECT_THROW(bethe::detail::take_observable_result(*result, scheduler), std::runtime_error);
    else
    {
      auto value = bethe::detail::take_observable_result(*result, scheduler);
      EXPECT_TRUE(value.has_value());
      EXPECT_EQ(live, 1);
    }
    result.reset();
    EXPECT_EQ(live, 0);
  }
}

TEST(ObservableGrid, WaitObserverFailureDrainsAcceptedWork)
{
  DebugScheduler scheduler;
  Counts counts;
  std::size_t observed = 0, collected = 0;
  EXPECT_THROW(bethe::evaluate_observable_grid(
                   100, [&](auto i) { return std::make_unique<Payload>(counts, i); }, [&](auto, auto) { ++collected; },
                   {&scheduler, 128, 3},
                   [&](auto, auto) {
                     if (++observed == 2) throw std::runtime_error("observer failed");
                     return 0;
                   }),
               std::runtime_error);
  EXPECT_EQ(collected, 1u);
  EXPECT_EQ(counts.live, 0);
  scheduler.run_all();
  EXPECT_EQ(counts.live, 0);
}

TEST(ObservableGrid, PartialAdmissionFailureAndPreflight)
{
  struct Reject : DebugScheduler
  {
      int calls = 0, reject = 1;
      void schedule(AsyncTask&& task) override
      {
        if (++calls == reject) throw std::bad_alloc{};
        DebugScheduler::schedule(std::move(task));
      }
  };
  for (bool async : {false, true})
    for (int reject : {1, 2, 4})
    {
      Reject scheduler;
      scheduler.reject = reject;
      Counts counts;
      EXPECT_THROW(run(
                       async, 100, [&](auto i) { return std::make_unique<Payload>(counts, i); }, [](auto, auto) {},
                       {&scheduler, 128, 3}),
                   std::bad_alloc);
      EXPECT_EQ(counts.live, 0);
      scheduler.run_all();
      EXPECT_EQ(counts.live, 0);
    }
  DebugScheduler scheduler;
  auto eval = [](auto i) { return i; };
  auto collect = [](auto, auto) {};
  EXPECT_THROW(bethe::evaluate_observable_grid(0, eval, collect, {&scheduler, 0, 3}), std::invalid_argument);
  EXPECT_THROW(bethe::evaluate_observable_grid(0, eval, collect, {&scheduler, 3, 0}), std::invalid_argument);
  EXPECT_THROW((void)bethe::evaluate_observable_grid_async(0, eval, collect, {&scheduler, 3, 0}),
               std::invalid_argument);
}

TEST(ObservableGrid, AsyncDependenciesOwnedCallablesAndDroppedHandle)
{
  DebugScheduler first, other;
  std::size_t calls = 0;
  auto owner = std::make_shared<int>(7);
  std::weak_ptr<int> weak = owner;
  auto pending = [&] {
    ScopedScheduler use(&first);
    return bethe::evaluate_observable_grid_async(
        10, [owner](auto i) { return i + *owner; },
        [&calls](auto i, auto v) {
          EXPECT_EQ(v, i + 7);
          ++calls;
        },
        {nullptr, 128, 3});
  }();
  owner.reset();
  EXPECT_FALSE(weak.expired());
  Async<std::size_t> result;
  first.schedule(downstream(pending.read(), result.write()));
  ScopedScheduler use_other(&other);
  EXPECT_EQ(result.get_wait(first), 11u);
  EXPECT_EQ(calls, 10u);
  first.run_all();
  EXPECT_TRUE(weak.expired());
  EXPECT_EQ(get_global_scheduler(), &other);
  {
    auto discarded = bethe::evaluate_observable_grid_async(
        10, [](auto i) { return i; }, [&calls](auto, auto) { ++calls; }, {&first, 128, 3});
  }
  first.run_all();
  EXPECT_EQ(calls, 20u);
}

namespace
{
template <typename Point> void equal_point(Point const& a, Point const& b)
{
  EXPECT_EQ(a.branch, b.branch);
  EXPECT_EQ(a.momentum, b.momentum);
  EXPECT_EQ(a.energy, b.energy);
  EXPECT_EQ(a.energy_error, b.energy_error);
  EXPECT_EQ(a.momentum_error, b.momentum_error);
  EXPECT_EQ(a.iterations, b.iterations);
  EXPECT_EQ(a.converged, b.converged);
  EXPECT_EQ(a.status, b.status);
  if constexpr (requires { a.parameter; }) EXPECT_EQ(a.parameter, b.parameter);
  if constexpr (requires { a.rapidity; }) EXPECT_EQ(a.rapidity, b.rapidity);
  if constexpr (requires { a.edge_distance; }) EXPECT_EQ(a.edge_distance, b.edge_distance);
  if constexpr (requires { a.evaluations; }) EXPECT_EQ(a.evaluations, b.evaluations);
  if constexpr (requires { a.symmetric_energy; }) EXPECT_EQ(a.symmetric_energy, b.symmetric_energy);
  if constexpr (requires { a.interaction; }) EXPECT_EQ(a.interaction, b.interaction);
  if constexpr (requires { a.energy_offset; }) EXPECT_EQ(a.energy_offset, b.energy_offset);
  if constexpr (requires { a.fermi_energy; }) EXPECT_EQ(a.fermi_energy, b.fermi_energy);
  if constexpr (requires { a.delta_particles; }) EXPECT_EQ(a.delta_particles, b.delta_particles);
  if constexpr (requires { a.spin; }) EXPECT_EQ(a.spin, b.spin);
  if constexpr (requires { a.convention; }) EXPECT_EQ(a.convention, b.convention);
  if constexpr (requires { a.reference; }) EXPECT_EQ(a.reference, b.reference);
}

template <typename Evaluate> void compare_grid(std::size_t n, Evaluate evaluate)
{
  using Point = std::invoke_result_t<Evaluate, std::size_t>;
  std::vector<Point> reference;
  for (std::size_t i = 0; i < n; ++i)
    reference.push_back(evaluate(i));
  DebugScheduler serial({.order = DebugSchedulerOrder::reverse});
  TbbScheduler parallel(4);
  for (auto* scheduler : std::vector<IAsyncScheduler*>{&serial, &parallel})
    for (bool async : {false, true})
    {
      std::size_t next = 0;
      run(async, n, evaluate,
          [&](auto i, auto const& value) {
            EXPECT_EQ(i, next++);
            equal_point(value, reference[i]);
          },
          {scheduler, 128, 3});
      EXPECT_EQ(next, n);
    }
}
} // namespace

template <typename T> class ObservableGridModels : public ::testing::Test {};
TYPED_TEST_SUITE(ObservableGridModels, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(ObservableGridModels, LiebLinigerPreparedBackground)
{
  using R = TypeParam;
  namespace model = bethe::lieb_liniger::thermo;
  R const pi = R{4} * std::atan(R{1});
  for (std::size_t iterations : {160, 0})
  {
    model::Options<R> options;
    options.max_momentum_iterations = iterations;
    model::Solver<R> const solver(R{4}, R{1}, options);
    ASSERT_TRUE(solver.background().converged);
    compare_grid(10, [&](std::size_t i) {
      return solver.at_momentum(i / 5 == 0 ? model::Branch::type_i : model::Branch::type_ii,
                                R{2} * pi * R(i % 5) / R{4});
    });
  }
}

TYPED_TEST(ObservableGridModels, HubbardDopedPreparedBackground)
{
  using R = TypeParam;
  namespace model = bethe::hubbard::thermo;
  for (std::size_t iterations : {160, 0})
  {
    model::DopedOptions<R> options;
    options.max_iterations = iterations;
    model::DopedSolver<R> const solver(R{4}, R{0.5}, options);
    ASSERT_TRUE(solver.background().converged);
    compare_grid(15, [&](std::size_t i) {
      auto b =
          std::array{model::DopedBranch::spinon, model::DopedBranch::holon, model::DopedBranch::charge_particle}[i / 5];
      auto [lo, hi] = solver.momentum_range(b);
      return solver.at_momentum(b, i % 5 == 4 ? hi : lo + (hi - lo) * R(i % 5) / R{4});
    });
  }
}

TYPED_TEST(ObservableGridModels, HubbardHalfFilledAndFailureValues)
{
  using R = TypeParam;
  namespace model = bethe::hubbard::thermo;
  R const pi = R{4} * std::atan(R{1});
  for (std::size_t budget : {1000000, 0})
    for (auto convention : {model::Convention::symmetric, model::Convention::unshifted})
    {
      model::Options<R> options;
      options.max_evaluations = budget;
      compare_grid(9, [&](std::size_t i) {
        auto b = std::array{model::Branch::spinon, model::Branch::holon, model::Branch::antiholon}[i / 3];
        return model::dispersion(b, R{4}, pi * R(i % 3) / R{2}, convention, options);
      });
    }
}
