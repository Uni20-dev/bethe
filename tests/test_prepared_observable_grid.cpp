// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include <bethe/hubbard_doped.hpp>
#include <bethe/lieb_liniger_thermo.hpp>
#include <bethe/prepared_observable_grid.hpp>
#include <uni20/async/tbb_scheduler.hpp>
#if UNI20_HAS_MPREAL
#include <uni20/core/mpreal.hpp>
#endif
#include <array>
#include <atomic>
#include <chrono>
#include <future>
#include <memory>
#include <optional>

namespace
{
using namespace uni20::async;
struct Counts
{
    std::atomic<int> prepared{0}, copies{0}, live{0}, queries{0}, writes{0};
};
struct Background
{
    Counts* counts;
    std::unique_ptr<int> value;
    explicit Background(Counts& c) : counts(&c), value(std::make_unique<int>(7)) { ++counts->live; }
    Background(Background const& b) : counts(b.counts), value(std::make_unique<int>(*b.value))
    {
      ++counts->copies;
      ++counts->live;
    }
    Background(Background&&) noexcept = default;
    ~Background()
    {
      if (value) --counts->live;
    }
};
AsyncTask overwrite(WriteBuffer<Background> out, Counts& counts, WriteBuffer<bool> completed)
{
  auto done = co_await completed;
  auto value = co_await out;
  EXPECT_EQ(counts.queries, 24);
  *value->value = 99;
  ++counts.writes;
  done = true;
}
} // namespace

TEST(PreparedGrid, TwoReadersRetainOneBackgroundWithoutCopies)
{
  for (auto order : {DebugSchedulerOrder::fifo, DebugSchedulerOrder::reverse, DebugSchedulerOrder::random})
  {
    DebugScheduler scheduler({.order = order, .random_seed = 731});
    Counts counts;
    std::array<std::vector<std::size_t>, 2> rows;
    Async<bool> written;
    auto results = [&] {
      auto background = bethe::prepare_background_async(
          [&] {
            ++counts.prepared;
            return Background(counts);
          },
          {&scheduler});
      auto submit = [&](std::size_t branch) {
        return bethe::evaluate_observable_grid_async(
            background.read(), 12,
            [&](Background const& b, std::size_t i) {
              EXPECT_EQ(counts.writes, 0);
              ++counts.queries;
              return i + *b.value;
            },
            [&, branch](std::size_t i, auto value) {
              EXPECT_EQ(i, rows[branch].size());
              rows[branch].push_back(value);
            },
            {&scheduler, 128, 3});
      };
      auto a = submit(0), b = submit(1);
      scheduler.schedule(overwrite(background.write(), counts, written.write()));
      EXPECT_EQ(counts.prepared, 0); // submission does not wait for preparation
      EXPECT_EQ(counts.queries, 0);
      return std::pair{std::move(a), std::move(b)};
    }(); // neither reader borrows the submission-local handle
    EXPECT_EQ(results.first.get_wait(scheduler), 12u);
    EXPECT_EQ(results.second.get_wait(scheduler), 12u);
    EXPECT_TRUE(written.get_wait(scheduler));
    EXPECT_EQ(counts.prepared, 1);
    EXPECT_EQ(counts.copies, 0);
    EXPECT_EQ(counts.live, 0);
    EXPECT_EQ(counts.writes, 1);
    EXPECT_EQ(rows[0], rows[1]);
    for (std::size_t i = 0; i < rows[0].size(); ++i)
      EXPECT_EQ(rows[0][i], i + 7);
  }
}

TEST(PreparedGrid, OwnedMoveOnlyCallbacksAndDiscardedHandles)
{
  DebugScheduler scheduler, other;
  Counts counts;
  auto owner = std::make_shared<int>(9);
  std::weak_ptr<int> weak = owner;
  std::size_t collected = 0;
  {
    ScopedScheduler active(&scheduler);
    auto background = bethe::prepare_background_async([&, owner, move_only = std::make_unique<int>(9)] {
      EXPECT_EQ(*owner, *move_only);
      ++counts.prepared;
      return Background(counts);
    });
    auto dropped = bethe::evaluate_observable_grid_async(
        background.read(), 20,
        [p = std::make_unique<int>(1)](Background const& b, auto i) { return i + *b.value + *p; },
        [&collected](auto i, auto v) {
          EXPECT_EQ(v, i + 8);
          ++collected;
        });
  }
  owner.reset();
  EXPECT_FALSE(weak.expired());
  ScopedScheduler active(&other);
  // Deliberately test discarded handles with the debug scheduler. Normal
  // callers retain and wait on their operation's result, not the whole queue.
  scheduler.run_all();
  EXPECT_EQ(collected, 20u);
  EXPECT_EQ(counts.prepared, 1);
  EXPECT_EQ(counts.live, 0);
  EXPECT_TRUE(weak.expired());
  EXPECT_EQ(get_global_scheduler(), &other);
}

TEST(PreparedGrid, PreparationExceptionsPropagateToBothReadersIncludingEmptyGrid)
{
  DebugScheduler scheduler;
  auto background =
      bethe::prepare_background_async([]() -> int { throw std::domain_error("preparation"); }, {&scheduler});
  int evaluated = 0, collected = 0;
  auto submit = [&](std::size_t n) {
    return bethe::evaluate_observable_grid_async(
        background.read(), n,
        [&](int, auto) {
          ++evaluated;
          return 0;
        },
        [&](auto, auto) { ++collected; }, {&scheduler});
  };
  auto a = submit(10), b = submit(0);
  EXPECT_THROW((void)a.get_wait(scheduler), std::domain_error);
  EXPECT_THROW((void)b.get_wait(scheduler), std::domain_error);
  EXPECT_EQ(evaluated, 0);
  EXPECT_EQ(collected, 0);
}

TEST(PreparedGrid, CancellationDoesNotCancelSharedPreparationOrSibling)
{
  DebugScheduler scheduler;
  std::stop_source stop;
  int prepared = 0, cancelled_queries = 0, other_queries = 0;
  auto background = bethe::prepare_background_async(
      [&] {
        ++prepared;
        return 7;
      },
      {&scheduler});
  auto a = bethe::evaluate_observable_grid_async(
      background.read(), 8,
      [&](int, auto) {
        ++cancelled_queries;
        return 0;
      },
      [](auto, auto) {}, {.scheduler = &scheduler, .stop = stop.get_token()});
  auto b = bethe::evaluate_observable_grid_async(
      background.read(), 8,
      [&](int, auto) {
        ++other_queries;
        return 0;
      },
      [](auto, auto) {}, {&scheduler});
  stop.request_stop();
  EXPECT_THROW((void)a.get_wait(scheduler), bethe::ExecutionCancelled);
  EXPECT_EQ(b.get_wait(scheduler), 8u);
  EXPECT_EQ(prepared, 1);
  EXPECT_EQ(cancelled_queries, 0);
  EXPECT_EQ(other_queries, 8);
}

TEST(PreparedGrid, CollectorFailureAndCancellationPreserveSibling)
{
  for (bool cancel : {false, true})
  {
    DebugScheduler scheduler;
    std::stop_source stop;
    Counts counts;
    auto background = bethe::prepare_background_async([&] { return Background(counts); }, {&scheduler});
    auto a = bethe::evaluate_observable_grid_async(
        background.read(), 20, [](Background const& b, auto i) { return i + *b.value; },
        [&](auto, auto) {
          if (!cancel) throw std::domain_error("collector");
          stop.request_stop();
        },
        {.scheduler = &scheduler, .max_in_flight = 3, .stop = stop.get_token()});
    std::size_t other = 0;
    auto b = bethe::evaluate_observable_grid_async(
        background.read(), 20, [](Background const& b, auto i) { return i + *b.value; },
        [&](auto i, auto value) {
          EXPECT_EQ(value, i + 7);
          ++other;
        },
        {&scheduler, 128, 3});
    if (cancel)
      EXPECT_THROW((void)a.get_wait(scheduler), bethe::ExecutionCancelled);
    else
      EXPECT_THROW((void)a.get_wait(scheduler), std::domain_error);
    EXPECT_EQ(b.get_wait(scheduler), 20u);
    EXPECT_EQ(other, 20u);
    EXPECT_EQ(counts.copies, 0);
  }
}

TEST(PreparedGrid, EmptySuccessAndPreflightValidation)
{
  DebugScheduler scheduler;
  int preparations = 0, queries = 0;
  auto prepare = [&] {
    ++preparations;
    return 3;
  };
  EXPECT_THROW((void)bethe::prepare_background_async(prepare, {&scheduler, 128, 0}), std::invalid_argument);
  auto background = bethe::prepare_background_async(prepare, {&scheduler});
  auto evaluate = [&](int, auto) {
    ++queries;
    return 3;
  };
  auto collect = [](auto, auto) {};
  EXPECT_THROW(
      (void)bethe::evaluate_observable_grid_async(background.read(), 1, evaluate, collect, {&scheduler, 128, 0}),
      std::invalid_argument);
  auto empty = bethe::evaluate_observable_grid_async(background.read(), 0, evaluate, collect, {&scheduler});
  EXPECT_EQ(empty.get_wait(scheduler), 0u);
  EXPECT_EQ(preparations, 1);
  EXPECT_EQ(queries, 0);
}

TEST(PreparedGrid, PreparationCancellationBeforeAndAfterKernel)
{
  for (bool before : {false, true})
  {
    DebugScheduler scheduler;
    std::stop_source stop;
    if (before) stop.request_stop();
    Counts counts;
    auto background = bethe::prepare_background_async(
        [&] {
          ++counts.prepared;
          stop.request_stop();
          return Background(counts);
        },
        {.scheduler = &scheduler, .stop = stop.get_token()});
    auto grid = bethe::evaluate_observable_grid_async(
        background.read(), 10,
        [&](Background const&, auto) {
          ++counts.queries;
          return 0;
        },
        [](auto, auto) {}, {&scheduler});
    EXPECT_THROW((void)grid.get_wait(scheduler), bethe::ExecutionCancelled);
    EXPECT_THROW((void)background.get_wait(scheduler), bethe::ExecutionCancelled);
    EXPECT_EQ(counts.prepared, before ? 0 : 1);
    EXPECT_EQ(counts.queries, 0);
  }
}

TEST(PreparedGrid, AdmissionAndQueryFailuresDrainBeforeBackgroundRelease)
{
  struct Reject : DebugScheduler
  {
      int calls = 0, reject = 0;
      void schedule(AsyncTask&& task) override
      {
        if (++calls == reject) throw std::bad_alloc{};
        DebugScheduler::schedule(std::move(task));
      }
  };
  for (int reject : {1, 2, 3, 5, 0})
  {
    Reject scheduler;
    scheduler.reject = reject;
    Counts counts;
    std::optional<Async<Background>> background;
    bool failed = false;
    try
    {
      background.emplace(bethe::prepare_background_async([&] { return Background(counts); }, {&scheduler}));
      auto grid = bethe::evaluate_observable_grid_async(
          background->read(), 20,
          [&](Background const&, auto i) {
            ++counts.queries;
            if (!reject && i == 1) throw std::runtime_error("query");
            return i;
          },
          [](auto, auto) {}, {&scheduler, 128, 3});
      (void)grid.get_wait(scheduler);
      ADD_FAILURE() << "expected admission or query failure";
    }
    catch (std::bad_alloc const&)
    {
      EXPECT_NE(reject, 0);
      failed = true;
    }
    catch (std::runtime_error const& error)
    {
      EXPECT_EQ(reject, 0);
      EXPECT_STREQ(error.what(), "query");
      failed = true;
    }
    EXPECT_TRUE(failed);
    // Preparation can be admitted even when grid submission fails. Observe
    // that specific result; a queue-wide drain could hide orphaned queries.
    if (background) EXPECT_NO_THROW((void)background->get_wait(scheduler));
    background.reset();
    EXPECT_EQ(counts.live, 0);
  }
}

namespace
{
enum class GridFailure
{
  query,
  collector,
  cancellation
};
class PreparedGridDrain : public ::testing::TestWithParam<GridFailure> {};

AsyncTask overwrite_after_queries(WriteBuffer<int> out, std::atomic<bool>& finished, std::promise<void>& overwritten)
{
  auto value = co_await out;
  EXPECT_TRUE(finished.load());
  value = 99;
  overwritten.set_value();
}
} // namespace

TEST_P(PreparedGridDrain, CompletionWaitsForHeldQuery)
{
  using namespace std::chrono_literals;
  TbbScheduler scheduler(4);
  std::promise<void> entered, release, triggered, overwritten;
  auto entered_future = entered.get_future().share();
  auto release_future = release.get_future().share();
  auto triggered_future = triggered.get_future();
  auto overwritten_future = overwritten.get_future();
  std::atomic<bool> finished{false};
  std::atomic<int> sibling_queries{0};
  std::stop_source stop;
  auto const failure = GetParam();
  auto background = bethe::prepare_background_async([] { return 7; }, {&scheduler});
  auto failed = bethe::evaluate_observable_grid_async(
      background.read(), 2,
      [&](int const& value, std::size_t index) {
        if (index == 1)
        {
          entered.set_value();
          if (release_future.wait_for(5s) != std::future_status::ready)
            throw std::runtime_error("held query timed out");
          EXPECT_EQ(value, 7);
          finished = true;
        }
        else
        {
          if (entered_future.wait_for(5s) != std::future_status::ready)
            throw std::runtime_error("second query did not enter");
          if (failure == GridFailure::query)
          {
            triggered.set_value();
            throw std::domain_error("query");
          }
        }
        return value;
      },
      [&](auto, auto) {
        triggered.set_value();
        if (failure == GridFailure::collector) throw std::domain_error("collector");
        stop.request_stop();
      },
      {.scheduler = &scheduler, .max_in_flight = 2, .stop = stop.get_token()});
  auto sibling = bethe::evaluate_observable_grid_async(
      background.read(), 12,
      [&](int const& value, auto) {
        ++sibling_queries;
        return value;
      },
      [](auto, auto value) { EXPECT_EQ(value, 7); }, {.scheduler = &scheduler, .max_in_flight = 2});
  scheduler.schedule(overwrite_after_queries(background.write(), finished, overwritten));
  auto ready = failed.read();

  // Check while the query is still held, before any final wait can conceal
  // premature completion. Other workers can finish the sibling and propagate
  // the failure, but neither the failed result nor the writer may become ready.
  // Keep these assertions nonfatal so we always release and join accepted work.
  EXPECT_EQ(triggered_future.wait_for(5s), std::future_status::ready);
  EXPECT_NO_THROW(EXPECT_EQ(sibling.get_wait(scheduler), 12u));
  EXPECT_EQ(sibling_queries.load(), 12);
  EXPECT_EQ(overwritten_future.wait_for(20ms), std::future_status::timeout);
  EXPECT_FALSE(ready.await_ready());
  EXPECT_FALSE(finished.load());
  release.set_value();
  try
  {
    (void)failed.get_wait(scheduler);
    ADD_FAILURE() << "expected grid failure";
  }
  catch (bethe::ExecutionCancelled const&)
  {
    EXPECT_EQ(failure, GridFailure::cancellation);
  }
  catch (std::domain_error const& error)
  {
    EXPECT_NE(failure, GridFailure::cancellation);
    EXPECT_STREQ(error.what(), failure == GridFailure::query ? "query" : "collector");
  }
  catch (std::exception const& error)
  {
    ADD_FAILURE() << "unexpected execution error: " << error.what();
  }
  EXPECT_TRUE(finished.load());
  EXPECT_EQ(background.get_wait(scheduler), 99); // joins the dependent writer
}

INSTANTIATE_TEST_SUITE_P(PreparedGrid, PreparedGridDrain,
                         ::testing::Values(GridFailure::query, GridFailure::collector, GridFailure::cancellation),
                         [](auto const& info) {
                           switch (info.param)
                           {
                             case GridFailure::query:
                               return "QueryFailure";
                             case GridFailure::collector:
                               return "CollectorFailure";
                             case GridFailure::cancellation:
                               return "Cancellation";
                           }
                           return "Unknown";
                         });

TEST(PreparedGrid, SingleWorkerComposesBackgroundAndGrid)
{
  TbbScheduler scheduler(1);
  auto background = bethe::prepare_background_async([] { return 7; }, {&scheduler});
  auto grid = bethe::evaluate_observable_grid_async(
      background.read(), 20, [](int value, auto i) { return value + i; },
      [](auto i, auto value) { EXPECT_EQ(value, i + 7); }, {&scheduler});
  EXPECT_EQ(grid.get_wait(scheduler), 20u);
}

TEST(PreparedGrid, BoundsConcurrentQueriesAndRetainedResultsAcrossTwoGrids)
{
  TbbScheduler scheduler(4);
  struct Live
  {
      std::atomic<int> current{0}, peak{0};
  } live;
  struct Row
  {
      Live& live;
      explicit Row(Live& l) : live(l)
      {
        int n = ++live.current, old = live.peak.load();
        while (old < n && !live.peak.compare_exchange_weak(old, n))
        {}
      }
      ~Row() { --live.current; }
  };
  Counts counts;
  auto background = bethe::prepare_background_async(
      [&] {
        ++counts.prepared;
        return Background(counts);
      },
      {&scheduler});
  auto submit = [&] {
    return bethe::evaluate_observable_grid_async(
        background.read(), 200,
        [&](Background const& b, auto) {
          EXPECT_EQ(*b.value, 7);
          ++counts.queries;
          return std::make_unique<Row>(live);
        },
        [](auto, auto) {}, {&scheduler, 128, 3});
  };
  auto a = submit(), b = submit();
  EXPECT_EQ(a.get_wait(scheduler), 200u);
  EXPECT_EQ(b.get_wait(scheduler), 200u);
  EXPECT_EQ(counts.prepared, 1);
  EXPECT_EQ(counts.copies, 0);
  EXPECT_EQ(counts.queries, 400);
  EXPECT_LE(live.peak, 6);
  EXPECT_EQ(live.current, 0);
}

template <typename T> class PreparedGridModels : public ::testing::Test {};
TYPED_TEST_SUITE(PreparedGridModels, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(PreparedGridModels, LiebLinigerBothBranchesAndFailureDiagnostics)
{
  using R = TypeParam;
  namespace model = bethe::lieb_liniger::thermo;
  R const endpoint = R{8} * std::atan(R{1});
  for (int failure : {0, 1, 2})
  {
    model::Options<R> options;
    if (failure == 1) options.max_iterations = 0;
    if (failure == 2) options.max_momentum_iterations = 0;
    model::Solver<R> const reference(R{4}, R{1}, options);
    DebugScheduler scheduler;
    int preparations = 0;
    auto background = bethe::prepare_background_async(
        [&, options] {
          ++preparations;
          return model::Solver<R>(R{4}, R{1}, options);
        },
        {&scheduler});
    std::array<std::vector<model::Point<R>>, 2> results;
    auto submit = [&](std::size_t branch) {
      auto kind = branch ? model::Branch::type_ii : model::Branch::type_i;
      return bethe::evaluate_observable_grid_async(
          background.read(), 5,
          [=](model::Solver<R> const& solver, auto i) { return solver.at_momentum(kind, endpoint * R(i) / R{4}); },
          [&, branch](auto, auto point) { results[branch].push_back(std::move(point)); }, {&scheduler, 128, 2});
    };
    auto a = submit(0), b = submit(1);
    EXPECT_EQ(preparations, 0);
    EXPECT_EQ(a.get_wait(scheduler), 5u);
    EXPECT_EQ(b.get_wait(scheduler), 5u);
    EXPECT_EQ(preparations, 1);
    for (std::size_t branch = 0; branch < 2; ++branch)
    {
      ASSERT_EQ(results[branch].size(), 5u);
      for (std::size_t i = 0; i < 5; ++i)
      {
        auto expected =
            reference.at_momentum(branch ? model::Branch::type_ii : model::Branch::type_i, endpoint * R(i) / R{4});
        auto const& got = results[branch][i];
        EXPECT_EQ(got.branch, expected.branch);
        EXPECT_EQ(got.momentum, expected.momentum);
        EXPECT_EQ(got.energy, expected.energy);
        EXPECT_EQ(got.energy_error, expected.energy_error);
        EXPECT_EQ(got.momentum_error, expected.momentum_error);
        EXPECT_EQ(got.rapidity, expected.rapidity);
        EXPECT_EQ(got.edge_distance, expected.edge_distance);
        EXPECT_EQ(got.iterations, expected.iterations);
        EXPECT_EQ(got.status, expected.status);
        EXPECT_EQ(got.converged, expected.converged);
        if (failure == 1)
        {
          EXPECT_FALSE(got.converged);
          EXPECT_FALSE(got.energy);
        }
        if (!failure && i == 0)
        {
          EXPECT_TRUE(got.converged);
          EXPECT_EQ(got.energy, R{0});
        }
      }
    }
  }
}

TYPED_TEST(PreparedGridModels, HubbardBranchesShareBackgroundAndPreserveFailures)
{
  using R = TypeParam;
  namespace model = bethe::hubbard::thermo;
  std::array const branches{model::DopedBranch::spinon, model::DopedBranch::holon, model::DopedBranch::charge_particle};
  // Exercise both queued dependency submission and actual simultaneous queries.
  auto check = [&](auto& scheduler, int failure) {
    model::DopedOptions<R> options;
    if (failure == 1) options.max_background_iterations = 0;
    if (failure == 2) options.max_iterations = 0;
    model::DopedSolver<R> const reference(R{4}, R{1} / R{2}, options);
    auto preparations = std::make_shared<std::atomic<int>>(0);
    auto results = std::make_shared<std::array<std::vector<model::DopedPoint<R>>, 3>>();
    auto background = bethe::prepare_background_async(
        [options, preparations] {
          ++*preparations;
          return model::DopedSolver<R>(R{4}, R{1} / R{2}, options);
        },
        {&scheduler});
    auto submit = [&](std::size_t b) {
      return bethe::evaluate_observable_grid_async(
          background.read(), 5,
          [branch = branches[b], b](model::DopedSolver<R> const& solver, std::size_t i) {
            auto const [low, high] = solver.momentum_range(branch);
            auto const p = i == 4 ? high : low + (high - low) * R(i) / R{4};
            return solver.at_momentum(branch, p, model::Convention::unshifted,
                                      b == 1 ? model::EnergyReference::fermi : model::EnergyReference::hamiltonian);
          },
          [results, b](std::size_t i, auto point) {
            EXPECT_EQ(i, (*results)[b].size());
            (*results)[b].push_back(std::move(point));
          },
          {.scheduler = &scheduler, .max_in_flight = 2});
    };
    auto a = submit(0), b = submit(1), c = submit(2);
    if constexpr (std::is_same_v<std::decay_t<decltype(scheduler)>, DebugScheduler>) EXPECT_EQ(*preparations, 0);
    EXPECT_EQ(a.get_wait(scheduler), 5u);
    EXPECT_EQ(b.get_wait(scheduler), 5u);
    EXPECT_EQ(c.get_wait(scheduler), 5u);
    EXPECT_EQ(*preparations, 1);
    auto const& prepared = background.get_wait(scheduler).background();
    EXPECT_EQ(prepared.status, reference.background().status);
    EXPECT_EQ(prepared.fermi_rapidity, reference.background().fermi_rapidity);
    EXPECT_EQ(prepared.mu_unshifted, reference.background().mu_unshifted);
    EXPECT_EQ(prepared.energy_per_site_unshifted, reference.background().energy_per_site_unshifted);
    std::size_t failed = 0;
    for (std::size_t branch = 0; branch < branches.size(); ++branch)
    {
      ASSERT_EQ((*results)[branch].size(), 5u);
      auto const [low, high] = reference.momentum_range(branches[branch]);
      for (std::size_t i = 0; i < 5; ++i)
      {
        auto const p = i == 4 ? high : low + (high - low) * R(i) / R{4};
        auto const expected =
            reference.at_momentum(branches[branch], p, model::Convention::unshifted,
                                  branch == 1 ? model::EnergyReference::fermi : model::EnergyReference::hamiltonian);
        auto const& got = (*results)[branch][i];
        EXPECT_EQ(got.branch, expected.branch);
        EXPECT_EQ(got.convention, expected.convention);
        EXPECT_EQ(got.reference, expected.reference);
        EXPECT_EQ(got.momentum, expected.momentum);
        EXPECT_EQ(got.parameter, expected.parameter);
        EXPECT_EQ(got.energy, expected.energy);
        EXPECT_EQ(got.symmetric_energy, expected.symmetric_energy);
        EXPECT_EQ(got.fermi_energy, expected.fermi_energy);
        EXPECT_EQ(got.delta_particles, expected.delta_particles);
        EXPECT_EQ(got.spin, expected.spin);
        EXPECT_EQ(got.energy_error, expected.energy_error);
        EXPECT_EQ(got.momentum_error, expected.momentum_error);
        EXPECT_EQ(got.iterations, expected.iterations);
        EXPECT_EQ(got.converged, expected.converged);
        EXPECT_EQ(got.status, expected.status);
        failed += !got.converged;
        if (failure == 1) EXPECT_EQ(got.status, model::DopedStatus::density_limit);
      }
    }
    if (failure)
      EXPECT_GT(failed, 0u);
    else
      EXPECT_EQ(failed, 0u);
  };
  for (int failure : {0, 1, 2})
  {
    DebugScheduler queued;
    check(queued, failure);
    TbbScheduler threaded(3);
    check(threaded, failure);
  }
}

#if UNI20_HAS_MPREAL
TEST(PreparedGrid, IndependentMpfrPrecisionsOnSharedScheduler)
{
  TbbScheduler scheduler(4);
  using R = uni20::mpreal;
  auto submit = [&](int bits) {
    auto precision = uni20::Precision::bits(bits);
    auto background = bethe::prepare_background_async([precision] { return R(2, precision); }, {&scheduler});
    return bethe::evaluate_observable_grid_async(
        background.read(), 20, [](R const& x, auto) { return sqrt(x); },
        [precision](auto, R const& value) {
          EXPECT_EQ(value.precision(), precision);
          EXPECT_EQ(value, sqrt(R(2, precision)));
        },
        {&scheduler, 128, 3});
  };
  auto a = submit(80), b = submit(333);
  EXPECT_EQ(a.get_wait(scheduler), 20u);
  EXPECT_EQ(b.get_wait(scheduler), 20u);
}
#endif
