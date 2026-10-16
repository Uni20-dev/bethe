// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <array>
#include <atomic>
#include <bethe/hubbard_doped.hpp>
#include <bethe/lieb_liniger_thermo.hpp>
#include <bethe/observable_grid.hpp>
#include <memory>
#include <numeric>
#include <thread>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
class RecordingScheduler : public uni20::async::DebugScheduler {
  public:
    std::vector<std::size_t> batches;

  private:
    void execute_batch_impl(uni20::async::LightweightTaskBatch const& batch) override
    {
      batches.push_back(batch.size());
      for (std::size_t i = batch.size(); i > 0; --i)
        batch(i - 1);
    }
};

struct BatchScope
{
    std::atomic<bool>& active;
    explicit BatchScope(std::atomic<bool>& a) : active(a) { active = true; }
    BatchScope(BatchScope const&) = delete;
    ~BatchScope() { active = false; }
};

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
  if constexpr (requires { a.fermi_energy; }) EXPECT_EQ(a.fermi_energy, b.fermi_energy);
  if constexpr (requires { a.delta_particles; }) EXPECT_EQ(a.delta_particles, b.delta_particles);
  if constexpr (requires { a.spin; }) EXPECT_EQ(a.spin, b.spin);
  if constexpr (requires { a.convention; }) EXPECT_EQ(a.convention, b.convention);
  if constexpr (requires { a.reference; }) EXPECT_EQ(a.reference, b.reference);
}

template <typename Evaluate> void check_schedulers(std::size_t count, Evaluate evaluate, bool complete)
{
  using Point = std::invoke_result_t<Evaluate, std::size_t>;
  std::vector<Point> reference;
  bool all_converged = true;
  for (std::size_t i = 0; i < count; ++i)
  {
    reference.push_back(evaluate(i)); // Independent scalar path, not the adapter.
    all_converged = all_converged && reference.back().converged;
  }
  EXPECT_EQ(all_converged, complete);
  RecordingScheduler reverse;
  uni20::async::TbbScheduler parallel(4);
  for (auto* scheduler : std::vector<uni20::async::IAsyncScheduler*>{&reverse, &parallel})
  {
    std::size_t next = 0;
    bethe::evaluate_observable_grid(count, evaluate,
                                    [&](std::size_t i, Point const& point) {
                                      ASSERT_EQ(i, next++);
                                      equal_point(reference[i], point);
                                    },
                                    {scheduler, 3});
    EXPECT_EQ(next, count);
  }
  EXPECT_EQ(std::accumulate(reverse.batches.begin(), reverse.batches.end(), std::size_t{0}), count);
}
} // namespace

TEST(ObservableGrid, OrderedBoundedMoveOnlyResultsAndScopedObservation)
{
  RecordingScheduler scheduler;
  auto const caller = std::this_thread::get_id();
  auto* original = uni20::async::get_global_scheduler();
  std::size_t live = 0, peak = 0, collected = 0, observed = 0;
  std::atomic<bool> active = false;
  struct Payload
  {
      std::size_t id;
      std::size_t& live;
      Payload(std::size_t i, std::size_t& n) : id(i), live(n) { ++live; }
      ~Payload() { --live; }
  };
  bethe::evaluate_observable_grid(
      10,
      [&](std::size_t i) {
        EXPECT_TRUE(active);
        auto value = std::make_unique<Payload>(i, live);
        peak = std::max(peak, live); // This recording scheduler is serial/reversed.
        return value;
      },
      [&](std::size_t i, auto value) {
        EXPECT_EQ(std::this_thread::get_id(), caller);
        EXPECT_FALSE(active);
        EXPECT_EQ(i, collected++);
        EXPECT_EQ(value->id, i);
      },
      {&scheduler, 3},
      [&](std::size_t offset, std::size_t size) {
        EXPECT_EQ(std::this_thread::get_id(), caller);
        EXPECT_FALSE(active);
        EXPECT_EQ(offset, observed);
        EXPECT_EQ(live, 0U);
        observed += size;
        return BatchScope(active);
      });
  EXPECT_EQ(scheduler.batches, (std::vector<std::size_t>{3, 3, 3, 1}));
  EXPECT_EQ(collected, 10U);
  EXPECT_EQ(observed, 10U);
  EXPECT_EQ(peak, 3U);
  EXPECT_EQ(live, 0U);
  EXPECT_EQ(uni20::async::get_global_scheduler(), original);
}

TEST(ObservableGrid, EmptySingletonActiveSchedulerAndNumericalFailureValues)
{
  RecordingScheduler active, explicit_scheduler;
  uni20::async::ScopedScheduler scope(&active);
  std::size_t evaluated = 0, collected = 0;
  auto evaluate = [&](std::size_t i) -> std::optional<std::size_t> {
    ++evaluated;
    return i % 2 ? std::nullopt : std::optional(i);
  };
  auto collect = [&](std::size_t i, auto result) {
    EXPECT_EQ(i, collected++);
    EXPECT_EQ(result, i % 2 ? std::nullopt : std::optional(i));
  };
  bethe::evaluate_observable_grid(0, evaluate, collect);
  EXPECT_EQ(evaluated, 0U);
  EXPECT_TRUE(active.batches.empty());
  EXPECT_THROW(bethe::evaluate_observable_grid(0, evaluate, collect, {nullptr, 0}), std::invalid_argument);
  bethe::evaluate_observable_grid(1, evaluate, collect);
  EXPECT_EQ(active.batches, (std::vector<std::size_t>{1}));
  collected = 0;
  bethe::evaluate_observable_grid(7, evaluate, collect, {&explicit_scheduler, 3});
  EXPECT_EQ(evaluated, 8U);
  EXPECT_EQ(collected, 7U); // Failed numerical points are not filtered out.
  EXPECT_EQ(active.batches.size(), 1U);
  EXPECT_EQ(uni20::async::get_global_scheduler(), &active);
}

TEST(ObservableGrid, WorkerAndCollectorExceptionsJoinBeforeReturning)
{
  uni20::async::TbbScheduler scheduler(4);
  auto const caller = std::this_thread::get_id();
  std::atomic<std::size_t> running = 0;
  std::atomic<bool> observing = false;
  std::vector<std::size_t> collected;
  auto observe = [&](std::size_t, std::size_t) { return BatchScope(observing); };
  auto evaluate = [&](std::size_t i) {
    struct Guard
    {
        std::atomic<std::size_t>& n;
        explicit Guard(std::atomic<std::size_t>& count) : n(count) { ++n; }
        ~Guard() { --n; }
    } guard(running);
    EXPECT_TRUE(observing);
    if (i == 4) throw std::runtime_error("injected evaluation failure");
    return i;
  };
  auto collect = [&](std::size_t i, std::size_t value) {
    EXPECT_EQ(std::this_thread::get_id(), caller);
    EXPECT_FALSE(observing);
    EXPECT_EQ(running, 0U);
    EXPECT_EQ(value, i);
    collected.push_back(i);
  };
  EXPECT_THROW(bethe::evaluate_observable_grid(9, evaluate, collect, {&scheduler, 3}, observe), std::runtime_error);
  EXPECT_EQ(collected, (std::vector<std::size_t>{0, 1, 2}));
  EXPECT_EQ(running, 0U);
  EXPECT_FALSE(observing);
  auto fail_collect = [&](std::size_t i, std::size_t value) {
    collect(i, value);
    throw std::runtime_error("injected output failure");
  };
  EXPECT_THROW(bethe::evaluate_observable_grid(9, evaluate, fail_collect, {&scheduler, 3}, observe),
               std::runtime_error);
  EXPECT_EQ(collected, (std::vector<std::size_t>{0, 1, 2, 0}));
  EXPECT_EQ(running, 0U);
  EXPECT_FALSE(observing);
}

template <typename T> class ObservableGridModels : public ::testing::Test {};
TYPED_TEST_SUITE(ObservableGridModels, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(ObservableGridModels, LiebLinigerPreparedBackground)
{
  using Real = TypeParam;
  namespace model = bethe::lieb_liniger::thermo;
  Real const pi = Real{4} * std::atan(Real{1});
  for (std::size_t iterations : {160, 0})
  {
    model::Options<Real> options;
    options.max_momentum_iterations = iterations;
    model::Solver<Real> const solver(Real{4}, Real{1}, options);
    ASSERT_TRUE(solver.background().converged);
    check_schedulers(
        10,
        [&](std::size_t i) {
          auto const branch = i / 5 == 0 ? model::Branch::type_i : model::Branch::type_ii;
          return solver.at_momentum(branch, Real{2} * pi * Real(i % 5) / Real{4});
        },
        iterations != 0);
  }
}

TYPED_TEST(ObservableGridModels, DopedHubbardPreparedBackground)
{
  using Real = TypeParam;
  namespace model = bethe::hubbard::thermo;
  for (std::size_t iterations : {160, 0})
  {
    model::DopedOptions<Real> options;
    options.max_iterations = iterations;
    model::DopedSolver<Real> const solver(Real{4}, Real{0.5}, options);
    ASSERT_TRUE(solver.background().converged);
    check_schedulers(
        15,
        [&](std::size_t i) {
          auto const branch = std::array{model::DopedBranch::spinon, model::DopedBranch::holon,
                                         model::DopedBranch::charge_particle}[i / 5];
          auto const [lo, hi] = solver.momentum_range(branch);
          Real const p = i % 5 == 4 ? hi : lo + (hi - lo) * Real(i % 5) / Real{4};
          return solver.at_momentum(branch, p);
        },
        iterations != 0);
  }
}

TYPED_TEST(ObservableGridModels, HalfFilledHubbardQuadrature)
{
  using Real = TypeParam;
  namespace model = bethe::hubbard::thermo;
  Real const pi = Real{4} * std::atan(Real{1});
  for (std::size_t budget : {1000000, 0})
  {
    model::Options<Real> options;
    options.max_evaluations = budget;
    check_schedulers(
        9,
        [&](std::size_t i) {
          auto const branch = std::array{model::Branch::spinon, model::Branch::holon, model::Branch::antiholon}[i / 3];
          return model::dispersion(branch, Real{4}, pi * Real(i % 3) / Real{2}, model::Convention::symmetric, options);
        },
        budget != 0);
  }
}
