// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/ladder.hpp>
#include <bethe/prepared_observable_grid.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
namespace model = bethe::ladder;
template <typename R> auto fields(model::State<R> const& s)
{
  auto const& b = s.highest_weight;
  return std::tie(s.rungs, s.singlets, s.iterations, s.branches, s.tableaux, s.populations, s.rung_coupling,
                  s.magnetic_field, s.magnetization, s.energy, s.converged, s.descendant, s.analytic, s.status, b.shape,
                  b.labels, b.rapidities, b.energy, b.residual, b.momentum_index, b.iterations, b.converged, b.status);
}
template <typename R> void compare(model::SectorScan<R> const& actual, model::SectorScan<R> const& expected)
{
  EXPECT_EQ(std::tie(actual.iterations, actual.branches, actual.tableaux, actual.complete, actual.status),
            std::tie(expected.iterations, expected.branches, expected.tableaux, expected.complete, expected.status));
  ASSERT_EQ(actual.sectors.size(), expected.sectors.size());
  for (std::size_t i = 0; i != actual.sectors.size(); ++i)
  {
    SCOPED_TRACE(i);
    EXPECT_TRUE(fields(actual.sectors[i]) == fields(expected.sectors[i]));
  }
}
template <typename R> class PreparedLadder : public ::testing::Test {};
TYPED_TEST_SUITE(PreparedLadder, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(PreparedLadder, PreservesCompleteAndPartialScansAndAnalyticBypasses)
{
  using R = TypeParam;
  std::vector<model::SolverOptions<R>> controls(7);
  controls[1].max_iterations = 0;
  controls[2].max_iterations = 1;
  controls[3].max_branches = 0;
  controls[4].max_branches = 1;
  controls[5].max_branches = 3;
  controls[6].max_branches = controls[6].max_iterations = 0;
  R const huge = R{1} / (uni20::numeric_limits<R>::epsilon() * uni20::numeric_limits<R>::epsilon());
  std::array<std::pair<R, R>, 8> parameters{{{R{0}, R{0}},
                                             {R{1}, R{1} / R{8}},
                                             {R{-1}, R{-1} / R{8}},
                                             {R{5}, R{1}},
                                             {R{-2}, R{-9}},
                                             {R{3}, R{9}},
                                             {huge, huge},
                                             {R{0}, huge}}};
  for (std::size_t rungs : {2u, 5u, 6u})
    for (auto const& options : controls)
    {
      model::PreparedSectors<R> prepared(rungs, options);
      auto const reference = model::sector_ground_states<R>(rungs, R{0}, R{0}, options);
      EXPECT_EQ(prepared.rungs(), rungs);
      EXPECT_EQ(prepared.branches(), reference.branches);
      EXPECT_EQ(prepared.iterations(), reference.iterations);
      EXPECT_EQ(prepared.tableaux(), reference.tableaux);
      EXPECT_EQ(prepared.complete(), reference.complete);
      EXPECT_EQ(prepared.status(), reference.status);
      for (auto const& [rung, field] : parameters)
      {
        auto const state = prepared.ground_state(rung, field);
        EXPECT_TRUE(fields(state) == fields(model::ground_state(rungs, rung, field, options)));
        compare(prepared.sector_ground_states(rung, field), model::sector_ground_states(rungs, rung, field, options));
      }
    }
}

TYPED_TEST(PreparedLadder, OwnsControlsAndRootsAndRejectsInvalidQueries)
{
  using R = TypeParam;
  model::SolverOptions<R> options;
  options.max_branches = 3;
  auto prepared = [&] {
    model::PreparedSectors<R> original(6, options);
    auto copy = original;
    options.max_branches = 0;
    return copy;
  }(); // the original object and its constructor's options are no longer authoritative
  EXPECT_EQ(prepared.options().max_branches, 3u);
  auto state = prepared.ground_state(R{1});
  ASSERT_TRUE(state.energy);
  state.highest_weight.rapidities.clear();
  state.highest_weight.labels.clear();
  model::SolverOptions<R> reference_options;
  reference_options.max_branches = 3;
  EXPECT_TRUE(fields(prepared.ground_state(R{1})) == fields(model::ground_state<R>(6, R{1}, reference_options)));
  R const nan = uni20::numeric_limits<R>::quiet_NaN(), inf = uni20::numeric_limits<R>::infinity();
  EXPECT_THROW((void)prepared.ground_state(nan), std::invalid_argument);
  EXPECT_THROW((void)prepared.ground_state(R{5}, inf), std::invalid_argument);
  EXPECT_THROW((void)prepared.sector_ground_states(inf), std::invalid_argument);
  EXPECT_THROW((void)prepared.sector_ground_states(R{0}, nan), std::invalid_argument);
  EXPECT_THROW((void)model::PreparedSectors<R>(1), std::invalid_argument);
  options.residual_tolerance = R{0};
  EXPECT_THROW((void)model::PreparedSectors<R>(6, options), std::invalid_argument);
}

TYPED_TEST(PreparedLadder, TwoDependentGridsOutlivePreparationHandle)
{
  using R = TypeParam;
  auto check = [](uni20::async::IAsyncScheduler& scheduler) {
    std::array<std::vector<model::State<R>>, 2> states;
    std::size_t preparations = 0;
    auto done = [&] {
      model::SolverOptions<R> options;
      auto background = bethe::prepare_background_async(
          [options, &preparations] {
            ++preparations;
            return model::PreparedSectors<R>(6, options);
          },
          {.scheduler = &scheduler});
      auto submit = [&](std::size_t curve) {
        return bethe::evaluate_observable_grid_async(
            background.read(), 9,
            [curve](model::PreparedSectors<R> const& p, std::size_t i) {
              R const field = (R(i) - R{4}) / R{8};
              return p.ground_state(R(curve), field);
            },
            [&, curve](std::size_t i, auto state) {
              EXPECT_EQ(i, states[curve].size());
              states[curve].push_back(std::move(state));
            },
            {.scheduler = &scheduler, .max_in_flight = 3});
      };
      return std::pair{submit(0), submit(1)};
    }();
    EXPECT_EQ(done.first.get_wait(scheduler), 9u);
    EXPECT_EQ(done.second.get_wait(scheduler), 9u);
    EXPECT_EQ(preparations, 1u);
    for (std::size_t curve = 0; curve != 2; ++curve)
      for (std::size_t i = 0; i != 9; ++i)
        EXPECT_TRUE(fields(states[curve][i]) == fields(model::ground_state<R>(6, R(curve), (R(i) - R{4}) / R{8})));
  };
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  check(deferred);
  for (int workers : {1, 4})
  {
    uni20::async::TbbScheduler scheduler(workers);
    check(scheduler);
  }
}

TYPED_TEST(PreparedLadder, InvalidQueryDoesNotPoisonSharedPreparation)
{
  using R = TypeParam;
  uni20::async::TbbScheduler scheduler(4);
  auto background = bethe::prepare_background_async([] { return model::PreparedSectors<R>(6); }, {&scheduler});
  auto bad = bethe::evaluate_observable_grid_async(
      background.read(), 3,
      [](auto const& p, std::size_t i) {
        return p.ground_state(i == 1 ? uni20::numeric_limits<R>::quiet_NaN() : R{1});
      },
      [](auto, auto) {}, {.scheduler = &scheduler, .max_in_flight = 2});
  std::vector<model::State<R>> states;
  auto good = bethe::evaluate_observable_grid_async(
      background.read(), 3, [](auto const& p, std::size_t i) { return p.ground_state(R(i)); },
      [&states](auto, auto state) { states.push_back(std::move(state)); },
      {.scheduler = &scheduler, .max_in_flight = 2});
  EXPECT_THROW((void)bad.get_wait(scheduler), std::invalid_argument);
  EXPECT_EQ(good.get_wait(scheduler), 3u);
  ASSERT_EQ(states.size(), 3u);
  for (std::size_t i = 0; i != 3; ++i)
    EXPECT_TRUE(fields(states[i]) == fields(model::ground_state<R>(6, R(i))));
}
} // namespace
