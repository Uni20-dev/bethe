// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/su_fermions.hpp>
#include <bethe/takhtajan_babujian.hpp>
#include <bethe/tj.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
template <typename R> auto fields(bethe::takhtajan_babujian::State<R> const& s)
{
  return std::tie(s.sites, s.iterations, s.string_quantum_numbers, s.centers, s.deviations, s.rapidities, s.energy,
                  s.momentum, s.residual_norm, s.phase_residual, s.modulus_residual, s.momentum_index, s.converged,
                  s.status);
}
template <typename R> auto fields(bethe::su_fermions::State<R> const& s)
{
  return std::tie(s.populations, s.component_order, s.particles, s.iterations, s.stages, s.length, s.interaction,
                  s.residual_norm, s.target_residual_norm, s.momentum, s.momentum_index, s.quantum_numbers,
                  s.rapidities, s.level_residuals, s.free_modes, s.energy, s.reached_interaction, s.free, s.converged,
                  s.status);
}
template <typename R> auto fields(bethe::tj::State<R> const& s)
{
  return std::tie(s.sites, s.num_up, s.num_down, s.holes, s.branch, s.quantum_numbers, s.rapidities, s.free_modes,
                  s.energy, s.momentum, s.residual_norm, s.level_residuals, s.momentum_index, s.iterations, s.converged,
                  s.spin_reversed, s.status);
}
template <typename R> struct TBRequest
{
    std::size_t sites;
    bethe::takhtajan_babujian::SolverOptions<R> options{};
    auto operator()() const { return bethe::takhtajan_babujian::ground_state<R>(sites, options); }
};
template <typename R> struct SURequest
{
    std::vector<std::size_t> populations;
    R length, c;
    bethe::su_fermions::SolverOptions<R> options{};
    auto operator()() const { return bethe::su_fermions::ground_state<R>(populations, length, c, options); }
};
template <typename R> struct TJRequest
{
    std::size_t sites, up, down;
    bethe::tj::SolverOptions<R> options{};
    auto operator()() const { return bethe::tj::ground_state<R>(sites, up, down, options); }
};
template <typename Request, typename State>
auto submit(std::vector<Request> requests, std::vector<State>& actual, uni20::async::IAsyncScheduler& scheduler)
{
  auto const count = requests.size();
  return bethe::evaluate_observable_grid_async(
      count, [requests = std::move(requests)](std::size_t i) { return requests[i](); },
      [&actual](std::size_t i, State state) {
        EXPECT_EQ(i, actual.size());
        actual.push_back(std::move(state));
      },
      {.scheduler = &scheduler, .max_in_flight = 2});
}
template <typename Request> auto direct(std::vector<Request> const& requests)
{
  std::vector<std::invoke_result_t<Request>> out;
  for (auto const& request : requests)
    out.push_back(request());
  return out;
}
template <typename State> void compare(std::vector<State> const& actual, std::vector<State> const& expected)
{
  ASSERT_EQ(actual.size(), expected.size());
  for (std::size_t i = 0; i != expected.size(); ++i)
  {
    SCOPED_TRACE(i);
    EXPECT_TRUE(fields(actual[i]) == fields(expected[i]));
  }
}
template <typename R> class AsyncNestedSweeps : public ::testing::Test {};
TYPED_TEST_SUITE(AsyncNestedSweeps, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(AsyncNestedSweeps, OwnsNestedRequestsAndPreservesPartialStates)
{
  using R = TypeParam;
  auto check = [](uni20::async::IAsyncScheduler& scheduler) {
    std::vector<TBRequest<R>> tb{{4}, {6}, {10}, {8, {.max_iterations = 0}}};
    std::vector<SURequest<R>> su{{{0, 0}, R{1}, R{1}},    {{2, 1}, R{2}, R{0}},
                                 {{3, 1, 1}, R{1}, R{1}}, {{1, 0, 3, 1}, R{2}, R{1} / R{10}},
                                 {{3, 1, 1}, R{1}, R{1}}, {{3, 1, 1}, R{1}, R{1}}};
    su[4].options.max_iterations = 0;
    su[5].options.max_stages = 1;
    std::vector<TJRequest<R>> tj{{8, 0, 0}, {7, 4, 0}, {6, 3, 3},
                                 {8, 3, 3}, {8, 1, 3}, {8, 3, 3, {.max_iterations = 0}}};
    auto const t_expected = direct(tb);
    auto const s_expected = direct(su);
    auto const j_expected = direct(tj);
    for (std::size_t i = 0; i != 3; ++i)
      ASSERT_TRUE(t_expected[i].converged);
    EXPECT_FALSE(t_expected.back().converged);
    for (std::size_t i = 0; i != 4; ++i)
      ASSERT_TRUE(s_expected[i].converged);
    EXPECT_FALSE(s_expected[4].energy);
    ASSERT_FALSE(s_expected[5].converged);
    ASSERT_TRUE(s_expected[5].energy);
    ASSERT_TRUE(s_expected[5].reached_interaction);
    EXPECT_GT(*s_expected[5].reached_interaction, s_expected[5].interaction);
    for (std::size_t i = 0; i != 5; ++i)
      ASSERT_TRUE(j_expected[i].converged);
    EXPECT_FALSE(j_expected.back().converged);
    EXPECT_TRUE(j_expected[4].spin_reversed);
    decltype(direct(tb)) t_actual;
    decltype(direct(su)) s_actual;
    decltype(direct(tj)) j_actual;
    // All three grids are admitted before waiting. In the deferred case even
    // the population vectors and their enclosing submission frames are gone.
    auto t_done = submit(std::move(tb), t_actual, scheduler);
    auto s_done = submit(std::move(su), s_actual, scheduler);
    auto j_done = submit(std::move(tj), j_actual, scheduler);
    EXPECT_EQ(t_done.get_wait(scheduler), t_expected.size());
    EXPECT_EQ(s_done.get_wait(scheduler), s_expected.size());
    EXPECT_EQ(j_done.get_wait(scheduler), j_expected.size());
    compare(t_actual, t_expected);
    compare(s_actual, s_expected);
    compare(j_actual, j_expected);
  };
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  check(deferred);
  for (int workers : {1, 4})
  {
    uni20::async::TbbScheduler scheduler(workers);
    check(scheduler);
  }
}

TYPED_TEST(AsyncNestedSweeps, UnsupportedShellDoesNotCancelOtherGrids)
{
  using R = TypeParam;
  auto check = [](uni20::async::IAsyncScheduler& scheduler) {
    std::vector<SURequest<R>> su{{{3, 1, 1}, R{1}, R{1}}, {{2, 1, 1}, R{1}, R{1}}};
    std::vector<TBRequest<R>> tb{{4}, {6}, {8}};
    std::vector<TJRequest<R>> tj{{8, 3, 3}, {7, 1, 3}, {6, 3, 3}};
    decltype(direct(su)) s_rows;
    decltype(direct(tb)) t_rows;
    decltype(direct(tj)) j_rows;
    auto s_done = submit(su, s_rows, scheduler);
    auto t_done = submit(tb, t_rows, scheduler);
    auto j_done = submit(tj, j_rows, scheduler);
    EXPECT_THROW((void)s_done.get_wait(scheduler), std::invalid_argument);
    EXPECT_EQ(t_done.get_wait(scheduler), tb.size());
    EXPECT_EQ(j_done.get_wait(scheduler), tj.size());
    compare(t_rows, direct(tb));
    compare(j_rows, direct(tj));
  };
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  check(deferred);
  for (int workers : {1, 4})
  {
    uni20::async::TbbScheduler scheduler(workers);
    check(scheduler);
  }
}
} // namespace
