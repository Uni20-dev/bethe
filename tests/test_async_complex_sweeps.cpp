// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/tasep.hpp>
#include <bethe/xyz.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
template <typename R> auto fields(bethe::tasep::GapState<R> const& s)
{
  return std::tie(s.sites, s.particles, s.effective_particles, s.iterations, s.seed_iterations, s.rate, s.residual_norm,
                  s.roots, s.log_y, s.eigenvalue, s.gap, s.converged, s.status);
}
template <typename R> auto fields(bethe::xyz::GroundState<R> const& s)
{
  return std::tie(s.sites, s.iterations, s.eta, s.t, s.residual_norm, s.momentum, s.exchange.x, s.exchange.y,
                  s.exchange.z, s.quantum_numbers, s.roots, s.energy, s.converged, s.status);
}
template <typename R> struct TasepRequest
{
    std::size_t sites, particles;
    R rate;
    bethe::tasep::Options<R> options{};
    auto operator()() const { return bethe::tasep::relaxation_gap(sites, particles, rate, options); }
};
template <typename R> struct XyzRequest
{
    std::size_t sites;
    R eta, t;
    bethe::xyz::SolverOptions<R> options{};
    auto operator()() const { return bethe::xyz::ground_state(sites, eta, t, options); }
};

template <typename Request, typename State>
auto submit(std::vector<Request> requests, std::vector<State>& rows, uni20::async::IAsyncScheduler& scheduler)
{
  auto const count = requests.size();
  // The submission-local request vector disappears before deferred work runs.
  return bethe::evaluate_observable_grid_async(
      count, [requests = std::move(requests)](std::size_t i) { return requests[i](); },
      [&rows](std::size_t i, State state) {
        EXPECT_EQ(i, rows.size());
        rows.push_back(std::move(state));
      },
      {.scheduler = &scheduler, .max_in_flight = 2});
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

template <typename R> class AsyncComplexSweeps : public ::testing::Test {};
TYPED_TEST_SUITE(AsyncComplexSweeps, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(AsyncComplexSweeps, OwnsRequestsAndPreservesEveryStateField)
{
  using R = TypeParam;
  auto check = [](uni20::async::IAsyncScheduler& scheduler) {
    std::vector<TasepRequest<R>> tasep{{8, 0, R{1}}, {8, 8, R{1}},  {5, 1, R{3}}, {5, 2, R{1}},
                                       {8, 5, R{2}}, {12, 6, R{1}}, {8, 3, R{1}}, {8, 3, R{1}}};
    tasep[6].options.max_seed_iterations = 0;
    tasep[7].options.max_iterations = 0;
    std::vector<XyzRequest<R>> xyz{{2, R{1} / R{5}, R{1} / R{2}},
                                   {8, R{1} / R{5}, R{2}},
                                   {12, R{1} / R{2}, R{3} / R{4}},
                                   {8, R{4} / R{5}, R{1} / R{2}},
                                   {8, R{1} / R{5}, R{2}}};
    xyz.back().options.max_iterations = 0;
    std::vector<bethe::tasep::GapState<R>> t_expected, t_actual;
    std::vector<bethe::xyz::GroundState<R>> x_expected, x_actual;
    for (auto const& request : tasep)
      t_expected.push_back(request());
    for (auto const& request : xyz)
      x_expected.push_back(request());
    for (std::size_t i = 0; i != 6; ++i)
      ASSERT_TRUE(t_expected[i].converged);
    EXPECT_EQ(t_expected[6].status, bethe::tasep::Status::seed_limit);
    EXPECT_EQ(t_expected[7].status, bethe::tasep::Status::iteration_limit);
    for (std::size_t i = 0; i != 4; ++i)
      ASSERT_TRUE(x_expected[i].converged);
    EXPECT_FALSE(x_expected.back().converged);
    // Both grids share one scheduler. Captured controls and inputs are owning;
    // no ambient precision, worker reservation or scheduler-wide drain is used.
    auto t_done = submit(std::move(tasep), t_actual, scheduler);
    auto x_done = submit(std::move(xyz), x_actual, scheduler);
    EXPECT_EQ(t_done.get_wait(scheduler), t_expected.size());
    EXPECT_EQ(x_done.get_wait(scheduler), x_expected.size());
    compare(t_actual, t_expected);
    compare(x_actual, x_expected);
  };
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  check(deferred);
  for (int workers : {1, 4})
  {
    uni20::async::TbbScheduler scheduler(workers);
    check(scheduler);
  }
}

TYPED_TEST(AsyncComplexSweeps, InvalidRequestDoesNotCancelOtherModel)
{
  using R = TypeParam;
  auto check = [](uni20::async::IAsyncScheduler& scheduler, bool fail_tasep) {
    std::vector<TasepRequest<R>> tasep{{8, 3, R{1}}, {8, fail_tasep ? 9u : 4u, R{1}}, {6, 2, R{1}}};
    std::vector<XyzRequest<R>> xyz{{4, R{1} / R{3}, R{2}}, {fail_tasep ? 6u : 5u, R{1} / R{2}, R{2}}};
    std::vector<bethe::tasep::GapState<R>> t_rows;
    std::vector<bethe::xyz::GroundState<R>> x_rows;
    auto t_done = submit(tasep, t_rows, scheduler), x_done = submit(xyz, x_rows, scheduler);
    if (fail_tasep)
    {
      EXPECT_THROW((void)t_done.get_wait(scheduler), std::invalid_argument);
      EXPECT_EQ(x_done.get_wait(scheduler), xyz.size());
      ASSERT_EQ(x_rows.size(), xyz.size());
      for (std::size_t i = 0; i != xyz.size(); ++i)
        EXPECT_TRUE(fields(x_rows[i]) == fields(xyz[i]()));
    }
    else
    {
      EXPECT_THROW((void)x_done.get_wait(scheduler), std::invalid_argument);
      EXPECT_EQ(t_done.get_wait(scheduler), tasep.size());
      ASSERT_EQ(t_rows.size(), tasep.size());
      for (std::size_t i = 0; i != tasep.size(); ++i)
        EXPECT_TRUE(fields(t_rows[i]) == fields(tasep[i]()));
    }
  };
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  check(deferred, false);
  check(deferred, true);
  for (int workers : {1, 4})
  {
    uni20::async::TbbScheduler scheduler(workers);
    check(scheduler, false);
    check(scheduler, true);
  }
}
} // namespace
