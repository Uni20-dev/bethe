// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <array>
#include <atomic>
#include <bethe/asep.hpp>
#include <bethe/prepared_observable_grid.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
template <typename R> auto fields(bethe::asep::GapState<R> const& s)
{
  return std::tie(s.sites, s.particles, s.effective_particles, s.wave_index, s.iterations, s.continuation_steps,
                  s.seed_iterations, s.seed_newton_iterations, s.right_rate, s.left_rate, s.reached_ratio,
                  s.residual_norm, s.wave_base, s.scaled_roots, s.eigenvalue, s.gap, s.converged, s.analytic, s.status);
}
template <typename R> auto fields(bethe::tasep::GapState<R> const& s)
{
  return std::tie(s.sites, s.particles, s.effective_particles, s.iterations, s.seed_iterations, s.rate, s.residual_norm,
                  s.roots, s.log_y, s.eigenvalue, s.gap, s.converged, s.status);
}
template <typename Real> class ASEPPrepared : public ::testing::Test {};
TYPED_TEST_SUITE(ASEPPrepared, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(ASEPPrepared, EveryFieldMatchesIndependentSolves)
{
  using R = TypeParam;
  using namespace bethe::asep;
  for (std::size_t particles : {0u, 1u, 3u, 4u, 7u, 8u})
  {
    PreparedGap<R> const prepared(8, particles);
    auto const seed = prepared.seed_state();
    for (R q : {R{0}, R{0.3}, R{0.9}, R{1} - R{8} * uni20::numeric_limits<R>::epsilon(), R{1}})
      for (bool reflect : {false, true})
      {
        R const right = reflect ? R{3} * q : R{3}, left = reflect ? R{3} : R{3} * q;
        auto const direct = relaxation_gap(8, particles, right, left);
        auto const query = prepared.at_rates(right, left);
        ASSERT_TRUE(direct.converged) << particles << ' ' << uni20::format_real(q);
        EXPECT_TRUE(fields(query) == fields(direct));
      }
    EXPECT_TRUE(fields(prepared.seed_state()) == fields(seed));
  }
}

TYPED_TEST(ASEPPrepared, NumericalFailuresAndAnalyticBypass)
{
  using R = TypeParam;
  using namespace bethe::asep;
  for (int mode = 0; mode != 4; ++mode)
  {
    Options<R> options;
    if (mode == 0) options.seed_options.max_seed_iterations = 0;
    if (mode == 1) options.seed_options.max_iterations = 0;
    if (mode == 2) options.max_iterations = 0;
    if (mode == 3) options.max_continuation_steps = 0;
    PreparedGap<R> const prepared(8, 3, options);
    for (R q : {R{0}, R{0.5}, R{1}})
    {
      auto const actual = prepared.at_rates(R{1}, q);
      auto const expected = relaxation_gap(8, 3, R{1}, q, options);
      EXPECT_TRUE(fields(actual) == fields(expected)) << mode;
      if (q == R{0.5}) EXPECT_FALSE(actual.converged);
      if (q == R{1}) EXPECT_TRUE(actual.converged); // independent analytic endpoint
    }
  }
  EXPECT_THROW((PreparedGap<R>(8, 9)), std::invalid_argument);
  EXPECT_THROW((PreparedGap<R>(8, 3, {.tolerance = R{0}})), std::invalid_argument);
  Options<R> small;
  small.seed_options.max_sites = 7;
  EXPECT_THROW((PreparedGap<R>(8, 3, small)), std::length_error);
  PreparedGap<R> const prepared(8, 3);
  EXPECT_THROW((void)prepared.at_rates(R{0}, R{0}), std::invalid_argument);
  EXPECT_THROW((void)prepared.at_rates(R{-1}, R{1}), std::invalid_argument);
  EXPECT_THROW((void)prepared.at_rates(uni20::numeric_limits<R>::infinity(), R{1}), std::invalid_argument);
}

TYPED_TEST(ASEPPrepared, TwoAsyncReadersShareOwnedPreparation)
{
  using R = TypeParam;
  using State = bethe::asep::GapState<R>;
  using Prepared = bethe::asep::PreparedGap<R>;
  using namespace uni20::async;
  auto check = [&](IAsyncScheduler& scheduler, bool fail_seed) {
    std::array<std::vector<State>, 2> actual, expected;
    std::atomic<unsigned> preparations{0};
    std::atomic<bethe::tasep::GapState<R> const*> observed_seed{nullptr};
    bethe::asep::Options<R> options;
    if (fail_seed) options.seed_options.max_seed_iterations = 0;
    constexpr std::size_t count = 5;
    for (std::size_t b = 0; b != 2; ++b)
      for (std::size_t i = 0; i != count; ++i)
        expected[b].push_back(bethe::asep::relaxation_gap(8, 3, R(b + 1), R(b + 1) * R(i) / R(count - 1), options));
    auto completions = [&] {
      // Preparation and both consumers are submitted before observing any
      // result. All submission-local handles and controls may then disappear.
      auto background = bethe::prepare_background_async(
          [options, &preparations] {
            ++preparations;
            return Prepared(8, 3, options);
          },
          {.scheduler = &scheduler});
      auto submit = [&](std::size_t b) {
        return bethe::evaluate_observable_grid_async(
            background.read(), count,
            [b, count, &observed_seed](Prepared const& prepared, std::size_t i) {
              auto* seen = static_cast<bethe::tasep::GapState<R> const*>(nullptr);
              if (!observed_seed.compare_exchange_strong(seen, &prepared.seed_state()))
                EXPECT_EQ(seen, &prepared.seed_state());
              return prepared.at_rates(R(b + 1), R(b + 1) * R(i) / R(count - 1));
            },
            [&, b](std::size_t i, State state) {
              EXPECT_EQ(i, actual[b].size());
              actual[b].push_back(std::move(state));
            },
            {.scheduler = &scheduler, .max_in_flight = 2});
      };
      auto first = submit(0), second = submit(1);
      return std::pair{std::move(first), std::move(second)};
    }();
    EXPECT_EQ(completions.first.get_wait(scheduler), count);
    EXPECT_EQ(completions.second.get_wait(scheduler), count);
    EXPECT_EQ(preparations.load(), 1u);
    for (std::size_t b = 0; b != 2; ++b)
    {
      ASSERT_EQ(actual[b].size(), count);
      for (std::size_t i = 0; i != count; ++i)
        EXPECT_TRUE(fields(actual[b][i]) == fields(expected[b][i]));
    }
  };
  DebugScheduler deferred({.order = DebugSchedulerOrder::reverse});
  check(deferred, false);
  check(deferred, true);
  for (int workers : {1, 4})
  {
    TbbScheduler scheduler(workers);
    check(scheduler, false);
    check(scheduler, true);
  }
}
} // namespace
