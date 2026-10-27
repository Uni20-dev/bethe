// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/heisenberg_excitations.hpp>
#include <gtest/gtest.h>
#include <uni20/async/tbb_scheduler.hpp>
#if UNI20_HAS_MPLAPACK_MPREAL
#include "two_stage_dsf_reference.hpp"
#endif

namespace
{
using R = uni20::mpreal;
using uni20::Precision;
namespace model = bethe::heisenberg;
using namespace uni20::async;

template <typename State> void check_state(State const& state, Precision p)
{
  EXPECT_TRUE(state.converged);
  EXPECT_EQ(state.energy.precision(), p);
  EXPECT_EQ(state.residual_norm.precision(), p);
  for (auto const& root : state.rapidities)
    EXPECT_EQ(root.precision(), p);
  if constexpr (requires { state.momentum; }) EXPECT_EQ(state.momentum.precision(), p);
}
} // namespace

TEST(AsyncMpfr, ConcurrentPeriodicAndOpenScansOwnTheirPrecision)
{
  TbbScheduler scheduler(4);
  auto const ambient = mpfr_get_default_prec();
  auto check = [&]<bool open>() {
    using State = std::conditional_t<open, model::open::RealState<R>, model::RealState<R>>;
    using Result = model::RealExcitationScan<State>;
    std::vector<Async<Result>> pending;
    for (int bits : {128, 333})
    {
      bethe::SolverOptions<R> options(Precision::bits(bits));
      bethe::RealExcitationOptions scan{100, 100, {&scheduler, 128, 3}};
      if constexpr (open)
        pending.push_back(model::open::real_excitations_async<R>(8, uni20::half_int{1}, scan, options));
      else
        pending.push_back(model::real_excitations_async<R>(8, uni20::half_int{1}, scan, options));
    } // options are gone; both operations were submitted before waiting
    for (std::size_t i = 0; i < pending.size(); ++i)
    {
      auto p = Precision::bits(i == 0 ? 128 : 333);
      auto const& result = pending[i].get_wait(scheduler);
      ASSERT_TRUE(result.converged());
      ASSERT_EQ(result.levels.size(), 10u);
      check_state(result.ground_state, p);
      for (auto const& level : result.levels)
      {
        check_state(level.state, p);
        ASSERT_TRUE(level.gap);
        EXPECT_EQ(level.gap->precision(), p);
        auto expected = [&] {
          bethe::SolverOptions<R> options(p);
          if constexpr (open)
            return model::open::solve_real<R>(8, level.state.quantum_numbers, options);
          else
            return model::solve_real<R>(8, level.state.quantum_numbers, options);
        }();
        EXPECT_EQ(level.state.rapidities, expected.rapidities);
        EXPECT_EQ(level.state.energy, expected.energy);
      }
    }
  };
  check.template operator()<false>();
  check.template operator()<true>();
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

#if UNI20_HAS_MPLAPACK_MPREAL
TEST(AsyncMpfr, FusedSpectraMatchPreparedTwoStageAtConcurrentPrecisions)
{
  TbbScheduler scheduler(4);
  DebugScheduler serial;
  auto const ambient = mpfr_get_default_prec();
  std::vector<Async<model::TwoSpinonStructureFactor<R>>> pending;
  for (int bits : {128, 333})
    pending.push_back(model::two_spinon_structure_factor_async<R>(8, bethe::SolverOptions<R>(Precision::bits(bits)),
                                                                  100, {&scheduler, 128, 3}));
  for (std::size_t i = 0; i < pending.size(); ++i)
  {
    auto p = Precision::bits(i == 0 ? 128 : 333);
    auto const& result = pending[i].get_wait(scheduler);
    auto expected = dsf_reference::two_spinon_structure_factor<R>(8, bethe::SolverOptions<R>(p), 100, {&serial, 3});
    ASSERT_TRUE(result.converged());
    ASSERT_TRUE(expected.converged());
    ASSERT_EQ(result.lines.size(), expected.lines.size());
    check_state(result.scan.ground_state, p);
    for (std::size_t j = 0; j < result.lines.size(); ++j)
    {
      check_state(result.scan.levels[j].state, p);
      auto const& a = result.lines[j];
      auto const& b = expected.lines[j];
      EXPECT_EQ(a.momentum_index, b.momentum_index);
      EXPECT_EQ(a.momentum, b.momentum);
      EXPECT_EQ(a.gap, b.gap);
      EXPECT_EQ(a.weight, b.weight);
      EXPECT_EQ(a.momentum.precision(), p);
      EXPECT_EQ(a.gap.precision(), p);
      EXPECT_EQ(a.weight.precision(), p);
      auto const& factor = result.form_factors[j];
      ASSERT_TRUE(factor.log_weight && factor.weight);
      EXPECT_EQ(factor.log_weight->precision(), p);
      EXPECT_EQ(factor.pivot_margin.precision(), p);
      EXPECT_EQ(factor.log_weight, expected.form_factors[j].log_weight);
      EXPECT_EQ(factor.pivot_margin, expected.form_factors[j].pivot_margin);
    }
    for (std::size_t q = 0; q < result.moments.size(); ++q)
    {
      EXPECT_EQ(result.moments[q].weight.precision(), p);
      EXPECT_EQ(result.moments[q].first_moment.precision(), p);
      EXPECT_EQ(result.moments[q].weight, expected.moments[q].weight);
      EXPECT_EQ(result.moments[q].first_moment, expected.moments[q].first_moment);
    }
    EXPECT_EQ(result.moments[0].weight, R(0, p));
    EXPECT_EQ(result.integrated_weight.precision(), p);
    EXPECT_EQ(result.weight_fraction, expected.weight_fraction);
    EXPECT_EQ(result.first_moment_fraction, expected.first_moment_fraction);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(AsyncMpfr, FailedReferencesRetainFinitePrecisionZerosAndRejectMixedOptions)
{
  DebugScheduler scheduler;
  auto p = Precision::bits(333);
  bethe::SolverOptions<R> options(p);
  options.max_iterations = 0;
  auto pending = model::two_spinon_structure_factor_async<R>(8, options, 100, {&scheduler});
  auto const& result = pending.get_wait(scheduler);
  EXPECT_FALSE(result.converged());
  EXPECT_TRUE(result.lines.empty());
  EXPECT_EQ(result.weight_fraction.precision(), p);
  for (auto const& m : result.moments)
  {
    EXPECT_EQ(m.weight, R(0, p));
    EXPECT_EQ(m.weight.precision(), p);
    EXPECT_EQ(m.first_moment.precision(), p);
  }
  options.residual_tolerance = options.residual_tolerance.at(Precision::bits(128));
  EXPECT_THROW((void)model::two_spinon_structure_factor_async<R>(8, options, 100, {&scheduler}), std::invalid_argument);
}
#endif
