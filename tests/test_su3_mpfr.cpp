// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/heisenberg.hpp>
#include <bethe/observable_grid.hpp>
#include <bethe/su3.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using Real = uni20::mpreal;
namespace model = bethe::su3;
namespace math = uni20::math;

void check_precision(Real const& value, uni20::Precision p)
{
  ASSERT_TRUE(value.initialized());
  EXPECT_FALSE(value.is_exact());
  EXPECT_EQ(value.precision(), p);
  EXPECT_TRUE(uni20::isfinite(value));
}
#if UNI20_HAS_MPLAPACK_MPREAL
void check_state(model::State<Real> const& state, uni20::Precision p)
{
  check_precision(state.energy, p);
  check_precision(state.momentum, p);
  check_precision(state.residual_norm, p);
  for (std::size_t a = 0; a < 2; ++a)
  {
    check_precision(state.level_residuals[a], p);
    auto const& roots = state.rapidities[a];
    for (std::size_t j = 0; j < roots.size(); ++j)
    {
      check_precision(roots[j], p);
      EXPECT_EQ(roots[j], -roots[roots.size() - 1 - j]);
      if (j) EXPECT_GT(roots[j], roots[j - 1]);
    }
  }
}
#endif

class SU3Mpfr : public ::testing::TestWithParam<int> {};
INSTANTIATE_TEST_SUITE_P(WorkingPrecisions, SU3Mpfr, ::testing::Values(128, 256, 333));

TEST_P(SU3Mpfr, FinitePrecisionFailuresAndInvalidControls)
{
  auto const p = uni20::Precision::bits(GetParam());
  EXPECT_THROW((void)model::SolverOptions<Real>(uni20::Precision::exact()), std::invalid_argument);
  model::SolverOptions<Real> options(p);
  // Scalar-only builds can form the equations, but cannot instantiate the
  // MPFR LU solve. The CLI must reject that unavailable provider explicitly.
  model::detail::GroundSystem<Real> const system(9, {6, 3}, p);
  auto const x = system.seed();
  ASSERT_TRUE(system.physical(x));
  for (auto const& family : system.expand(x))
    for (auto const& root : family)
      check_precision(root, p);
  auto const evaluation = system.evaluate(x);
  for (auto const& residual : evaluation.residual)
    check_precision(residual, p);
  for (auto const& residual : evaluation.levels)
    check_precision(residual, p);
#if UNI20_HAS_MPLAPACK_MPREAL
  options.max_iterations = 0;
  auto const seed = model::ground_state<Real>(9, options);
  check_state(seed, p);
  EXPECT_FALSE(seed.converged);
  EXPECT_EQ(seed.status, model::SolveStatus::iteration_limit);
  EXPECT_EQ(seed.iterations, 0u);
  EXPECT_EQ(seed.momentum, 0);
  EXPECT_EQ(seed.rapidities[1][1], 0);
#endif
  for (auto const& invalid :
       {Real{1}, Real(1, uni20::Precision::bits(GetParam() + 32)), Real(uni20::uninitialized), Real(0, p), Real(-1, p)})
  {
    options.residual_tolerance = invalid;
    EXPECT_THROW(model::detail::check_options(6, options), std::invalid_argument);
  }
}

#if UNI20_HAS_MPLAPACK_MPREAL
TEST_P(SU3Mpfr, AnalyticEnergiesAndNestedRoots)
{
  auto const p = uni20::Precision::bits(GetParam());
  model::SolverOptions<Real> options(p);
  auto const eps = uni20::numeric_limits<Real>::epsilon(p);
  auto const three = model::ground_state<Real>(3, options);
  ASSERT_TRUE(three.converged);
  check_state(three, p);
  EXPECT_REAL_NEAR(three.energy, Real(-3, p), 256 * eps);
  EXPECT_REAL_NEAR(three.rapidities[0][1], 1 / math::sqrt(Real(12, p)), 64 * eps);
  EXPECT_EQ(three.rapidities[1][0], 0);
  auto const six = model::ground_state<Real>(6, options);
  ASSERT_TRUE(six.converged);
  check_state(six, p);
  EXPECT_REAL_NEAR(six.energy, -1 - math::sqrt(Real(13, p)), 1024 * eps);
  for (std::size_t n : {9, 24})
  {
    auto const state = model::ground_state<Real>(n, options);
    ASSERT_TRUE(state.converged);
    check_state(state, p);
    EXPECT_LE(state.residual_norm, options.residual_tolerance);
  }
}

TEST(SU3MpfrShared, AnalyticReferenceResolvesAccuracyBeyondBinary128)
{
  // The independent five-dimensional singlet Hamiltonian has characteristic
  // polynomial E*(E+2)^2*(E^2+2E-12), so its lowest energy is -1-sqrt(13).
  // Compare in 512-bit arithmetic to avoid rounding the reference to each solve.
  auto const high = uni20::Precision::bits(512);
  Real const exact = -1 - math::sqrt(Real(13, high));
  Real previous(1, high);
  for (int bits : {128, 256, 333})
  {
    auto const p = uni20::Precision::bits(bits);
    auto const state = model::ground_state<Real>(6, model::SolverOptions<Real>(p));
    ASSERT_TRUE(state.converged);
    Real const error = math::abs(state.energy.at(high) - exact);
    EXPECT_LT(error, 1024 * uni20::numeric_limits<Real>::epsilon(p).at(high));
    EXPECT_LT(error, previous * Real("1e-18", high));
    previous = error;
  }
}

TEST_P(SU3Mpfr, SU2ReductionAgainstIndependentXXX)
{
  auto const p = uni20::Precision::bits(GetParam());
  model::SolverOptions<Real> options(p);
  auto const eps = uni20::numeric_limits<Real>::epsilon(p);
  model::detail::GroundSystem<Real> const system(12, {6, 0}, p);
  auto const state = model::detail::solve(system, options);
  auto const xxx = bethe::heisenberg::ground_state<Real>(12, options);
  ASSERT_TRUE(state.converged && xxx.converged);
  check_state(state, p);
  EXPECT_REAL_NEAR(state.energy, 2 * xxx.energy + 6, 4096 * eps);
  for (std::size_t j = 0; j < 6; ++j)
    EXPECT_REAL_NEAR(2 * state.rapidities[0][j], xxx.rapidities[j], 4096 * eps);
  EXPECT_TRUE(state.rapidities[1].empty());
  EXPECT_EQ(state.level_residuals[1], 0);
}

TEST(SU3MpfrShared, OwningAsyncSweepKeepsIndependentPrecisionsAndFailures)
{
  auto const ambient = mpfr_get_default_prec();
  uni20::async::TbbScheduler scheduler(2);
  std::vector<model::State<Real>> actual;
  auto evaluate = [bits = std::array{128, 256, 333, 128}](std::size_t i) {
    model::SolverOptions<Real> options(uni20::Precision::bits(bits[i]));
    if (i == 3) options.max_iterations = 0;
    return model::ground_state<Real>(6 + 3 * i, options);
  };
  auto done = bethe::evaluate_observable_grid_async(4, evaluate,
                                                    [&](std::size_t i, auto const& state) {
                                                      EXPECT_EQ(i, actual.size());
                                                      actual.push_back(state);
                                                    },
                                                    {.scheduler = &scheduler, .max_in_flight = 2});
  EXPECT_EQ(done.get_wait(scheduler), 4u);
  ASSERT_EQ(actual.size(), 4u);
  for (std::size_t i = 0; i < actual.size(); ++i)
  {
    auto const expected = evaluate(i);
    check_state(actual[i], expected.energy.precision());
    EXPECT_EQ(actual[i].energy, expected.energy);
    EXPECT_EQ(actual[i].rapidities, expected.rapidities);
    EXPECT_EQ(actual[i].level_residuals, expected.level_residuals);
    EXPECT_EQ(actual[i].status, expected.status);
    EXPECT_EQ(actual[i].iterations, expected.iterations);
    EXPECT_EQ(actual[i].converged, i != 3);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
#endif
} // namespace
