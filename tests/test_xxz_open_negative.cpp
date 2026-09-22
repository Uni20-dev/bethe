// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "exact_spectrum.hpp"
#include "test_support.hpp"
#include <bethe/xxz_excitations.hpp>

namespace
{
namespace model = bethe::xxz::open;
using Convention = model::GroundResidualConvention;
using Status = model::GroundSolveStatus;
using uni20::half_int;
template <typename Real> class XXZOpenNegative : public ::testing::Test {};
TYPED_TEST_SUITE(XXZOpenNegative, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXZOpenNegative, PublicSectorsAndGroundAgainstED)
{
  using Real = TypeParam;
  for (unsigned n = 2; n <= 9; ++n)
    for (Real d : {-Real{1} / Real{100}, -Real{1} / Real{2}, -Real{9} / Real{10}, -Real{99} / Real{100}})
    {
      SCOPED_TRACE(::testing::Message() << "N=" << n << " Delta=" << uni20::format_scalar(d));
      auto const sectors = model::sector_ground_states<Real>(n, d);
      auto const ground = model::ground_state<Real>(n, d);
      ASSERT_EQ(sectors.size(), n + 1);
      ASSERT_TRUE(ground.converged);
      EXPECT_EQ(ground.sz, uni20::from_twice(std::int64_t(n % 2)));
      EXPECT_EQ(ground.energy, sectors[(n + 1) / 2].energy);
      for (unsigned index = 0; index <= n; ++index)
      {
        auto const& state = sectors[index];
        auto const sz = uni20::from_twice(2 * std::int64_t(index) - std::int64_t(n));
        auto const m = std::min(index, n - index);
        ASSERT_TRUE(state.converged);
        EXPECT_EQ(state.status, Status::converged);
        EXPECT_EQ(state.residual_convention, Convention::negative_rank_scaled);
        EXPECT_EQ(state.delta, d);
        EXPECT_EQ(state.root_delta, d);
        EXPECT_EQ(state.sz, sz);
        EXPECT_EQ(state.spin_reversed, sz.twice() < 0);
        EXPECT_FALSE(state.boundary_root);
        EXPECT_EQ(state.rapidities.size(), m);
        EXPECT_EQ(state.log_rapidities.size(), m);
        EXPECT_EQ(state.quantum_numbers, model::sector_ground_quantum_numbers(n, sz));
        EXPECT_EQ(state.rapidities, sectors[n - index].rapidities);
        EXPECT_EQ(state.log_rapidities, sectors[n - index].log_rapidities);
        EXPECT_EQ(state.energy, sectors[n - index].energy);
        auto const direct = model::sector_ground_state<Real>(n, d, sz);
        EXPECT_EQ(state.energy, direct.energy);
        auto const ed = test_support::exact_spectrum(n, m, 0, false, static_cast<double>(d));
        EXPECT_REAL_NEAR(static_cast<double>(state.energy), ed.front(), 3e-11);
        EXPECT_LE(ground.energy, state.energy + Real{256} * Real(n) * uni20::numeric_limits<Real>::epsilon());
      }
    }
}

// Evaluate the reported equations independently of NegativeGroundSystem,
// directly from the public lambda array, including failed iterates.
template <typename Real> void check_scaled_state(std::size_t n, model::GroundState<Real> const& state)
{
  Real const d = state.delta, s = std::sqrt((Real{1} + d) / (Real{1} - d));
  Real const b = std::sqrt(Real{1} + d) * std::sqrt(Real{1} - d);
  Real norm = Real{0}, energy = Real(n - 1) * d / Real{4};
  auto phase = [&](Real x) { return std::copysign(std::atan2(b, -d * std::abs(std::tanh(x))) / s, x); };
  for (std::size_t i = 0; i < state.log_rapidities.size(); ++i)
  {
    Real const x = state.log_rapidities[i], v = std::tanh(x), z = s * v, sech = Real{1} / std::cosh(x);
    EXPECT_EQ(state.rapidities[i], z);
    EXPECT_GT(x, Real{0});
    if (i) EXPECT_GT(x, state.log_rapidities[i - 1]);
    Real residual = Real(2 * n) * std::atan(z) / s - Real{2} * std::atan(s / v) / s;
    for (std::size_t j = 0; j < state.log_rapidities.size(); ++j)
      if (i != j) residual -= phase(x - state.log_rapidities[j]) + phase(x + state.log_rapidities[j]);
    norm = std::max(norm, std::abs(residual) / Real(n));
    energy -= (Real{1} + d) * sech * sech / (Real{1} + z * z);
  }
  Real const allowance = Real{64} * Real(n) * uni20::numeric_limits<Real>::epsilon();
  EXPECT_REAL_NEAR(norm, state.residual_norm, allowance);
  EXPECT_REAL_NEAR(energy, state.energy, allowance);
}

TYPED_TEST(XXZOpenNegative, NativePrecisionAndEndpoint)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real d : {-std::sqrt(eps), -Real{1} / Real{3}, -Real{999} / Real{1000}, -Real{1} + Real{128} * eps})
  {
    for (std::size_t n : {2, 3, 32, 33})
    {
      SCOPED_TRACE(::testing::Message() << "N=" << n << " Delta=" << uni20::format_scalar(d));
      auto const state = model::ground_state<Real>(n, d);
      ASSERT_TRUE(state.converged);
      EXPECT_LE(state.residual_norm, bethe::SolverOptions<Real>{}.residual_tolerance);
      check_scaled_state(n, state);
      if (n == 2) EXPECT_REAL_NEAR(state.energy, -Real{1} / Real{2} - d / Real{4}, Real{64} * eps);
      if (n == 3) EXPECT_REAL_NEAR(state.energy, -(d + std::sqrt(d * d + Real{8})) / Real{4}, Real{128} * eps);
      if (d == -std::sqrt(eps))
      {
        auto const free = model::ground_state<Real>(n, Real{0});
        ASSERT_TRUE(free.converged);
        EXPECT_LE(std::abs(state.energy - free.energy),
                  Real(n - 1) * std::abs(d) / Real{4} + Real{256} * Real(n) * eps);
      }
    }
  }
}

TYPED_TEST(XXZOpenNegative, BudgetsAndConventionTransitions)
{
  using Real = TypeParam;
  Real const d = -Real{9} / Real{10};
  for (std::size_t budget : {0, 1})
  {
    auto const state = model::ground_state<Real>(7, d, {.max_iterations = budget});
    EXPECT_FALSE(state.converged);
    EXPECT_EQ(state.status, Status::iteration_limit);
    EXPECT_EQ(state.iterations, budget);
    EXPECT_EQ(state.root_delta, d);
    EXPECT_EQ(state.residual_convention, Convention::negative_rank_scaled);
    EXPECT_GT(state.residual_norm, bethe::SolverOptions<Real>{}.residual_tolerance);
    check_scaled_state(7, state);
  }
  auto const vacuum =
      model::sector_ground_state<Real>(7, d, uni20::from_twice(std::int64_t{-7}), {.max_iterations = 0});
  EXPECT_TRUE(vacuum.converged);
  EXPECT_EQ(vacuum.energy, Real{6} * d / Real{4});
  EXPECT_EQ(vacuum.iterations, 0);
  EXPECT_TRUE(vacuum.log_rapidities.empty());
  EXPECT_EQ(vacuum.residual_convention, Convention::negative_rank_scaled);
  for (Real nonnegative : {Real{0}, Real{1}, Real{3}})
  {
    auto const state = model::ground_state<Real>(4, nonnegative);
    ASSERT_TRUE(state.converged);
    EXPECT_TRUE(state.log_rapidities.empty());
    EXPECT_EQ(state.residual_convention,
              nonnegative > Real{1} ? Convention::massive_regularized : Convention::logarithmic_phase);
  }
}

TYPED_TEST(XXZOpenNegative, DomainAndExcitationSeparation)
{
  using Real = TypeParam;
  Real const d = -Real{1} / Real{2};
  for (Real invalid :
       {-Real{1}, -Real{2}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    EXPECT_THROW((void)model::ground_state<Real>(4, invalid), std::invalid_argument);
    EXPECT_THROW((void)model::sector_ground_state<Real>(4, invalid, half_int{0}), std::invalid_argument);
    EXPECT_THROW((void)model::sector_ground_states<Real>(4, invalid), std::invalid_argument);
  }
  for (Real tol : {Real{0}, -Real{1}, uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW((void)model::ground_state<Real>(4, d, {.residual_tolerance = tol}), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state<Real>(1, d), std::invalid_argument);
  EXPECT_THROW((void)model::sector_ground_state<Real>(5, d, half_int{0}), std::invalid_argument);
  EXPECT_THROW((void)model::sector_ground_state<Real>(4, d, half_int{3}), std::invalid_argument);
  EXPECT_THROW((void)model::real_excitations<Real>(4, d, half_int{1}), std::invalid_argument);
  EXPECT_THROW((void)model::solve_real<Real>(4, d, {}), std::invalid_argument);
  EXPECT_THROW((void)model::real_quantum_number_window(4, d, half_int{1}), std::invalid_argument);
}
} // namespace
