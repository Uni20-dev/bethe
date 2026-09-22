// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "exact_spectrum.hpp"
#include "test_support.hpp"
#include <bethe/xxz_excitations.hpp>
#include <type_traits>

namespace
{
namespace model = bethe::xxz;
using Convention = model::GroundResidualConvention;
using Status = model::GroundSolveStatus;
using uni20::half_int;
static_assert(!std::is_convertible_v<model::GroundState<double>, model::RealState<double>>);
static_assert(!std::is_convertible_v<model::RealState<double>, model::GroundState<double>>);
static_assert(std::is_same_v<decltype(model::ground_state(4, 0.5)), model::GroundState<double>>);
static_assert(
    std::is_same_v<decltype(model::real_excitations(4, 0.5, half_int{1}).ground_state), model::RealState<double>>);
template <typename Real> class XXZPeriodicNegative : public ::testing::Test {};
TYPED_TEST_SUITE(XXZPeriodicNegative, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXZPeriodicNegative, PublicSectorsAndMomentumAgainstED)
{
  using Real = TypeParam;
  for (unsigned n : {2, 4, 6, 8})
    for (Real d : {-Real{1} / Real{100}, -Real{1} / Real{2}, -Real{9} / Real{10}, -Real{99} / Real{100}})
    {
      SCOPED_TRACE(::testing::Message() << "N=" << n << " Delta=" << uni20::format_scalar(d));
      auto const sectors = model::sector_ground_states<Real>(n, d);
      auto const ground = model::ground_state<Real>(n, d);
      ASSERT_TRUE(ground.converged);
      ASSERT_EQ(sectors.size(), n + 1);
      EXPECT_EQ(ground.sz, half_int{0});
      EXPECT_EQ(ground.energy, sectors[n / 2].energy);
      for (unsigned i = 0; i <= n; ++i)
      {
        auto const& state = sectors[i];
        auto const sz = half_int(std::int64_t(i) - std::int64_t(n / 2));
        auto const m = std::min(i, n - i);
        ASSERT_TRUE(state.converged);
        EXPECT_EQ(state.status, Status::converged);
        EXPECT_EQ(state.residual_convention, Convention::negative_rank_scaled);
        EXPECT_EQ(state.delta, d);
        EXPECT_EQ(state.sz, sz);
        EXPECT_EQ(state.spin_reversed, sz.twice() < 0);
        EXPECT_EQ(state.rapidities.size(), m);
        EXPECT_EQ(state.log_rapidities.size(), m);
        EXPECT_EQ(state.quantum_numbers, model::sector_ground_quantum_numbers(n, sz));
        EXPECT_EQ(state.rapidities, sectors[n - i].rapidities);
        EXPECT_EQ(state.log_rapidities, sectors[n - i].log_rapidities);
        EXPECT_EQ(state.energy, sectors[n - i].energy);
        EXPECT_EQ(state.momentum_index, m % 2 ? n / 2 : 0);
        EXPECT_EQ(state.momentum, sectors[n - i].momentum);
        auto const direct = model::sector_ground_state<Real>(n, d, sz);
        EXPECT_EQ(state.energy, direct.energy);
        auto const ed = test_support::exact_spectrum(n, m, 0, true, static_cast<double>(d));
        EXPECT_REAL_NEAR(static_cast<double>(state.energy), ed.front(), 3e-11);
        // Joint H + a*(T+T^-1)/2 spectrum detects a wrong translation sector.
        auto const joint = test_support::exact_spectrum(n, m, .371, true, static_cast<double>(d));
        double const target = static_cast<double>(state.energy + Real{371} / Real{1000} * std::cos(state.momentum));
        EXPECT_TRUE(std::any_of(joint.begin(), joint.end(), [&](double e) { return std::abs(e - target) < 3e-11; }));
        EXPECT_LE(ground.energy, state.energy + Real{256} * Real(n) * uni20::numeric_limits<Real>::epsilon());
      }
    }
}

template <typename Real> void check_scaled_state(std::size_t n, model::GroundState<Real> const& state)
{
  Real const d = state.delta, s = std::sqrt((Real{1} + d) / (Real{1} - d));
  Real const b = std::sqrt(Real{1} + d) * std::sqrt(Real{1} - d);
  Real norm = Real{0}, energy = Real(n) * d / Real{4}, momentum = Real{0};
  Real const pi = Real{4} * std::atan(Real{1});
  for (std::size_t i = 0; i < state.log_rapidities.size(); ++i)
  {
    Real const x = state.log_rapidities[i], z = s * std::tanh(x), sech = Real{1} / std::cosh(x);
    EXPECT_EQ(state.rapidities[i], z);
    if (i) EXPECT_GT(x, state.log_rapidities[i - 1]);
    Real f = Real(n) * std::atan(z) / s;
    for (std::size_t j = 0; j < state.log_rapidities.size(); ++j)
      if (i != j)
      {
        Real const difference = x - state.log_rapidities[j];
        f -= std::copysign(std::atan2(b, -d * std::abs(std::tanh(difference))) / s, difference);
      }
    norm = std::max(norm, std::abs(f) / Real(n));
    energy -= (Real{1} + d) * sech * sech / (Real{1} + z * z);
    momentum += pi - Real{2} * std::atan(z);
  }
  Real const allowance = Real{64} * Real(n) * uni20::numeric_limits<Real>::epsilon();
  EXPECT_REAL_NEAR(norm, state.residual_norm, allowance);
  EXPECT_REAL_NEAR(energy, state.energy, allowance);
  if (state.converged)
  {
    EXPECT_REAL_NEAR(std::cos(momentum), std::cos(state.momentum), allowance);
    EXPECT_REAL_NEAR(std::sin(momentum), std::sin(state.momentum), allowance);
  }
}

TYPED_TEST(XXZPeriodicNegative, NativeEndpointsAndFailedIterates)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real d : {-std::sqrt(eps), -Real{1} / Real{3}, -Real{99} / Real{100}, -Real{1} + Real{128} * eps})
    for (std::size_t n : {2, 4, 32, 64})
    {
      SCOPED_TRACE(::testing::Message() << "N=" << n << " Delta=" << uni20::format_scalar(d));
      auto const state = model::ground_state<Real>(n, d);
      ASSERT_TRUE(state.converged);
      EXPECT_LE(state.residual_norm, bethe::SolverOptions<Real>{}.residual_tolerance);
      check_scaled_state(n, state);
      if (n == 2) EXPECT_REAL_NEAR(state.energy, -Real{1} - d / Real{2}, Real{128} * eps);
      if (n == 4) EXPECT_REAL_NEAR(state.energy, -(d + std::sqrt(d * d + Real{8})) / Real{2}, Real{256} * eps);
      if (n > 4)
      {
        auto const seed = model::ground_state<Real>(n, d, {.max_iterations = 0});
        EXPECT_FALSE(seed.converged);
        EXPECT_EQ(seed.status, Status::iteration_limit);
        check_scaled_state(n, seed);
      }
    }
  for (std::size_t budget : {0, 1})
  {
    auto const state = model::ground_state<Real>(8, -Real{1} / Real{2}, {.max_iterations = budget});
    EXPECT_FALSE(state.converged);
    EXPECT_EQ(state.status, Status::iteration_limit);
    EXPECT_EQ(state.iterations, budget);
  }
  auto const vacuum = model::sector_ground_state<Real>(8, -Real{1} / Real{2}, half_int{-4}, {.max_iterations = 0});
  EXPECT_TRUE(vacuum.converged);
  EXPECT_EQ(vacuum.energy, -Real{1});
  EXPECT_EQ(vacuum.momentum_index, 0);
  EXPECT_EQ(vacuum.residual_convention, Convention::negative_rank_scaled);
}

TYPED_TEST(XXZPeriodicNegative, NonnegativeGroundResultPreservesSolverValues)
{
  using Real = TypeParam;
  for (std::size_t n : {5, 6})
    for (Real d : {Real{0}, Real{1} / Real{2}, Real{1}, Real{3}})
      for (std::size_t budget : {std::size_t{0}, std::size_t{10000}})
      {
        auto const sz = uni20::from_twice(std::int64_t(n % 2));
        auto const numbers = model::sector_ground_quantum_numbers(n, sz);
        auto const original = model::detail::solve_validated<Real>(n, d, numbers, {.max_iterations = budget});
        auto const ground = model::ground_state<Real>(n, d, {.max_iterations = budget});
        EXPECT_EQ(ground.energy, original.energy);
        EXPECT_EQ(ground.rapidities, original.rapidities);
        EXPECT_EQ(ground.quantum_numbers, original.quantum_numbers);
        EXPECT_EQ(ground.momentum_index, original.momentum_index);
        EXPECT_EQ(ground.momentum, original.momentum);
        EXPECT_EQ(ground.residual_norm, original.residual_norm);
        EXPECT_EQ(ground.iterations, original.iterations);
        EXPECT_EQ(ground.converged, original.converged);
        EXPECT_EQ(ground.status, original.converged ? Status::converged : Status::iteration_limit);
        EXPECT_EQ(ground.residual_convention, Convention::logarithmic_phase);
        EXPECT_TRUE(ground.log_rapidities.empty());
      }
}

TYPED_TEST(XXZPeriodicNegative, DomainAndExcitationSeparation)
{
  using Real = TypeParam;
  Real const d = -Real{1} / Real{2};
  for (Real invalid :
       {-Real{1}, -Real{2}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    EXPECT_THROW((void)model::ground_state<Real>(4, invalid), std::invalid_argument);
    EXPECT_THROW((void)model::sector_ground_states<Real>(4, invalid), std::invalid_argument);
  }
  for (Real tol : {Real{0}, -Real{1}, uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW((void)model::ground_state<Real>(4, d, {.residual_tolerance = tol}), std::invalid_argument);
  // Reject even an empty odd-ring sector rather than advertise a partial
  // odd-ring contract which a sector scan/global selection cannot fulfill.
  EXPECT_THROW((void)model::ground_state<Real>(5, d), std::invalid_argument);
  EXPECT_THROW((void)model::sector_ground_states<Real>(5, d), std::invalid_argument);
  EXPECT_THROW((void)model::sector_ground_state<Real>(5, d, uni20::from_twice(std::int64_t{5})), std::invalid_argument);
  EXPECT_THROW((void)model::sector_ground_state<Real>(4, d, half_int{3}), std::invalid_argument);
  EXPECT_THROW((void)model::real_excitations<Real>(4, d, half_int{1}), std::invalid_argument);
  EXPECT_THROW((void)model::solve_real<Real>(4, d, {}), std::invalid_argument);
  EXPECT_THROW((void)model::real_excitation_count(4, d, half_int{1}), std::invalid_argument);
}
} // namespace
