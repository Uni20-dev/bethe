// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "tj_ed.hpp"
#include <bethe/tj.hpp>
#include <complex>

namespace
{
namespace model = bethe::tj;
template <typename Real> class TJ : public ::testing::Test {};
TYPED_TEST_SUITE(TJ, test_support::RealTypes, test_support::PrecisionNames);

template <typename Real> void check_equations(model::State<Real> const& state)
{
  using C = std::complex<Real>;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  auto const e = [](Real x, Real w) { return C{x, w / Real{2}} / C{x, -w / Real{2}}; };
  std::array<Real, 2> residual{};
  auto const& roots = state.rapidities;
  Real energy = Real{2} * Real(state.holes);
  for (std::size_t a = 0; a < 2; ++a)
    for (std::size_t j = 0; j < roots[a].size(); ++j)
    {
      Real const x = roots[a][j];
      if (j) EXPECT_GT(x, roots[a][j - 1]);
      Real f = -pi * Real(state.quantum_numbers[a][j].twice());
      C lhs{1}, rhs{1};
      if (a == 0)
      {
        f += Real(state.sites) * Real{2} * std::atan(Real{2} * x);
        for (std::size_t k = 0; k < state.sites; ++k)
          lhs *= e(x, Real{1});
        energy -= Real{1} / (x * x + Real{1} / Real{4});
        for (std::size_t k = 0; k < roots[0].size(); ++k)
          if (k != j)
          {
            f -= Real{2} * std::atan(x - roots[0][k]);
            rhs *= e(x - roots[0][k], Real{2});
          }
        for (Real other : roots[1])
        {
          f += Real{2} * std::atan(Real{2} * (x - other));
          rhs /= e(x - other, Real{1});
        }
      }
      else
        for (Real other : roots[0])
        {
          f += Real{2} * std::atan(Real{2} * (x - other));
          rhs *= e(x - other, Real{1});
        }
      residual[a] = std::max(residual[a], std::abs(f) / Real(state.sites));
      if (state.converged) EXPECT_LT(std::abs(lhs - rhs), Real{512} * Real(state.sites) * eps);
    }
  for (std::size_t a = 0; a < 2; ++a)
    EXPECT_REAL_NEAR(state.level_residuals[a], residual[a], Real{128} * eps);
  EXPECT_REAL_NEAR(state.energy, energy, Real{128} * Real(state.sites) * eps);
}

TYPED_TEST(TJ, FiniteSizeGroundStatesAndOriginalEquations)
{
  using Real = TypeParam;
  for (auto counts :
       {std::array<std::size_t, 3>{3, 1, 1}, {7, 3, 3}, {8, 3, 3}, {12, 5, 3}, {24, 7, 5}, {48, 13, 11}, {96, 25, 23}})
  {
    auto const [n, up, down] = counts;
    SCOPED_TRACE(n);
    auto const state = model::ground_state<Real>(n, up, down);
    ASSERT_TRUE(state.converged) << int(state.status) << " residual=" << uni20::format_scalar(state.residual_norm);
    EXPECT_EQ(state.branch, model::Branch::sutherland);
    EXPECT_EQ(state.rapidities[0].size(), n - up);
    EXPECT_EQ(state.rapidities[1].size(), n - up - down);
    EXPECT_EQ(state.momentum_index, 0U);
    ASSERT_NO_FATAL_FAILURE(check_equations(state));
    auto const reversed = model::ground_state<Real>(n, down, up);
    EXPECT_EQ(state.energy, reversed.energy);
    EXPECT_EQ(state.rapidities, reversed.rapidities);
  }
}

TYPED_TEST(TJ, ProjectedFockSpaceGroundStatesAndMomentum)
{
  using Real = TypeParam;
  for (unsigned n = 3; n <= 8; ++n)
    for (unsigned up = 0; up <= n; ++up)
      for (unsigned down = 0; down <= std::min(up, n - up); ++down)
      {
        if (down && up + down < n && (up % 2 == 0 || down % 2 == 0)) continue;
        SCOPED_TRACE(::testing::Message() << "L=" << n << " up=" << up << " down=" << down);
        auto const state = model::ground_state<Real>(n, up, down);
        ASSERT_TRUE(state.converged);
        auto const exact = bethe::test::tj_exact_ground(n, up, down);
        EXPECT_REAL_NEAR(static_cast<double>(state.energy), exact.energy, 3e-10);
        EXPECT_GT(exact.momentum_weights[state.momentum_index], 0.9);
      }
}

TYPED_TEST(TJ, AnalyticLimitsAtNativePrecision)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {3, 4, 5, 6, 9, 16, 32})
  {
    auto const pair = model::ground_state<Real>(n, 1, 1);
    ASSERT_TRUE(pair.converged);
    EXPECT_REAL_NEAR(pair.energy, -Real{4}, Real{512} * Real(n) * eps);
    EXPECT_EQ(model::ground_state<Real>(n, 0, 0).energy, Real{0});
    EXPECT_EQ(model::ground_state<Real>(n, n, 0).energy, Real{0});
    EXPECT_REAL_NEAR(model::ground_state<Real>(n, 1, 0).energy, -Real{2}, Real{8} * eps);
  }
  auto const four = model::ground_state<Real>(5, 3, 1);
  EXPECT_REAL_NEAR(four.energy, -Real{3} - std::sqrt(Real{5}), Real{512} * eps);
  auto const triangle = model::ground_state<Real>(3, 1, 1);
  EXPECT_REAL_NEAR(triangle.rapidities[0].back(), Real{1} / std::sqrt(Real{12}), Real{64} * eps);
  for (std::size_t n : {5, 6, 9, 10})
  {
    auto const down = n / 2, up = n - down;
    auto const state = model::ground_state<Real>(n, up, down);
    auto const xxx = bethe::heisenberg::sector_ground_state<Real>(n, uni20::from_twice(std::int64_t(up - down)));
    EXPECT_EQ(state.branch, model::Branch::no_holes_xxx);
    EXPECT_REAL_NEAR(state.energy, Real{2} * xxx.energy - Real(n) / Real{2}, Real{128} * Real(n) * eps);
    ASSERT_NO_FATAL_FAILURE(check_equations(state));
  }
}

TYPED_TEST(TJ, NoHoleEnergyDoesNotCancelExtensiveConstants)
{
  using Real = TypeParam;
  // One overturned spin on an even ring has lambda=0 and E_tJ=-4.
  // Computing 2*E_XXX-L/2 loses that energy entirely in fp64 here.
  std::size_t const n = 1000000000000000000ULL;
  auto const state = model::ground_state<Real>(n, n - 1, 1, {.max_iterations = 0});
  ASSERT_TRUE(state.converged);
  EXPECT_EQ(state.energy, -Real{4});
  EXPECT_EQ(state.rapidities[0].size(), 1U);
  EXPECT_EQ(state.momentum_index, 0U);
}

TYPED_TEST(TJ, AnalyticJacobian)
{
  using Real = TypeParam;
  Real const h = std::cbrt(uni20::numeric_limits<Real>::epsilon());
  model::detail::GroundSystem<Real> const system(8, 3, 3);
  auto x = system.seed();
  std::vector<Real> jac;
  (void)system.evaluate(x, &jac);
  for (std::size_t j = 0; j < x.size(); ++j)
  {
    auto a = x, b = x;
    a[j] += h;
    b[j] -= h;
    auto const fa = system.evaluate(a), fb = system.evaluate(b);
    for (std::size_t i = 0; i < x.size(); ++i)
      EXPECT_REAL_NEAR(jac[i * x.size() + j], (fa.residual[i] - fb.residual[i]) / (Real{2} * h), Real{2048} * h * h);
  }
}

TYPED_TEST(TJ, BudgetsAndRecoverableLinearFailures)
{
  using Real = TypeParam;
  for (std::size_t budget : {0, 1, 2})
  {
    auto const state = model::ground_state<Real>(8, 3, 3, {.max_iterations = budget});
    EXPECT_FALSE(state.converged);
    EXPECT_EQ(state.iterations, budget);
    EXPECT_EQ(state.status, model::SolveStatus::iteration_limit);
    ASSERT_NO_FATAL_FAILURE(check_equations(state));
  }
  std::vector<Real> rhs{1, 1};
  EXPECT_FALSE(model::detail::newton_step<Real>({0, 0, 0, 0}, rhs));
  EXPECT_FALSE(model::detail::newton_step<Real>({1, 1, 1, 1}, rhs));
  rhs = {1, 2};
  EXPECT_TRUE(model::detail::newton_step<Real>({0, 1, 1, 0}, rhs));
  EXPECT_EQ(rhs[0], Real{2});
  EXPECT_EQ(rhs[1], Real{1});
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  auto const fine = model::ground_state<Real>(8, 3, 3, {.residual_tolerance = eps * eps, .max_iterations = 40});
  EXPECT_FALSE(fine.converged);
}

TYPED_TEST(TJ, InvalidAndUnsupportedInputs)
{
  using Real = TypeParam;
  EXPECT_THROW((void)model::ground_state<Real>(2, 1, 1), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state<Real>(5, 4, 2), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state<Real>(5, std::numeric_limits<std::size_t>::max(), 1), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state<Real>(8, 2, 2), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state<Real>(8, 3, 2), std::invalid_argument);
  for (Real tol :
       {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW((void)model::ground_state<Real>(8, 3, 3, {.residual_tolerance = tol}), std::invalid_argument);
}
} // namespace
