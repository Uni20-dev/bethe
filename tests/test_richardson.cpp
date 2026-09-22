// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "richardson_ed.hpp"
#include "test_support.hpp"
#include <bethe/richardson.hpp>
#include <complex>

namespace
{
namespace model = bethe::richardson;
template <typename Real> class Richardson : public ::testing::Test {};
TYPED_TEST_SUITE(Richardson, test_support::RealTypes, test_support::PrecisionNames);

template <uni20::Real Real> void check_equations(model::State<Real> const& s)
{
  Real const eps = uni20::numeric_limits<Real>::epsilon(), g = s.reached_coupling;
  auto const& y = s.eigenvalue_variables;
  Real sum{}, energy{};
  for (auto i : s.blocked)
    energy += s.levels[i];
  for (std::size_t i = 0; i < y.size(); ++i)
  {
    Real f = y[i] * (y[i] - Real{1}), scale = Real{1} + std::abs(f);
    for (std::size_t j = 0; j < y.size(); ++j)
      if (i != j)
      {
        Real const term = g * (y[i] - y[j]) / (Real{2} * (s.levels[s.active[i]] - s.levels[s.active[j]]));
        f -= term;
        scale += std::abs(g / (Real{2} * (s.levels[s.active[i]] - s.levels[s.active[j]]))) *
                 (std::abs(y[i]) + std::abs(y[j]));
      }
    EXPECT_REAL_NEAR(f / scale, Real{0}, Real{256} * Real(y.size() + 1) * eps);
    sum += y[i];
    energy += Real{2} * s.levels[s.active[i]] * y[i];
  }
  energy -= g * Real(s.pairs) * Real(y.size() - s.pairs + 1);
  EXPECT_REAL_NEAR(sum, Real(s.pairs), Real{64} * Real(y.size() + 1) * eps);
  EXPECT_REAL_NEAR(energy, s.energy, Real{512} * (Real{1} + std::abs(energy)) * Real(y.size() + 1) * eps);
}

TYPED_TEST(Richardson, AnalyticOnePairAndExactLimits)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  std::vector<Real> levels{Real{0}, Real{1}};
  for (Real g : {Real{0}, Real{1} / Real{100}, Real{1}, Real{10}, Real{100}})
  {
    auto const s = model::ground_state<Real>(levels, 1, g);
    ASSERT_TRUE(s.converged) << int(s.status) << " stages=" << s.stages;
    Real const exact = Real{1} - g - std::sqrt(Real{1} + g * g);
    EXPECT_REAL_NEAR(s.energy, exact, Real{1024} * (Real{1} + g) * eps);
    ASSERT_NO_FATAL_FAILURE(check_equations(s));
  }
  // With epsilon_0=0, the weak energy remains O(g). Eliminating an EMPTY
  // y_i as 1-sum(y) would erase it and return an incorrect leading coefficient.
  Real const weak = eps * eps;
  auto const s = model::ground_state<Real>(levels, 1, weak);
  ASSERT_TRUE(s.converged);
  EXPECT_REAL_NEAR(s.energy / weak, -Real{1}, Real{256} * eps);
  model::SolverOptions<Real> zero;
  zero.max_iterations = 0;
  for (std::size_t m : {0, 1, 2})
  {
    auto const free = model::ground_state<Real>(levels, m, Real{0}, {}, zero);
    EXPECT_TRUE(free.converged);
    EXPECT_EQ(free.energy, m == 2 ? Real{2} : Real{0});
    if (m != 1)
    {
      auto const exact = model::ground_state<Real>(levels, m, Real{3}, {}, zero);
      EXPECT_TRUE(exact.converged);
      EXPECT_EQ(exact.energy, m ? -Real{4} : Real{0});
      EXPECT_EQ(exact.iterations, 0u);
    }
  }
}

TYPED_TEST(Richardson, AllSmallPairSectorsAndBlockedLevels)
{
  using Real = TypeParam;
  for (std::vector<double> input :
       {std::vector<double>{-2, -.8, .1, .9, 2, 3.5}, std::vector<double>{0, .01, .013, .2, 1}})
  {
    std::vector<Real> levels(input.begin(), input.end());
    for (std::vector<std::size_t> blocked : {std::vector<std::size_t>{}, std::vector<std::size_t>{1, 3}})
      for (unsigned m = 0; m <= levels.size() - blocked.size(); ++m)
        for (Real g : {Real{1} / Real{100}, Real{1} / Real{2}, Real{2}, Real{20}})
        {
          SCOPED_TRACE(::testing::Message() << m << "," << double(g) << "," << blocked.size());
          auto const s = model::ground_state<Real>(levels, m, g, blocked);
          ASSERT_TRUE(s.converged) << int(s.status) << " at " << double(s.reached_coupling) << " stages=" << s.stages;
          EXPECT_NEAR(double(s.energy), bethe::test::richardson_exact_ground(input, m, double(g), blocked), 2e-10);
          ASSERT_NO_FATAL_FAILURE(check_equations(s));
        }
  }
}

TYPED_TEST(Richardson, PairRootsAcrossACollision)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  std::vector<Real> levels{Real{0}, Real{1}, Real{2}, Real{3}};
  for (Real g : {Real{3} / Real{5}, Real{7} / Real{10}, Real{1}, Real{2}})
  {
    auto const s = model::ground_state<Real>(levels, 2, g);
    ASSERT_TRUE(s.converged);
    auto const& y = s.eigenvalue_variables;
    // P(z)=z^2-E*z+T; at e_0=0, y_0*T=-g*E.
    Real const t = -g * s.energy / y[0], discriminant = s.energy * s.energy - Real{4} * t;
    if (g == Real{3} / Real{5})
      EXPECT_GT(discriminant, Real{0});
    else
      EXPECT_LT(discriminant, Real{0});
    C const d = std::sqrt(C{discriminant, Real{0}});
    std::array<C, 2> roots{(C{s.energy, Real{0}} + d) / Real{2}, (C{s.energy, Real{0}} - d) / Real{2}};
    for (std::size_t j = 0; j < 2; ++j)
    {
      C f = Real{1} / g - Real{2} / (roots[j] - roots[1 - j]);
      for (Real e : levels)
        f += Real{1} / (roots[j] - Real{2} * e);
      EXPECT_REAL_NEAR(f.real(), Real{0}, Real{8192} * eps);
      EXPECT_REAL_NEAR(f.imag(), Real{0}, Real{8192} * eps);
    }
    for (std::size_t j = 0; j < levels.size(); ++j)
    {
      C const from_roots =
          g * (Real{1} / (Real{2} * levels[j] - roots[0]) + Real{1} / (Real{2} * levels[j] - roots[1]));
      EXPECT_REAL_NEAR(from_roots.real(), y[j], Real{8192} * eps);
      EXPECT_REAL_NEAR(from_roots.imag(), Real{0}, Real{8192} * eps);
    }
  }
}

TYPED_TEST(Richardson, ExactCollisionAndLinearization)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  std::vector<Real> levels{Real{0}, Real{1}, Real{2}, Real{3}};
  auto const collision = model::ground_state<Real>(levels, 2, Real{2} / Real{3});
  ASSERT_TRUE(collision.converged);
  EXPECT_REAL_NEAR(collision.energy, Real{0}, Real{2048} * eps);
  std::array<Real, 4> expected{Real{7} / Real{9}, Real{2} / Real{3}, Real{1} / Real{3}, Real{2} / Real{9}};
  for (std::size_t j = 0; j < 4; ++j)
    EXPECT_REAL_NEAR(collision.eigenvalue_variables[j], expected[j], Real{512} * eps);
  // At this point both pair roots equal the pole E_alpha=0; direct rational
  // equations are singular, while the EBVs and physical ground state are regular.
  EXPECT_NEAR(bethe::test::richardson_exact_ground(std::vector<double>{0, 1, 2, 3}, 2, 2.0 / 3.0), 0.0, 1e-12);

  model::detail::System<Real> const system({Real{0}, Real{2}, Real{4}, Real{6}}, 2);
  std::vector<Real> x{Real{9} / Real{10}, Real{1} / Real{20}, Real{1} / Real{20}}, jac;
  Real const g = Real{3} / Real{10}, h = std::cbrt(eps);
  (void)system.evaluate(x, g, &jac);
  auto raw = [&](std::vector<Real> const& values, std::size_t i) {
    auto const y = system.expand(values);
    Real f = y[i] * (y[i] - Real{1});
    for (std::size_t j = 0; j < 4; ++j)
      if (j != i) f -= g * (y[i] - y[j]) / (system.levels[i] - system.levels[j]);
    return f;
  };
  auto const y = system.expand(x);
  for (std::size_t i = 0; i < 4; ++i)
  {
    Real scale = Real{1} + std::abs(y[i] * (y[i] - Real{1}));
    for (std::size_t j = 0; j < 4; ++j)
      if (j != i) scale += std::abs(g / (system.levels[i] - system.levels[j])) * (std::abs(y[i]) + std::abs(y[j]));
    for (std::size_t j = 0; j < 3; ++j)
    {
      auto p = x, m = x;
      p[j] += h;
      m[j] -= h;
      EXPECT_REAL_NEAR((raw(p, i) - raw(m, i)) / (Real{2} * h * scale), jac[i * 3 + j], Real{2048} * h * h);
    }
  }
  // A genuinely rectangular system, with nonzero least-squares residual.
  std::vector<Real> rhs{Real{1}, Real{2}, Real{4}};
  ASSERT_TRUE(bethe::detail::least_squares_step(std::vector<Real>{Real{1}, Real{0}, Real{0}, Real{1}, Real{1}, Real{1}},
                                                rhs, 2));
  EXPECT_REAL_NEAR(rhs[0], Real{4} / Real{3}, Real{64} * eps);
  EXPECT_REAL_NEAR(rhs[1], Real{7} / Real{3}, Real{64} * eps);
  rhs.assign(3, Real{1});
  EXPECT_FALSE(bethe::detail::least_squares_step(std::vector<Real>(6, Real{0}), rhs, 2));
  EXPECT_FALSE(bethe::detail::least_squares_step(std::vector<Real>(1, Real{0}), rhs, 2));
}

TYPED_TEST(Richardson, ScalingShiftsAndLargerSystems)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  std::vector<Real> levels{Real{0}, Real{1}, Real{2}, Real{3}, Real{5}};
  std::vector<std::size_t> blocked{2};
  auto const base = model::ground_state<Real>(levels, 2, Real{1}, blocked);
  ASSERT_TRUE(base.converged);
  for (Real factor : {Real{1} / Real{100}, Real{10}})
  {
    auto shifted = levels;
    for (auto& e : shifted)
      e = factor * e + Real{7};
    auto const s = model::ground_state<Real>(shifted, 2, factor, blocked);
    ASSERT_TRUE(s.converged);
    EXPECT_REAL_NEAR(s.energy, factor * base.energy + Real{35}, Real{8192} * eps);
  }
  for (std::size_t n : {12, 24, 48})
  {
    levels.resize(n);
    for (std::size_t i = 0; i < n; ++i)
      levels[i] = Real(i) / Real(n);
    auto const s = model::ground_state<Real>(levels, n / 2, Real{1} / Real{5});
    ASSERT_TRUE(s.converged) << n << "," << int(s.status) << "," << s.stages;
    ASSERT_NO_FATAL_FAILURE(check_equations(s));
  }
}

TYPED_TEST(Richardson, BudgetsAndInputValidation)
{
  using Real = TypeParam;
  std::vector<Real> levels{Real{0}, Real{1}, Real{2}, Real{3}};
  for (std::size_t budget : {0, 1, 2})
  {
    model::SolverOptions<Real> options;
    options.max_iterations = budget;
    auto const s = model::ground_state<Real>(levels, 2, Real{1}, {}, options);
    EXPECT_FALSE(s.converged);
    EXPECT_EQ(s.status, model::SolveStatus::iteration_limit);
    EXPECT_EQ(s.iterations, budget);
    EXPECT_LT(s.reached_coupling, s.coupling);
    ASSERT_NO_FATAL_FAILURE(check_equations(s));
  }
  model::SolverOptions<Real> limited;
  limited.max_stages = 0;
  auto const no_stage = model::ground_state<Real>(levels, 2, Real{1}, {}, limited);
  EXPECT_FALSE(no_stage.converged);
  EXPECT_EQ(no_stage.status, model::SolveStatus::stage_limit);
  EXPECT_EQ(no_stage.stages, 0u);
  EXPECT_EQ(no_stage.energy, Real{2});
  EXPECT_THROW((void)model::ground_state<Real>(levels, 5, Real{1}), std::invalid_argument);
  std::vector<std::size_t> duplicates{1, 1}, outside{4};
  EXPECT_THROW((void)model::ground_state<Real>(levels, 1, Real{1}, duplicates), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state<Real>(levels, 1, Real{1}, outside), std::invalid_argument);
  for (Real g : {-Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW((void)model::ground_state<Real>(levels, 1, g), std::invalid_argument);
  for (std::vector<Real> bad :
       {std::vector<Real>{}, std::vector<Real>{Real{1}, Real{0}}, std::vector<Real>{Real{1}, Real{1}},
        std::vector<Real>{uni20::numeric_limits<Real>::infinity()}})
    EXPECT_THROW((void)model::ground_state<Real>(bad, 0, Real{0}), std::invalid_argument);
  for (Real tol :
       {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    model::SolverOptions<Real> o;
    o.residual_tolerance = tol;
    EXPECT_THROW((void)model::ground_state<Real>(levels, 2, Real{1}, {}, o), std::invalid_argument);
  }
}
} // namespace
