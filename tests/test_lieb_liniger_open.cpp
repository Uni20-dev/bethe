// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/lieb_liniger_open.hpp>
#include <set>

namespace
{
namespace model = bethe::lieb_liniger::open;
using Numbers = bethe::lieb_liniger::QuantumNumbers;
using uni20::half_int;
template <typename Real> class LiebLinigerOpen : public ::testing::Test {};
TYPED_TEST_SUITE(LiebLinigerOpen, test_support::RealTypes, test_support::PrecisionNames);

template <typename Real> void check_state(model::State<Real> const& state)
{
  using C = uni20::complex<Real>;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  Real energy{};
  for (std::size_t j = 0; j < state.particles; ++j)
  {
    Real const k = state.momenta[j];
    ASSERT_GT(k, Real{0});
    if (j)
    {
      ASSERT_GT(k, state.momenta[j - 1]);
    }
    energy += k * k;
    if (!state.converged) continue;
    Real logarithm = state.length * k;
    C product{1, 0};
    for (std::size_t l = 0; l < state.particles; ++l)
      if (j != l)
        for (Real x : {k - state.momenta[l], k + state.momenta[l]})
        {
          logarithm += std::atan(x / state.interaction);
          product *= C{x, state.interaction} / C{x, -state.interaction};
        }
    Real const target = pi * Real(state.quantum_numbers[j].twice()) / Real{2};
    EXPECT_REAL_NEAR(logarithm, target, Real{256} * Real(state.particles + 1) * eps * (Real{1} + target));
    C const difference = std::exp(C{0, Real{2} * state.length * k}) - product;
    EXPECT_LT(std::abs(difference), Real{1024} * Real(state.particles + 1) * eps * (Real{1} + target));
  }
  EXPECT_REAL_NEAR(energy, state.energy, Real{32} * Real(state.particles + 1) * eps * (Real{1} + energy));
  EXPECT_EQ(state.converged, state.status == bethe::lieb_liniger::SolveStatus::converged);
}

// For two particles, q2+q1 and q2-q1 independently obey x+2 atan(x/g)=pi*K.
// Bisection of the complementary equation has no many-body Newton/Jacobian.
template <typename Real> Real scalar_root(Real g, int label)
{
  Real const pi = Real{4} * std::atan(Real{1}), origin = pi * Real(label - 1);
  Real lo = origin, hi = pi * Real(label);
  for (int i = 0; i < uni20::numeric_limits<Real>::digits + 30; ++i)
  {
    Real const mid = lo + (hi - lo) / Real{2};
    if (mid == lo || mid == hi) break;
    if (mid - origin - Real{2} * std::atan2(g, mid) > Real{0})
      hi = mid;
    else
      lo = mid;
  }
  return lo + (hi - lo) / Real{2};
}

TYPED_TEST(LiebLinigerOpen, VacuumAndSingleParticle)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  auto const vacuum = model::ground_state(0, Real{2}, Real{3});
  ASSERT_TRUE(vacuum.converged);
  EXPECT_EQ(vacuum.energy, Real{0});
  for (int label : {1, 2, 7})
  {
    auto const s = model::solve_real(Real{3}, Real{4}, Numbers{half_int(label)}, {.max_iterations = 0});
    ASSERT_TRUE(s.converged);
    EXPECT_EQ(s.iterations, 0u);
    Real const exact = pi * pi * Real(label * label) / Real{9};
    EXPECT_REAL_NEAR(s.energy, exact, Real{8} * eps * exact);
    if constexpr (uni20::numeric_limits<Real>::digits > 53)
      if (label == 1) EXPECT_GT(std::abs(Real(double(exact)) - exact), Real{8} * eps * exact);
    check_state(s);
  }
}

TYPED_TEST(LiebLinigerOpen, IndependentTwoBodyEquations)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real g : {Real{1} / Real{100000000}, Real{1} / Real{10}, Real{1}, Real{4}, Real{100000000}})
    for (int i : {1, 2, 5})
      for (int j : {i + 1, i + 3})
      {
        SCOPED_TRACE(::testing::Message() << uni20::format_real(g) << " " << i << " " << j);
        auto const s = model::solve_real(Real{2}, g / Real{2}, Numbers{half_int(i), half_int(j)});
        ASSERT_TRUE(s.converged) << s.iterations;
        Real const sum = scalar_root(g, i + j), diff = scalar_root(g, j - i);
        EXPECT_REAL_NEAR(s.energy, (sum * sum + diff * diff) / Real{8}, Real{512} * eps * s.energy);
        EXPECT_REAL_NEAR(s.momenta[0], (sum - diff) / Real{4}, Real{512} * eps * (Real{1} + sum));
        EXPECT_REAL_NEAR(s.momenta[1], (sum + diff) / Real{4}, Real{512} * eps * (Real{1} + sum));
        check_state(s);
      }
}

TYPED_TEST(LiebLinigerOpen, SizesLimitsAndScaling)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {2, 3, 8, 32})
  {
    for (Real c : {Real{1} / Real{100000000}, Real{1} / Real{10}, Real{1}, Real{100000000}})
    {
      SCOPED_TRACE(::testing::Message() << n << " " << uni20::format_real(c));
      auto const s = model::ground_state(n, Real(n), c);
      ASSERT_TRUE(s.converged) << s.iterations;
      check_state(s);
      auto const scaled = model::ground_state(n, Real(2 * n), c / Real{2});
      ASSERT_TRUE(scaled.converged);
      EXPECT_REAL_NEAR(scaled.energy * Real{4}, s.energy, Real{256} * eps * s.energy);
    }
    Real const g = Real{1} / Real{1000000};
    auto const weak = model::ground_state(n, Real{1}, g);
    ASSERT_TRUE(weak.converged);
    EXPECT_REAL_NEAR(weak.energy, Real(n) * pi * pi + Real{1.5} * g * Real(n * (n - 1)),
                     Real{256} * eps * weak.energy + g * g * Real(n * n * n));
    auto const strong = model::ground_state(n, Real{1}, Real{100000000});
    ASSERT_TRUE(strong.converged);
    Real const tonks = pi * pi * Real(n) * Real(n + 1) * Real(2 * n + 1) / Real{6};
    EXPECT_LT(strong.energy, tonks);
    EXPECT_REAL_NEAR(strong.energy / tonks, Real{1}, Real(5 * n) / Real{100000000});
  }
}

TYPED_TEST(LiebLinigerOpen, JacobianAndFiniteWindow)
{
  using Real = TypeParam;
  Real const h = std::cbrt(uni20::numeric_limits<Real>::epsilon());
  auto const labels = model::ground_quantum_numbers(4);
  for (Real g : {Real{1} / Real{10}, Real{1}, Real{10}})
  {
    model::detail::System<Real> const system{labels, g};
    auto const q = system.seed();
    uni20::DenseMatrix<Real> jac(4, 4);
    (void)system.evaluate(q, &jac);
    for (std::size_t col = 0; col < 4; ++col)
    {
      auto plus = q, minus = q;
      plus[col] += h;
      minus[col] -= h;
      auto const a = system.evaluate(plus), b = system.evaluate(minus);
      for (std::size_t row = 0; row < 4; ++row)
      {
        Real const numerical = (a.residual[row] - b.residual[row]) / (Real{2} * h);
        EXPECT_REAL_NEAR((jac[row, col]), numerical, Real{2000} * h * h * (Real{1} + std::abs(numerical)));
        EXPECT_EQ((jac[row, col]), (jac[col, row]));
      }
    }
  }
  for (std::size_t n : {0, 1, 2, 3, 4})
    for (std::size_t padding : {0, 1, 2})
    {
      auto const scan =
          model::real_excitations(n, Real{4}, Real{2}, padding, {.count = std::numeric_limits<std::size_t>::max()});
      ASSERT_TRUE(scan.converged());
      EXPECT_EQ(scan.candidate_count, model::excitation_count(n, padding));
      EXPECT_EQ(scan.levels.size(), scan.candidate_count);
      ASSERT_FALSE(scan.levels.empty());
      EXPECT_EQ(scan.levels.front().state.quantum_numbers, model::ground_quantum_numbers(n));
      std::set<Numbers> labels_seen;
      for (std::size_t i = 0; i < scan.levels.size(); ++i)
      {
        auto const& level = scan.levels[i];
        ASSERT_EQ(level.state.quantum_numbers.size(), n);
        for (auto label : level.state.quantum_numbers)
        {
          EXPECT_GE(label, half_int{1});
          EXPECT_LE(label, half_int(std::int64_t(n + padding)));
        }
        ASSERT_TRUE(labels_seen.insert(level.state.quantum_numbers).second);
        ASSERT_TRUE(level.gap);
        EXPECT_GE(*level.gap, Real{0});
        if (i) EXPECT_LE(scan.levels[i - 1].state.energy, level.state.energy);
        check_state(level.state);
      }
    }
}

TYPED_TEST(LiebLinigerOpen, BudgetsRestartAndDomain)
{
  using Real = TypeParam;
  auto const labels = model::ground_quantum_numbers(4);
  for (std::size_t budget : {0, 1})
  {
    auto const failed = model::solve_real(Real{2}, Real{1}, labels, {.max_iterations = budget});
    EXPECT_FALSE(failed.converged);
    EXPECT_EQ(failed.iterations, budget);
    EXPECT_EQ(failed.status, bethe::lieb_liniger::SolveStatus::iteration_limit);
  }
  auto const s = model::ground_state(4, Real{2}, Real{1});
  ASSERT_TRUE(s.converged);
  auto const restart =
      model::solve_real(Real{2}, Real{1}, labels, {.max_iterations = 0}, std::span<Real const>(s.momenta));
  EXPECT_TRUE(restart.converged);
  EXPECT_EQ(restart.momenta, s.momenta);
  auto const failed = model::real_excitations(4, Real{2}, Real{1}, 1, {}, {.max_iterations = 0});
  EXPECT_FALSE(failed.converged());
  EXPECT_TRUE(failed.first_unconverged);
  EXPECT_TRUE(failed.levels.empty());
  for (Numbers const& bad : {Numbers{half_int{0}}, Numbers{half_int{-1}}, Numbers{half_int{2}, half_int{1}},
                             Numbers{half_int{1}, half_int{1}}, Numbers{uni20::from_twice(std::int64_t{1})}})
    EXPECT_THROW((model::solve_real(Real{1}, Real{1}, bad)), std::invalid_argument);
  EXPECT_THROW((model::ground_state(2, Real{0}, Real{1})), std::invalid_argument);
  EXPECT_THROW((model::ground_state(2, Real{1}, Real{0})), std::invalid_argument);
  EXPECT_THROW((model::ground_state(2, Real{1}, -Real{1})), std::invalid_argument);
  EXPECT_THROW((model::ground_state(2, Real{1}, Real{1}, {.residual_tolerance = Real{0}})), std::invalid_argument);
  EXPECT_THROW((model::excitation_count(4, std::numeric_limits<std::size_t>::max())), std::length_error);
  EXPECT_THROW((model::real_excitations(4, Real{1}, Real{1}, 1, {.max_candidates = 4})), std::length_error);
  EXPECT_THROW((model::real_excitations(4, Real{1}, Real{1}, 1, {.count = 0})), std::invalid_argument);
  EXPECT_EQ(model::excitation_count(4, 1), 5u);
  auto bad = s.momenta;
  bad[0] = Real{0};
  EXPECT_THROW((model::solve_real(Real{2}, Real{1}, labels, {}, std::span<Real const>(bad))), std::invalid_argument);
  EXPECT_THROW((model::ground_state(2, Real{1}, uni20::numeric_limits<Real>::min())), std::invalid_argument);
  EXPECT_THROW((model::ground_state(1, uni20::numeric_limits<Real>::min(), Real{16})), std::overflow_error);
  EXPECT_THROW((model::ground_state(1, uni20::numeric_limits<Real>::max(), Real{1})), std::underflow_error);
}
} // namespace
