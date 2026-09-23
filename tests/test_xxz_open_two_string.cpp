// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "tl_ed.hpp"
#include <bethe/biquadratic.hpp>
#include <bethe/biquadratic_two_string.hpp>
#include <complex>

namespace
{
namespace ts = bethe::xxz::quantum_group::two_string;
template <typename Real> class OpenTwoString : public ::testing::Test {};
TYPED_TEST_SUITE(OpenTwoString, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(OpenTwoString, FourSiteNativePrecision)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon();
  auto const state = bethe::biquadratic::two_string::singlet<R>(4);
  ASSERT_TRUE(state.reference.converged) << int(state.reference.status);
  R const exact = (-R{15} + std::sqrt(R{17})) / R{2};
  EXPECT_REAL_NEAR(state.energy, exact, R{512} * eps);
  if constexpr (uni20::numeric_limits<R>::digits > 53) EXPECT_GT(std::abs(R(double(exact)) - exact), R{512} * eps);
  EXPECT_EQ(state.multiplicity, 1);
  EXPECT_EQ(state.through_lines, 0);
}

TYPED_TEST(OpenTwoString, AnalyticJacobian)
{
  using R = TypeParam;
  R const h = std::cbrt(uni20::numeric_limits<R>::epsilon());
  for (R delta : {R{11} / R{10}, R{3} / R{2}, R{3}})
    for (bool ideal : {false, true})
    {
      ts::detail::System<R> system(8, delta);
      auto x = system.seed();
      x.back() = R{4}; // Resolve deviation derivatives as well as the ideal limit.
      std::vector<R> jac;
      system.evaluate(x, &jac, ideal);
      for (std::size_t j = 0; j < x.size(); ++j)
      {
        auto plus = x, minus = x;
        plus[j] += h;
        minus[j] -= h;
        auto const a = system.evaluate(plus, nullptr, ideal), b = system.evaluate(minus, nullptr, ideal);
        for (std::size_t i = 0; i < x.size(); ++i)
          EXPECT_REAL_NEAR((a.residual[i] - b.residual[i]) / (R{2} * h), jac[i * x.size() + j],
                           R{4096} * h * h * std::max(R{1}, std::abs(jac[i * x.size() + j])));
      }
    }
}

TYPED_TEST(OpenTwoString, OriginalComplexEquations)
{
  using R = TypeParam;
  using C = std::complex<R>;
  R const eps = uni20::numeric_limits<R>::epsilon();
  for (std::size_t n : {4, 6, 8, 16, 32})
  {
    auto const s = ts::singlet<R>(n, R{3} / R{2});
    ASSERT_TRUE(s.converged) << n << " " << int(s.status);
    R const eta = std::acosh(s.delta), d = std::exp(-s.log_deviation);
    std::vector<C> roots;
    for (R alpha : s.rapidities)
      roots.emplace_back(R{0}, alpha / R{2});
    roots.emplace_back((eta + d) / R{2}, s.center / R{2});
    roots.emplace_back((eta + d) / R{2}, -s.center / R{2});
    for (std::size_t i = 0; i < roots.size(); ++i)
    {
      C const drive = std::sinh(roots[i] + eta / R{2}) / std::sinh(roots[i] - eta / R{2});
      C lhs{1}, rhs{1};
      for (std::size_t k = 0; k < 2 * n; ++k)
        lhs *= drive;
      for (std::size_t j = 0; j < roots.size(); ++j)
        if (i != j)
          rhs *= std::sinh(roots[i] - roots[j] + eta) * std::sinh(roots[i] + roots[j] + eta) /
                 (std::sinh(roots[i] - roots[j] - eta) * std::sinh(roots[i] + roots[j] - eta));
      // Reconstructed roots lose relative O(eps/d) at the string pole.
      EXPECT_REAL_NEAR(std::abs(lhs / rhs - C{1}), R{0}, R{8192} * R(n) * eps / d);
    }
  }
}

TYPED_TEST(OpenTwoString, LargeChainsAndBudgets)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon();
  for (std::size_t n : {64, 128, 256, 512})
  {
    auto const s = ts::singlet<R>(n, R{3} / R{2});
    ASSERT_TRUE(s.converged) << n << " " << int(s.status) << " " << uni20::format_real(s.residual_norm);
    EXPECT_LE(s.residual_norm, R{32} * eps);
    EXPECT_GT(s.log_deviation, R{10});
    EXPECT_EQ(s.rapidities.size(), n / 2 - 2);
    EXPECT_LT(s.iterations, 30);
    // Independent original-complex-equation check, including real roots that
    // cross the string center. Evaluate only the singular reflected factor
    // from d itself, never from rounded (u_+ + u_- - eta).
    using C = std::complex<R>;
    R const eta = std::acosh(s.delta), d = std::exp(-s.log_deviation);
    std::vector<C> roots;
    for (R alpha : s.rapidities)
      roots.emplace_back(R{0}, alpha / R{2});
    auto const pair = roots.size();
    roots.emplace_back((eta + d) / R{2}, s.center / R{2});
    roots.emplace_back((eta + d) / R{2}, -s.center / R{2});
    for (std::size_t i : {std::size_t{0}, pair - 1, pair, pair + 1})
    {
      bethe::detail::CompensatedSum<R> magnitude, phase;
      auto factor = [&](C value, R weight) {
        magnitude.add(weight * std::log(std::abs(value)));
        phase.add(weight * std::arg(value));
      };
      factor(std::sinh(roots[i] + eta / R{2}), R{2} * R(n));
      factor(std::sinh(roots[i] - eta / R{2}), -R{2} * R(n));
      for (std::size_t j = 0; j < roots.size(); ++j)
        if (i != j)
        {
          factor(std::sinh(roots[i] - roots[j] + eta), -R{1});
          factor(std::sinh(roots[i] + roots[j] + eta), -R{1});
          factor(std::sinh(roots[i] - roots[j] - eta), R{1});
          factor(i >= pair && j >= pair ? C{std::sinh(d)} : std::sinh(roots[i] + roots[j] - eta), R{1});
        }
      EXPECT_REAL_NEAR(magnitude.value() / (R{2} * R(n)), R{0}, R{128} * eps);
      EXPECT_REAL_NEAR(std::sin(phase.value() / R{2}), R{0}, R{128} * R(n) * eps);
    }
    auto const ground = bethe::biquadratic::ground_state<R>(n);
    ASSERT_TRUE(ground.reference.converged);
    EXPECT_GT(s.energy, ground.reference.energy);
  }
  for (std::size_t budget : {0, 1})
  {
    auto const s = ts::singlet<R>(16, R{3} / R{2}, {.max_iterations = budget});
    EXPECT_FALSE(s.converged);
    EXPECT_EQ(s.iterations, budget);
    EXPECT_EQ(s.status, bethe::xxz::quantum_group::SolveStatus::iteration_limit);
    auto x = s.rapidities;
    x.push_back(s.center);
    x.push_back(s.log_deviation);
    ts::detail::System<R> system(16, s.delta);
    EXPECT_EQ(s.residual_norm, system.evaluate(x).norm);
    EXPECT_EQ(s.energy, system.energy(x, s.delta));
  }
  // L remains meaningful even when exp(-L) underflows to zero.
  ts::detail::System<R> system(8, R{3} / R{2});
  auto x = system.seed();
  x.back() = -R{2} * std::log(uni20::numeric_limits<R>::min());
  ASSERT_EQ(std::exp(-x.back()), R{0});
  auto const a = system.evaluate(x);
  x.back() += R{1};
  auto const b = system.evaluate(x);
  EXPECT_REAL_NEAR(b.residual.back() - a.residual.back(), -R{1} / R{16}, R{32} * eps * x.back());
}

TYPED_TEST(OpenTwoString, InvalidInputs)
{
  using R = TypeParam;
  for (std::size_t n : {0, 2, 3, 5})
    EXPECT_THROW((void)ts::singlet<R>(n, R{3} / R{2}), std::invalid_argument);
  for (R v : {R{0}, R{1}, -R{1}, uni20::numeric_limits<R>::infinity(), uni20::numeric_limits<R>::quiet_NaN()})
    EXPECT_THROW((void)ts::singlet<R>(4, v), std::invalid_argument);
  for (R v : {R{0}, -R{1}, uni20::numeric_limits<R>::infinity(), uni20::numeric_limits<R>::quiet_NaN()})
    EXPECT_THROW((void)ts::singlet<R>(4, R{3} / R{2}, {.residual_tolerance = v}), std::invalid_argument);
  auto const strict =
      ts::singlet<R>(16, R{3} / R{2},
                     {.residual_tolerance = uni20::numeric_limits<R>::epsilon() * uni20::numeric_limits<R>::epsilon(),
                      .max_iterations = 30});
  EXPECT_FALSE(strict.converged);
}

TEST(OpenTwoStringED, LowestExcitedSinglet)
{
  for (unsigned n : {4, 6, 8, 10})
    for (double delta : {1.1, 1.5, 2.0, 3.0})
    {
      auto const s = ts::singlet<double>(n, delta);
      ASSERT_TRUE(s.converged) << n << " Delta=" << delta;
      auto const exact = bethe::test::quantum_group_module_ed(n, 0, delta);
      ASSERT_GE(exact.size(), 2);
      EXPECT_NEAR(s.energy, exact[1], 1e-11) << n << " Delta=" << delta;
    }
}
} // namespace
