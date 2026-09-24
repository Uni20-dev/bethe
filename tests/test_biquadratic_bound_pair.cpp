// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "tl_ed.hpp"
#include <bethe/biquadratic_ferromagnetic.hpp>
#include <bethe/biquadratic_two_string.hpp>
#include <complex>

namespace
{
namespace ferro = bethe::biquadratic::ferromagnetic;
namespace ts = bethe::xxz::quantum_group::two_string;
template <typename Real> class BiquadraticBoundPair : public ::testing::Test {};
TYPED_TEST_SUITE(BiquadraticBoundPair, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(BiquadraticBoundPair, NativePrecisionAndFamily)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon();
  auto const s = ferro::bound_pair<R>(4);
  ASSERT_TRUE(s.reference.converged);
  R const exact = (R{15} - std::sqrt(R{17})) / R{2};
  EXPECT_REAL_NEAR(s.energy, exact, R{128} * eps);
  EXPECT_REAL_NEAR(s.tl_energy, exact - R{3}, R{128} * eps);
  if constexpr (uni20::numeric_limits<R>::digits > 53) EXPECT_GT(std::abs(R(double(exact)) - exact), R{128} * eps);
  auto const af = bethe::biquadratic::two_string::singlet<R>(4);
  EXPECT_REAL_NEAR(s.energy, -af.energy, R{128} * eps);
  EXPECT_EQ(s.multiplicity, 1);
  EXPECT_EQ(s.through_lines, 0);
  for (std::size_t n : {5, 6, 7, 8, 16, 32})
  {
    R previous{};
    for (std::size_t mode = 1; mode <= n - 3; ++mode)
    {
      auto const p = ferro::bound_pair<R>(n, mode);
      ASSERT_TRUE(p.reference.converged) << n << " mode=" << mode << " status=" << int(p.reference.status);
      EXPECT_EQ(p.through_lines, n - 4);
      EXPECT_EQ(p.reference.string_label, n - 2 - mode);
      EXPECT_EQ(p.reference.deviation_sign, mode % 2 ? 1 : -1);
      EXPECT_TRUE(p.reference.rapidities.empty());
      EXPECT_GT(p.tl_energy, previous);
      EXPECT_REAL_NEAR(p.energy, R(n - 1) + p.tl_energy, R{128} * R(n) * eps);
      EXPECT_REAL_NEAR(p.energy, R{7} * R(n - 1) / R{4} - R{2} * p.reference.energy, R{128} * R(n) * eps);
      previous = p.tl_energy;
    }
  }
}

TYPED_TEST(BiquadraticBoundPair, SignedJacobian)
{
  using R = TypeParam;
  R const h = std::cbrt(uni20::numeric_limits<R>::epsilon());
  for (int sign : {-1, 1})
    for (bool ideal : {false, true})
    {
      ts::detail::System<R> system(16, R{3} / R{2}, 2, sign > 0 ? 13 : 12, sign);
      std::vector<R> x{R{2}, R{4}}, jac;
      system.evaluate(x, &jac, ideal);
      for (std::size_t j = 0; j < 2; ++j)
      {
        auto plus = x, minus = x;
        plus[j] += h;
        minus[j] -= h;
        auto a = system.evaluate(plus, nullptr, ideal), b = system.evaluate(minus, nullptr, ideal);
        for (std::size_t i = 0; i < 2; ++i)
          EXPECT_REAL_NEAR((a.residual[i] - b.residual[i]) / (R{2} * h), jac[2 * i + j], R{4096} * h * h);
      }
    }
}

TYPED_TEST(BiquadraticBoundPair, OriginalComplexEquationsBothSigns)
{
  using R = TypeParam;
  using C = std::complex<R>;
  R const eps = uni20::numeric_limits<R>::epsilon(), pi = R{4} * std::atan(R{1});
  for (std::size_t n : {4, 5, 8, 16, 64, 129, 1024})
    for (std::size_t mode = 1; mode <= std::min(n - 3, std::size_t{4}); ++mode)
    {
      auto const s = ferro::bound_pair<R>(n, mode).reference;
      ASSERT_TRUE(s.converged) << n << " " << mode;
      R const eta = std::acosh(s.delta), d = R(s.deviation_sign) * std::exp(-s.log_deviation);
      C const u{(eta + d) / R{2}, s.center / R{2}};
      C const drive = std::sinh(u + eta / R{2}) / std::sinh(u - eta / R{2});
      C const scatter = std::sinh(C{eta, s.center}) / std::sinh(C{-eta, s.center});
      // Original complex equation, with only its singular reflected factor
      // evaluated via L. Includes its negative sign on even modes; omitting
      // that sign would give a phase residual of pi, not a small error.
      R const phase = R{2} * R(n) * std::arg(drive) - std::arg(scatter) - (s.deviation_sign < 0 ? pi : R{0});
      EXPECT_REAL_NEAR(std::sin(phase / R{2}), R{0}, R{256} * R(n) * eps);
      R const correction = std::abs(d) < std::sqrt(eps) ? d * d / R{6} : std::log(std::sinh(d) / d);
      R const magnitude = std::log(std::abs(drive)) -
                          (std::log(std::sinh(R{2} * eta + d)) + s.log_deviation - correction) / (R{2} * R(n));
      EXPECT_REAL_NEAR(magnitude, R{0}, R{128} * eps);
      if (n <= 8 && mode <= 2)
      {
        C const rhs = scatter * (std::sinh(R{2} * eta + d) / std::sinh(d));
        C lhs{1};
        for (std::size_t j = 0; j < 2 * n; ++j)
          lhs *= drive;
        EXPECT_REAL_NEAR(std::abs(lhs / rhs - C{1}), R{0}, R{4096} * eps);
      }
    }
}

TYPED_TEST(BiquadraticBoundPair, LongChainsAndBudgets)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon();
  for (std::size_t n : {64, 128, 1024, 100000})
  {
    auto const s = ferro::bound_pair<R>(n);
    ASSERT_TRUE(s.reference.converged) << n;
    EXPECT_GT(s.tl_energy, R{5} / R{3});
    EXPECT_LT(s.tl_energy - R{5} / R{3}, R{5} / (R(n) * R(n)));
    EXPECT_LE(s.reference.residual_norm, R{32} * eps);
    EXPECT_LT(s.reference.iterations, 30);
    EXPECT_FALSE(s.multiplicity);
    if (n == 100000) EXPECT_EQ(std::exp(-s.reference.log_deviation), R{0});
    // Other edge of the family: a->0, where direct subtraction of cosh
    // terms loses digits in the energy. The ideal limit is accurate here
    // because the log-deviation grows as N*log(N).
    auto const top = ferro::bound_pair<R>(n, n - 3);
    ASSERT_TRUE(top.reference.converged) << n << " status=" << int(top.reference.status)
                                         << " residual=" << uni20::format_real(top.reference.residual_norm);
    R const a = top.reference.center;
    R const ideal_gap = (R{15} / R{2}) / (R{7} / R{2} - std::cos(a));
    EXPECT_REAL_NEAR(top.tl_energy, ideal_gap, R{128} * eps);
  }
  for (std::size_t budget : {0, 1})
  {
    auto const s = ferro::bound_pair<R>(16, 2, {.max_iterations = budget});
    EXPECT_FALSE(s.reference.converged);
    EXPECT_EQ(s.reference.iterations, budget);
    EXPECT_EQ(s.reference.status, bethe::xxz::quantum_group::SolveStatus::iteration_limit);
    ts::detail::System<R> system(16, R{3} / R{2}, 2, 12, -1);
    std::vector<R> x{s.reference.center, s.reference.log_deviation};
    test_support::expect_exact(s.reference.residual_norm, system.evaluate(x).norm, "final finite-deviation residual");
    test_support::expect_exact(s.tl_energy, -R{2} * system.energy_shift(x, R{3} / R{2}), "gap estimate");
  }
  auto strict = ferro::bound_pair<R>(16, 2, {.residual_tolerance = eps * eps, .max_iterations = 30});
  // A rounded residual can be exactly zero; an overly strict request need
  // not fail, but must never be silently relaxed to the default tolerance.
  if (strict.reference.converged) EXPECT_LE(strict.reference.residual_norm, eps * eps);
  for (std::size_t n : {0, 1, 2, 3})
    EXPECT_THROW((void)ferro::bound_pair<R>(n), std::invalid_argument);
  for (std::size_t mode : {0, 6})
    EXPECT_THROW((void)ferro::bound_pair<R>(8, mode), std::invalid_argument);
  if constexpr (uni20::numeric_limits<R>::digits < std::numeric_limits<std::size_t>::digits)
    EXPECT_THROW((void)ferro::bound_pair<R>(std::numeric_limits<std::size_t>::max()), std::overflow_error);
  for (R bad : {R{0}, -R{1}, uni20::numeric_limits<R>::infinity(), uni20::numeric_limits<R>::quiet_NaN()})
  {
    EXPECT_THROW((void)ferro::bound_pair<R>(8, 1, {.residual_tolerance = bad}), std::invalid_argument);
    EXPECT_THROW((void)ts::bound_pair<R>(8, bad), std::invalid_argument);
  }
}

TEST(BiquadraticBoundPairED, EachModeAndModuleMinimumOddAndEven)
{
  for (unsigned n = 4; n <= 10; ++n)
  {
    auto const exact = bethe::test::quantum_group_module_ed(n, n - 4, 1.5);
    auto remaining = exact;
    for (std::size_t mode = 1; mode <= n - 3; ++mode)
    {
      auto const s = ferro::bound_pair<double>(n, mode);
      ASSERT_TRUE(s.reference.converged);
      auto match = std::find_if(remaining.begin(), remaining.end(),
                                [&](double e) { return std::abs(e - s.reference.energy) < 1e-11; });
      ASSERT_NE(match, remaining.end()) << n << " mode=" << mode;
      remaining.erase(match);
      if (mode == 1) EXPECT_NEAR(s.reference.energy, exact.back(), 1e-11);
    }
  }
}
} // namespace
