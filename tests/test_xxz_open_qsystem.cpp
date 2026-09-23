// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "tl_ed.hpp"
#include <bethe/biquadratic_qsystem.hpp>
#include <bethe/xxz_open_qsystem.hpp>

namespace
{
namespace qs = bethe::xxz::quantum_group::qsystem;
template <typename Real> class OpenQSystem : public ::testing::Test {};
TYPED_TEST_SUITE(OpenQSystem, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(OpenQSystem, MissingFourSiteSinglet)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon();
  auto const s = qs::solve<R>(4, R{3} / R{2}, std::vector<R>{R{3} / R{2}, -R{4} / R{3}});
  ASSERT_TRUE(s.converged) << int(s.status) << " residual=" << uni20::format_real(s.residual_norm)
                           << " Bethe=" << uni20::format_real(s.bethe_residual);
  ASSERT_TRUE(s.energy);
  R const physical = R{2} * *s.energy - R{21} / R{4};
  R const exact = (-R{15} + std::sqrt(R{17})) / R{2};
  EXPECT_REAL_NEAR(physical, exact, R{1024} * eps);
  if constexpr (uni20::numeric_limits<R>::digits > 53) EXPECT_GT(std::abs(R(double(exact)) - exact), R{1024} * eps);
  ASSERT_EQ(s.roots.roots.size(), 2);
  EXPECT_GT(std::abs(s.roots.roots[0].imag()), R{1});
  EXPECT_LE(s.residual_norm, R{32} * eps);
  EXPECT_LE(s.bethe_residual, R{4096} * eps);
  EXPECT_REAL_NEAR(s.coefficients[1], -R{4} / R{3}, R{128} * eps);
}

TYPED_TEST(OpenQSystem, AnalyticJacobian)
{
  using R = TypeParam;
  R const h = std::cbrt(uni20::numeric_limits<R>::epsilon());
  for (R delta : {R{11} / R{10}, R{3} / R{2}, R{3}})
    for (std::size_t m : {1, 2, 3})
    {
      qs::System<R> system(6, m, delta);
      std::vector<R> c(m, R{1} / R{4}), jac;
      (void)system.evaluate(c, &jac);
      for (std::size_t j = 0; j < m; ++j)
      {
        auto a = c, b = c;
        a[j] += h;
        b[j] -= h;
        auto plus = system.evaluate(a), minus = system.evaluate(b);
        for (std::size_t i = 0; i < m; ++i)
          EXPECT_REAL_NEAR((plus.residual[i] - minus.residual[i]) / (R{2} * h), jac[i * m + j],
                           R{8192} * h * h * std::max(R{1}, std::abs(jac[i * m + j])));
      }
    }
}

TYPED_TEST(OpenQSystem, ComplexBranchContinuationAndWronskian)
{
  using R = TypeParam;
  using C = std::complex<R>;
  R const eps = uni20::numeric_limits<R>::epsilon();
  std::vector<R> seed{R{3} / R{2}, -R{4} / R{3}};
  for (R delta : {R{3} / R{2}, R{7} / R{4}, R{2}, R{5} / R{2}})
  {
    auto const state = qs::solve<R>(4, delta, seed);
    ASSERT_TRUE(state.converged) << int(state.status);
    EXPECT_REAL_NEAR(*state.energy, -R{3} * delta / R{4} + std::sqrt(delta * delta + R{2}) / R{2}, R{8192} * eps);
    seed = state.coefficients;
    auto q = seed;
    q.push_back(R{1});
    auto at = [](std::vector<R> const& c, C x) {
      C value{};
      for (std::size_t k = c.size(); k-- > 0;)
        value = value * x + c[k];
      return value;
    };
    C const u{R{1} / R{5}, R{7} / R{10}};
    R const eta = std::acosh(delta);
    C const xp = std::cosh(R{2} * u + eta), xm = std::cosh(R{2} * u - eta), x = std::cosh(R{2} * u);
    C const lhs = at(state.partner, xp) * at(q, xm) - at(state.partner, xm) * at(q, xp);
    C rhs = R{2} * std::sinh(eta) * std::sinh(R{2} * u);
    for (std::size_t k = 0; k < 4; ++k)
      rhs *= x - C{1};
    EXPECT_REAL_NEAR(std::abs(lhs - rhs), R{0}, R{8192} * eps * std::max(R{1}, std::abs(rhs)));
  }
}

TYPED_TEST(OpenQSystem, InvalidInputsBudgetsAndAdmissibility)
{
  using R = TypeParam;
  std::vector<R> const seed{R{3} / R{2}, -R{4} / R{3}};
  for (std::size_t n : {0, 1, 3, 5})
    EXPECT_THROW((void)qs::solve<R>(n, R{3} / R{2}, seed), std::invalid_argument);
  EXPECT_THROW((void)qs::solve<R>(2, R{3} / R{2}, seed), std::invalid_argument);
  EXPECT_THROW((void)qs::spectrum<R>(4, R{3} / R{2}, 1), std::invalid_argument);
  for (R value : {R{0}, -R{1}, uni20::numeric_limits<R>::infinity(), uni20::numeric_limits<R>::quiet_NaN()})
  {
    EXPECT_THROW((void)qs::solve<R>(4, value, seed), std::invalid_argument);
    EXPECT_THROW((void)qs::solve<R>(4, R{3} / R{2}, seed, {.residual_tolerance = value}), std::invalid_argument);
  }
  EXPECT_THROW((void)qs::solve<R>(4, R{1}, seed), std::invalid_argument);
  EXPECT_THROW((void)qs::solve<R>(4, R{3} / R{2}, std::vector<R>{uni20::numeric_limits<R>::quiet_NaN()}),
               std::invalid_argument);
  for (std::size_t budget : {0, 1})
  {
    auto const state = qs::solve<R>(4, R{3} / R{2}, seed, {.max_iterations = budget});
    EXPECT_FALSE(state.converged);
    EXPECT_EQ(state.status, qs::Status::iteration_limit);
    EXPECT_EQ(state.iterations, budget);
    EXPECT_EQ(state.residual_norm, qs::System<R>(4, 2, R{3} / R{2}).evaluate(state.coefficients).norm);
  }
  auto const strict =
      qs::solve<R>(4, R{3} / R{2}, seed,
                   {.residual_tolerance = uni20::numeric_limits<R>::epsilon() * uni20::numeric_limits<R>::epsilon(),
                    .max_iterations = 20});
  EXPECT_FALSE(strict.converged);
  auto const loose = qs::solve<R>(4, R{3} / R{2}, std::vector<R>{R{0}}, {.residual_tolerance = R{1}});
  EXPECT_FALSE(loose.converged); // Coefficient tolerance cannot bypass physical validation.
  EXPECT_EQ(loose.status, qs::Status::inadmissible);
  for (R root : {R{-1}, R{1}, R{3} / R{2}})
  {
    qs::State<R> state;
    state.sites = 4;
    state.delta = R{3} / R{2};
    state.energy = R{0};
    state.coefficients = {-root};
    qs::detail::check_roots(state);
    EXPECT_FALSE(state.converged);
    EXPECT_EQ(state.status, qs::Status::inadmissible);
  }
  auto const pole = qs::solve<R>(4, R{3} / R{2}, std::vector<R>{-R{3} / R{2}}, {.max_iterations = 0});
  EXPECT_FALSE(pole.energy);
  auto const zero = qs::spectrum<R>(4, R{3} / R{2}, 0, {.max_attempts = 0});
  EXPECT_FALSE(zero.complete());
  EXPECT_EQ(zero.attempts, 0);
  EXPECT_TRUE(zero.states.empty());
  auto const limited = qs::spectrum<R>(4, R{3} / R{2}, 0, {.max_attempts = 2}, {.max_iterations = 0});
  EXPECT_FALSE(limited.complete());
  EXPECT_EQ(limited.attempts, 2);
  EXPECT_EQ(limited.failed_attempts, 2);
  auto const vacuum = qs::solve<R>(4, R{3} / R{2}, std::vector<R>{});
  ASSERT_TRUE(vacuum.converged);
  EXPECT_EQ(*vacuum.energy, R{9} / R{8});
  auto const larger = qs::solve<R>(34, R{3} / R{2}, std::vector<R>{});
  ASSERT_TRUE(larger.converged);
  EXPECT_EQ(*larger.energy, R{99} / R{8});
  auto const budgeted = qs::spectrum<R>(12, R{3} / R{2}, 0, {.max_attempts = 0});
  EXPECT_EQ(budgeted.expected_count, 132);
  EXPECT_FALSE(budgeted.complete());
  auto const overflow = qs::spectrum<R>(74, R{3} / R{2}, 0, {.max_attempts = 0});
  EXPECT_FALSE(overflow.expected_count);
  EXPECT_FALSE(overflow.complete());
  EXPECT_EQ(overflow.attempts, 0);
  EXPECT_THROW((void)qs::System<R>(1000000000000000000ULL, 0, R{3} / R{2}), std::length_error);
}

TEST(OpenQSystemSearch, ModuleDimensionWithoutBinomialOverflow)
{
  EXPECT_EQ(qs::module_dimension(10, 0), 42);
  EXPECT_EQ(qs::module_dimension(1000, 1000), 1);
  EXPECT_EQ(qs::module_dimension(1000, 998), 999);
  if constexpr (std::numeric_limits<std::size_t>::digits == 64)
  {
    EXPECT_EQ(qs::module_dimension(70, 0), 3116285494907301262ULL);
    EXPECT_EQ(qs::module_dimension(72, 0), 11959798385860453492ULL);
    EXPECT_FALSE(qs::module_dimension(74, 0));
  }
  auto const found = qs::spectrum<double>(10, 1.5, 8, {.max_attempts = 4000}, {.max_iterations = 100});
  ASSERT_TRUE(found.complete());
  EXPECT_EQ(found.states.size(), 9);
  auto const exact = bethe::test::quantum_group_module_ed(10, 8, 1.5);
  for (std::size_t i = 0; i < exact.size(); ++i)
    EXPECT_NEAR(*found.states[i].energy, exact[i], 1e-10);
}

TEST(OpenQSystemSearch, SmallModuleSpectra)
{
  for (unsigned n : {2, 4, 6, 8})
    for (unsigned ell = 0; ell <= n; ell += 2)
    {
      SCOPED_TRACE(::testing::Message() << "N=" << n << " ell=" << ell);
      auto const found = qs::spectrum<double>(n, 1.5, ell, {.max_attempts = 12000}, {.max_iterations = 100});
      auto const exact = bethe::test::quantum_group_module_ed(n, ell, 1.5);
      EXPECT_EQ(found.expected_count, exact.size());
      if (n <= 6) ASSERT_TRUE(found.complete()) << found.states.size() << "/" << found.expected_count.value_or(0);
      auto available = exact;
      for (auto const& state : found.states)
      {
        auto match = std::find_if(available.begin(), available.end(),
                                  [&](double e) { return std::abs(e - *state.energy) < 1e-8; });
        ASSERT_NE(match, available.end());
        available.erase(match);
      }
      if (n == 8 && uni20::numeric_limits<long double>::digits > 53)
      {
        auto const higher = qs::spectrum<long double>(n, 1.5L, ell, {.max_attempts = 12000}, {.max_iterations = 100});
        ASSERT_TRUE(higher.complete()) << higher.states.size() << "/" << higher.expected_count.value_or(0);
        for (std::size_t i = 0; i < exact.size(); ++i)
          EXPECT_NEAR(double(*higher.states[i].energy), exact[i], 1e-10);
      }
    }
}

TEST(OpenQSystemSearch, EntirePhysicalSpinOneSpectrum)
{
  for (unsigned n : {2, 4, 6})
  {
    std::vector<double> reconstructed;
    for (unsigned ell = 0; ell <= n; ell += 2)
    {
      auto const scan = bethe::biquadratic::qsystem::spectrum<double>(n, ell);
      ASSERT_TRUE(scan.complete());
      for (auto const& state : scan.states)
      {
        ASSERT_TRUE(state.multiplicity);
        for (std::uint64_t k = 0; k < *state.multiplicity; ++k)
          reconstructed.push_back(*state.energy);
      }
    }
    auto const exact = bethe::test::biquadratic_ed(n);
    ASSERT_EQ(reconstructed.size(), exact.size());
    std::sort(reconstructed.begin(), reconstructed.end());
    for (std::size_t i = 0; i < exact.size(); ++i)
      EXPECT_NEAR(reconstructed[i], exact[i], 1e-9);
  }
}
} // namespace
