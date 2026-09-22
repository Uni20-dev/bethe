// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "tb_ed.hpp"
#include "test_support.hpp"
#include <bethe/takhtajan_babujian.hpp>

namespace
{
namespace model = bethe::takhtajan_babujian;
template <typename Real> class TB : public ::testing::Test {};
TYPED_TEST_SUITE(TB, test_support::RealTypes, test_support::PrecisionNames);

template <uni20::Real Real> void check_original(model::State<Real> const& s)
{
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  C const i{Real{0}, Real{1}};
  C energy{}, momentum{Real{1}, Real{0}};
  ASSERT_EQ(s.rapidities.size(), s.sites);
  for (std::size_t j = 0; j < s.sites; ++j)
  {
    C const z = s.rapidities[j], q = (z + i) / (z - i);
    C product{Real{1}, Real{0}};
    for (std::size_t site = 0; site < s.sites; ++site)
      product *= q;
    for (std::size_t k = 0; k < s.sites; ++k)
      if (j != k) product *= (z - s.rapidities[k] - i) / (z - s.rapidities[k] + i);
    EXPECT_REAL_NEAR(product.real(), Real{1}, Real{512} * Real(s.sites) * eps);
    EXPECT_REAL_NEAR(product.imag(), Real{0}, Real{512} * Real(s.sites) * eps);
    energy -= Real{4} / (Real{1} + z * z);
    momentum *= q;
  }
  EXPECT_REAL_NEAR(energy.real(), s.energy, Real{128} * Real(s.sites) * eps);
  EXPECT_REAL_NEAR(energy.imag(), Real{0}, Real{128} * Real(s.sites) * eps);
  EXPECT_REAL_NEAR(momentum.real(), Real{1}, Real{128} * Real(s.sites) * eps);
  EXPECT_REAL_NEAR(momentum.imag(), Real{0}, Real{128} * Real(s.sites) * eps);
}

TYPED_TEST(TB, FiniteDeviationsAndOriginalEquations)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {4, 6, 8, 10, 12, 24, 48, 64, 96, 128})
  {
    SCOPED_TRACE(n);
    auto const s = model::ground_state<Real>(n);
    ASSERT_TRUE(s.converged) << int(s.status) << ", iterations=" << s.iterations;
    EXPECT_EQ(s.status, model::SolveStatus::converged);
    EXPECT_LE(s.residual_norm, Real{32} * eps);
    EXPECT_EQ(s.residual_norm, std::max(s.phase_residual, s.modulus_residual));
    EXPECT_EQ(s.momentum_index, 0u);
    EXPECT_EQ(s.momentum, Real{0});
    for (std::size_t j = 0; j < n / 2; ++j)
    {
      EXPECT_EQ(s.centers[j], -s.centers[n / 2 - 1 - j]);
      EXPECT_EQ(s.deviations[j], s.deviations[n / 2 - 1 - j]);
      EXPECT_GT(s.deviations[j], Real{0});
      EXPECT_LT(s.deviations[j], Real{1} / Real{2});
      EXPECT_EQ(s.string_quantum_numbers[j].twice(), 2 * std::int64_t(j) - std::int64_t(n / 2 - 1));
      if (j) EXPECT_GT(s.centers[j], s.centers[j - 1]);
      EXPECT_TRUE(s.rapidities[2 * j] == std::conj(s.rapidities[2 * j + 1]));
    }
    ASSERT_NO_FATAL_FAILURE(check_original(s));
  }
}

TYPED_TEST(TB, AnalyticFourSitesAndNonidealStrings)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  auto const s = model::ground_state<Real>(4);
  ASSERT_TRUE(s.converged);
  Real const exact = -Real{11} - std::sqrt(Real{41});
  EXPECT_REAL_NEAR(s.energy, exact, Real{512} * eps);
  // Finite deviations are essential, not an output decoration on an ideal-string calculation.
  using C = std::complex<Real>;
  Real ideal_energy{};
  for (Real x : s.centers)
    ideal_energy -= Real{8} * (Real{1} / (Real{1} + C{x, Real{1} / Real{2}} * C{x, Real{1} / Real{2}})).real();
  EXPECT_GT(std::abs(s.energy - ideal_energy), Real{1} / Real{10});
  if constexpr (uni20::numeric_limits<Real>::digits > 53)
    EXPECT_GT(std::abs(Real(double(exact)) - exact), Real{512} * eps);
  test_support::expect_exact(uni20::parse_real<Real>(uni20::format_real(s.energy)), s.energy, "TB energy I/O");
}

TYPED_TEST(TB, IndependentSpinBasisGroundAndTranslation)
{
  using Real = TypeParam;
  for (unsigned n : {4, 6, 8})
  {
    SCOPED_TRACE(n);
    auto const exact = bethe::test::tb_exact_ground(n);
    auto const s = model::ground_state<Real>(n);
    ASSERT_TRUE(s.converged);
    EXPECT_NEAR(double(s.energy), exact.energy, 2e-10);
    EXPECT_NEAR(exact.translation, 1.0, 1e-10);
  }
}

TYPED_TEST(TB, AnalyticJacobian)
{
  using Real = TypeParam;
  Real const h = std::cbrt(uni20::numeric_limits<Real>::epsilon());
  for (std::size_t n : {4, 6, 8})
  {
    model::detail::GroundSystem<Real> const system(n);
    auto x = system.seed();
    std::vector<Real> jac;
    auto const eval = system.evaluate(x, &jac);
    for (std::size_t j = 0; j < x.size(); ++j)
    {
      auto p = x, m = x;
      p[j] += h;
      m[j] -= h;
      auto const plus = system.evaluate(p), minus = system.evaluate(m);
      for (std::size_t k = 0; k < x.size(); ++k)
        EXPECT_REAL_NEAR((plus.residual[k] - minus.residual[k]) / (Real{2} * h), jac[k * x.size() + j],
                         Real{65536} * h * h);
    }
  }
}

TYPED_TEST(TB, BudgetsAndInvalidInputs)
{
  using Real = TypeParam;
  for (std::size_t budget : {0, 1, 2})
  {
    auto const s = model::ground_state<Real>(8, {.max_iterations = budget});
    EXPECT_FALSE(s.converged);
    EXPECT_EQ(s.status, model::SolveStatus::iteration_limit);
    EXPECT_EQ(s.iterations, budget);
    model::detail::GroundSystem<Real> const system(8);
    auto x = system.seed();
    for (std::size_t j = 2; j < 4; ++j)
    {
      x[j - 2] = s.centers[j];
      x[j] = s.deviations[j];
    }
    EXPECT_EQ(system.evaluate(x).norm(), s.residual_norm);
  }
  EXPECT_THROW((void)model::ground_state<Real>(0), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state<Real>(2), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state<Real>(5), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state<Real>(std::numeric_limits<std::size_t>::max()), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state<Real>(1000000000000000000ULL), std::length_error);
  for (Real tol :
       {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW((void)model::ground_state<Real>(4, {.residual_tolerance = tol}), std::invalid_argument);
  std::vector<Real> rhs(2, Real{1});
  EXPECT_FALSE(bethe::detail::newton_step(std::vector<Real>{Real{1}}, rhs));
  EXPECT_FALSE(bethe::detail::newton_step(std::vector<Real>(4, Real{0}), rhs));
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  auto const strict = model::ground_state<Real>(8, {.residual_tolerance = eps * eps, .max_iterations = 30});
  EXPECT_FALSE(strict.converged);
  EXPECT_NE(strict.status, model::SolveStatus::converged);
}
} // namespace
