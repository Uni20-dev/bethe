// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "tl_ed.hpp"
#include <bethe/biquadratic_ferromagnetic.hpp>
#include <complex>

namespace
{
namespace ferro = bethe::biquadratic::ferromagnetic;
namespace ts = bethe::xxz::quantum_group::three_string;
template <typename Real> class BiquadraticBoundTriple : public ::testing::Test {};
TYPED_TEST_SUITE(BiquadraticBoundTriple, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(BiquadraticBoundTriple, SixSiteNativePrecision)
{
  using R = TypeParam;
  // Independent 5-dimensional TL link-pattern matrix at loop weight 3:
  // charpoly=(g-6)(g-7)(g^3-17g^2+80g-106). Smallest root is the droplet.
  R exact = R{23} / R{10};
  for (int i = 0; i < 12; ++i)
    exact -= (((exact - R{17}) * exact + R{80}) * exact - R{106}) / ((R{3} * exact - R{34}) * exact + R{80});
  R const eps = uni20::numeric_limits<R>::epsilon();
  auto const s = ferro::bound_triple<R>(6);
  ASSERT_TRUE(s.reference.converged);
  EXPECT_REAL_NEAR(s.tl_energy, exact, R{128} * eps);
  EXPECT_REAL_NEAR(s.energy, R{5} + exact, R{128} * eps);
  EXPECT_EQ(s.multiplicity, 1);
  if constexpr (uni20::numeric_limits<R>::digits > 53) EXPECT_GT(std::abs(R(double(exact)) - exact), R{128} * eps);
}

TYPED_TEST(BiquadraticBoundTriple, FamilyOddAndEven)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon();
  for (std::size_t n : {6, 7, 8, 9, 10, 16, 32})
  {
    R previous{2};
    for (std::size_t mode = 1; mode <= n - 5; ++mode)
    {
      auto const s = ferro::bound_triple<R>(n, mode);
      ASSERT_TRUE(s.reference.converged) << n << " mode=" << mode << " status=" << int(s.reference.status)
                                         << " residual=" << uni20::format_real(s.reference.residual_norm);
      EXPECT_EQ(s.through_lines, n - 6);
      EXPECT_EQ(s.reference.string_label, n - 4 - mode);
      EXPECT_GT(s.tl_energy, previous);
      EXPECT_LT(s.tl_energy, R{5} / R{2});
      EXPECT_REAL_NEAR(s.energy, R(n - 1) + s.tl_energy, R{128} * R(n) * eps);
      EXPECT_REAL_NEAR(s.energy, R{7} * R(n - 1) / R{4} - R{2} * s.reference.energy, R{128} * R(n) * eps);
      previous = s.tl_energy;
    }
  }
}

TYPED_TEST(BiquadraticBoundTriple, ComplexDeviationJacobian)
{
  using R = TypeParam;
  R const h = std::cbrt(uni20::numeric_limits<R>::epsilon());
  for (R phi : {-R{2}, R{1}})
    for (bool ideal : {false, true})
    {
      ts::detail::System<R> system(16, R{3} / R{2}, 11);
      std::vector<R> x{R{2}, R{4}, phi}, jac;
      system.evaluate(x, &jac, ideal);
      for (std::size_t j = 0; j < 3; ++j)
      {
        auto plus = x, minus = x;
        plus[j] += h;
        minus[j] -= h;
        auto a = system.evaluate(plus, nullptr, ideal), b = system.evaluate(minus, nullptr, ideal);
        for (std::size_t i = 0; i < 3; ++i)
          EXPECT_REAL_NEAR((a.residual[i] - b.residual[i]) / (R{2} * h), jac[3 * i + j], R{4096} * h * h);
      }
    }
}

TYPED_TEST(BiquadraticBoundTriple, OriginalComplexEquations)
{
  using R = TypeParam;
  using C = std::complex<R>;
  R const eps = uni20::numeric_limits<R>::epsilon();
  for (std::size_t n : {6, 7, 8, 16, 64, 129, 1024})
    for (std::size_t mode = 1; mode <= std::min(n - 5, std::size_t{4}); ++mode)
    {
      auto const s = ferro::bound_triple<R>(n, mode).reference;
      ASSERT_TRUE(s.converged) << n << " " << mode;
      R const eta = std::acosh(s.delta), a = s.center;
      C const z = std::exp(-s.log_deviation) * C{std::cos(s.deviation_phase), std::sin(s.deviation_phase)};
      C const u = C{eta, a / R{2}} + z, v{0, a / R{2}};
      C const drive = std::sinh(u + eta / R{2}) / std::sinh(u - eta / R{2});
      // Original u+ equation. Only the singular denominator sinh(z) is
      // reconstructed logarithmically; all other scattering is direct.
      C const regular = std::sinh(R{2} * eta + z) * std::sinh(u + v + eta) / std::sinh(u + v - eta) *
                        std::sinh(u - std::conj(u) + eta) / std::sinh(u - std::conj(u) - eta) *
                        std::sinh(u + std::conj(u) + eta) / std::sinh(u + std::conj(u) - eta);
      C const corr = std::abs(z) < std::sqrt(eps) ? z * z / R{6} : std::log(std::sinh(z) / z);
      R const phase = R{2} * R(n) * std::arg(drive) - std::arg(regular) + s.deviation_phase + corr.imag();
      EXPECT_REAL_NEAR(std::sin(phase / R{2}), R{0}, R{256} * R(n) * eps);
      R const magnitude =
          std::log(std::abs(drive)) - (std::log(std::abs(regular)) + s.log_deviation - corr.real()) / (R{2} * R(n));
      EXPECT_REAL_NEAR(magnitude, R{0}, R{128} * eps);
      // Central-root equation has two singular factors which cancel. Use
      // the exact z/-z ratio, avoiding subtraction of rounded string roots.
      C const cdrive = std::sinh(v + eta / R{2}) / std::sinh(v - eta / R{2});
      C const b = std::sinh(u + v + eta) / std::sinh(u + v - eta);
      R const cphase = R{2} * R(n) * std::arg(cdrive) -
                       R{2} * (s.deviation_phase + corr.imag() - std::arg(std::sinh(R{2} * eta + z)) + std::arg(b));
      EXPECT_REAL_NEAR(std::sin(cphase / R{2}), R{0}, R{512} * R(n) * eps);
      // Direct all-root products are independently useful while z is well
      // resolved. Below that scale they are deliberately not an oracle.
      if (n <= 8 && mode == 1)
      {
        std::array<C, 3> roots{v, u, std::conj(u)};
        for (std::size_t j = 0; j < 3; ++j)
        {
          C lhs{1}, rhs{1};
          C const d = std::sinh(roots[j] + eta / R{2}) / std::sinh(roots[j] - eta / R{2});
          for (std::size_t k = 0; k < 2 * n; ++k)
            lhs *= d;
          for (std::size_t k = 0; k < 3; ++k)
            if (k != j)
              rhs *= std::sinh(roots[j] - roots[k] + eta) / std::sinh(roots[j] - roots[k] - eta) *
                     std::sinh(roots[j] + roots[k] + eta) / std::sinh(roots[j] + roots[k] - eta);
          EXPECT_REAL_NEAR(std::abs(lhs / rhs - C{1}), R{0}, R{4096} * eps / std::abs(z));
        }
      }
    }
}

TYPED_TEST(BiquadraticBoundTriple, LongChainsAndBudgets)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon();
  for (std::size_t n : {64, 128, 1024, 100000})
    for (std::size_t mode : {std::size_t{1}, n - 5})
    {
      auto const s = ferro::bound_triple<R>(n, mode);
      ASSERT_TRUE(s.reference.converged) << n << " " << mode << " status=" << int(s.reference.status);
      EXPECT_GT(s.tl_energy, R{2});
      if (mode == 1) EXPECT_LT(s.tl_energy - R{2}, R{2} / (R(n) * R(n)));
      EXPECT_LE(s.reference.residual_norm, R{32} * eps);
      EXPECT_LT(s.reference.iterations, 40);
      EXPECT_FALSE(s.multiplicity);
      if (n == 100000) EXPECT_EQ(std::exp(-s.reference.log_deviation), R{0});
      R const ideal_gap = R{20} / (R{9} - std::cos(s.reference.center));
      EXPECT_REAL_NEAR(s.tl_energy, ideal_gap, R{256} * eps);
    }
  for (std::size_t budget : {0, 1})
  {
    auto const s = ferro::bound_triple<R>(16, 2, {.max_iterations = budget});
    EXPECT_FALSE(s.reference.converged);
    EXPECT_EQ(s.reference.iterations, budget);
    EXPECT_EQ(s.reference.status, bethe::xxz::quantum_group::SolveStatus::iteration_limit);
    ts::detail::System<R> system(16, R{3} / R{2}, 10);
    std::vector<R> x{s.reference.center, s.reference.log_deviation, s.reference.deviation_phase};
    test_support::expect_exact(s.reference.residual_norm, system.evaluate(x).norm, "final finite-deviation residual");
    test_support::expect_exact(s.tl_energy, -R{2} * system.energy_shift(x), "gap estimate");
  }
  auto strict = ferro::bound_triple<R>(16, 2, {.residual_tolerance = eps * eps, .max_iterations = 30});
  if (strict.reference.converged) EXPECT_LE(strict.reference.residual_norm, eps * eps);
  for (std::size_t n : {0, 1, 2, 3, 4, 5})
    EXPECT_THROW((void)ferro::bound_triple<R>(n), std::invalid_argument);
  for (std::size_t mode : {0, 4})
    EXPECT_THROW((void)ferro::bound_triple<R>(8, mode), std::invalid_argument);
  if constexpr (uni20::numeric_limits<R>::digits < std::numeric_limits<std::size_t>::digits)
    EXPECT_THROW((void)ferro::bound_triple<R>(std::numeric_limits<std::size_t>::max()), std::overflow_error);
  for (R bad : {R{0}, -R{1}, uni20::numeric_limits<R>::infinity(), uni20::numeric_limits<R>::quiet_NaN()})
  {
    EXPECT_THROW((void)ferro::bound_triple<R>(8, 1, {.residual_tolerance = bad}), std::invalid_argument);
    EXPECT_THROW((void)ts::bound_triple<R>(8, bad), std::invalid_argument);
  }
}

TEST(BiquadraticBoundTripleED, EachModeAndModuleMinimumOddAndEven)
{
  for (double delta : {1.25, 1.5, 2.0, 3.0})
    for (unsigned n = 6; n <= 10; ++n)
    {
      auto const exact = bethe::test::quantum_group_module_ed(n, n - 6, delta);
      auto remaining = exact;
      for (std::size_t mode = 1; mode <= n - 5; ++mode)
      {
        auto const s = ts::bound_triple<double>(n, delta, mode);
        ASSERT_TRUE(s.converged) << n << " " << delta << " " << mode;
        auto match =
            std::find_if(remaining.begin(), remaining.end(), [&](double e) { return std::abs(e - s.energy) < 1e-11; });
        ASSERT_NE(match, remaining.end()) << n << " mode=" << mode << " energy=" << s.energy;
        remaining.erase(match);
        if (mode == 1) EXPECT_NEAR(s.energy, exact.back(), 1e-11);
      }
    }
}
} // namespace
