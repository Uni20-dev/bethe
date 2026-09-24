// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "tl_ed.hpp"
#include <bethe/biquadratic_ferromagnetic.hpp>
#include <complex>
#include <numeric>

namespace
{
namespace ferro = bethe::biquadratic::ferromagnetic;
namespace ts = bethe::xxz::quantum_group::two_string;
template <typename Real> class BiquadraticPairDefect : public ::testing::Test {};
TYPED_TEST_SUITE(BiquadraticPairDefect, test_support::RealTypes, test_support::PrecisionNames);

template <typename R> void expect_original_equations(ts::State<R> const& s)
{
  using C = std::complex<R>;
  auto const n = s.sites, sea = s.rapidities.size();
  R const eps = uni20::numeric_limits<R>::epsilon(), pi = R{4} * std::atan(R{1});
  R const eta = std::acosh(s.delta), d = R(s.deviation_sign) * std::exp(-s.log_deviation);
  std::vector<C> roots;
  for (R alpha : s.rapidities)
    roots.emplace_back(R{0}, alpha / R{2});
  roots.emplace_back((eta + d) / R{2}, s.center / R{2});
  roots.emplace_back((eta + d) / R{2}, -s.center / R{2});
  for (std::size_t k = 0; k < roots.size(); ++k)
  {
    bethe::detail::CompensatedSum<R> magnitude, phase;
    auto factor = [&](C v, R weight) {
      magnitude.add(weight * std::log(std::abs(v)));
      phase.add(weight * std::arg(v));
    };
    factor(std::sinh(roots[k] + eta / R{2}), R{2} * R(n));
    factor(std::sinh(roots[k] - eta / R{2}), -R{2} * R(n));
    for (std::size_t l = 0; l < roots.size(); ++l)
      if (k != l)
      {
        factor(std::sinh(roots[k] - roots[l] + eta), -R{1});
        factor(std::sinh(roots[k] + roots[l] + eta), -R{1});
        factor(std::sinh(roots[k] - roots[l] - eta), R{1});
        if (k >= sea && l >= sea)
        {
          // Only the singular reflected factor uses L rather than rounded roots.
          magnitude.add(-s.log_deviation + ts::detail::System<R>::sinh_correction(d).first);
          phase.add(s.deviation_sign < 0 ? pi : R{0});
        }
        else
          factor(std::sinh(roots[k] + roots[l] - eta), R{1});
      }
    EXPECT_REAL_NEAR(magnitude.value() / (R{2} * R(n)), R{0}, R{256} * eps);
    EXPECT_REAL_NEAR(std::sin(phase.value() / R{2}), R{0}, R{256} * R(n) * eps);
  }
}

TEST(BiquadraticPairDefectED, PairWithTwoRealRoots)
{
  for (unsigned n = 8; n <= 10; ++n)
  {
    auto ed = bethe::test::quantum_group_module_ed(n, n - 8, 1.5);
    for (std::size_t i = 1; i < n - 4; ++i)
      for (std::size_t k = i + 1; k <= n - 4; ++k)
        for (std::size_t j = 1; j <= n - 7; ++j)
        {
          std::vector<std::size_t> labels{i, k};
          auto const s = ferro::pair_with_real_roots<double>(n, labels, j);
          ASSERT_TRUE(s.reference.converged) << n << " " << i << " " << k << " " << j;
          auto nearest = std::min_element(ed.begin(), ed.end(), [&](double a, double b) {
            return std::abs(a - s.reference.energy) < std::abs(b - s.reference.energy);
          });
          ASSERT_NE(nearest, ed.end());
          EXPECT_NEAR(*nearest, s.reference.energy, 2e-11) << n << " " << i << " " << k << " " << j;
          ed.erase(nearest);
        }
  }
}

TEST(BiquadraticPairDefectED, ThreeFamiliesExhaustSmallModules)
{
  for (unsigned n = 6; n <= 10; ++n)
  {
    auto const ed = bethe::test::quantum_group_module_ed(n, n - 6, 1.5);
    std::vector<double> obtained;
    for (std::size_t i = 1; i <= n - 3; ++i)
      for (std::size_t j = 1; j <= n - 5; ++j)
      {
        auto const s = ferro::pair_defect<double>(n, i, j);
        ASSERT_TRUE(s.reference.converged) << n << " " << i << " " << j;
        EXPECT_EQ(s.through_lines, n - 6);
        obtained.push_back(s.reference.energy);
      }
    for (std::size_t mode = 1; mode <= n - 5; ++mode)
    {
      auto const s = ferro::bound_triple<double>(n, mode);
      ASSERT_TRUE(s.reference.converged);
      obtained.push_back(s.reference.energy);
    }
    auto const real = ferro::real_excitations<double>(n, n - 6, {.count = 1000});
    ASSERT_TRUE(real.family_converged());
    for (auto const& s : real.levels)
      obtained.push_back(s.state.reference.energy);
    std::sort(obtained.begin(), obtained.end());
    ASSERT_EQ(obtained.size(), ed.size()) << n;
    for (std::size_t i = 0; i < ed.size(); ++i)
      EXPECT_NEAR(obtained[i], ed[i], 2e-11) << n << " " << i;
  }
}

TEST(BiquadraticPairDefectED, PairWithThreeRealRoots)
{
  auto ed = bethe::test::quantum_group_module_ed(10, 0, 1.5);
  for (std::size_t i = 1; i <= 3; ++i)
    for (std::size_t j = i + 1; j <= 4; ++j)
      for (std::size_t k = j + 1; k <= 5; ++k)
      {
        std::vector<std::size_t> labels{i, j, k};
        auto const s = ferro::pair_with_real_roots<double>(10, labels, 1);
        ASSERT_TRUE(s.reference.converged);
        auto nearest = std::min_element(ed.begin(), ed.end(), [&](double a, double b) {
          return std::abs(a - s.reference.energy) < std::abs(b - s.reference.energy);
        });
        ASSERT_NE(nearest, ed.end());
        EXPECT_NEAR(*nearest, s.reference.energy, 2e-11);
        ed.erase(nearest);
      }
}

TYPED_TEST(BiquadraticPairDefect, GeneralSeaNativeEquationsAndLimits)
{
  using R = TypeParam;
  for (std::size_t sea : {2, 3, 4})
  {
    R previous = R{100};
    for (std::size_t n : {16, 65, 128, 1024, 100000})
    {
      std::vector<std::size_t> labels(sea);
      auto const m = sea + 2;
      for (std::size_t i = 0; i < sea; ++i)
        labels[i] = n - m - sea + 1 + i;
      for (std::size_t offset : {0, 1})
      {
        auto const s = ferro::pair_with_real_roots<R>(n, labels, n - 2 * m + 1 - offset);
        ASSERT_TRUE(s.reference.converged) << n << " " << sea << " " << offset;
        EXPECT_EQ(s.through_lines, n - 2 * m);
        EXPECT_EQ(s.reference.real_labels, labels);
        EXPECT_EQ(s.reference.deviation_sign, offset % 2 ? -1 : 1);
        // Use modest N for direct complex equations: roots separated by O(1/N)
        // become too close to resolve unscaled residuals on enormous chains.
        if (n <= 128) expect_original_equations(s.reference);
        if (offset == 0)
        {
          R const threshold = R{5} / R{3} + R(sea);
          EXPECT_GT(s.tl_energy, threshold);
          EXPECT_LT(s.tl_energy - threshold, R{10000} / (R(n) * R(n)));
          EXPECT_LT(s.tl_energy, previous);
          previous = s.tl_energy;
        }
      }
    }
  }
}

TYPED_TEST(BiquadraticPairDefect, GeneralSeaReductionsAndValidation)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon();
  for (std::size_t n : {6, 7, 16})
    for (std::size_t j = 1; j <= n - 3; ++j)
    {
      auto const generic = ferro::pair_with_real_roots<R>(n, {}, j);
      auto const pair = ferro::bound_pair<R>(n, n - 2 - j);
      ASSERT_TRUE(generic.reference.converged);
      test_support::expect_exact(generic.tl_energy, pair.tl_energy, "empty sea reproduces bound pair");
    }
  for (std::size_t n : {8, 16, 32})
  {
    std::vector<std::size_t> labels(n / 2 - 2);
    std::iota(labels.begin(), labels.end(), 1);
    auto const generic = ts::pair_with_real_roots<R>(n, R{3} / R{2}, labels, 1);
    auto const singlet = ts::singlet<R>(n, R{3} / R{2});
    ASSERT_TRUE(generic.converged);
    EXPECT_REAL_NEAR(generic.energy, singlet.energy, R{512} * R(n) * eps);
  }
  for (std::vector<std::size_t> labels : {std::vector<std::size_t>{0, 1}, {1, 1}, {2, 1}, {1, 5}, {1, 2, 3}})
    EXPECT_THROW((void)ferro::pair_with_real_roots<R>(8, labels, 1), std::invalid_argument);
  std::vector<std::size_t> labels{3, 4};
  for (std::size_t j : {0, 2})
    EXPECT_THROW((void)ferro::pair_with_real_roots<R>(8, labels, j), std::invalid_argument);
  for (std::size_t n : {0, 1, 2, 3})
    EXPECT_THROW((void)ferro::pair_with_real_roots<R>(n, {}, 1), std::invalid_argument);
  for (std::size_t budget : {0, 1})
  {
    auto const s = ts::pair_with_real_roots<R>(8, R{3} / R{2}, labels, 1, {.max_iterations = budget});
    EXPECT_FALSE(s.converged);
    EXPECT_EQ(s.iterations, budget);
    ts::detail::System<R> system(8, s.delta, labels, 1);
    auto x = s.rapidities;
    x.push_back(s.center);
    x.push_back(s.log_deviation);
    EXPECT_EQ(s.residual_norm, system.evaluate(x).norm);
  }
}

TYPED_TEST(BiquadraticPairDefect, OriginalEquationsAndBothSigns)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon();
  for (std::size_t n : {6, 7, 8, 9, 16, 65, 128})
    for (std::size_t i : {std::size_t{1}, n - 3})
      for (std::size_t j : {std::size_t{1}, n - 5, std::max(std::size_t{1}, n - 6)})
      {
        auto const state = ferro::pair_defect<R>(n, i, j);
        auto const& s = state.reference;
        ASSERT_TRUE(s.converged) << n << " " << i << " " << j << " status=" << int(s.status);
        EXPECT_EQ(s.real_labels, (std::vector<std::size_t>{i}));
        EXPECT_EQ(s.string_label, j);
        EXPECT_EQ(s.deviation_sign, (n - j) % 2 ? 1 : -1);
        EXPECT_REAL_NEAR(state.energy, R(n - 1) + state.tl_energy, R{128} * R(n) * eps);
        EXPECT_EQ(state.tl_energy, -R{2} * s.energy_shift);
        expect_original_equations(s);
      }
}

TYPED_TEST(BiquadraticPairDefect, LongChainThresholdAndNativePrecision)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon();
  R previous = R{100};
  for (std::size_t n : {16, 64, 128, 1024, 100000})
  {
    auto const s = ferro::pair_defect<R>(n, n - 3, n - 5);
    ASSERT_TRUE(s.reference.converged) << n;
    EXPECT_GT(s.tl_energy, R{8} / R{3});
    EXPECT_LT(s.tl_energy, previous);
    EXPECT_LT(s.tl_energy - R{8} / R{3}, R{100} / (R(n) * R(n)));
    EXPECT_LE(s.reference.residual_norm, R{32} * eps);
    if (n == 100000)
    {
      EXPECT_EQ(std::exp(-s.reference.log_deviation), R{0});
      EXPECT_TRUE(uni20::isfinite(s.reference.log_deviation));
    }
    previous = s.tl_energy;
  }
  // N=6 has exact mixed levels g=6,7 and the middle root of this cubic.
  auto const s = ferro::pair_defect<R>(6, 3, 1);
  ASSERT_TRUE(s.reference.converged);
  R const g = s.tl_energy;
  EXPECT_REAL_NEAR(((g - R{17}) * g + R{80}) * g, R{106}, R{4096} * eps);
  // Even this nearly representable root must not take an fp64 round trip.
  if constexpr (uni20::numeric_limits<R>::digits > 53) EXPECT_GT(std::abs(g - R(double(g))), eps);
}

TYPED_TEST(BiquadraticPairDefect, JacobianAndFailureContract)
{
  using R = TypeParam;
  R const h = std::cbrt(uni20::numeric_limits<R>::epsilon());
  for (std::size_t j : {10, 11})
    for (bool ideal : {false, true})
    {
      ts::detail::System<R> system(16, R{3} / R{2}, 13, j);
      std::vector<R> x{R{2}, R{5} / R{2}, R{4}}, jac;
      system.evaluate(x, &jac, ideal);
      for (std::size_t col = 0; col < 3; ++col)
      {
        auto a = x, b = x;
        a[col] += h;
        b[col] -= h;
        auto const fa = system.evaluate(a, nullptr, ideal), fb = system.evaluate(b, nullptr, ideal);
        for (std::size_t row = 0; row < 3; ++row)
          EXPECT_REAL_NEAR((fa.residual[row] - fb.residual[row]) / (R{2} * h), jac[3 * row + col], R{4096} * h * h);
      }
    }
  for (std::size_t budget : {0, 1})
  {
    auto const s = ts::pair_defect<R>(16, R{3} / R{2}, 13, 11, {.max_iterations = budget});
    EXPECT_FALSE(s.converged);
    EXPECT_EQ(s.iterations, budget);
    EXPECT_EQ(s.status, bethe::xxz::quantum_group::SolveStatus::iteration_limit);
    ts::detail::System<R> system(16, s.delta, 13, 11);
    std::vector<R> x{s.rapidities[0], s.center, s.log_deviation};
    EXPECT_EQ(s.residual_norm, system.evaluate(x).norm);
    EXPECT_EQ(s.energy_shift, system.energy_shift(x, s.delta));
  }
  for (std::size_t n : {0, 1, 4, 5})
    EXPECT_THROW((void)ferro::pair_defect<R>(n, 1, 1), std::invalid_argument);
  for (std::size_t i : {0, 6})
    EXPECT_THROW((void)ferro::pair_defect<R>(8, i, 1), std::invalid_argument);
  for (std::size_t j : {0, 4})
    EXPECT_THROW((void)ferro::pair_defect<R>(8, 1, j), std::invalid_argument);
  for (R bad : {R{0}, R{1}, uni20::numeric_limits<R>::infinity(), uni20::numeric_limits<R>::quiet_NaN()})
    EXPECT_THROW((void)ts::pair_defect<R>(8, bad, 1, 1), std::invalid_argument);
  EXPECT_THROW((void)ferro::pair_defect<R>(std::numeric_limits<std::size_t>::max(), 1, 1), std::invalid_argument);
  EXPECT_THROW((void)ferro::pair_defect<R>(8, 1, 1, {.residual_tolerance = R{0}}), std::invalid_argument);
  // Storing the explicit sea labels must not allocate before the shared
  // quadratic-matrix preflight for the pre-existing AF singlet API.
  EXPECT_THROW((void)ts::singlet<R>(std::size_t{1} << 40, R{3} / R{2}), std::length_error);
}
} // namespace
