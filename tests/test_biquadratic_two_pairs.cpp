// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "tl_ed.hpp"
#include <bethe/biquadratic_ferromagnetic.hpp>
#include <complex>

namespace
{
namespace tp = bethe::xxz::quantum_group::two_pairs;
namespace ferro = bethe::biquadratic::ferromagnetic;
template <typename R> class TwoPairs : public ::testing::Test {};
TYPED_TEST_SUITE(TwoPairs, test_support::RealTypes, test_support::PrecisionNames);
TEST(TwoPairsED, SmallModules)
{
  for (double delta : {1.25, 1.5, 2.0, 3.0})
    for (unsigned n = 8; n <= 10; ++n)
    {
      auto ed = bethe::test::quantum_group_module_ed(n, n - 8, delta);
      for (std::size_t i = 1; i < n - 6; ++i)
        for (std::size_t j = i + 1; j <= n - 6; ++j)
        {
          auto const s = tp::solve<double>(n, delta, {i, j});
          ASSERT_TRUE(s.converged) << n << " " << i << " " << j << " " << int(s.status);
          auto nearest = std::min_element(ed.begin(), ed.end(), [&](double a, double b) {
            return std::abs(a - s.energy) < std::abs(b - s.energy);
          });
          ASSERT_NE(nearest, ed.end());
          EXPECT_NEAR(*nearest, s.energy, 2e-11) << n << " " << i << " " << j;
          ed.erase(nearest);
        }
    }
}

TYPED_TEST(TwoPairs, OriginalComplexEquationsAndSigns)
{
  using R = TypeParam;
  using C = std::complex<R>;
  R const eps = uni20::numeric_limits<R>::epsilon(), pi = R{4} * std::atan(R{1});
  for (std::size_t n : {8, 9, 10, 16, 65, 128, 1024})
    for (auto labels : {std::array<std::size_t, 2>{1, 2}, {n - 7, n - 6}, {1, n - 6}})
    {
      auto const s = tp::solve<R>(n, R{3} / R{2}, labels);
      ASSERT_TRUE(s.converged) << n << " " << labels[0] << " " << labels[1] << " " << int(s.status);
      EXPECT_LT(s.centers[0], s.centers[1]);
      EXPECT_EQ(s.deviation_signs[0], (n - labels[0]) % 2 ? 1 : -1);
      EXPECT_EQ(s.deviation_signs[1], (n - labels[1]) % 2 ? -1 : 1);
      R const eta = std::acosh(s.delta);
      std::array<R, 2> d;
      std::vector<C> roots;
      for (std::size_t i = 0; i < 2; ++i)
      {
        d[i] = R(s.deviation_signs[i]) * std::exp(-s.log_deviations[i]);
        roots.emplace_back((eta + d[i]) / R{2}, s.centers[i] / R{2});
        roots.emplace_back((eta + d[i]) / R{2}, -s.centers[i] / R{2});
      }
      for (std::size_t k = 0; k < 4; ++k)
      {
        bethe::detail::CompensatedSum<R> magnitude, phase;
        auto factor = [&](C v, R weight) {
          magnitude.add(weight * std::log(std::abs(v)));
          phase.add(weight * std::arg(v));
        };
        factor(std::sinh(roots[k] + eta / R{2}), R{2} * R(n));
        factor(std::sinh(roots[k] - eta / R{2}), -R{2} * R(n));
        for (std::size_t l = 0; l < 4; ++l)
          if (k != l)
          {
            factor(std::sinh(roots[k] - roots[l] + eta), -R{1});
            factor(std::sinh(roots[k] + roots[l] + eta), -R{1});
            factor(std::sinh(roots[k] - roots[l] - eta), R{1});
            if (k / 2 == l / 2)
            {
              auto const p = k / 2;
              magnitude.add(-s.log_deviations[p] + tp::detail::System<R>::Pair::sinh_correction(d[p]).first);
              phase.add(s.deviation_signs[p] < 0 ? pi : R{0});
            }
            else
              factor(std::sinh(roots[k] + roots[l] - eta), R{1});
          }
        EXPECT_REAL_NEAR(magnitude.value() / (R{2} * R(n)), R{0}, R{512} * eps);
        EXPECT_REAL_NEAR(std::sin(phase.value() / R{2}), R{0}, R{512} * R(n) * eps);
      }
    }
}

TEST(TwoPairsED, KnownFourDefectFamiliesAreDistinctSubsets)
{
  for (unsigned n = 8; n <= 10; ++n)
  {
    auto ed = bethe::test::quantum_group_module_ed(n, n - 8, 1.5);
    auto remove = [&](double energy) {
      auto nearest = std::min_element(ed.begin(), ed.end(),
                                      [&](double a, double b) { return std::abs(a - energy) < std::abs(b - energy); });
      ASSERT_NE(nearest, ed.end());
      ASSERT_NEAR(*nearest, energy, 2e-11);
      ed.erase(nearest);
    };
    auto const real = ferro::real_excitations<double>(n, n - 8, {.count = 1000});
    ASSERT_TRUE(real.converged());
    for (auto const& s : real.levels)
      ASSERT_NO_FATAL_FAILURE(remove(s.state.reference.energy));
    for (std::size_t i = 1; i < n - 4; ++i)
      for (std::size_t k = i + 1; k <= n - 4; ++k)
        for (std::size_t j = 1; j <= n - 7; ++j)
        {
          auto const s = ferro::pair_with_real_roots<double>(n, std::array<std::size_t, 2>{i, k}, j);
          ASSERT_TRUE(s.reference.converged);
          ASSERT_NO_FATAL_FAILURE(remove(s.reference.energy));
        }
    for (std::size_t i = 1; i < n - 6; ++i)
      for (std::size_t j = i + 1; j <= n - 6; ++j)
      {
        auto const s = ferro::two_bound_pairs<double>(n, {i, j});
        ASSERT_TRUE(s.reference.converged);
        ASSERT_NO_FATAL_FAILURE(remove(s.reference.energy));
      }
    EXPECT_EQ(ed.size(), (n - 2) * (n - 7)); // Other string topologies remain missing.
  }
}

TYPED_TEST(TwoPairs, CoupledJacobian)
{
  using R = TypeParam;
  R const h = std::cbrt(uni20::numeric_limits<R>::epsilon());
  for (auto labels : {std::array<std::size_t, 2>{1, 2}, {1, 3}, {2, 4}, {2, 3}})
    for (bool ideal : {false, true})
    {
      tp::detail::System<R> system(16, R{3} / R{2}, labels);
      std::vector<R> x{R{1}, R{3}, R{2}, R{4}}, jac;
      system.evaluate(x, &jac, ideal);
      for (std::size_t j = 0; j < 4; ++j)
      {
        auto a = x, b = x;
        a[j] += h;
        b[j] -= h;
        auto const fa = system.evaluate(a, nullptr, ideal), fb = system.evaluate(b, nullptr, ideal);
        for (std::size_t i = 0; i < 4; ++i)
          EXPECT_REAL_NEAR((fa.residual[i] - fb.residual[i]) / (R{2} * h), jac[4 * i + j], R{4096} * h * h);
      }
    }
}

TYPED_TEST(TwoPairs, LongChainThresholdAndFailureContract)
{
  using R = TypeParam;
  R previous{100};
  for (std::size_t n : {8, 16, 65, 128, 1024, 100000})
  {
    auto const s = ferro::two_bound_pairs<R>(n, {n - 7, n - 6});
    ASSERT_TRUE(s.reference.converged) << n;
    EXPECT_EQ(s.through_lines, n - 8);
    EXPECT_GT(s.tl_energy, R{10} / R{3});
    EXPECT_LT(s.tl_energy, previous);
    // An asymptotic scaling check, not a bound on the shortest N=8 chain.
    if (n >= 16) EXPECT_LT(s.tl_energy - R{10} / R{3}, R{100} / (R(n) * R(n)));
    EXPECT_EQ(s.tl_energy, -R{2} * s.reference.energy_shift);
    previous = s.tl_energy;
    if (n == 100000)
      for (R L : s.reference.log_deviations)
      {
        EXPECT_EQ(std::exp(-L), R{0});
        EXPECT_TRUE(uni20::isfinite(L));
      }
  }
  for (std::size_t budget : {0, 1})
  {
    auto const s = tp::solve<R>(16, R{3} / R{2}, {9, 10}, {.max_iterations = budget});
    EXPECT_FALSE(s.converged);
    EXPECT_EQ(s.iterations, budget);
    EXPECT_EQ(s.status, bethe::xxz::quantum_group::SolveStatus::iteration_limit);
    tp::detail::System<R> system(16, s.delta, s.string_labels);
    std::vector<R> x{s.centers[0], s.log_deviations[0], s.centers[1], s.log_deviations[1]};
    EXPECT_EQ(s.residual_norm, system.evaluate(x).norm);
    EXPECT_EQ(s.energy_shift, system.energy_shift(x, s.delta));
  }
  for (std::size_t n : {0, 2, 6, 7})
    EXPECT_THROW((void)ferro::two_bound_pairs<R>(n, {1, 2}), std::invalid_argument);
  for (auto labels : {std::array<std::size_t, 2>{0, 1}, {1, 1}, {2, 1}, {1, 3}})
    EXPECT_THROW((void)ferro::two_bound_pairs<R>(8, labels), std::invalid_argument);
  for (R bad : {R{0}, R{1}, uni20::numeric_limits<R>::infinity(), uni20::numeric_limits<R>::quiet_NaN()})
    EXPECT_THROW((void)tp::solve<R>(8, bad, {1, 2}), std::invalid_argument);
  EXPECT_THROW((void)ferro::two_bound_pairs<R>(8, {1, 2}, {.residual_tolerance = R{0}}), std::invalid_argument);
  EXPECT_THROW((void)ferro::two_bound_pairs<R>(std::numeric_limits<std::size_t>::max(), {1, 2}), std::invalid_argument);
  if constexpr (uni20::numeric_limits<R>::digits <= 53)
  {
    std::size_t const n = 1000000000000000000ULL;
    EXPECT_THROW((void)ferro::two_bound_pairs<R>(n, {n - 7, n - 6}), std::overflow_error);
  }
}
} // namespace
