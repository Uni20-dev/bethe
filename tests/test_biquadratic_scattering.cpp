// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "tl_ed.hpp"
#include <bethe/biquadratic_ferromagnetic.hpp>
#include <complex>

namespace
{
namespace ferro = bethe::biquadratic::ferromagnetic;
namespace qg = bethe::xxz::quantum_group;
template <typename Real> class BiquadraticScattering : public ::testing::Test {};
TYPED_TEST_SUITE(BiquadraticScattering, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(BiquadraticScattering, WindowsAgreeWithFullFamily)
{
  using R = TypeParam;
  for (std::size_t n : {6, 7, 8, 9})
    for (std::size_t m : {1, 2, 3})
    {
      auto const ell = n - 2 * m;
      auto const full = ferro::real_excitations<R>(n, ell, {.count = 1000});
      auto const entire = ferro::real_excitations_window<R>(n, ell, n - m, {.count = 1000});
      ASSERT_TRUE(full.converged());
      ASSERT_TRUE(entire.converged());
      ASSERT_EQ(full.levels.size(), entire.levels.size());
      for (std::size_t i = 0; i < full.levels.size(); ++i)
      {
        EXPECT_EQ(full.levels[i].state.reference.quantum_numbers, entire.levels[i].state.reference.quantum_numbers);
        test_support::expect_exact(*full.levels[i].gap, *entire.levels[i].gap, "full window preserves gap");
      }
      auto const width = std::min(m + 1, n - m);
      auto const window = ferro::real_excitations_window<R>(n, ell, width, {.count = 1000});
      auto const one = ferro::real_excitations_window<R>(n, ell, width, {.count = 1});
      ASSERT_TRUE(window.converged());
      ASSERT_TRUE(one.converged());
      EXPECT_EQ(window.candidate_count, bethe::detail::bounded_binomial(width, m, 1000));
      ASSERT_EQ(one.levels.size(), 1);
      test_support::expect_exact(*one.levels[0].gap, *window.levels[0].gap, "selection precedes truncation");
      std::size_t found{};
      for (auto const& s : full.levels)
        if (s.state.reference.quantum_numbers.front() >= uni20::half_int(std::int64_t(n - m - width + 1)))
        {
          ASSERT_LT(found, window.levels.size());
          EXPECT_EQ(s.state.reference.quantum_numbers, window.levels[found].state.reference.quantum_numbers);
          test_support::expect_exact(*s.gap, *window.levels[found++].gap, "window is a subset, not new equations");
        }
      EXPECT_EQ(found, window.levels.size());
    }
}

TYPED_TEST(BiquadraticScattering, NativeGapAndOrderingSurviveExtensiveOffset)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon();
  // Native four-site unbound singlet, distinct from the bound-pair root.
  auto const small = ferro::real_excitations_window<R>(4, 0, 2);
  ASSERT_TRUE(small.converged());
  R const exact = (R{9} + std::sqrt(R{17})) / R{2};
  EXPECT_REAL_NEAR(*small.levels[0].gap, exact, R{256} * eps);
  // Even the long-double total energy loses the spacing here, while the
  // native gap resolves it. Ranking by total energy would choose lowest I.
  std::size_t const n = 100000000;
  auto const band = ferro::real_excitations_window<R>(n, n - 2, 4, {.count = 4});
  ASSERT_TRUE(band.converged());
  ASSERT_EQ(band.levels.size(), 4);
  for (std::size_t i = 0; i < 4; ++i)
  {
    auto const& s = band.levels[i];
    EXPECT_EQ(s.state.reference.quantum_numbers.front(), uni20::half_int(std::int64_t(n - 1 - i)));
    EXPECT_REAL_NEAR(*s.gap, ferro::one_defect_level<R>(n, n - 1 - i).gap, R{32} * eps);
    test_support::expect_exact(*s.gap, s.state.tl_energy, "direct gap matches TL energy");
    if (i) EXPECT_GT(*s.gap, *band.levels[i - 1].gap);
  }
  if constexpr (uni20::numeric_limits<R>::digits <= 64)
    EXPECT_EQ(band.levels.front().state.energy, band.levels.back().state.energy);
  for (std::size_t length : {127, 128, 100000})
    for (std::size_t m : {2, 3})
    {
      auto const scan = ferro::real_excitations_window<R>(length, length - 2 * m, m + 2, {.count = 3});
      ASSERT_TRUE(scan.converged()) << length << " " << m;
      ASSERT_EQ(scan.levels.size(), 3);
      EXPECT_GT(*scan.levels[0].gap, R(m));
      EXPECT_LT(*scan.levels[0].gap - R(m), R{256} / (R(length) * R(length)));
      for (auto const& s : scan.levels)
        test_support::expect_exact(*s.gap, s.state.tl_energy, "gap does not subtract extensive energies");
    }
}

TYPED_TEST(BiquadraticScattering, OriginalEquationsOddAndEven)
{
  using R = TypeParam;
  using C = std::complex<R>;
  R const eps = uni20::numeric_limits<R>::epsilon(), eta = std::acosh(R{3} / R{2});
  for (std::size_t n : {5, 6, 9, 128})
  {
    auto const scan = ferro::real_excitations_window<R>(n, n - 4, 3, {.count = 3});
    ASSERT_TRUE(scan.converged());
    for (auto const& level : scan.levels)
    {
      auto const& s = level.state.reference;
      for (std::size_t j = 0; j < s.rapidities.size(); ++j)
      {
        C const u{0, s.rapidities[j] / R{2}};
        C lhs{1}, rhs{1};
        C const drive = std::sinh(u + eta / R{2}) / std::sinh(u - eta / R{2});
        for (std::size_t k = 0; k < 2 * n; ++k)
          lhs *= drive;
        for (std::size_t k = 0; k < s.rapidities.size(); ++k)
          if (k != j)
          {
            C const v{0, s.rapidities[k] / R{2}};
            rhs *= std::sinh(u - v + eta) / std::sinh(u - v - eta) * std::sinh(u + v + eta) / std::sinh(u + v - eta);
          }
        EXPECT_REAL_NEAR(std::abs(lhs / rhs - C{1}), R{0}, R{256} * R(n) * eps);
      }
    }
  }
}

TYPED_TEST(BiquadraticScattering, BudgetsParityAndDomains)
{
  using R = TypeParam;
  for (std::size_t width : {0, 1, 7})
    EXPECT_THROW((void)ferro::real_excitations_window<R>(8, 4, width), std::invalid_argument);
  EXPECT_THROW((void)ferro::real_excitations_window<R>(8, 8, 2), std::invalid_argument);
  EXPECT_THROW((void)ferro::real_excitations_window<R>(7, 2, 3), std::invalid_argument);
  EXPECT_THROW((void)ferro::real_excitations_window<R>(128, 124, 8, {.max_candidates = 27}), std::length_error);
  EXPECT_THROW((void)ferro::real_excitations_window<R>(128, 124, 8, {.count = 0}), std::invalid_argument);
  EXPECT_THROW((void)ferro::real_excitations_window<R>(1000000000000000000ULL, 0, 500000000000000000ULL),
               std::length_error);
  auto const empty = ferro::real_excitations_window<R>(129, 125, 8, {}, {.max_iterations = 0});
  EXPECT_FALSE(empty.converged());
  EXPECT_EQ(empty.candidate_count, 28);
  EXPECT_EQ(empty.converged_count, 0);
  EXPECT_TRUE(empty.levels.empty());
  EXPECT_TRUE(empty.first_unconverged);
  EXPECT_TRUE(empty.ground_state.reference.converged);
  EXPECT_EQ(empty.ground_state.energy, R{128});
  // Ground-sector and Q-system restrictions are deliberately unchanged.
  EXPECT_THROW((void)bethe::biquadratic::ground_state<R>(7), std::invalid_argument);
  EXPECT_THROW((void)ferro::qsystem::spectrum<R>(7, 3), std::invalid_argument);
  if constexpr (uni20::numeric_limits<R>::digits < 60)
  {
    std::size_t const n = std::size_t{1} << 60;
    EXPECT_THROW((void)ferro::solve_real<R>(n, std::vector<uni20::half_int>{uni20::half_int(std::int64_t(n - 1))}),
                 std::overflow_error);
  }
}

TEST(BiquadraticScatteringED, TwoDefectScatteringAndPairsExhaustSmallModules)
{
  for (unsigned n = 4; n <= 10; ++n)
  {
    auto expected = bethe::test::quantum_group_module_ed(n, n - 4, 1.5);
    auto const scattering = ferro::real_excitations<double>(n, n - 4, {.count = 1000});
    ASSERT_TRUE(scattering.converged());
    std::vector<double> obtained;
    for (auto const& s : scattering.levels)
      obtained.push_back(s.state.reference.energy);
    for (std::size_t mode = 1; mode <= n - 3; ++mode)
      obtained.push_back(ferro::bound_pair<double>(n, mode).reference.energy);
    std::sort(obtained.begin(), obtained.end());
    ASSERT_EQ(obtained.size(), expected.size());
    for (std::size_t i = 0; i < obtained.size(); ++i)
      EXPECT_NEAR(obtained[i], expected[i], 1e-11);
  }
}

TEST(BiquadraticScatteringED, OddSelectedRealRootsInEveryModule)
{
  for (unsigned n : {3, 5, 7, 9})
    for (unsigned ell = 1; ell <= n; ell += 2)
    {
      auto expected = bethe::test::quantum_group_module_ed(n, ell, 1.5);
      auto const scan = ferro::real_excitations<double>(n, ell, {.count = 1000});
      ASSERT_TRUE(scan.converged()) << n << " " << ell;
      for (auto const& level : scan.levels)
      {
        auto match = std::find_if(expected.begin(), expected.end(),
                                  [&](double e) { return std::abs(e - level.state.reference.energy) < 1e-11; });
        ASSERT_NE(match, expected.end()) << n << " " << ell;
        expected.erase(match);
      }
    }
}
} // namespace
