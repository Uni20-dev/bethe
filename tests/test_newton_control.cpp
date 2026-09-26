// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/detail/newton_backtracking.hpp>

namespace
{
template <typename Real> class NewtonControl : public ::testing::Test {};
TYPED_TEST_SUITE(NewtonControl, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(NewtonControl, DomainRejectionHalvesFromOriginalIterate)
{
  using Real = TypeParam;
  std::vector<Real> x{Real{1}};
  unsigned trials = 0, evaluations = 0;
  EXPECT_TRUE(bethe::detail::backtrack_newton(
      x, [](std::size_t) { return Real{-2}; },
      [&](auto const& trial, Real damping) {
        ++trials;
        if (trial[0] <= Real{0}) return false;
        ++evaluations;
        return bethe::detail::newton_decreases(trial[0], Real{1}, damping, Real{0});
      }));
  EXPECT_EQ(trials, 3u);
  EXPECT_EQ(evaluations, 1u);
  test_support::expect_exact(x[0], Real{0.5}, "accepted quarter step");
}

TYPED_TEST(NewtonControl, ExhaustionKeepsIterateAndNativeHalvingBudget)
{
  using Real = TypeParam;
  std::vector<Real> x{Real{1}};
  unsigned trials = 0;
  Real last{};
  EXPECT_FALSE(bethe::detail::backtrack_newton(
      x, [](std::size_t) { return Real{1}; },
      [&](auto const&, Real damping) {
        ++trials;
        last = damping;
        return false;
      }));
  EXPECT_EQ(trials, unsigned(uni20::numeric_limits<Real>::digits + 1));
  Real expected{1};
  for (int i = 0; i < uni20::numeric_limits<Real>::digits; ++i)
    expected /= Real{2};
  test_support::expect_exact(last, expected, "last damping");
  test_support::expect_exact(x[0], Real{1}, "rejected iterate");
}

TYPED_TEST(NewtonControl, NormalizationRunsOnFreshTrials)
{
  using Real = TypeParam;
  std::vector<Real> x{Real{1}};
  unsigned trials = 0;
  EXPECT_TRUE(bethe::detail::backtrack_newton(
      x, [](std::size_t) { return Real{2}; },
      [&](auto const& trial, Real damping) {
        ++trials;
        test_support::expect_exact(trial[0], Real{11} + Real{2} * damping, "normalized trial");
        return trials == 2;
      },
      [](auto& trial) { trial[0] += Real{10}; }));
  test_support::expect_exact(x[0], Real{12}, "normalized accepted iterate");
}

TYPED_TEST(NewtonControl, ToleranceShortcutAndNativeResolution)
{
  using Real = TypeParam;
  EXPECT_TRUE(bethe::detail::newton_decreases(Real{1}, Real{1}, Real{1}, Real{1}));
  EXPECT_FALSE(bethe::detail::newton_decreases(Real{1}, Real{1}, Real{1}, Real{0}));
  EXPECT_FALSE(bethe::detail::newton_decreases(uni20::numeric_limits<Real>::quiet_NaN(), Real{1}, Real{1}, Real{1}));
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  std::vector<Real> x{Real{1}};
  EXPECT_TRUE(bethe::detail::backtrack_newton(
      x, [&](std::size_t) { return eps; }, [](auto const& trial, Real) { return trial[0] > Real{1}; }));
  test_support::expect_exact(x[0], Real{1} + eps, "native correction");
}
} // namespace
