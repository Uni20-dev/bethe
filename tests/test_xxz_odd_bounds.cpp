// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/xxz_odd_bounds.hpp>

#include "exact_spectrum.hpp"
#include "test_support.hpp"
#include <complex>
#include <numeric>

namespace
{
namespace engine = bethe::xxz::detail;
template <typename Real> class XXZOddBounds : public ::testing::Test {};
TYPED_TEST_SUITE(XXZOddBounds, test_support::RealTypes, test_support::PrecisionNames);

// Independent factorial Slater determinant and direct spin-basis action.
// No Bethe roots, continuation coefficients, or wavefunction library helpers.
template <typename Real> std::pair<Real, Real> trial_expectations(unsigned n, unsigned m, Real delta)
{
  using C = std::complex<Real>;
  Real const pi = Real{4} * std::atan(Real{1});
  std::vector<C> sea(1U << n), helix(1U << n);
  for (unsigned bits = 0; bits < (1U << n); ++bits)
    if (std::popcount(bits) == int(m))
    {
      std::vector<unsigned> occupied, permutation(m);
      std::iota(permutation.begin(), permutation.end(), 0U);
      unsigned sum = 0;
      for (unsigned j = 0; j < n; ++j)
        if (bits & (1U << j))
        {
          occupied.push_back(j);
          sum += j;
        }
      Real const pitch = pi + pi / Real(n), helix_angle = pitch * Real(sum);
      helix[bits] = {std::cos(helix_angle), std::sin(helix_angle)};
      do
      {
        Real angle{};
        unsigned inversions = 0;
        for (unsigned a = 0; a < m; ++a)
        {
          angle += (pi + Real{2} * pi * (Real(permutation[a]) - Real(m) / Real{2}) / Real(n)) * Real(occupied[a]);
          for (unsigned b = a + 1; b < m; ++b)
            inversions += permutation[a] > permutation[b];
        }
        sea[bits] += Real(inversions % 2 ? -1 : 1) * C{std::cos(angle), std::sin(angle)};
      }
      while (std::next_permutation(permutation.begin(), permutation.end()));
    }
  auto rayleigh = [&](std::vector<C> const& wave) {
    Real norm{}, expectation{};
    for (unsigned bits = 0; bits < (1U << n); ++bits)
      if (std::popcount(bits) == int(m))
      {
        norm += std::norm(wave[bits]);
        C action{};
        for (unsigned j = 0; j < n; ++j)
        {
          unsigned const k = (j + 1) % n;
          bool const unlike = ((bits >> j) & 1U) != ((bits >> k) & 1U);
          action += delta * (unlike ? -Real{1} : Real{1}) * wave[bits] / Real{4};
          if (unlike) action += wave[bits ^ (1U << j) ^ (1U << k)] / Real{2};
        }
        expectation += (std::conj(wave[bits]) * action).real();
      }
    return expectation / norm;
  };
  return {rayleigh(sea), rayleigh(helix)};
}

TYPED_TEST(XXZOddBounds, ExplicitTrialVectorsAndSectorMinima)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (unsigned n : {3, 5, 7, 9})
    for (unsigned m = 0; m <= n / 2; ++m)
      for (Real d : {Real{0}, -Real{1} / Real{5}, -Real{7} / Real{10}, -Real{99} / Real{100}})
      {
        SCOPED_TRACE(::testing::Message() << "N=" << n << " M=" << m << " Delta=" << uni20::format_scalar(d));
        auto const result = engine::OddSectorVariationalBounds<Real>(n, m).evaluate(d);
        auto const [sea, helix] = trial_expectations(n, m, d);
        EXPECT_REAL_NEAR(result.free_sea, sea, Real{2048} * Real(n) * eps);
        EXPECT_REAL_NEAR(result.helix, helix, Real{2048} * Real(n) * eps);
        auto const ed = test_support::exact_spectrum(n, m, 0, true, static_cast<double>(d));
        EXPECT_LE(ed.front(), static_cast<double>(result.upper_bound()) + 3e-11);
        if (d == Real{0}) EXPECT_REAL_NEAR(static_cast<double>(result.free_sea), ed.front(), 3e-11);
      }
}

TYPED_TEST(XXZOddBounds, ExcludesKnownWrongBranchIndependentlyOfHistory)
{
  using Real = TypeParam;
  // N=17, M=8, Delta=-0.97: the former wrong branch passed residual,
  // momentum, conditioning and Lipschitz guards. No preceding energies
  // or exact diagonalization are needed to exclude it variationally.
  auto const trial = engine::OddSectorVariationalBounds<Real>(17, 8).evaluate(-Real{97} / Real{100});
  Real const wrong = -Real{402462082373559LL} / Real{100000000000000LL};
  EXPECT_GT(wrong - trial.upper_bound(), Real{15} / Real{100});
  EXPECT_LT(trial.helix, trial.free_sea);
  // Away from the helix coupling the other trial can be stronger.
  auto const free = engine::OddSectorVariationalBounds<Real>(17, 8).evaluate(Real{0});
  EXPECT_LT(free.free_sea, free.helix);
}

TYPED_TEST(XXZOddBounds, AnalyticLimitsAndInvalidInput)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  for (std::size_t n : {3, 5, 17, 101})
  {
    Real const d = -std::cos(pi / Real(n));
    for (std::size_t m = 0; m <= n / 2; ++m)
    {
      auto const collision = engine::OddSectorVariationalBounds<Real>(n, m).evaluate(d);
      EXPECT_REAL_NEAR(collision.helix, Real(n) * d / Real{4}, Real{16} * Real(n) * eps);
    }
    Real const delta = -Real{1} + Real{128} * eps;
    auto const vacuum = engine::OddSectorVariationalBounds<Real>(n, 0).evaluate(delta);
    EXPECT_EQ(vacuum.free_sea, Real(n) * delta / Real{4});
    EXPECT_EQ(vacuum.helix, vacuum.free_sea);
    auto const one = engine::OddSectorVariationalBounds<Real>(n, 1).evaluate(delta);
    Real const expected = d + (Real(n) / Real{4} - Real{1}) * delta;
    EXPECT_REAL_NEAR(one.free_sea, expected, Real{16} * Real(n) * eps);
    EXPECT_REAL_NEAR(one.helix, expected, Real{16} * Real(n) * eps);
  }
  EXPECT_THROW((engine::OddSectorVariationalBounds<Real>{1, 0}), std::invalid_argument);
  EXPECT_THROW((engine::OddSectorVariationalBounds<Real>{6, 2}), std::invalid_argument);
  EXPECT_THROW((engine::OddSectorVariationalBounds<Real>{7, 4}), std::invalid_argument);
  for (Real d : {-Real{1}, Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW((engine::OddSectorVariationalBounds<Real>{7, 3}.evaluate(d)), std::invalid_argument);
}
} // namespace
