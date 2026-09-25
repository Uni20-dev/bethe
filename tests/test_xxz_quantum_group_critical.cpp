// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxz_quantum_group_critical.hpp>

namespace
{
namespace model = bethe::xxz::quantum_group::critical;
template <typename Real> class CriticalQuantumGroup : public ::testing::Test {};
TYPED_TEST_SUITE(CriticalQuantumGroup, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(CriticalQuantumGroup, NativeSmallChains)
{
  using Real = TypeParam;
  auto const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real d : {Real{1} / Real{100}, Real{1} / Real{4}, Real{3} / Real{5}, Real{9} / Real{10}, Real{99} / Real{100}})
  {
    auto const two = model::sea_state<Real>(2, d);
    ASSERT_TRUE(two.converged);
    EXPECT_REAL_NEAR(*two.energy, -Real{3} * d / Real{4}, Real{128} * eps);
    auto const four = model::sea_state<Real>(4, d);
    ASSERT_TRUE(four.converged);
    Real const exact = -Real{3} * d / Real{4} - std::sqrt(d * d + Real{2}) / Real{2};
    EXPECT_REAL_NEAR(*four.energy, exact, Real{512} * eps);
    auto const three = model::sea_state<Real>(3, d, 1);
    ASSERT_TRUE(three.converged);
    EXPECT_REAL_NEAR(*three.energy, -(Real{1} + d) / Real{2}, Real{128} * eps);
    test_support::expect_exact(uni20::parse_real<Real>(uni20::format_real(*four.energy)), *four.energy,
                               "critical QG I/O");
  }
}

TYPED_TEST(CriticalQuantumGroup, OriginalComplexBetheEquations)
{
  using Real = TypeParam;
  using C = uni20::complex<Real>;
  auto const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {4, 7, 16, 32})
    for (Real d : {Real{1} / Real{4}, Real{3} / Real{5}, Real{9} / Real{10}})
      for (bool excited : {false, true})
      {
        auto labels = bethe::xxz::quantum_group::detail::consecutive(n / 2 - (excited ? 1 : 0));
        if (excited) labels.back() += uni20::half_int{1};
        auto const s = model::solve_real<Real>(n, d, labels);
        ASSERT_TRUE(s.converged) << n << ' ' << uni20::format_real(d);
        Real const gamma = std::acos(d);
        auto ratio = [&](Real x, Real width) { return std::sinh(C{x, width}) / std::sinh(C{x, -width}); };
        Real energy = Real(n - 1) * d / Real{4};
        for (std::size_t i = 0; i < s.rapidities.size(); ++i)
        {
          auto const lambda = s.rapidities[i];
          ASSERT_TRUE(uni20::isfinite(lambda));
          EXPECT_GT(lambda, Real{0});
          C lhs{1, 0}, rhs{1, 0};
          auto const drive = ratio(lambda, gamma / Real{2});
          for (std::size_t k = 0; k < 2 * n; ++k)
            lhs *= drive;
          for (std::size_t j = 0; j < s.rapidities.size(); ++j)
            if (i != j) rhs *= ratio(lambda - s.rapidities[j], gamma) * ratio(lambda + s.rapidities[j], gamma);
          EXPECT_REAL_NEAR(lhs.real(), rhs.real(), Real{1024} * Real(n) * eps);
          EXPECT_REAL_NEAR(lhs.imag(), rhs.imag(), Real{1024} * Real(n) * eps);
          energy -= (Real{1} - d * d) / (std::cosh(Real{2} * lambda) - d);
        }
        EXPECT_REAL_NEAR(*s.energy, energy, Real{512} * Real(n) * eps);
      }
}

TYPED_TEST(CriticalQuantumGroup, InputsBudgetsAndExplicitLabels)
{
  using Real = TypeParam;
  Real const d = Real{1} / Real{4};
  auto const failed = model::sea_state<Real>(8, d, 0, {.max_iterations = 0});
  EXPECT_FALSE(failed.converged);
  EXPECT_EQ(failed.status, model::Status::iteration_limit);
  EXPECT_FALSE(failed.energy);
  EXPECT_FALSE(failed.energy_shift);
  auto const vacuum = model::sea_state<Real>(8, d, 8, {.max_iterations = 0});
  ASSERT_TRUE(vacuum.converged);
  EXPECT_EQ(*vacuum.energy, Real{7} * d / Real{4});
  std::vector<uni20::half_int> labels{uni20::half_int{2}};
  auto const one = model::solve_real<Real>(8, d, labels);
  ASSERT_TRUE(one.converged);
  Real const pi = Real{4} * std::atan(Real{1});
  EXPECT_REAL_NEAR(*one.energy, Real{7} * d / Real{4} - d - std::cos(pi / Real{4}),
                   Real{128} * uni20::numeric_limits<Real>::epsilon());
  labels[0] = uni20::half_int{7};
  EXPECT_THROW(model::solve_real<Real>(8, d, labels), std::invalid_argument);
  for (Real invalid : {Real{0}, Real{1}, -d, uni20::numeric_limits<Real>::infinity()})
    EXPECT_THROW(model::sea_state<Real>(4, invalid), std::invalid_argument);
  EXPECT_THROW(model::sea_state<Real>(5, d), std::invalid_argument);
  EXPECT_THROW(model::sea_state<Real>(4, d, 1), std::invalid_argument);
  EXPECT_THROW(model::sea_state<Real>(4, d, 0, {.residual_tolerance = Real{0}}), std::invalid_argument);
}

TYPED_TEST(CriticalQuantumGroup, IndependentNonHermitianSpinBasis)
{
  using Real = TypeParam;
  // Direct complex spin-basis eigensolver, not a Hermitian diagonalization
  // or the logarithmic equations. See scripts/reference_xxz_nonhermitian.py.
  struct Reference
  {
      std::size_t n;
      int numerator;
      double energy;
  };
  for (auto const& r : {Reference{5, 25, -1.3016421401931328},
                        {5, 60, -1.582686814421488},
                        {5, 90, -1.8397334893284407},
                        {6, 25, -1.605519747797113},
                        {6, 60, -1.9906369970263755},
                        {6, 90, -2.362192244915388},
                        {8, 25, -2.3007582273337226},
                        {8, 60, -2.7645096966847658},
                        {8, 90, -3.2148081413167406}})
  {
    auto const s = model::sea_state<Real>(r.n, Real(r.numerator) / Real{100}, r.n % 2);
    ASSERT_TRUE(s.converged);
    EXPECT_REAL_NEAR(*s.energy, Real(r.energy), Real{1e-11});
  }
}
} // namespace
