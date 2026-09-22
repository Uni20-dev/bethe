// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "exact_spectrum.hpp"
#include "test_support.hpp"
#include <bethe/xxz_negative.hpp>

namespace
{
namespace engine = bethe::xxz::detail;
template <typename Real> class XXZNegative : public ::testing::Test {};
TYPED_TEST_SUITE(XXZNegative, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXZNegative, BipartiteSectorsAgainstED)
{
  using Real = TypeParam;
  for (bool open : {false, true})
    for (unsigned n = 2; n <= 9; ++n)
    {
      if (!open && n % 2) continue;
      for (unsigned m = 0; m <= n / 2; ++m)
        for (Real d : {-Real{1} / Real{100}, -Real{1} / Real{2}, -Real{9} / Real{10}, -Real{99} / Real{100}})
        {
          SCOPED_TRACE(::testing::Message()
                       << "open=" << open << " N=" << n << " M=" << m << " Delta=" << uni20::format_scalar(d));
          auto const sz = uni20::from_twice(std::int64_t(n - 2 * m));
          auto const state = engine::negative_ground_roots(n, d, sz, open);
          ASSERT_TRUE(state.converged) << uni20::format_scalar(state.residual_norm);
          EXPECT_EQ(state.rapidities.size(), m);
          auto const ed = test_support::exact_spectrum(n, m, 0, !open, static_cast<double>(d));
          EXPECT_REAL_NEAR(static_cast<double>(state.energy), ed.front(), 3e-11);
          if (!open)
          {
            auto const p_index = engine::momentum_index(n, state.quantum_numbers);
            double const p = 8 * std::atan(1.0) * double(p_index) / double(n);
            auto const resolved = test_support::exact_spectrum(n, m, .371, true, static_cast<double>(d));
            double const target = static_cast<double>(state.energy) + .371 * std::cos(p);
            EXPECT_TRUE(
                std::any_of(resolved.begin(), resolved.end(), [&](double e) { return std::abs(e - target) < 3e-11; }));
          }
          auto const reversed = engine::negative_ground_roots(n, d, -sz, open);
          EXPECT_EQ(state.energy, reversed.energy);
          EXPECT_EQ(state.rapidities, reversed.rapidities);
        }
    }
}

TYPED_TEST(XXZNegative, OriginalEquationsAndJacobian)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (bool open : {false, true})
    for (Real d : {-Real{1} / Real{10}, -Real{1} / Real{2}, -Real{9} / Real{10}})
    {
      std::size_t const n = open ? 7 : 8, m = 3;
      auto const sz = uni20::from_twice(std::int64_t(n - 2 * m));
      auto const state = engine::negative_ground_roots(n, d, sz, open);
      ASSERT_TRUE(state.converged);
      for (std::size_t i = 0; i < m; ++i)
      {
        Real const z = state.rapidities[i];
        Real f = Real(open ? 2 * n : n) * std::atan(z) - pi * Real(state.quantum_numbers[i].twice()) / Real{2};
        if (open) f += Real{2} * std::atan((Real{1} - d) * z / (Real{1} + d));
        for (std::size_t j = 0; j < m; ++j)
          if (i != j)
          {
            Real const w = state.rapidities[j];
            f -= std::atan(d * (z - w) / (Real{1} + d - (Real{1} - d) * z * w));
            if (open) f -= std::atan(d * (z + w) / (Real{1} + d + (Real{1} - d) * z * w));
          }
        EXPECT_REAL_NEAR(f, Real{0}, Real{256} * Real(n) * eps);
      }
      engine::NegativeGroundSystem<Real> const system(n, m, d, open);
      uni20::DenseMatrix<Real> jac(m, m);
      auto const x = system.seed();
      (void)system.evaluate(x, &jac);
      Real const h = std::cbrt(eps);
      for (std::size_t j = 0; j < m; ++j)
      {
        auto plus = x, minus = x;
        plus[j] += h;
        minus[j] -= h;
        auto const fp = system.evaluate(plus), fm = system.evaluate(minus);
        for (std::size_t i = 0; i < m; ++i)
        {
          Real const finite = (fp.residual[i] - fm.residual[i]) / (Real{2} * h);
          EXPECT_REAL_NEAR((jac[i, j]), finite, Real{4096} * h * h);
          EXPECT_REAL_NEAR((jac[i, j]), (jac[j, i]), Real{32} * eps);
        }
      }
    }
}

TYPED_TEST(XXZNegative, NativePrecisionAndFerromagneticLimit)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real d : {-eps, -Real{1} / Real{2}, -Real{1} + Real{128} * eps})
  {
    auto const two = engine::negative_ground_roots(2, d, uni20::half_int{0}, true);
    ASSERT_TRUE(two.converged);
    EXPECT_REAL_NEAR(two.energy, -Real{1} / Real{2} - d / Real{4}, Real{128} * eps);
    Real const scaled = two.rapidities.front() / std::sqrt((Real{1} + d) / (Real{1} - d));
    EXPECT_REAL_NEAR(scaled, std::sqrt((Real{1} - d) / (Real{3} - d)), Real{256} * eps);
    auto const three = engine::negative_ground_roots(3, d, uni20::from_twice(std::int64_t{1}), true);
    ASSERT_TRUE(three.converged);
    EXPECT_REAL_NEAR(three.energy, -(d + std::sqrt(d * d + Real{8})) / Real{4}, Real{256} * eps);
  }
  for (bool open : {false, true})
    for (std::size_t n : {16, 64})
    {
      Real const d = -Real{1} + Real{128} * eps;
      auto const state = engine::negative_ground_roots(n, d, uni20::half_int{0}, open);
      ASSERT_TRUE(state.converged) << uni20::format_scalar(state.residual_norm);
      EXPECT_LE(state.residual_norm, bethe::SolverOptions<Real>{}.residual_tolerance);
      Real const limit = -Real(open ? n - 1 : n) / Real{4};
      EXPECT_REAL_NEAR(state.energy, limit, Real{128} * Real(n) * eps);
      auto const seed = engine::negative_ground_roots(n, d, uni20::half_int{0}, open, {.max_iterations = 0});
      EXPECT_FALSE(seed.converged);
      EXPECT_GT(seed.residual_norm, Real{1} / Real{100});
    }
}

TYPED_TEST(XXZNegative, FerromagneticEnergySlope)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), h = std::sqrt(eps);
  // At Delta=-1 a staggered rotation gives the isotropic ferromagnet.
  // Its symmetric fixed-Sz state has <Sz_i Sz_j> = ((N-2M)^2-N)/(4N(N-1)).
  // This first-order slope distinguishes the sector minimum from a spurious
  // polarized-limit answer even when both energies approach -bonds/4.
  for (bool open : {false, true})
    for (std::size_t n : {4, 8, 16})
      for (std::size_t m : {std::size_t{1}, n / 2})
      {
        SCOPED_TRACE(::testing::Message() << "open=" << open << " N=" << n << " M=" << m);
        auto const sz = uni20::from_twice(std::int64_t(n - 2 * m));
        auto const state = engine::negative_ground_roots(n, -Real{1} + h, sz, open);
        ASSERT_TRUE(state.converged);
        Real const bonds = Real(open ? n - 1 : n), magnetization = Real(n - 2 * m);
        Real const slope = bonds * (magnetization * magnetization - Real(n)) / (Real{4} * Real(n) * Real(n - 1));
        EXPECT_REAL_NEAR(state.energy, -bonds / Real{4} + h * slope, Real{32} * Real(n * n) * eps);
      }
}

TYPED_TEST(XXZNegative, BudgetsAndDomain)
{
  using Real = TypeParam;
  for (bool open : {false, true})
    for (std::size_t budget : {0, 1})
    {
      Real const d = -Real{1} / Real{2};
      auto const state = engine::negative_ground_roots(8, d, uni20::half_int{0}, open, {.max_iterations = budget});
      EXPECT_FALSE(state.converged);
      EXPECT_EQ(state.iterations, budget);
      EXPECT_EQ(state.status, engine::NegativeSolveStatus::iteration_limit);
      engine::NegativeGroundSystem<Real> const system(8, 4, d, open);
      EXPECT_EQ(state.residual_norm, system.evaluate(state.log_rapidities).norm);
      Real energy = Real(open ? 7 : 8) * d / Real{4};
      for (Real z : state.rapidities)
        energy -= (Real{1} + d - (Real{1} - d) * z * z) / (Real{1} + z * z);
      EXPECT_REAL_NEAR(energy, state.energy, Real{128} * uni20::numeric_limits<Real>::epsilon());
    }
  EXPECT_THROW((void)engine::negative_ground_roots(5, -Real{1} / Real{2}, uni20::from_twice(std::int64_t{1}), false),
               std::invalid_argument);
  for (Real bad : {-Real{2}, -Real{1}, Real{0}, Real{1}, uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW((void)engine::negative_ground_roots(4, bad, uni20::half_int{0}, true), std::invalid_argument);
  EXPECT_THROW((void)engine::negative_ground_roots(4, -Real{1} / Real{2}, uni20::half_int{0}, true,
                                                   {.residual_tolerance = Real{0}}),
               std::invalid_argument);
}

TEST(XXZNegativeAudit, OddRingGlobalGroundStateIsNotAlwaysSmallestSz)
{
  // Independent counterexample to extending the old ground_state wrapper by
  // changing only its Delta check. Here the polarized state is lower.
  auto const polarized = test_support::exact_spectrum(5, 0, 0, true, -.9);
  auto const smallest_sz = test_support::exact_spectrum(5, 2, 0, true, -.9);
  EXPECT_LT(polarized.front(), smallest_sz.front());
  EXPECT_REAL_NEAR(polarized.front(), -1.125, 1e-14);
}
} // namespace
