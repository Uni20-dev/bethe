// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "exact_spectrum.hpp"
#include "test_support.hpp"
#include <bethe/xxz_open.hpp>
#include <bethe/xxz_open_massive.hpp>

namespace
{
namespace model = bethe::xxz::open::massive;
using uni20::half_int;
template <typename Real> class XXZOpenMassive : public ::testing::Test {};
TYPED_TEST_SUITE(XXZOpenMassive, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXZOpenMassive, AllSectorsAgainstED)
{
  using Real = TypeParam;
  for (unsigned n = 2; n <= 9; ++n)
    for (unsigned m = 0; m <= n / 2; ++m)
      for (Real d : {Real{101} / Real{100}, Real{2}, Real{3}, Real{10}, Real{100}})
      {
        SCOPED_TRACE(::testing::Message() << "N=" << n << " M=" << m << " Delta=" << uni20::format_scalar(d));
        auto const sz = uni20::from_twice(std::int64_t(n - 2 * m));
        auto const state = model::sector_ground_state(n, d, sz);
        ASSERT_TRUE(state.converged) << uni20::format_scalar(state.residual_norm);
        EXPECT_EQ(state.status, model::SolveStatus::converged);
        EXPECT_EQ(state.delta, d);
        EXPECT_EQ(state.root_delta, d);
        EXPECT_EQ(state.sz, sz);
        EXPECT_EQ(state.boundary_root.has_value(), 2 * m == n);
        EXPECT_EQ(state.rapidities.size() + state.boundary_root.has_value(), m);
        ASSERT_EQ(state.rapidities.size(), state.quantum_numbers.size());
        for (std::size_t i = 0; i < state.rapidities.size(); ++i)
        {
          EXPECT_GT(state.rapidities[i], Real{0});
          if (i) EXPECT_LT(state.rapidities[i - 1], state.rapidities[i]);
          EXPECT_EQ(state.quantum_numbers[i], half_int(std::int64_t(i + 1)));
        }
        if (state.boundary_root) EXPECT_EQ(state.boundary_root->quantum_number, half_int(m));
        auto const ed = test_support::exact_spectrum(n, m, 0, false, static_cast<double>(d));
        EXPECT_REAL_NEAR(static_cast<double>(state.energy), ed.front(), 5e-11);
        auto const reversed = model::sector_ground_state(n, d, -sz);
        EXPECT_EQ(reversed.energy, state.energy);
        EXPECT_EQ(reversed.rapidities, state.rapidities);
        EXPECT_EQ(reversed.spin_reversed, sz.twice() != 0);
      }
}

TYPED_TEST(XXZOpenMassive, BoundaryCrossingAndNativePrecision)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real d : {Real{1} + Real{128} * eps, Real{2}, Real{299} / Real{100}, Real{3}, Real{301} / Real{100}, Real{10},
                 Real{100}, Real{10000}})
  {
    SCOPED_TRACE(uni20::format_scalar(d));
    auto const two = model::ground_state(2, d);
    ASSERT_TRUE(two.converged);
    ASSERT_TRUE(two.boundary_root);
    EXPECT_TRUE(two.rapidities.empty());
    Real const exact_y = (Real{3} - d) / (Real{1} + d);
    EXPECT_REAL_NEAR(two.boundary_root->inverse_square, exact_y, Real{2048} * eps);
    EXPECT_REAL_NEAR(two.energy, -Real{1} / Real{2} - d / Real{4}, Real{2048} * d * eps);
    auto const three = model::ground_state(3, d);
    ASSERT_TRUE(three.converged);
    EXPECT_FALSE(three.boundary_root);
    Real const exact = -(d + std::sqrt(d * d + Real{8})) / Real{4};
    EXPECT_REAL_NEAR(three.energy, exact, Real{512} * d * eps);
  }
  // Infinity must not be a spurious solution at every anisotropy: removing
  // the vanishing 1/z_B factor leaves a finite, nonzero equation off Delta=3.
  model::detail::GroundSystem<Real> const system(2, 1);
  std::vector<Real> infinity{Real{0}};
  EXPECT_GT(system.evaluate(infinity, Real{2}).norm, Real{1} / Real{10});
  EXPECT_EQ(system.evaluate(infinity, Real{3}).norm, Real{0});
  EXPECT_GT(system.evaluate(infinity, Real{4}).norm, Real{1} / Real{10});
  // The regularized equation must retain a linear signal on either side of
  // infinity, even much closer than a finite-difference Jacobian step.
  Real const tiny = std::sqrt(eps);
  std::vector<Real> left{-tiny}, right{tiny};
  EXPECT_LT(system.evaluate(left, Real{3}).residual[0], -tiny / Real{8});
  EXPECT_GT(system.evaluate(right, Real{3}).residual[0], tiny / Real{8});
}

TYPED_TEST(XXZOpenMassive, IndependentComplexReflectionEquations)
{
  using Real = TypeParam;
  using Complex = uni20::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Complex const one{Real{1}, Real{0}}, ii{Real{0}, Real{1}};
  auto ratio = [&](Complex z) { return (one + ii * z) / (one - ii * z); };
  for (std::size_t n : {4, 5, 6})
    for (Real d : {Real{6} / Real{5}, Real{2}, Real{4}})
    {
      auto const state = model::ground_state(n, d);
      ASSERT_TRUE(state.converged);
      Real const r = (d - Real{1}) / (d + Real{1}), p = Real{1} + d, q = Real{1} - d;
      std::vector<Complex> roots;
      for (Real z : state.rapidities)
        roots.emplace_back(z, Real{0});
      if (state.boundary_root)
      {
        Real const y = state.boundary_root->inverse_square;
        ASSERT_NE(y, Real{0});
        roots.push_back(y > Real{0} ? Complex{Real{1} / std::sqrt(y), Real{0}}
                                    : Complex{Real{0}, Real{1} / std::sqrt(-y)});
      }
      Complex energy{Real(n - 1) * d / Real{4}, Real{0}};
      for (std::size_t i = 0; i < roots.size(); ++i)
      {
        auto const z = roots[i];
        auto lhs = one, rhs = one;
        for (std::size_t j = 0; j < 2 * n; ++j)
          lhs *= ratio(z);
        lhs /= ratio(r * z) * ratio(r * z);
        for (std::size_t j = 0; j < roots.size(); ++j)
          if (i != j)
          {
            rhs *= ratio(d * (z - roots[j]) / (p - q * z * roots[j]));
            rhs *= ratio(d * (z + roots[j]) / (p + q * z * roots[j]));
          }
        Real const scale = Real{1} + std::abs(rhs.real()) + std::abs(rhs.imag());
        EXPECT_REAL_NEAR(lhs.real() / scale, rhs.real() / scale, Real{16384} * eps);
        EXPECT_REAL_NEAR(lhs.imag() / scale, rhs.imag() / scale, Real{16384} * eps);
        energy -= (p - q * z * z) / (one + z * z);
      }
      EXPECT_REAL_NEAR(energy.real(), state.energy, Real{16384} * eps);
      EXPECT_REAL_NEAR(energy.imag(), Real{0}, Real{16384} * eps);
    }
}

TYPED_TEST(XXZOpenMassive, ContinuityToXXX)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const d = Real{1} + Real{128} * eps;
  for (std::size_t n : {7, 8, 16})
    for (std::size_t m = 1; m <= n / 2; ++m)
    {
      auto const sz = uni20::from_twice(std::int64_t(n - 2 * m));
      auto const state = model::sector_ground_state(n, d, sz);
      auto const xxx = bethe::xxz::open::sector_ground_state(n, Real{1}, sz);
      ASSERT_TRUE(state.converged);
      ASSERT_TRUE(xxx.converged);
      EXPECT_EQ(state.delta, d); // Never snap the requested anisotropy to one.
      EXPECT_REAL_NEAR(state.energy, xxx.energy, Real{512} * Real(n) * eps);
      for (std::size_t i = 0; i < state.rapidities.size(); ++i)
        EXPECT_REAL_NEAR(state.rapidities[i], xxx.rapidities[i], Real{4096} * Real(n) * eps);
      if (state.boundary_root)
      {
        Real const z = Real{1} / std::sqrt(state.boundary_root->inverse_square);
        EXPECT_REAL_NEAR(z, xxx.rapidities.back(), Real{4096} * Real(n) * eps);
      }
    }
}

TYPED_TEST(XXZOpenMassive, LargeChainsAndExponentialBoundaryDeviation)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {16, 33, 64})
    for (Real d : {Real{2}, Real{10}, Real{100}})
    {
      SCOPED_TRACE(::testing::Message() << "N=" << n << " Delta=" << uni20::format_scalar(d));
      auto const state = model::ground_state(n, d);
      ASSERT_TRUE(state.converged) << uni20::format_scalar(state.residual_norm);
      EXPECT_LE(state.residual_norm, bethe::SolverOptions<Real>{}.residual_tolerance);
      // Neel product-state variational upper bound; lower bound follows by
      // subtracting at most 1/2 per XY bond from the minimal Ising energy.
      Real const ising = -Real(n - 1) * d / Real{4};
      EXPECT_LE(state.energy, ising + Real{64} * Real(n) * d * eps);
      EXPECT_GE(state.energy, ising - Real(n - 1) / Real{2});
      if (n == 64 && d == Real{100})
      {
        ASSERT_TRUE(state.boundary_root);
        EXPECT_GT(state.boundary_root->log_distance, -std::log(eps));
        EXPECT_TRUE(uni20::isfinite(state.boundary_root->log_distance));
      }
    }
}

TYPED_TEST(XXZOpenMassive, BudgetAndRequestedHamiltonianDiagnostics)
{
  using Real = TypeParam;
  for (std::size_t n : {6, 7})
    for (std::size_t budget : {0, 1, 7})
    {
      auto const state = model::ground_state(n, Real{10}, {.max_iterations = budget});
      EXPECT_FALSE(state.converged);
      EXPECT_EQ(state.status, model::SolveStatus::iteration_limit);
      EXPECT_EQ(state.iterations, budget);
      EXPECT_EQ(state.delta, Real{10});
      EXPECT_LT(state.root_delta, state.delta);
      std::vector<Real> x;
      Real energy = (Real(n - 1) / Real{4} - Real(n / 2)) * state.delta;
      for (Real z : state.rapidities)
      {
        x.push_back(std::atan(z));
        energy += (z * z - Real{1}) / (z * z + Real{1});
      }
      if (state.boundary_root)
      {
        x.push_back(state.boundary_root->log_distance);
        Real const y = state.boundary_root->inverse_square;
        energy += (Real{1} - y) / (Real{1} + y);
      }
      Real const allowance = Real{256} * Real(n) * uni20::numeric_limits<Real>::epsilon();
      EXPECT_REAL_NEAR(energy, state.energy, allowance);
      model::detail::GroundSystem<Real> const system(n, n / 2);
      EXPECT_REAL_NEAR(system.evaluate(x, state.delta).norm, state.residual_norm, allowance);
    }
}

TYPED_TEST(XXZOpenMassive, StrongAnisotropyRootLimit)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  for (std::size_t n : {7, 8, 16})
    for (Real d : {Real{1000}, Real{1000000}, Real{1000000000}})
    {
      SCOPED_TRACE(::testing::Message() << "N=" << n << " Delta=" << uni20::format_scalar(d));
      auto const state = model::ground_state(n, d);
      ASSERT_TRUE(state.converged) << uni20::format_scalar(state.residual_norm);
      for (std::size_t i = 0; i < state.rapidities.size(); ++i)
      {
        Real const limiting_angle = pi * Real(i + 1) / Real(2 * (n - n / 2));
        EXPECT_REAL_NEAR(std::atan(state.rapidities[i]), limiting_angle, Real{4} / d + Real{128} * eps);
      }
      Real const limiting_energy = -Real(n - 1) / Real{4};
      EXPECT_REAL_NEAR(state.energy / d, limiting_energy, Real(n) / d + Real{128} * eps);
    }
}

TYPED_TEST(XXZOpenMassive, Validation)
{
  using Real = TypeParam;
  for (Real d :
       {Real{-1}, Real{0}, Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW((void)model::ground_state(4, d), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state(1, Real{2}), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state(std::size_t{1} << 40, Real{2}), std::length_error);
  EXPECT_THROW((void)model::sector_ground_state(4, Real{2}, uni20::from_twice(std::int64_t{1})), std::invalid_argument);
  EXPECT_THROW((void)model::sector_ground_state(4, Real{2}, half_int{3}), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state(4, Real{2}, {.residual_tolerance = Real{0}}), std::invalid_argument);
  auto const polarized = model::sector_ground_state(4, Real{2}, half_int{-2}, {.max_iterations = 0});
  EXPECT_TRUE(polarized.converged);
  EXPECT_EQ(polarized.iterations, 0);
  EXPECT_EQ(polarized.energy, Real{3} / Real{2});
}
} // namespace
