// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "exact_spectrum.hpp"
#include "test_support.hpp"
#include <bethe/xxz_odd_continuation.hpp>

namespace
{
namespace engine = bethe::xxz::detail;
template <typename Real> class XXZOddContinuation : public ::testing::Test {};
TYPED_TEST_SUITE(XXZOddContinuation, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXZOddContinuation, SmallSectorsAgainstED)
{
  using Real = TypeParam;
  for (unsigned n : {3, 5, 7, 9})
    for (unsigned m = 0; m <= n / 2; ++m)
      for (Real d : {-Real{1} / Real{5}, -Real{7} / Real{10}, -Real{9} / Real{10}, -Real{999} / Real{1000}})
      {
        SCOPED_TRACE(::testing::Message() << "N=" << n << " M=" << m << " Delta=" << uni20::format_scalar(d));
        auto const sz = uni20::from_twice(std::int64_t(n - 2 * m));
        auto const state = engine::continue_odd_polynomial(n, d, sz);
        ASSERT_TRUE(state.equations_converged)
            << "status=" << int(state.status) << " residual=" << uni20::format_scalar(state.residual_norm)
            << " root_delta=" << uni20::format_scalar(state.root_delta)
            << " condition=" << uni20::format_scalar(state.reciprocal_condition);
        EXPECT_EQ(state.root_delta, d);
        EXPECT_LE(state.residual_norm, bethe::SolverOptions<Real>{}.residual_tolerance);
        EXPECT_LE(state.momentum_error, bethe::SolverOptions<Real>{}.residual_tolerance);
        auto const upper = engine::OddSectorVariationalBounds<Real>(n, m).evaluate(d).upper_bound();
        EXPECT_EQ(state.variational_upper_bound, upper);
        EXPECT_LE(state.energy, upper + Real{1024} * Real(n) * uni20::numeric_limits<Real>::epsilon());
        auto const ed = test_support::exact_spectrum(n, m, 0, true, static_cast<double>(d));
        EXPECT_REAL_NEAR(static_cast<double>(state.energy), ed.front(), 3e-10);
        engine::PolynomialBetheSystem<Real> const system(n, m, state.center, state.coordinate_scale);
        auto const joint = test_support::exact_spectrum(n, m, .371, true, static_cast<double>(d));
        double const target = static_cast<double>(state.energy + Real{371} / Real{1000} *
                                                                     system.momentum_phase(state.coefficients).real());
        EXPECT_TRUE(std::any_of(joint.begin(), joint.end(), [&](double e) { return std::abs(e - target) < 3e-10; }));
      }
}

TYPED_TEST(XXZOddContinuation, CriticalCollisions)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {5, 7, 9})
    for (std::size_t m = 2; m <= n / 2; ++m)
    {
      SCOPED_TRACE(::testing::Message() << "N=" << n << " M=" << m);
      Real const d = -std::cos(pi / Real(n));
      auto const state = engine::continue_odd_polynomial(n, d, uni20::from_twice(std::int64_t(n - 2 * m)));
      ASSERT_TRUE(state.equations_converged)
          << "status=" << int(state.status) << " residual=" << uni20::format_scalar(state.residual_norm)
          << " root_delta=" << uni20::format_scalar(state.root_delta);
      EXPECT_REAL_NEAR(state.energy, Real(n) * d / Real{4}, Real{4096} * Real(n) * eps);
    }
}

TYPED_TEST(XXZOddContinuation, BudgetsAndRequestedDiagnostics)
{
  using Real = TypeParam;
  Real const d = -Real{9} / Real{10};
  for (std::size_t budget : {0, 1, 3})
  {
    auto const state =
        engine::continue_odd_polynomial(7, d, uni20::from_twice(std::int64_t{1}), {.max_iterations = budget});
    EXPECT_FALSE(state.equations_converged);
    EXPECT_EQ(state.iterations, budget);
    EXPECT_EQ(state.status, engine::PolynomialContinuationStatus::iteration_limit);
    EXPECT_NE(state.root_delta, d);
    EXPECT_EQ(state.variational_upper_bound, engine::OddSectorVariationalBounds<Real>(7, 3).evaluate(d).upper_bound());
    engine::PolynomialBetheSystem<Real> const system(7, 3, state.center, state.coordinate_scale);
    EXPECT_EQ(state.energy, system.energy(state.coefficients, d));
    EXPECT_EQ(state.residual_norm, system.evaluate(state.coefficients, d).norm);
  }
  for (auto sz : {uni20::from_twice(std::int64_t{7}), uni20::from_twice(std::int64_t{5})})
  {
    auto const state = engine::continue_odd_polynomial(7, d, sz, {.max_iterations = 0});
    EXPECT_TRUE(state.equations_converged);
    EXPECT_EQ(state.iterations, 0);
  }
}

TYPED_TEST(XXZOddContinuation, CoordinatesAndJacobian)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), h = std::cbrt(eps);
  engine::OddPolynomialCoordinates<Real> const coordinates(9, 4);
  for (Real d : {-Real{1} / Real{2}, -Real{99} / Real{100}})
  {
    auto c = coordinates.seed();
    auto const old = c;
    Real const s = coordinates.scale(d), s0 = coordinates.scale(Real{0});
    coordinates.rescale(c, s0, s);
    engine::PolynomialBetheSystem<Real> const original(9, 4, coordinates.center, s0),
        system(9, 4, coordinates.center, s);
    EXPECT_REAL_NEAR(system.energy(c, d), original.energy(old, d), Real{512} * eps);
    EXPECT_REAL_NEAR(std::abs(system.momentum_phase(c) - original.momentum_phase(old)), Real{0}, Real{512} * eps);
    coordinates.impose_momentum(c, coordinates.weights(s));
    EXPECT_LT(std::abs(system.momentum_phase(c) - coordinates.target_phase), Real{512} * eps);
    uni20::DenseMatrix<Real> jac(4, 4);
    (void)system.evaluate(c, d, &jac);
    for (std::size_t j = 0; j < 4; ++j)
    {
      auto plus = c, minus = c;
      plus[j] += h;
      minus[j] -= h;
      auto const fp = system.evaluate(plus, d), fm = system.evaluate(minus, d);
      for (std::size_t i = 0; i < 4; ++i)
        EXPECT_REAL_NEAR((jac[i, j]), (fp.residual[i] - fm.residual[i]) / (Real{2} * h),
                         Real{8192} * h * h * (Real{1} + std::abs(jac[i, j])));
    }
  }
}
TYPED_TEST(XXZOddContinuation, LargerChains)
{
  using Real = TypeParam;
  // Independent sparse spin-basis diagonalization, regenerated by
  // reference_xxz_odd_ed.py. This oracle is double precision; the analytic
  // five-site test below separately checks native-precision accuracy.
  double const expected[3][3] = {{-3.4642927862816908, -3.1557977107952437, -3.1498333742645181},
                                 {-4.5868880833802352, -4.1813111949074635, -4.1743772174219691},
                                 {-5.7001667598194876, -5.1974076892523398, -5.1892895735101501}};
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j)
    {
      std::size_t const n = 13 + 4 * i;
      Real const deltas[] = {-Real{1} / Real{2}, -Real{97} / Real{100}, -Real{999} / Real{1000}};
      Real const d = deltas[j];
      auto const state =
          engine::continue_odd_polynomial(n, d, uni20::from_twice(std::int64_t{1}), {.max_iterations = 500});
      ASSERT_TRUE(state.equations_converged)
          << "N=" << n << " delta=" << uni20::format_scalar(d) << " status=" << int(state.status)
          << " residual=" << uni20::format_scalar(state.residual_norm)
          << " root_delta=" << uni20::format_scalar(state.root_delta) << " iterations=" << state.iterations
          << " rcond=" << uni20::format_scalar(state.reciprocal_condition);
      EXPECT_REAL_NEAR(static_cast<double>(state.energy), expected[i][j], 3e-10);
      EXPECT_LE(state.residual_norm, bethe::SolverOptions<Real>{}.residual_tolerance);
      EXPECT_LE(state.momentum_error, bethe::SolverOptions<Real>{}.residual_tolerance);
    }
}

TYPED_TEST(XXZOddContinuation, NativeFiveSiteEnergyAndSpinReversal)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const C = std::cos(Real{4} * std::atan(Real{1}) / Real{5});
  for (Real d : {-Real{1} / Real{2}, -Real{9} / Real{10}, -Real{1} + Real{128} * eps})
  {
    Real const energy =
        Real{5} * d / Real{4} + (-Real{3} * d - C - std::sqrt((d + C) * (d + C) + Real{4} * C * C)) / Real{2};
    auto const up = engine::continue_odd_polynomial(5, d, uni20::from_twice(std::int64_t{1}));
    auto const down = engine::continue_odd_polynomial(5, d, uni20::from_twice(std::int64_t{-1}));
    ASSERT_TRUE(up.equations_converged);
    ASSERT_TRUE(down.equations_converged);
    EXPECT_REAL_NEAR(up.energy, energy, Real{1024} * eps);
    EXPECT_EQ(up.energy, down.energy);
    EXPECT_EQ(up.coefficients, down.coefficients);
  }
}

TYPED_TEST(XXZOddContinuation, InvalidInputAndPrecisionFloor)
{
  using Real = TypeParam;
  auto const sz = uni20::from_twice(std::int64_t{1});
  for (Real scale :
       {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW((engine::PolynomialBetheSystem<Real>{7, 3, Real{0}, scale}), std::invalid_argument);
  EXPECT_THROW((engine::PolynomialBetheSystem<Real>{7, 3, uni20::numeric_limits<Real>::quiet_NaN(), Real{1}}),
               std::invalid_argument);
  for (Real d : {-Real{1}, Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW(engine::continue_odd_polynomial(7, d, sz), std::invalid_argument);
  EXPECT_THROW(engine::continue_odd_polynomial(6, -Real{1} / Real{2}, uni20::half_int{0}), std::invalid_argument);
  EXPECT_THROW(engine::continue_odd_polynomial(7, -Real{1} / Real{2}, uni20::half_int{0}), std::invalid_argument);
  EXPECT_THROW(engine::continue_odd_polynomial(7, -Real{1} / Real{2}, uni20::half_int{4}), std::invalid_argument);
  for (Real tol :
       {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW(engine::continue_odd_polynomial(7, -Real{1} / Real{2}, sz, {.residual_tolerance = tol}),
                 std::invalid_argument);
  // A requested precision below the arithmetic floor must not silently use
  // a looser stopping tolerance or claim a certified eigenstate.
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  auto const state = engine::continue_odd_polynomial(7, -Real{9} / Real{10}, sz,
                                                     {.residual_tolerance = eps * eps, .max_iterations = 100});
  EXPECT_FALSE(state.equations_converged);
  EXPECT_NE(state.status, engine::PolynomialContinuationStatus::equations_converged);
  EXPECT_LE(state.iterations, 100);
}
} // namespace
