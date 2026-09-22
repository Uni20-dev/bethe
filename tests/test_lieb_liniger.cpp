// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/lieb_liniger.hpp>
#include <bit>
#include <set>

namespace
{
namespace model = bethe::lieb_liniger;
using uni20::half_int;
using Numbers = model::QuantumNumbers;
template <typename Real> class LiebLiniger : public ::testing::Test {};
TYPED_TEST_SUITE(LiebLiniger, test_support::RealTypes, test_support::PrecisionNames);

template <uni20::Real Real> void check_state(model::State<Real> const& state)
{
  using std::abs;
  using std::atan;
  Real const pi = Real{4} * atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  ASSERT_EQ(state.momenta.size(), state.particles);
  ASSERT_EQ(state.quantum_numbers.size(), state.particles);
  EXPECT_EQ(state.converged, state.status == model::SolveStatus::converged);
  Real energy = Real{0}, momentum = Real{0};
  for (std::size_t j = 0; j < state.particles; ++j)
  {
    SCOPED_TRACE(j);
    Real const k = state.momenta[j];
    if (j)
    {
      ASSERT_GT(k, state.momenta[j - 1]);
    }
    energy += k * k;
    momentum += k;
    if (state.converged)
    {
      // Original logarithmic equation in PHYSICAL k, independently of the
      // solver's dimensionless variables and complementary weak-coupling phase.
      Real phase = state.length * k;
      for (std::size_t l = 0; l < state.particles; ++l)
        if (l != j) phase += Real{2} * atan((k - state.momenta[l]) / state.interaction);
      Real const target = pi * Real(state.quantum_numbers[j].twice());
      EXPECT_REAL_NEAR(phase, target, Real{256} * Real(state.particles + 1) * eps * (Real{1} + abs(target)));
    }
  }
  EXPECT_REAL_NEAR(energy, state.energy, Real{32} * Real(state.particles + 1) * eps * (Real{1} + energy));
  EXPECT_REAL_NEAR(momentum, state.momentum, Real{256} * Real(state.particles + 1) * eps * (Real{1} + abs(momentum)));
  EXPECT_REAL_NEAR(state.momentum, Real{2} * pi * Real(state.momentum_index) / state.length, Real{8} * eps);
}

// Independent two-body even-relative-wave boundary condition:
// 2*x*tan((x-pi*(I2-I1-1))/2)=c*length, q1,2=pi*(I1+I2) +/- x.
template <uni20::Real Real> Real two_body_x(Real g, std::int64_t difference)
{
  using std::atan;
  using std::sqrt;
  using std::tan;
  Real const pi = Real{4} * atan(Real{1});
  Real lo = pi * Real(difference - 1), hi = pi * Real(difference);
  if (difference == 1) hi = std::min(hi, sqrt(g));
  Real const origin = lo;
  for (int iteration = 0; iteration < uni20::numeric_limits<Real>::digits + 20; ++iteration)
  {
    Real const mid = lo + (hi - lo) / Real{2};
    if (mid == lo || mid == hi) break;
    if (Real{2} * mid * tan((mid - origin) / Real{2}) > g)
      hi = mid;
    else
      lo = mid;
  }
  return lo + (hi - lo) / Real{2};
}

TYPED_TEST(LiebLiniger, VacuumOneParticleAndParity)
{
  using Real = TypeParam;
  auto const vacuum = model::ground_state(0, Real{3}, Real{1});
  EXPECT_TRUE(vacuum.converged);
  EXPECT_EQ(vacuum.energy, Real{0});
  EXPECT_TRUE(vacuum.momenta.empty());
  auto const one = model::solve_real(Real{3}, Real{2}, Numbers{half_int{-3}}, {.max_iterations = 0});
  ASSERT_TRUE(one.converged);
  EXPECT_EQ(one.iterations, 0u);
  EXPECT_EQ(one.momentum_index, -3);
  ASSERT_NO_FATAL_FAILURE(check_state(one));
  EXPECT_EQ(model::ground_quantum_numbers(3), (Numbers{half_int{-1}, half_int{0}, half_int{1}}));
  EXPECT_EQ(model::ground_quantum_numbers(2),
            (Numbers{uni20::from_twice(std::int64_t{-1}), uni20::from_twice(std::int64_t{1})}));
}

TYPED_TEST(LiebLiniger, AnalyticNativePrecision)
{
  using Real = TypeParam;
  using std::atan;
  Real const pi = Real{4} * atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (Real length : {Real{1}, Real{7} / Real{3}})
  {
    SCOPED_TRACE(uni20::format_scalar(length));
    auto const two = model::ground_state(2, length, pi / length);
    ASSERT_TRUE(two.converged);
    Real const exact = pi * pi / (Real{2} * length * length);
    EXPECT_REAL_NEAR(two.energy, exact, Real{128} * eps * exact);
    if constexpr (uni20::numeric_limits<Real>::digits > 53)
      if (length == Real{1}) ASSERT_GT(std::abs(Real(static_cast<double>(exact)) - exact), Real{128} * eps * exact);
    EXPECT_REAL_NEAR(two.momenta[1], pi / (Real{2} * length), Real{128} * eps);
    test_support::expect_exact(uni20::parse_real<Real>(uni20::format_real(two.energy)), two.energy, "energy I/O");
    ASSERT_NO_FATAL_FAILURE(check_state(two));
  }
  Real const a = Real{3} * pi / Real{2} - Real{2} * atan(Real{2});
  auto const three = model::ground_state(3, Real{1}, a);
  ASSERT_TRUE(three.converged);
  EXPECT_REAL_NEAR(three.energy, Real{2} * a * a, Real{256} * eps);
  EXPECT_REAL_NEAR(three.momenta[2], a, Real{64} * eps);
}

TYPED_TEST(LiebLiniger, TwoBodyBoundaryCondition)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real g : {uni20::parse_real<Real>("1e-40"), uni20::parse_real<Real>("1e-8"), Real{1} / Real{10}, Real{1},
                 uni20::parse_real<Real>("4.12345678901234567890123456789"), Real{100}, Real{1000000}})
  {
    SCOPED_TRACE(uni20::format_scalar(g));
    auto const state = model::ground_state(2, Real{1}, g);
    ASSERT_TRUE(state.converged) << state.iterations;
    Real const x = two_body_x(g, 1);
    EXPECT_REAL_NEAR(state.energy, Real{2} * x * x, Real{128} * eps * state.energy);
    ASSERT_NO_FATAL_FAILURE(check_state(state));
  }
  for (int difference : {1, 2, 3, 4})
  {
    SCOPED_TRACE(difference);
    Numbers const labels{uni20::from_twice(std::int64_t{-1}), uni20::from_twice(std::int64_t(-1 + 2 * difference))};
    auto const state = model::solve_real(Real{1}, Real{2}, labels);
    ASSERT_TRUE(state.converged);
    Real const pi = Real{4} * std::atan(Real{1}), center = pi * Real(difference - 1);
    Real const x = two_body_x(Real{2}, difference);
    EXPECT_REAL_NEAR(state.energy, Real{2} * (center * center + x * x), Real{1024} * eps * state.energy);
    ASSERT_NO_FATAL_FAILURE(check_state(state));
  }
}

TYPED_TEST(LiebLiniger, CouplingsSizesAndLimits)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  for (std::size_t n : {2, 3, 4, 7, 16, 32, 64})
  {
    SCOPED_TRACE(n);
    for (Real c : {Real{1} / Real{1000000}, Real{1} / Real{10}, Real{1}, Real{10}, Real{1000000}})
    {
      SCOPED_TRACE(uni20::format_scalar(c));
      auto const state = model::ground_state(n, Real(n), c);
      ASSERT_TRUE(state.converged) << state.iterations;
      EXPECT_LE(state.residual_norm, Real{32} * eps);
      EXPECT_EQ(state.momentum_index, 0);
      ASSERT_NO_FATAL_FAILURE(check_state(state));
      for (std::size_t j = 0; j < n; ++j)
        EXPECT_REAL_NEAR(state.momenta[j], -state.momenta[n - 1 - j], Real{128} * eps);
    }
    Real const small = uni20::parse_real<Real>("1e-30");
    auto const weak = model::ground_state(n, Real{1}, small);
    ASSERT_TRUE(weak.converged);
    EXPECT_REAL_NEAR(weak.energy / (small * Real(n) * Real(n - 1)), Real{1}, Real{256} * eps + small * Real(n));
    auto const hard = model::ground_state(n, Real(n), Real{1000000});
    Real const tonks = pi * pi * Real(n * n - 1) / (Real{3} * Real(n));
    EXPECT_LT(hard.energy, tonks);
    EXPECT_REAL_NEAR(hard.energy / tonks, Real{1}, Real{5} / Real{1000000});
  }
}

TYPED_TEST(LiebLiniger, BoostReflectionAndLengthScaling)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  Numbers const labels{half_int{-3}, half_int{-1}, half_int{2}};
  auto const state = model::solve_real(Real{5}, Real{3}, labels);
  ASSERT_TRUE(state.converged);
  Numbers reflected, boosted;
  for (auto i = labels.rbegin(); i != labels.rend(); ++i)
    reflected.push_back(-*i);
  for (auto i : labels)
    boosted.push_back(i + half_int{2});
  auto const mirror = model::solve_real(Real{5}, Real{3}, reflected);
  auto const boost = model::solve_real(Real{5}, Real{3}, boosted);
  ASSERT_TRUE(mirror.converged && boost.converged);
  EXPECT_REAL_NEAR(mirror.energy, state.energy, Real{256} * eps);
  EXPECT_EQ(mirror.momentum_index, -state.momentum_index);
  Real const shift = Real{4} * pi / Real{5};
  EXPECT_REAL_NEAR(boost.energy, state.energy + Real{2} * shift * state.momentum + Real{3} * shift * shift,
                   Real{1024} * eps);
  EXPECT_EQ(boost.momentum_index, state.momentum_index + 6);
  for (std::size_t j = 0; j < 3; ++j)
    EXPECT_REAL_NEAR(boost.momenta[j], state.momenta[j] + shift, Real{128} * eps);
  auto const scaled = model::solve_real(Real{10}, Real{3} / Real{2}, labels);
  ASSERT_TRUE(scaled.converged);
  EXPECT_REAL_NEAR(scaled.energy * Real{4}, state.energy, Real{128} * eps);
  ASSERT_NO_FATAL_FAILURE(check_state(boost));
}

TYPED_TEST(LiebLiniger, Jacobian)
{
  using Real = TypeParam;
  using std::abs;
  Real const h = std::cbrt(uni20::numeric_limits<Real>::epsilon());
  auto const labels = model::ground_quantum_numbers(5);
  for (Real g : {Real{1} / Real{10}, Real{1}, Real{10}})
  {
    SCOPED_TRACE(uni20::format_scalar(g));
    model::detail::System<Real> const system{labels, g};
    auto q = system.seed();
    uni20::DenseMatrix<Real> jacobian(q.size(), q.size());
    (void)system.evaluate(q, &jacobian);
    for (std::size_t col = 0; col < q.size(); ++col)
    {
      auto plus = q, minus = q;
      plus[col] += h;
      minus[col] -= h;
      auto const fp = system.evaluate(plus), fm = system.evaluate(minus);
      for (std::size_t row = 0; row < q.size(); ++row)
      {
        SCOPED_TRACE(::testing::Message() << row << "," << col);
        Real const numerical = (fp.residual[row] - fm.residual[row]) / (Real{2} * h);
        EXPECT_REAL_NEAR((jacobian[row, col]), numerical, Real{1000} * h * h * (Real{1} + abs(numerical)));
        EXPECT_FLOATING_EQ((jacobian[row, col]), (jacobian[col, row]), 0);
      }
    }
  }
}

TYPED_TEST(LiebLiniger, FiniteWindowAndRanking)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {0, 1, 2, 3, 4})
    for (std::size_t padding : {0, 1, 2})
    {
      SCOPED_TRACE(::testing::Message() << n << "," << padding);
      auto const all =
          model::real_excitations(n, Real{4}, Real{2}, padding, {.count = std::numeric_limits<std::size_t>::max()});
      ASSERT_TRUE(all.converged());
      EXPECT_EQ(all.candidate_count, model::excitation_count(n, padding));
      ASSERT_EQ(all.levels.size(), all.candidate_count);
      std::set<Numbers> expected;
      unsigned const slots = static_cast<unsigned>(n + 2 * padding);
      for (unsigned mask = 0; mask < (1U << slots); ++mask)
        if (std::popcount(mask) == static_cast<int>(n))
        {
          Numbers q;
          for (unsigned j = 0; j < slots; ++j)
            if (mask & (1U << j))
              q.push_back(uni20::from_twice(2 * std::int64_t(j) - std::int64_t(n) - 2 * std::int64_t(padding) + 1));
          expected.insert(q);
        }
      for (std::size_t j = 0; j < all.levels.size(); ++j)
      {
        auto const& level = all.levels[j];
        ASSERT_EQ(expected.erase(level.state.quantum_numbers), 1u);
        ASSERT_TRUE(level.gap.has_value());
        EXPECT_EQ(*level.gap, level.state.energy - all.ground_state.energy);
        EXPECT_GE(*level.gap, -Real{128} * eps);
        if (j) EXPECT_LE(all.levels[j - 1].state.energy, level.state.energy);
        ASSERT_NO_FATAL_FAILURE(check_state(level.state));
      }
      EXPECT_TRUE(expected.empty());
      auto const few = model::real_excitations(n, Real{4}, Real{2}, padding, {.count = 2});
      ASSERT_EQ(few.levels.size(), std::min(std::size_t{2}, all.levels.size()));
      for (std::size_t j = 0; j < few.levels.size(); ++j)
        EXPECT_EQ(few.levels[j].state.quantum_numbers, all.levels[j].state.quantum_numbers);
    }
}

TYPED_TEST(LiebLiniger, BudgetsRestartAndInvalidInputs)
{
  using Real = TypeParam;
  auto const labels = model::ground_quantum_numbers(4);
  for (std::size_t budget : {0, 1})
  {
    auto const state = model::solve_real(Real{2}, Real{1}, labels, {.max_iterations = budget});
    EXPECT_FALSE(state.converged);
    EXPECT_EQ(state.iterations, budget);
    EXPECT_EQ(state.status, model::SolveStatus::iteration_limit);
    ASSERT_NO_FATAL_FAILURE(check_state(state));
  }
  auto const state = model::solve_real(Real{2}, Real{1}, labels);
  ASSERT_TRUE(state.converged);
  auto const restart =
      model::solve_real(Real{2}, Real{1}, labels, {.max_iterations = 0}, std::span<Real const>(state.momenta));
  EXPECT_TRUE(restart.converged);
  EXPECT_EQ(restart.iterations, 0u);
  EXPECT_EQ(restart.momenta, state.momenta);
  auto const failed = model::real_excitations(4, Real{2}, Real{1}, 1, {}, {.max_iterations = 0});
  EXPECT_FALSE(failed.converged());
  EXPECT_TRUE(failed.first_unconverged.has_value());
  EXPECT_TRUE(failed.levels.empty());
  for (Real bad :
       {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    EXPECT_THROW((model::ground_state(2, bad, Real{1})), std::invalid_argument);
    EXPECT_THROW((model::ground_state(2, Real{1}, bad)), std::invalid_argument);
    EXPECT_THROW((model::ground_state(2, Real{1}, Real{1}, {.residual_tolerance = bad})), std::invalid_argument);
  }
  for (Numbers const& bad : {Numbers{half_int{0}, half_int{1}}, Numbers{half_int{0}, half_int{0}, half_int{1}},
                             Numbers{half_int{1}, half_int{0}, half_int{-1}}})
    EXPECT_THROW((model::solve_real(Real{1}, Real{1}, bad)), std::invalid_argument);
  EXPECT_THROW((model::ground_state(std::numeric_limits<std::size_t>::max(), Real{1}, Real{1})), std::invalid_argument);
  EXPECT_THROW((model::excitation_count(1000, 1000)), std::length_error);
  EXPECT_THROW((model::real_excitations(4, Real{1}, Real{1}, 1, {.max_candidates = 14})), std::length_error);
  EXPECT_EQ(model::excitation_count(4, 1, 15), 15u);
  EXPECT_THROW((model::real_excitations(4, Real{1}, Real{1}, 1, {.count = 0})), std::invalid_argument);
  EXPECT_THROW((model::real_excitations(4, Real{1}, Real{1}, 1, {.max_candidates = 0})), std::invalid_argument);
  std::vector<Real> bad_roots(4, Real{0});
  EXPECT_THROW((model::solve_real(Real{1}, Real{1}, labels, {}, std::span<Real const>(bad_roots))),
               std::invalid_argument);
  bad_roots.resize(3);
  EXPECT_THROW((model::solve_real(Real{1}, Real{1}, labels, {}, std::span<Real const>(bad_roots))),
               std::invalid_argument);
  Real const largest = uni20::numeric_limits<Real>::max();
  EXPECT_THROW((model::ground_state(2, largest, Real{2})), std::invalid_argument);
  EXPECT_THROW((model::ground_state(2, uni20::numeric_limits<Real>::min(), uni20::numeric_limits<Real>::min())),
               std::invalid_argument);
  EXPECT_THROW((model::ground_state(2, uni20::numeric_limits<Real>::min(), Real{16})), std::overflow_error);
  EXPECT_THROW((model::ground_state(2, largest, Real{1})), std::underflow_error);
  auto const big = std::numeric_limits<std::int64_t>::max() / 2;
  EXPECT_THROW((model::solve_real(Real{1}, Real{1}, Numbers{half_int{big - 2}, half_int{big - 1}, half_int{big}})),
               std::overflow_error);
}

TYPED_TEST(LiebLiniger, NegligibleRootDoesNotCauseEnergyUnderflow)
{
  using Real = TypeParam;
  Real const g = Real{1024} * uni20::numeric_limits<Real>::min();
  Real const outer = std::sqrt(Real{3} * g);
  // The exact central root is zero. A negligible perturbation whose square
  // underflows must not invalidate the representable total energy.
  std::vector<Real> const roots{-outer, g, outer};
  EXPECT_EQ(g * g, Real{0});
  auto const state = model::solve_real(Real{1}, g, model::ground_quantum_numbers(3), {.max_iterations = 0},
                                       std::span<Real const>(roots));
  ASSERT_TRUE(state.converged);
  EXPECT_REAL_NEAR(state.energy / (Real{6} * g), Real{1}, Real{128} * uni20::numeric_limits<Real>::epsilon());
}

TEST(LiebLinigerFailures, UnattainableToleranceIsNotConvergence)
{
  auto const state = model::ground_state(4, 4.0, 1.0, {.residual_tolerance = 1e-30});
  EXPECT_FALSE(state.converged);
  EXPECT_EQ(state.status, model::SolveStatus::stalled);
  EXPECT_GT(state.residual_norm, 1e-30);
  EXPECT_LT(state.iterations, 10000u);
}

// Independent thermodynamic oracle, not a finite-root solve. At unit density
// solve rho(k)-integral K(k-k')rho(k')/(2*pi) dk'=1/(2*pi) on [-Q,Q],
// adjusting Q until integral rho=1. Composite Simpson quadrature is refined
// separately from the finite-N comparison. See Lieb I; Essler/de Klerk (32)-(33).
double bulk_energy(double c, std::size_t intervals)
{
  double const pi = 4 * std::atan(1.0);
  double lo = 0, hi = pi, energy = 0;
  for (int iteration = 0; iteration < 44; ++iteration)
  {
    double const q = (lo + hi) / 2, h = 2 * q / intervals;
    std::vector<double> k(intervals + 1), weights(intervals + 1);
    for (std::size_t j = 0; j <= intervals; ++j)
    {
      k[j] = -q + h * j;
      weights[j] = h / 3 * (j == 0 || j == intervals ? 1 : (j % 2 ? 4 : 2));
    }
    uni20::DenseMatrix<double> matrix(intervals + 1, intervals + 1), rho(intervals + 1, 1);
    for (std::size_t j = 0; j <= intervals; ++j)
    {
      rho[j, 0] = 1 / (2 * pi);
      for (std::size_t l = 0; l <= intervals; ++l)
        matrix[j, l] = double(j == l) - weights[l] * c / (pi * (c * c + (k[j] - k[l]) * (k[j] - k[l])));
    }
    uni20::linalg::solve_inplace(matrix, rho);
    double density = 0;
    energy = 0;
    for (std::size_t j = 0; j <= intervals; ++j)
    {
      density += weights[j] * rho[j, 0];
      energy += weights[j] * k[j] * k[j] * rho[j, 0];
    }
    if (density > 1)
      hi = q;
    else
      lo = q;
  }
  return energy;
}

TEST(LiebLinigerBulk, FiniteRingsApproachIntegralEquation)
{
  for (double c : {1.0, 4.0})
  {
    SCOPED_TRACE(c);
    double const coarse = bulk_energy(c, 64), fine = bulk_energy(c, 128);
    ASSERT_NEAR(coarse, fine, 2e-7);
    double previous_error = 1;
    for (std::size_t n : {16, 32, 64, 128})
    {
      auto const state = model::ground_state(n, double(n), c);
      ASSERT_TRUE(state.converged);
      double const error = std::abs(state.energy / n - fine);
      EXPECT_LT(error, previous_error / 3); // leading finite-size correction is O(N^-2)
      previous_error = error;
    }
    EXPECT_LT(previous_error, 2e-4);
  }
}
} // namespace
