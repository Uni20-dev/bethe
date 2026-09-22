// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/gaudin_yang.hpp>
#include <bethe/heisenberg.hpp>
#include <bethe/hubbard.hpp>
#include <bethe/lieb_liniger.hpp>

namespace
{
namespace model = bethe::gaudin_yang;
template <typename Real> class GaudinYang : public ::testing::Test {};
TYPED_TEST_SUITE(GaudinYang, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(GaudinYang, TwoParticles)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (Real c : {Real{1} / Real{10000}, Real{1}, Real{2}, Real{100}, Real{1000000}})
  {
    auto const state = model::ground_state(1, 1, Real{1}, c);
    ASSERT_TRUE(state.converged) << state.iterations;
    auto const bosons = bethe::lieb_liniger::ground_state(2, Real{1}, c);
    ASSERT_TRUE(bosons.converged);
    EXPECT_REAL_NEAR(state.energy, bosons.energy, Real{512} * eps * state.energy);
    EXPECT_EQ(state.spin_rapidities[0], Real{0});
    EXPECT_EQ(state.momentum, Real{0});
  }
  // Jump condition q*tan(q/2)=c/2 has q=pi/2 at c=pi.
  auto const exact = model::ground_state(1, 1, Real{1}, pi);
  ASSERT_TRUE(exact.converged);
  EXPECT_REAL_NEAR(exact.energy, pi * pi / Real{2}, Real{512} * eps);
  if constexpr (uni20::numeric_limits<Real>::digits > 53)
    EXPECT_GT(std::abs(Real(double(pi * pi / Real{2})) - pi * pi / Real{2}), Real{512} * eps);
  Real const tiny = uni20::parse_real<Real>("1e-40");
  auto const weak = model::ground_state(1, 1, Real{1}, tiny);
  ASSERT_TRUE(weak.converged);
  EXPECT_REAL_NEAR(weak.energy / (Real{2} * tiny), Real{1}, Real{1024} * eps);
}

TYPED_TEST(GaudinYang, CouplingSweep)
{
  using Real = TypeParam;
  for (std::size_t n : {2, 4, 6, 8, 10, 18, 32})
    for (std::size_t m = 1; m <= n / 2; m += 2)
      for (Real c : {Real{1} / Real{1000}, Real{1}, Real{100}, Real{100000}})
      {
        SCOPED_TRACE(::testing::Message() << n << ',' << m << ',' << double(c));
        auto const state = model::ground_state(n - m, m, Real(n), c);
        ASSERT_TRUE(state.converged) << state.iterations << " residual=" << double(state.residual_norm)
                                     << " root_c=" << double(state.root_interaction);
        EXPECT_EQ(state.momenta.size(), n);
        EXPECT_EQ(state.spin_rapidities.size(), m);
        EXPECT_LT(state.residual_norm, Real{33} * uni20::numeric_limits<Real>::epsilon());
      }
}

TYPED_TEST(GaudinYang, FreeLimits)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t up = 0; up < 5; ++up)
    for (std::size_t down = 0; down < 5; ++down)
    {
      auto const state = model::ground_state(up, down, Real{2}, Real{0});
      EXPECT_TRUE(state.converged && state.free);
      EXPECT_EQ(state.iterations, 0u);
      EXPECT_TRUE(state.spin_rapidities.empty());
      Real energy = Real{0};
      for (auto const& modes : state.free_modes)
        for (auto j : modes)
          energy += pi * pi * Real(j * j);
      EXPECT_REAL_NEAR(state.energy, energy, Real{128} * eps * (Real{1} + energy));
      Real exact = Real{0};
      for (auto count : {up, down})
        exact += pi * pi * Real(count) * (Real(count) * Real(count) + (count % 2 ? -Real{1} : Real{2})) / Real{12};
      EXPECT_REAL_NEAR(state.energy, exact, Real{128} * eps * (Real{1} + exact));
      EXPECT_EQ(state.momentum_index, std::int64_t((up % 2 ? 0 : up / 2) + (down % 2 ? 0 : down / 2)));
    }
  auto const up = model::ground_state(6, 0, Real{2}, Real{100});
  auto const down = model::ground_state(0, 6, Real{2}, Real{100});
  EXPECT_EQ(up.energy, down.energy);
  EXPECT_EQ(up.energy, model::ground_state(6, 0, Real{2}, Real{0}).energy);
}

template <uni20::Real Real> void check_equations(model::State<Real> const& state)
{
  using Complex = uni20::complex<Real>;
  auto const n = state.up + state.down, m = std::min(state.up, state.down);
  Real const pi = Real{4} * std::atan(Real{1}), g = state.interaction * state.length,
             eps = uni20::numeric_limits<Real>::epsilon();
  auto const ratio = [](Real d, Real w) { return Complex{d, w} / Complex{d, -w}; };
  std::vector<Real> q, l;
  for (Real k : state.momenta)
    q.push_back(k * state.length);
  for (Real lambda : state.spin_rapidities)
    l.push_back(lambda * state.length);
  ASSERT_EQ(q.size(), n);
  ASSERT_EQ(l.size(), m);
  Real charge = Real{0}, spin = Real{0};
  for (std::size_t j = 0; j < n; ++j)
  {
    Real f = q[j] - pi * Real(state.quantum_numbers.charge[j].twice());
    Complex product{Real{1}, Real{0}};
    for (Real lambda : l)
    {
      f += Real{2} * std::atan(Real{2} * (q[j] - lambda) / g);
      product *= ratio(q[j] - lambda, g / Real{2});
    }
    Real const norm = g < Real{1} ? std::max(std::sqrt(g), std::abs(q[j])) : Real(n);
    charge = std::max(charge, std::abs(f) / norm);
    if (state.converged)
    {
      EXPECT_REAL_NEAR(product.real(), std::cos(q[j]), Real{512} * Real(n) * eps);
      EXPECT_REAL_NEAR(product.imag(), std::sin(q[j]), Real{512} * Real(n) * eps);
    }
    EXPECT_EQ(q[j], -q[n - 1 - j]);
    if (j) EXPECT_GT(q[j], q[j - 1]);
  }
  for (std::size_t a = 0; a < m; ++a)
  {
    Real f = -pi * Real(state.quantum_numbers.spin[a].twice());
    Complex lhs{Real{1}, Real{0}}, rhs = lhs;
    for (Real k : q)
    {
      f += Real{2} * std::atan(Real{2} * (l[a] - k) / g);
      lhs *= ratio(l[a] - k, g / Real{2});
    }
    for (std::size_t b = 0; b < m; ++b)
      if (a != b)
      {
        f -= Real{2} * std::atan((l[a] - l[b]) / g);
        rhs *= ratio(l[a] - l[b], g);
      }
    spin = std::max(spin, std::abs(f) / Real(n));
    if (state.converged)
    {
      EXPECT_REAL_NEAR(lhs.real(), rhs.real(), Real{512} * Real(n) * eps);
      EXPECT_REAL_NEAR(lhs.imag(), rhs.imag(), Real{512} * Real(n) * eps);
    }
    EXPECT_EQ(l[a], -l[m - 1 - a]);
    if (a) EXPECT_GT(l[a], l[a - 1]);
  }
  EXPECT_REAL_NEAR(charge, state.charge_residual, Real{2048} * Real(n) * eps);
  EXPECT_REAL_NEAR(spin, state.spin_residual, Real{2048} * Real(n) * eps);
  EXPECT_EQ(state.residual_norm, std::max(state.charge_residual, state.spin_residual));
}

TYPED_TEST(GaudinYang, OriginalEquationsAndDiagnostics)
{
  using Real = TypeParam;
  for (std::size_t n : {2, 6, 10, 18})
    for (std::size_t m = 1; m <= n / 2; m += 2)
      for (Real c : {Real{1} / Real{1000}, Real{1}, Real{1000}})
      {
        auto const state = model::ground_state(n - m, m, Real{1}, c);
        ASSERT_TRUE(state.converged);
        ASSERT_NO_FATAL_FAILURE(check_equations(state));
      }
}

TYPED_TEST(GaudinYang, Jacobian)
{
  using Real = TypeParam;
  Real const base = std::cbrt(uni20::numeric_limits<Real>::epsilon());
  for (std::size_t n : {6, 10})
    for (Real g : {Real{1} / Real{100}, Real{1}, Real{100}})
    {
      auto const state = model::ground_state(n / 2, n / 2, Real{1}, g);
      ASSERT_TRUE(state.converged);
      model::detail::GroundSystem<Real> system(n, n / 2);
      auto x = system.seed();
      for (std::size_t j = 0; j < system.nk; ++j)
        x[j] = state.momenta[system.nk + j];
      for (std::size_t a = 0; a < system.ns; ++a)
        x[system.nk + a] = state.spin_rapidities[system.ns + 1 + a] / system.scale(g);
      uni20::DenseMatrix<Real> jacobian(system.order, system.order);
      (void)system.evaluate(x, g, &jacobian);
      for (std::size_t col = 0; col < system.order; ++col)
      {
        Real const h = base * (Real{1} + std::abs(x[col]));
        auto plus = x, minus = x;
        plus[col] += h;
        minus[col] -= h;
        auto const fp = system.evaluate(plus, g), fm = system.evaluate(minus, g);
        for (std::size_t row = 0; row < system.order; ++row)
        {
          Real const numerical = (fp.residual[row] - fm.residual[row]) / (Real{2} * h);
          EXPECT_REAL_NEAR((jacobian[row, col]), numerical, Real{1000} * h * h * (Real{1} + std::abs(numerical)));
        }
      }
    }
}

TYPED_TEST(GaudinYang, WeakAndStrongLimits)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1});
  for (std::size_t n : {6, 10, 18})
    for (std::size_t m : {std::size_t{1}, n / 2})
    {
      auto const free = model::ground_state(n - m, m, Real{1}, Real{0});
      Real previous = Real{1};
      for (Real c : {Real{1} / Real{100}, Real{1} / Real{1000}, Real{1} / Real{10000}})
      {
        auto const weak = model::ground_state(n - m, m, Real{1}, c);
        ASSERT_TRUE(weak.converged);
        Real const error = std::abs((weak.energy - free.energy) / (Real{2} * c * Real(n - m) * Real(m)) - Real{1});
        EXPECT_LT(error, previous / Real{4});
        previous = error;
      }
      auto const xxx = bethe::heisenberg::sector_ground_state<Real>(n, uni20::from_twice(std::int64_t(n - 2 * m)));
      ASSERT_TRUE(xxx.converged);
      Real const omega = Real(n) / Real{2} - Real{2} * xxx.energy;
      Real const limiting = pi * pi * Real(n) * (Real(n) * Real(n) - Real{1}) / Real{3};
      previous = Real{1};
      for (Real c : {Real{1000} * Real(n), Real{10000} * Real(n), Real{100000} * Real(n)})
      {
        auto const strong = model::ground_state(n - m, m, Real{1}, c);
        ASSERT_TRUE(strong.converged);
        Real const error = std::abs((Real{1} - strong.energy / limiting) * c / (Real{2} * omega) - Real{1});
        EXPECT_LT(error, previous / Real{4});
        previous = error;
        for (std::size_t a = 0; a < m; ++a)
          EXPECT_REAL_NEAR(strong.spin_rapidities[a] / c, xxx.rapidities[a] / Real{2},
                           Real{100} * Real(n) * Real(n) / (c * c));
      }
    }
}

TYPED_TEST(GaudinYang, DiluteHubbardLimit)
{
  using Real = TypeParam;
  // Lattice spacing a: t=1/a^2, U/t=2*c*a, E_cont=(E_Hubbard+2*N)/a^2.
  for (std::size_t m : {1, 3})
  {
    std::size_t const n = 6;
    Real const length = Real{6}, c = Real{1};
    auto const continuum = model::ground_state(n - m, m, length, c);
    ASSERT_TRUE(continuum.converged);
    Real previous = Real{100};
    for (std::size_t sites : {24, 48, 96, 192})
    {
      Real const a = length / Real(sites);
      auto const lattice =
          bethe::hubbard::sector_ground_state(sites, n, uni20::from_twice(std::int64_t(n - 2 * m)), Real{2} * c * a);
      ASSERT_TRUE(lattice.converged);
      Real const energy = (lattice.energy + Real{2} * Real(n)) / (a * a);
      Real const error = std::abs(energy - continuum.energy);
      EXPECT_LT(error, previous / Real{3});
      previous = error;
    }
    EXPECT_LT(previous, Real{1} / Real{100});
  }
}

TYPED_TEST(GaudinYang, ScalingAndSpinReversal)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real c : {Real{1} / Real{100}, Real{1}, Real{100}})
  {
    auto const state = model::ground_state(5, 3, Real{1}, c);
    auto const reversed = model::ground_state(3, 5, Real{1}, c);
    auto const scaled = model::ground_state(5, 3, Real{2}, c / Real{2});
    ASSERT_TRUE(state.converged && reversed.converged && scaled.converged);
    EXPECT_FALSE(state.spin_reversed);
    EXPECT_TRUE(reversed.spin_reversed);
    EXPECT_EQ(state.energy, reversed.energy);
    EXPECT_REAL_NEAR(state.energy, Real{4} * scaled.energy, Real{128} * eps * state.energy);
    for (std::size_t j = 0; j < state.momenta.size(); ++j)
      EXPECT_EQ(state.momenta[j], Real{2} * scaled.momenta[j]);
    for (std::size_t a = 0; a < state.spin_rapidities.size(); ++a)
      EXPECT_EQ(state.spin_rapidities[a], Real{2} * scaled.spin_rapidities[a]);
  }
}

TYPED_TEST(GaudinYang, BudgetsAndInvalidInputs)
{
  using Real = TypeParam;
  for (std::size_t budget : {0, 1, 10})
  {
    auto const state = model::ground_state(3, 3, Real{1}, Real{1}, {.max_iterations = budget});
    EXPECT_FALSE(state.converged);
    EXPECT_EQ(state.iterations, budget);
    EXPECT_EQ(state.status, model::SolveStatus::iteration_limit);
    EXPECT_GT(state.root_interaction, state.interaction);
    ASSERT_NO_FATAL_FAILURE(check_equations(state));
  }
  for (auto counts : {std::array<std::size_t, 2>{2, 2}, {2, 1}, {1, 2}, {3, 2}})
    EXPECT_THROW((model::ground_state(counts[0], counts[1], Real{1}, Real{1})), std::invalid_argument);
  for (Real bad : {-Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    EXPECT_THROW((model::ground_state(1, 1, Real{1}, bad)), std::invalid_argument);
    EXPECT_THROW((model::ground_state(1, 1, bad, Real{1})), std::invalid_argument);
    EXPECT_THROW((model::ground_state(1, 1, Real{1}, Real{1}, {.residual_tolerance = bad})), std::invalid_argument);
  }
  EXPECT_THROW((model::ground_state(1, 1, Real{0}, Real{1})), std::invalid_argument);
  EXPECT_THROW((model::ground_state(1, 1, Real{1}, Real{1}, {.residual_tolerance = Real{0}})), std::invalid_argument);
  EXPECT_THROW((model::ground_state(std::numeric_limits<std::size_t>::max(), 1, Real{1}, Real{1})),
               std::invalid_argument);
  EXPECT_THROW((model::ground_state(6000000001ULL, 1, Real{1}, Real{1})), std::length_error);
  EXPECT_THROW((model::ground_state(1, 1, uni20::numeric_limits<Real>::max(), Real{2})), std::invalid_argument);
  Real const huge = uni20::numeric_limits<Real>::max();
  EXPECT_THROW((model::ground_state(2, 0, huge, Real{0})), std::underflow_error);
  EXPECT_THROW((model::ground_state(1, 1, huge, Real{1} / huge)), std::underflow_error);
  EXPECT_THROW((model::ground_state(2, 0, uni20::numeric_limits<Real>::min(), Real{0})), std::overflow_error);
  EXPECT_THROW((model::detail::physical_root(Real{1} / huge, huge)), std::underflow_error);
}
} // namespace
