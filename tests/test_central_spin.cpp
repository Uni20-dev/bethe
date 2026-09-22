// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "central_spin_ed.hpp"
#include "test_support.hpp"
#include <bethe/central_spin.hpp>
#include <complex>

namespace
{
namespace model = bethe::central_spin;
template <typename Real> class CentralSpin : public ::testing::Test {};
TYPED_TEST_SUITE(CentralSpin, test_support::RealTypes, test_support::PrecisionNames);

template <uni20::Real Real> void check_equations(model::State<Real> const& s)
{
  auto const& v = s.eigenvalue_variables;
  if (v.empty()) return;
  ASSERT_TRUE(s.reached_field.has_value());
  Real a{};
  for (Real value : s.couplings)
    a = std::max(a, std::abs(value));
  Real const h = std::abs(*s.reached_field) / a, t = Real{0.5} / (h + Real{0.5}), p = h / (h + Real{0.5});
  std::vector<Real> q(v.size());
  for (std::size_t j = 1; j < v.size(); ++j)
    q[j] = a / s.couplings[j - 1];
  Real sum{};
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t i = 0; i < v.size(); ++i)
  {
    Real f = v[i] * (v[i] - p), norm = Real{1} + std::abs(f);
    for (std::size_t j = 0; j < v.size(); ++j)
      if (i != j)
      {
        f -= t * (v[i] - v[j]) / (q[i] - q[j]);
        norm += std::abs(t / (q[i] - q[j])) * (std::abs(v[i]) + std::abs(v[j]));
      }
    EXPECT_REAL_NEAR(f / norm, Real{0}, Real{256} * Real(v.size()) * eps);
    sum += v[i];
  }
  auto const m = s.spin_reversed ? v.size() - s.up_spins : s.up_spins;
  EXPECT_REAL_NEAR(sum, Real(m) * p, Real{256} * Real(v.size()) * eps);
}

TYPED_TEST(CentralSpin, AllSmallSectorsAgainstIndependentED)
{
  using Real = TypeParam;
  for (auto const& input : {std::vector<double>{1}, std::vector<double>{1, .7, .3}, std::vector<double>{1, -.7, .3},
                            std::vector<double>{-1, -.7, -.3}, std::vector<double>{1, -2, .3, -.1}})
  {
    std::vector<Real> a(input.begin(), input.end());
    for (unsigned m = 0; m <= a.size() + 1; ++m)
      for (Real b : {Real{10}, Real{1}, Real{0.01}, Real{0}, Real{-1}})
      {
        SCOPED_TRACE(::testing::Message() << a.size() << "," << m << "," << double(b) << "," << input.front());
        auto const s =
            model::sector_ground_state<Real>(a, b, uni20::from_twice(std::int64_t(2 * m) - std::int64_t(a.size() + 1)));
        ASSERT_TRUE(s.converged) << int(s.status) << " t=" << double(s.continuation_parameter)
                                 << " stages=" << s.stages;
        ASSERT_TRUE(s.energy.has_value());
        EXPECT_NEAR(double(*s.energy), bethe::test::central_spin_exact_ground(input, double(b), m), 2e-10);
        ASSERT_NO_FATAL_FAILURE(check_equations(s));
      }
  }
}

TYPED_TEST(CentralSpin, TwoSpinAnalyticAndFieldReversal)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real a : {Real{1}, Real{-1}})
    for (Real b : {Real{0}, Real{1}, Real{10}, eps * eps})
    {
      auto const s = model::sector_ground_state<Real>(std::vector<Real>{a}, b, uni20::half_int{0});
      ASSERT_TRUE(s.converged);
      Real const exact = -a / Real{4} - std::hypot(a, b) / Real{2};
      EXPECT_REAL_NEAR(*s.energy, exact, Real{4096} * (Real{1} + std::abs(exact)) * eps);
      EXPECT_EQ(*s.reached_field, b);
      auto const r = model::sector_ground_state<Real>(std::vector<Real>{a}, -b, uni20::half_int{0});
      ASSERT_TRUE(r.converged);
      EXPECT_REAL_NEAR(*r.energy, *s.energy, Real{64} * eps);
    }
}

TYPED_TEST(CentralSpin, OneRootOriginalRationalEquation)
{
  using Real = TypeParam;
  std::vector<Real> a{Real{1}, Real{0.7}, Real{-0.3}};
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real b : {Real{1}, Real{0.1}, Real{0}})
  {
    auto const s = model::sector_ground_state<Real>(a, b, uni20::half_int{-1});
    ASSERT_TRUE(s.converged);
    Real const t = s.continuation_parameter;
    // v0=t*sum 1/(0-z), z=-Astar*lambda, Astar=1.
    Real const lambda = t / s.eigenvalue_variables[0];
    Real f = -Real{2} * b + Real{1} / lambda;
    for (Real aj : a)
      f += Real{1} / (lambda + Real{1} / aj);
    EXPECT_REAL_NEAR(f, Real{0}, Real{8192} * eps);
    Real energy = Real{0.5} / lambda - b / Real{2};
    for (Real aj : a)
      energy += aj / Real{4};
    EXPECT_REAL_NEAR(*s.energy, energy, Real{512} * eps);
  }
}

TYPED_TEST(CentralSpin, TwoRootsOriginalRationalEquations)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  std::vector<Real> a{Real{1}, Real{7} / Real{10}, Real{3} / Real{10}};
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real b : {Real{2}, Real{1}, Real{1} / Real{4}})
  {
    auto const state = model::sector_ground_state<Real>(a, b, uni20::half_int{0});
    ASSERT_TRUE(state.converged);
    Real const t = state.continuation_parameter;
    auto const& v = state.eigenvalue_variables;
    // P(z)=z^2-s*z+r. v_i=t*P'(q_i)/P(q_i) fixes s,r from q0=0,q1=1.
    Real const s = (Real{2} * t - v[1]) / (t - v[1] * (Real{1} + t / v[0]));
    Real const r = -t * s / v[0];
    C const d = std::sqrt(C{s * s - Real{4} * r, Real{0}});
    std::array<C, 2> lambda{-(C{s, Real{0}} + d) / Real{2}, -(C{s, Real{0}} - d) / Real{2}};
    for (std::size_t j = 0; j < 2; ++j)
    {
      C f = -Real{2} * b + Real{1} / lambda[j] - Real{2} / (lambda[j] - lambda[1 - j]);
      for (Real aj : a)
        f += Real{1} / (lambda[j] + Real{1} / aj);
      EXPECT_REAL_NEAR(f.real(), Real{0}, Real{65536} * eps);
      EXPECT_REAL_NEAR(f.imag(), Real{0}, Real{65536} * eps);
    }
    C energy = (Real{1} / lambda[0] + Real{1} / lambda[1]) / Real{2} - b / Real{2};
    for (Real aj : a)
      energy += aj / Real{4};
    EXPECT_REAL_NEAR(energy.real(), *state.energy, Real{65536} * eps);
    EXPECT_REAL_NEAR(energy.imag(), Real{0}, Real{65536} * eps);
  }
}

TYPED_TEST(CentralSpin, BudgetsAndExactLimits)
{
  using Real = TypeParam;
  std::vector<Real> a{Real{1}, Real{7} / Real{10}, Real{3} / Real{10}};
  model::SolverOptions<Real> options;
  options.max_iterations = 0;
  auto const seed = model::sector_ground_state<Real>(a, Real{1}, uni20::half_int{0}, options);
  EXPECT_FALSE(seed.converged);
  EXPECT_FALSE(seed.energy);
  EXPECT_FALSE(seed.reached_field);
  EXPECT_EQ(seed.iterations, 0U);
  EXPECT_EQ(seed.stages, 0U);
  EXPECT_EQ(seed.status, model::SolveStatus::iteration_limit);
  options.max_iterations = 10000;
  options.max_stages = 1;
  auto const partial = model::sector_ground_state<Real>(a, Real{-1}, uni20::half_int{0}, options);
  EXPECT_FALSE(partial.converged);
  ASSERT_TRUE(partial.energy);
  ASSERT_TRUE(partial.reached_field);
  EXPECT_LT(*partial.reached_field, Real{-1});
  EXPECT_EQ(partial.stages, 1U);
  auto const complete = model::sector_ground_state<Real>(a, *partial.reached_field, uni20::half_int{0});
  ASSERT_TRUE(complete.converged);
  EXPECT_REAL_NEAR(*partial.energy, *complete.energy, Real{8192} * uni20::numeric_limits<Real>::epsilon());
  options.max_stages = 0;
  auto const stage_seed = model::sector_ground_state<Real>(a, Real{1}, uni20::half_int{0}, options);
  EXPECT_EQ(stage_seed.status, model::SolveStatus::stage_limit);
  EXPECT_FALSE(stage_seed.energy);
  auto const polarized = model::sector_ground_state<Real>(a, Real{-2}, uni20::half_int{2}, options);
  ASSERT_TRUE(polarized.converged);
  EXPECT_REAL_NEAR(*polarized.energy, Real{-0.5}, Real{16} * uni20::numeric_limits<Real>::epsilon());
  auto const isolated = model::sector_ground_state<Real>({}, Real{3}, uni20::from_twice(std::int64_t{-1}), options);
  ASSERT_TRUE(isolated.converged);
  EXPECT_EQ(*isolated.energy, Real{-1.5});
}

TYPED_TEST(CentralSpin, Validation)
{
  using Real = TypeParam;
  auto const sz = uni20::half_int{0};
  EXPECT_THROW((void)model::sector_ground_state<Real>(std::vector<Real>{Real{1}, Real{1}, Real{2}}, Real{1}, sz),
               std::invalid_argument);
  EXPECT_THROW((void)model::sector_ground_state<Real>(std::vector<Real>{Real{0}}, Real{1}, sz), std::invalid_argument);
  EXPECT_THROW((void)model::sector_ground_state<Real>(std::vector<Real>{Real{1}}, Real{1}, uni20::half_int{2}),
               std::invalid_argument);
  EXPECT_THROW((void)model::sector_ground_state<Real>(std::vector<Real>{Real{1}, Real{2}}, Real{1}, sz),
               std::invalid_argument);
  Real const nan = uni20::numeric_limits<Real>::quiet_NaN();
  EXPECT_THROW((void)model::sector_ground_state<Real>(std::vector<Real>{nan}, Real{1}, sz), std::invalid_argument);
  EXPECT_THROW((void)model::sector_ground_state<Real>(std::vector<Real>{Real{1}}, nan, sz), std::invalid_argument);
  model::SolverOptions<Real> options;
  options.residual_tolerance = Real{0};
  EXPECT_THROW((void)model::sector_ground_state<Real>(std::vector<Real>{Real{1}}, Real{1}, sz, options),
               std::invalid_argument);
}

TYPED_TEST(CentralSpin, JacobianAndTangent)
{
  using Real = TypeParam;
  model::detail::System<Real> const system({Real{0}, Real{1}, Real{2}, Real{-3}}, 2);
  std::vector<Real> x{Real{-1} / Real{10}, Real{1} / Real{5}, Real{1} / Real{20}}, jac;
  Real const t = Real{1} / Real{4}, p = Real{1} - t;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), h = std::cbrt(eps);
  (void)system.evaluate(x, t, p, &jac);
  auto raw = [&](std::vector<Real> const& values, Real time, std::size_t i) {
    auto const v = system.expand(values, Real{1} - time);
    Real f = v[i] * (v[i] - (Real{1} - time));
    for (std::size_t j = 0; j < 4; ++j)
      if (j != i) f -= time * (v[i] - v[j]) / (system.q[i] - system.q[j]);
    return f;
  };
  auto const v = system.expand(x, p);
  for (std::size_t i = 0; i < 4; ++i)
  {
    Real scale = Real{1} + std::abs(v[i] * (v[i] - p));
    for (std::size_t j = 0; j < 4; ++j)
      if (j != i) scale += std::abs(t / (system.q[i] - system.q[j])) * (std::abs(v[i]) + std::abs(v[j]));
    for (std::size_t j = 0; j < 3; ++j)
    {
      auto plus = x, minus = x;
      plus[j] += h;
      minus[j] -= h;
      EXPECT_REAL_NEAR((raw(plus, t, i) - raw(minus, t, i)) / (Real{2} * h * scale), jac[i * 3 + j],
                       Real{2048} * h * h);
    }
  }
  // At an actual solution the analytic tangent must annihilate dF/dt.
  std::vector<Real> a{Real{1}, Real{1} / Real{2}, Real{-1} / Real{3}};
  auto const state = model::sector_ground_state<Real>(a, Real{1.5}, uni20::half_int{0});
  ASSERT_TRUE(state.converged);
  x = {state.eigenvalue_variables[0], state.eigenvalue_variables[2], state.eigenvalue_variables[3]};
  std::vector<Real> slope;
  ASSERT_TRUE(system.tangent(x, t, p, slope));
  auto plus = x, minus = x;
  for (std::size_t j = 0; j < 3; ++j)
  {
    plus[j] += h * slope[j];
    minus[j] -= h * slope[j];
  }
  for (std::size_t i = 0; i < 4; ++i)
    EXPECT_REAL_NEAR((raw(plus, t + h, i) - raw(minus, t - h, i)) / (Real{2} * h), Real{0}, Real{4096} * h * h);
}

TYPED_TEST(CentralSpin, PermutationScaleAndSpinReversal)
{
  using Real = TypeParam;
  std::vector<Real> a{Real{1}, Real{-2}, Real{1} / Real{3}, Real{1} / Real{7}};
  auto const sz = uni20::from_twice(std::int64_t{-1});
  auto const base = model::sector_ground_state<Real>(a, Real{1}, sz);
  ASSERT_TRUE(base.converged);
  auto const reverse = model::sector_ground_state<Real>(a, Real{-1}, -sz);
  ASSERT_TRUE(reverse.converged);
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  EXPECT_REAL_NEAR(*reverse.energy, *base.energy, Real{64} * eps);
  std::reverse(a.begin(), a.end());
  auto const perm = model::sector_ground_state<Real>(a, Real{1}, sz);
  ASSERT_TRUE(perm.converged);
  EXPECT_REAL_NEAR(*perm.energy, *base.energy, Real{64} * eps);
  for (std::size_t j = 0; j < a.size(); ++j)
    EXPECT_REAL_NEAR(perm.eigenvalue_variables[j + 1], base.eigenvalue_variables[a.size() - j], Real{64} * eps);
  Real const factor = eps * eps;
  for (auto& v : a)
    v *= factor;
  auto const scaled = model::sector_ground_state<Real>(a, factor, sz);
  ASSERT_TRUE(scaled.converged);
  EXPECT_REAL_NEAR(*scaled.energy / factor, *base.energy, Real{16384} * eps);
}

TYPED_TEST(CentralSpin, WeakInverseFieldAndZeroFieldMultiplet)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const large = Real{1} / (eps * eps);
  auto const high = model::sector_ground_state<Real>(std::vector<Real>{Real{1}}, large, uni20::half_int{0});
  ASSERT_TRUE(high.converged);
  EXPECT_REAL_NEAR(high.eigenvalue_variables[0] / high.continuation_parameter, Real{-1}, Real{64} * eps);
  std::vector<Real> a{Real{1}, Real{7} / Real{10}, Real{3} / Real{10}};
  auto const left = model::sector_ground_state<Real>(a, Real{0}, uni20::half_int{-1});
  auto const middle = model::sector_ground_state<Real>(a, Real{0}, uni20::half_int{0});
  auto const right = model::sector_ground_state<Real>(a, Real{0}, uni20::half_int{1});
  ASSERT_TRUE(left.converged && middle.converged && right.converged);
  EXPECT_REAL_NEAR(*left.energy, *middle.energy, Real{8192} * eps);
  EXPECT_REAL_NEAR(*right.energy, *middle.energy, Real{8192} * eps);
}

TYPED_TEST(CentralSpin, ClusteredCouplingsAndLargerBath)
{
  using Real = TypeParam;
  for (auto const& input :
       {std::vector<double>{1, 1.001, .25, -.7, -.7005}, std::vector<double>{-2, -1, .125, .5, .75, 1, 1.5}})
  {
    std::vector<Real> a(input.begin(), input.end());
    for (unsigned m = 1; m <= a.size(); ++m)
      for (Real b : {Real{0}, Real{0.01}, Real{1}})
      {
        SCOPED_TRACE(::testing::Message() << a.size() << "," << m << "," << double(b));
        auto const s =
            model::sector_ground_state<Real>(a, b, uni20::from_twice(std::int64_t(2 * m) - std::int64_t(a.size() + 1)));
        ASSERT_TRUE(s.converged) << int(s.status) << " stages=" << s.stages
                                 << " t=" << double(s.continuation_parameter);
        EXPECT_NEAR(double(*s.energy), bethe::test::central_spin_exact_ground(input, double(b), m), 2e-9);
        ASSERT_NO_FATAL_FAILURE(check_equations(s));
      }
  }
  std::vector<Real> a(23);
  for (std::size_t j = 0; j < a.size(); ++j)
    a[j] = Real{1} / Real(j + 1);
  for (Real b : {Real{0}, Real{1}})
  {
    auto const s = model::sector_ground_state<Real>(a, b, uni20::half_int{0});
    ASSERT_TRUE(s.converged) << int(s.status) << " stages=" << s.stages;
    ASSERT_NO_FATAL_FAILURE(check_equations(s));
  }
}
} // namespace
