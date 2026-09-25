// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/lieb_liniger.hpp>
#include <bethe/lieb_liniger_thermo.hpp>

namespace
{
namespace model = bethe::lieb_liniger::thermo;
template <typename Real> class LiebLinigerThermo : public ::testing::Test {};
TYPED_TEST_SUITE(LiebLinigerThermo, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(LiebLinigerThermo, SharedQuadraturePolynomialMoments)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {1u, 3u, 16u, 31u})
  {
    auto const rule = bethe::detail::gauss_legendre<Real>(n);
    for (std::size_t power = 0; power < 2 * n; ++power)
    {
      bethe::detail::CompensatedSum<Real> sum;
      for (std::size_t j = 0; j < n; ++j)
      {
        Real monomial{1};
        for (std::size_t p = 0; p < power; ++p)
          monomial *= rule.x[j];
        sum.add(rule.w[j] * monomial);
      }
      Real const expected = power % 2 ? Real{0} : Real{2} / Real(power + 1);
      EXPECT_REAL_NEAR(sum.value(), expected, Real{256} * eps);
    }
  }
}

TYPED_TEST(LiebLinigerThermo, NativePrecisionScalingAndFiniteRings)
{
  using Real = TypeParam;
  auto const s = model::ground_state(Real{4}, Real{1});
  ASSERT_TRUE(s.converged) << int(s.status) << " nodes=" << s.nodes << " iterations=" << s.iterations;
  auto const scaled = model::ground_state(Real{8}, Real{2});
  ASSERT_TRUE(scaled.converged);
  EXPECT_EQ(*scaled.fermi_rapidity, Real{2} * *s.fermi_rapidity);
  EXPECT_EQ(*scaled.energy_per_length, Real{8} * *s.energy_per_length);
  EXPECT_EQ(*scaled.chemical_potential, Real{4} * *s.chemical_potential);
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  if constexpr (uni20::numeric_limits<Real>::digits > 53)
    EXPECT_GT(std::abs(Real(double(*s.energy_per_length)) - *s.energy_per_length), Real{8} * eps);
  auto const small = bethe::lieb_liniger::ground_state(24, Real{24}, Real{4});
  auto const large = bethe::lieb_liniger::ground_state(48, Real{48}, Real{4});
  ASSERT_TRUE(small.converged);
  ASSERT_TRUE(large.converged);
  Real const e24 = small.energy / Real{24}, e48 = large.energy / Real{48};
  EXPECT_LT(e24, *s.energy_per_length);
  EXPECT_LT(e48, *s.energy_per_length);
  EXPECT_LT(std::abs(e48 - *s.energy_per_length), std::abs(e24 - *s.energy_per_length) / Real{3});
  EXPECT_REAL_NEAR((Real{4} * e48 - e24) / Real{3}, *s.energy_per_length, Real{1} / Real{1000000});
}

TYPED_TEST(LiebLinigerThermo, ChemicalPotentialAndKernelEquations)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  auto const s = model::ground_state(Real{4}, Real{1});
  ASSERT_TRUE(s.converged);
  auto const m = model::detail::mesh(Real{4}, *s.fermi_rapidity, bethe::detail::gauss_legendre<Real>(s.nodes));
  ASSERT_TRUE(m);
  for (std::size_t i = 0; i < m->k.size(); ++i)
  {
    Real rho = Real{1} / (Real{2} * pi), e = m->k[i] * m->k[i] - m->mu;
    for (std::size_t j = 0; j < m->k.size(); ++j)
    {
      Real const x = m->k[i] - m->k[j], kernel = Real{4} / (pi * (Real{16} + x * x));
      rho += m->w[j] * kernel * m->rho[j];
      e += m->w[j] * kernel * m->epsilon[j];
    }
    EXPECT_REAL_NEAR(rho, m->rho[i], Real{512} * eps);
    EXPECT_REAL_NEAR(e, m->epsilon[i], Real{4096} * eps);
    EXPECT_LT(m->epsilon[i], Real{0});
  }
  Real const step = std::sqrt(std::sqrt(eps));
  auto const plus = model::ground_state(Real{4}, Real{1} + step);
  auto const minus = model::ground_state(Real{4}, Real{1} - step);
  ASSERT_TRUE(plus.converged);
  ASSERT_TRUE(minus.converged);
  EXPECT_REAL_NEAR((*plus.energy_per_length - *minus.energy_per_length) / (Real{2} * step), *s.chemical_potential,
                   Real{100} * step * step);
}

TYPED_TEST(LiebLinigerThermo, StrongAndWeakCouplingLimits)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  Real const gamma = Real{1} / std::sqrt(std::sqrt(eps));
  auto const strong = model::ground_state(gamma, Real{1});
  ASSERT_TRUE(strong.converged);
  EXPECT_REAL_NEAR(*strong.energy_per_length, pi * pi / Real{3} * (Real{1} - Real{4} / gamma),
                   Real{100} / (gamma * gamma));
  EXPECT_REAL_NEAR(*strong.chemical_potential, pi * pi * (Real{1} - Real{16} / (Real{3} * gamma)),
                   Real{400} / (gamma * gamma));
  model::Options<Real> options;
  options.tolerance = Real{1} / Real{100000000};
  options.max_nodes = 512;
  auto const weak = model::ground_state(Real{1} / Real{100}, Real{1}, options);
  ASSERT_TRUE(weak.converged) << int(weak.status) << " " << weak.nodes;
  Real const expected = Real{1} / Real{100} - Real{4} / (Real{3000} * pi);
  EXPECT_REAL_NEAR(*weak.energy_per_length, expected, Real{1} / Real{10000});
}

TYPED_TEST(LiebLinigerThermo, FailureContracts)
{
  using Real = TypeParam;
  EXPECT_THROW(model::ground_state(Real{0}, Real{1}), std::invalid_argument);
  EXPECT_THROW(model::ground_state(Real{1}, Real{-1}), std::invalid_argument);
  EXPECT_THROW(model::ground_state(Real{1}, Real{1}, {.initial_nodes = 3}), std::invalid_argument);
  auto const budget = model::ground_state(Real{4}, Real{1}, {.max_iterations = 0});
  EXPECT_EQ(budget.status, model::Status::density_limit);
  EXPECT_FALSE(budget.energy_per_length);
  auto const mesh = model::ground_state(Real{4}, Real{1}, {.initial_nodes = 16, .max_nodes = 16});
  EXPECT_EQ(mesh.status, model::Status::mesh_limit);
  EXPECT_FALSE(mesh.converged);
  EXPECT_FALSE(mesh.energy_per_length);
  auto const weak = model::ground_state(Real{1} / Real{100000000}, Real{1}, {.initial_nodes = 4, .max_nodes = 8});
  EXPECT_FALSE(weak.converged);
  EXPECT_FALSE(weak.energy_per_length);
  // A tiny energy is not evidence of mesh convergence at a tiny coupling.
  auto const unresolved =
      model::ground_state(Real{1} / Real{1000000000000000000LL}, Real{1},
                          {.tolerance = Real{1} / Real{1000000}, .initial_nodes = 4, .max_nodes = 16});
  EXPECT_FALSE(unresolved.converged);
  EXPECT_FALSE(unresolved.chemical_potential);
  Real const huge = uni20::numeric_limits<Real>::max() / Real{4};
  auto const overflow = model::ground_state(huge, huge);
  EXPECT_EQ(overflow.status, model::Status::precision_limit);
  EXPECT_FALSE(overflow.energy_per_length);
}
TYPED_TEST(LiebLinigerThermo, DispersionsBackflowAndReflection)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1});
  for (Real gamma : {Real{1}, Real{4}, Real{20}})
  {
    model::Solver<Real> solver(gamma, Real{1});
    ASSERT_TRUE(solver.background().converged);
    auto const& bg = solver.background();
    auto const mesh = model::detail::mesh(gamma, *bg.fermi_rapidity, bethe::detail::gauss_legendre<Real>(bg.nodes));
    ASSERT_TRUE(mesh);
    auto const n = mesh->k.size();
    for (auto branch : {model::Branch::type_i, model::Branch::type_ii})
      for (Real p : {pi / Real{4}, pi / Real{2}, pi, Real{2} * pi})
      {
        auto const point = solver.at_momentum(branch, p);
        ASSERT_TRUE(point.converged) << int(point.status) << " gamma=" << uni20::format_real(gamma)
                                     << " p=" << uni20::format_real(p) << " branch=" << int(branch);
        bool const hole = branch == model::Branch::type_ii;
        if (hole && p == Real{2} * pi)
        {
          EXPECT_EQ(*point.energy, Real{0});
          continue;
        }
        EXPECT_GT(*point.energy, Real{0});
        Real const lambda = *point.rapidity, q = *bg.fermi_rapidity;
        std::vector<Real> a(n * n), d(n);
        for (std::size_t i = 0; i < n; ++i)
        {
          d[i] = (std::atan((mesh->k[i] - q) / gamma) - std::atan((mesh->k[i] - lambda) / gamma)) / pi;
          for (std::size_t j = 0; j < n; ++j)
          {
            Real const x = mesh->k[i] - mesh->k[j];
            a[i * n + j] = Real(i == j) - mesh->w[j] * gamma / (pi * (gamma * gamma + x * x));
          }
        }
        ASSERT_TRUE(bethe::detail::newton_step(a, d));
        Real momentum = lambda - q, energy = lambda * lambda - q * q;
        for (std::size_t j = 0; j < n; ++j)
        {
          momentum += mesh->w[j] * d[j];
          energy += Real{2} * mesh->w[j] * mesh->k[j] * d[j];
        }
        Real const tolerance = Real{32768} * uni20::numeric_limits<Real>::epsilon();
        EXPECT_REAL_NEAR((hole ? -momentum : momentum), p, tolerance * (Real{1} + p));
        EXPECT_REAL_NEAR((hole ? -energy : energy), *point.energy, tolerance * (Real{1} + *point.energy));
        if (hole)
        {
          auto const mirror = solver.at_momentum(branch, Real{2} * pi - p);
          ASSERT_TRUE(mirror.converged);
          EXPECT_REAL_NEAR(*mirror.energy, *point.energy, tolerance * *point.energy);
        }
      }
  }
}

TYPED_TEST(LiebLinigerThermo, DispersionLimitsScalingAndTinyGaps)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  model::Solver<Real> solver(Real{4}, Real{1}), scaled(Real{8}, Real{2});
  Real const tiny = eps * eps;
  for (auto branch : {model::Branch::type_i, model::Branch::type_ii})
  {
    auto const a = solver.at_momentum(branch, pi / Real{2});
    auto const b = scaled.at_momentum(branch, pi);
    ASSERT_TRUE(a.converged);
    ASSERT_TRUE(b.converged);
    EXPECT_EQ(*b.energy, Real{4} * *a.energy);
    auto const near = solver.at_momentum(branch, tiny);
    ASSERT_TRUE(near.converged) << int(near.status);
    EXPECT_GT(*near.energy, Real{0});
    EXPECT_GT(*near.edge_distance, Real{0});
    auto const zero = solver.at_momentum(branch, Real{0});
    ASSERT_TRUE(zero.converged);
    EXPECT_EQ(*zero.energy, Real{0});
  }
  auto const particle = solver.at_momentum(model::Branch::type_i, tiny);
  auto const hole = solver.at_momentum(model::Branch::type_ii, tiny);
  EXPECT_REAL_NEAR(*particle.energy / tiny, *hole.energy / tiny, Real{8192} * eps);
  Real const gamma = Real{1} / std::sqrt(std::sqrt(eps));
  model::Solver<Real> strong(gamma, Real{1});
  for (auto branch : {model::Branch::type_i, model::Branch::type_ii})
    for (Real p : {pi / Real{4}, pi, Real{2} * pi})
    {
      auto const a = strong.at_momentum(branch, p);
      ASSERT_TRUE(a.converged) << int(a.status);
      Real const exact = p * (Real{2} * pi + (branch == model::Branch::type_i ? p : -p));
      EXPECT_REAL_NEAR(*a.energy, exact, Real{1000} / gamma);
    }
}

TYPED_TEST(LiebLinigerThermo, DispersionsAgainstFiniteRings)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1});
  model::Solver<Real> solver(Real{4}, Real{1});
  for (auto branch : {model::Branch::type_i, model::Branch::type_ii})
  {
    auto const point = solver.at_momentum(branch, pi / Real{2});
    ASSERT_TRUE(point.converged);
    Real previous_error{};
    for (std::size_t n : {24u, 48u, 96u})
    {
      auto const ground = bethe::lieb_liniger::ground_state(n, Real(n), Real{4});
      ASSERT_TRUE(ground.converged);
      auto numbers = ground.quantum_numbers;
      if (branch == model::Branch::type_i)
        numbers.back() += uni20::half_int(n / 4);
      else
        for (std::size_t j = n - n / 4; j < n; ++j)
          numbers[j] += uni20::half_int(1);
      auto const state = bethe::lieb_liniger::solve_real(Real(n), Real{4}, numbers);
      ASSERT_TRUE(state.converged);
      Real const error = std::abs(state.energy - ground.energy - *point.energy);
      if (n > 24) EXPECT_LT(error, previous_error * Real{3} / Real{5});
      previous_error = error;
      EXPECT_REAL_NEAR(state.momentum, pi / Real{2}, Real{64} * uni20::numeric_limits<Real>::epsilon());
    }
    EXPECT_LT(previous_error, Real{1} / Real{5});
  }
}

TYPED_TEST(LiebLinigerThermo, DispersionFailureContracts)
{
  using Real = TypeParam;
  model::Solver<Real> solver(Real{4}, Real{1});
  EXPECT_THROW(solver.at_momentum(model::Branch::type_i, Real{-1}), std::invalid_argument);
  EXPECT_THROW(solver.at_momentum(model::Branch::type_ii, Real{7}), std::invalid_argument);
  EXPECT_THROW(solver.at_momentum(static_cast<model::Branch>(123), Real{0}), std::invalid_argument);
  EXPECT_THROW(solver.at_momentum(model::Branch::type_i, uni20::numeric_limits<Real>::infinity()),
               std::invalid_argument);
  model::Solver<Real> failed(Real{4}, Real{1}, {.max_iterations = 0});
  auto const no_background = failed.at_momentum(model::Branch::type_i, Real{0});
  EXPECT_FALSE(no_background.energy);
  EXPECT_EQ(no_background.status, model::Status::density_limit);
  model::Solver<Real> budget(Real{4}, Real{1}, {.max_momentum_iterations = 0});
  auto const limited = budget.at_momentum(model::Branch::type_i, Real{1});
  EXPECT_FALSE(limited.energy);
  EXPECT_EQ(limited.status, model::Status::momentum_limit);
  auto const huge = solver.at_momentum(model::Branch::type_i, uni20::numeric_limits<Real>::max() / Real{4});
  EXPECT_FALSE(huge.energy);
  EXPECT_EQ(huge.status, model::Status::precision_limit);
  // A nonzero subnormal input need not resolve its requested relative error.
  Real const subnormal = uni20::numeric_limits<Real>::min() * uni20::numeric_limits<Real>::epsilon();
  ASSERT_GT(subnormal, Real{0});
  auto const underresolved = solver.at_momentum(model::Branch::type_i, subnormal);
  EXPECT_FALSE(underresolved.energy);
  EXPECT_EQ(underresolved.status, model::Status::precision_limit);
}
} // namespace
