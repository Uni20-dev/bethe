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
} // namespace
