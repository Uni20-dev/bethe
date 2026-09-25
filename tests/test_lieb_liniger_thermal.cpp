// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/lieb_liniger_thermal.hpp>
#include <bethe/lieb_liniger_thermo.hpp>

namespace
{
namespace model = bethe::lieb_liniger::thermal;
template <typename Real> class LiebLinigerThermal : public ::testing::Test {};
TYPED_TEST_SUITE(LiebLinigerThermal, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(LiebLinigerThermal, StableThermalFactors)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real e : {Real{0}, Real{1}, Real{20}, uni20::numeric_limits<Real>::max()})
  {
    auto const a = bethe::detail::thermal_factors(e, Real{1});
    auto const b = bethe::detail::thermal_factors(-e, Real{1});
    EXPECT_TRUE(uni20::isfinite(a.log_weight));
    EXPECT_TRUE(uni20::isfinite(b.log_weight));
    EXPECT_TRUE(uni20::isfinite(a.entropy));
    EXPECT_REAL_NEAR(a.filling + b.filling, Real{1}, Real{8} * eps);
    EXPECT_EQ(a.entropy, b.entropy);
    EXPECT_REAL_NEAR(b.log_weight - a.log_weight, e, Real{8} * eps * std::max(Real{1}, e));
  }
  auto const zero = bethe::detail::thermal_factors(Real{0}, Real{2});
  EXPECT_EQ(zero.filling, Real{1} / Real{2});
  EXPECT_REAL_NEAR(zero.entropy, std::log(Real{2}), Real{4} * eps);
  auto const extreme =
      bethe::detail::thermal_factors(uni20::numeric_limits<Real>::max(), uni20::numeric_limits<Real>::min());
  EXPECT_EQ(extreme.entropy, Real{0});
  EXPECT_EQ(extreme.filling, Real{0});
}

TYPED_TEST(LiebLinigerThermal, TonksGasNativeFugacitySeries)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  auto const state = model::equilibrium(Real{1} / eps, Real{1}, Real{-1});
  ASSERT_TRUE(state.converged) << int(state.status) << " nodes=" << state.nodes << " cutoff=" << state.cutoffs;
  // Gaussian integrals of the alternating fugacity series, no quadrature.
  bethe::detail::CompensatedSum<Real> pressure, density;
  Real const z = std::exp(Real{-1});
  Real power = z;
  for (int l = 1; l < 4 * uni20::numeric_limits<Real>::digits; ++l)
  {
    Real const term = power / std::sqrt(Real(l));
    density.add(term);
    pressure.add(term / Real(l));
    power *= -z;
  }
  Real const p = pressure.value() / (Real{2} * std::sqrt(pi)), n = density.value() / (Real{2} * std::sqrt(pi));
  Real const tol = Real{32768} * eps;
  EXPECT_REAL_NEAR(*state.pressure, p, tol * p);
  EXPECT_REAL_NEAR(*state.density, n, tol * n);
  EXPECT_REAL_NEAR(*state.energy_per_length, p / Real{2}, tol * p);
  EXPECT_REAL_NEAR(*state.entropy_per_length, Real{3} * p / Real{2} + n, tol);
  if constexpr (uni20::numeric_limits<Real>::digits > 53)
    EXPECT_GT(std::abs(Real(double(*state.pressure)) - *state.pressure), Real{4} * eps * p);
}

TYPED_TEST(LiebLinigerThermal, ScalingThermodynamicIdentitiesAndDerivatives)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  auto const a = model::equilibrium(Real{4}, Real{1}, Real{-1});
  ASSERT_TRUE(a.converged) << int(a.status) << " nodes=" << a.nodes
                           << " residual=" << uni20::format_real(a.nonlinear_residual);
  auto const b = model::equilibrium(Real{8}, Real{4}, Real{-4});
  ASSERT_TRUE(b.converged);
  EXPECT_EQ(*b.pressure, Real{8} * *a.pressure);
  EXPECT_EQ(*b.density, Real{2} * *a.density);
  EXPECT_EQ(*b.entropy_per_length, Real{2} * *a.entropy_per_length);
  EXPECT_REAL_NEAR(*a.pressure + *a.energy_per_length + *a.density, *a.entropy_per_length, Real{65536} * eps);
  Real const step = std::sqrt(std::sqrt(eps));
  auto const mup = model::equilibrium(Real{4}, Real{1}, Real{-1} + step);
  auto const mum = model::equilibrium(Real{4}, Real{1}, Real{-1} - step);
  auto const tp = model::equilibrium(Real{4}, Real{1} + step, Real{-1});
  auto const tm = model::equilibrium(Real{4}, Real{1} - step, Real{-1});
  ASSERT_TRUE(mup.converged);
  ASSERT_TRUE(mum.converged);
  ASSERT_TRUE(tp.converged);
  ASSERT_TRUE(tm.converged);
  EXPECT_REAL_NEAR((*mup.pressure - *mum.pressure) / (Real{2} * step), *a.density, Real{100} * step * step);
  EXPECT_REAL_NEAR((*tp.pressure - *tm.pressure) / (Real{2} * step), *a.entropy_per_length, Real{100} * step * step);
}

TYPED_TEST(LiebLinigerThermal, OriginalEquationsAndZeroTemperatureLimit)
{
  using Real = TypeParam;
  model::Options<Real> options;
  options.tolerance = Real{1} / Real{100000000};
  options.max_nodes = 512;
  std::size_t iterations = 0;
  auto const mesh = model::detail::solve_mesh(Real{4}, Real{1}, Real{-1}, Real{8},
                                              bethe::detail::gauss_legendre<Real>(64), options, iterations);
  ASSERT_TRUE(mesh.converged);
  Real const pi = Real{4} * std::atan(Real{1});
  for (std::size_t i = 0; i < mesh.k.size(); ++i)
  {
    Real epsilon = mesh.k[i] * mesh.k[i] + Real{1}, rho = Real{1} / (Real{2} * pi);
    for (std::size_t j = 0; j < mesh.k.size(); ++j)
    {
      Real const minus = mesh.k[i] - mesh.k[j], plus = mesh.k[i] + mesh.k[j];
      Real const kernel = Real{4} / pi * (Real{1} / (Real{16} + minus * minus) + Real{1} / (Real{16} + plus * plus));
      epsilon -= mesh.w[j] * kernel * std::log1p(std::exp(-mesh.epsilon[j]));
      rho += mesh.w[j] * kernel * mesh.rho_total[j] / (Real{1} + std::exp(mesh.epsilon[j]));
    }
    EXPECT_REAL_NEAR(epsilon, mesh.epsilon[i], options.tolerance * (Real{1} + epsilon));
    EXPECT_REAL_NEAR(rho, mesh.rho_total[i], options.tolerance);
  }
  auto const ground = bethe::lieb_liniger::thermo::ground_state(Real{4}, Real{1});
  ASSERT_TRUE(ground.converged);
  Real const p0 = *ground.chemical_potential - *ground.energy_per_length;
  auto const warm = model::equilibrium(Real{4}, Real{1} / Real{2}, *ground.chemical_potential, options);
  auto const cold = model::equilibrium(Real{4}, Real{1} / Real{4}, *ground.chemical_potential, options);
  ASSERT_TRUE(warm.converged) << int(warm.status) << " " << warm.nodes;
  ASSERT_TRUE(cold.converged) << int(cold.status) << " " << cold.nodes;
  EXPECT_GT(*warm.pressure, p0);
  EXPECT_GT(*cold.pressure, p0);
  EXPECT_LT(*cold.pressure - p0, (*warm.pressure - p0) / Real{3});
  EXPECT_LT(*cold.entropy_per_length, *warm.entropy_per_length);
  EXPECT_REAL_NEAR(*cold.density, Real{1}, Real{1} / Real{100});
}

TYPED_TEST(LiebLinigerThermal, FixedDensityRoundTripAndScaling)
{
  using Real = TypeParam;
  auto const reference = model::equilibrium(Real{4}, Real{1}, Real{-1});
  ASSERT_TRUE(reference.converged);
  auto const canonical = model::at_density(Real{4}, Real{1}, *reference.density);
  ASSERT_TRUE(canonical.converged) << int(canonical.status) << " evaluations=" << canonical.evaluations;
  ASSERT_TRUE(canonical.state);
  Real const tol = Real{262144} * uni20::numeric_limits<Real>::epsilon();
  EXPECT_REAL_NEAR(canonical.state->chemical_potential, Real{-1}, tol);
  EXPECT_REAL_NEAR(*canonical.state->pressure, *reference.pressure, tol * *reference.pressure);
  auto const scaled = model::at_density(Real{8}, Real{4}, Real{2} * *reference.density);
  ASSERT_TRUE(scaled.converged) << int(scaled.status);
  EXPECT_REAL_NEAR(scaled.state->chemical_potential, Real{-4}, Real{4} * tol);
  EXPECT_REAL_NEAR(*scaled.state->energy_per_length, Real{8} * *canonical.state->energy_per_length,
                   Real{8} * tol * *canonical.state->energy_per_length);
  EXPECT_LT(canonical.evaluations, 32u);
  EXPECT_LE(canonical.density_error, model::DensityOptions<Real>{}.tolerance / Real{2});
}

TYPED_TEST(LiebLinigerThermal, FixedDensityFailures)
{
  using Real = TypeParam;
  EXPECT_THROW(model::at_density(Real{4}, Real{1}, Real{0}), std::invalid_argument);
  EXPECT_THROW(model::at_density(Real{4}, Real{1}, Real{1}, {.tolerance = Real{0}}), std::invalid_argument);
  EXPECT_THROW(model::at_density(Real{4}, Real{1}, Real{1}, {.equilibrium = {.max_nodes = 513}, .max_evaluations = 0}),
               std::invalid_argument);
  for (std::size_t budget : {0u, 1u})
  {
    auto const failed = model::at_density(Real{4}, Real{1}, Real{1} / Real{10}, {.max_evaluations = budget});
    EXPECT_FALSE(failed.converged);
    EXPECT_FALSE(failed.state);
    EXPECT_EQ(failed.status, model::Status::density_limit);
    EXPECT_EQ(failed.evaluations, budget);
  }
  auto const inner = model::at_density(Real{4}, Real{1}, Real{1}, {.equilibrium = {.max_cutoffs = 0}});
  EXPECT_EQ(inner.status, model::Status::cutoff_limit);
  EXPECT_FALSE(inner.state);
  EXPECT_EQ(inner.evaluations, 1u);
}

TYPED_TEST(LiebLinigerThermal, FixedDensityDiluteAndDegenerate)
{
  using Real = TypeParam;
  model::DensityOptions<Real> options;
  options.tolerance = Real{1} / Real{100000000};
  options.equilibrium.tolerance = options.tolerance / Real{8};
  for (Real density : {Real{1} / Real{100}, Real{1}})
  {
    auto const state = model::at_density(Real{4}, Real{1}, density, options);
    ASSERT_TRUE(state.converged) << int(state.status) << " evaluations=" << state.evaluations;
    EXPECT_REAL_NEAR(*state.state->density, density, options.tolerance * density);
    EXPECT_REAL_NEAR(*state.state->energy_per_length + *state.state->pressure -
                         state.state->chemical_potential * *state.state->density,
                     *state.state->entropy_per_length, options.tolerance);
    if (density == Real{1})
      EXPECT_GT(state.state->chemical_potential, Real{0});
    else
      EXPECT_LT(state.state->chemical_potential, Real{0});
  }
}

TYPED_TEST(LiebLinigerThermal, IndependentBudgetsAndValidation)
{
  using Real = TypeParam;
  EXPECT_THROW(model::equilibrium(Real{0}, Real{1}, Real{0}), std::invalid_argument);
  EXPECT_THROW(model::equilibrium(Real{1}, Real{0}, Real{0}), std::invalid_argument);
  EXPECT_THROW(model::equilibrium(Real{1}, Real{1}, uni20::numeric_limits<Real>::infinity()), std::invalid_argument);
  EXPECT_THROW(model::equilibrium(Real{1}, Real{1}, Real{0}, {.max_nodes = 513}), std::invalid_argument);
  EXPECT_THROW(model::equilibrium(Real{1}, Real{1}, Real{0}, {.initial_cutoff = Real{0}}), std::invalid_argument);
  auto const iter = model::equilibrium(Real{4}, Real{1}, Real{-1}, {.max_iterations = 0});
  EXPECT_EQ(iter.status, model::Status::iteration_limit);
  EXPECT_FALSE(iter.pressure);
  auto const mesh = model::equilibrium(Real{4}, Real{1}, Real{-1}, {.initial_nodes = 16, .max_nodes = 16});
  EXPECT_EQ(mesh.status, model::Status::mesh_limit);
  EXPECT_FALSE(mesh.density);
  auto const cutoff = model::equilibrium(Real{4}, Real{1}, Real{-1}, {.max_cutoffs = 1});
  EXPECT_EQ(cutoff.status, model::Status::cutoff_limit);
  EXPECT_FALSE(cutoff.energy_per_length);
  auto const empty = model::equilibrium(Real{4}, Real{1}, Real{-1}, {.max_cutoffs = 0});
  EXPECT_EQ(empty.status, model::Status::cutoff_limit);
  EXPECT_EQ(empty.iterations, 0u);
  auto const short_domain =
      model::equilibrium(Real{4}, Real{1}, Real{-1}, {.max_cutoffs = 2, .initial_cutoff = Real{1}});
  EXPECT_EQ(short_domain.status, model::Status::cutoff_limit);
  EXPECT_FALSE(short_domain.pressure);
  Real const large = uni20::numeric_limits<Real>::max() / Real{16};
  auto const overflow = model::equilibrium(Real{4} * std::sqrt(large), large, -large);
  EXPECT_EQ(overflow.status, model::Status::precision_limit);
  EXPECT_FALSE(overflow.pressure);
  auto const dilute =
      model::equilibrium(Real{4}, Real{1}, -Real(4 * uni20::numeric_limits<Real>::max_exponent), {.max_nodes = 32});
  EXPECT_EQ(dilute.status, model::Status::precision_limit);
  EXPECT_FALSE(dilute.density);
  // Nonzero subnormal thermal weights are not necessarily accurate enough.
  std::size_t updates = 0;
  auto const subnormal =
      model::detail::solve_mesh(Real{4}, Real{1}, std::log(uni20::numeric_limits<Real>::min()) - Real{32}, Real{4},
                                bethe::detail::gauss_legendre<Real>(16), model::Options<Real>{}, updates);
  EXPECT_FALSE(subnormal.converged);
  EXPECT_EQ(subnormal.status, model::Status::precision_limit);
}
} // namespace
