// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/hubbard_continuum.hpp>

namespace
{
namespace model = bethe::hubbard::thermo;
template <typename Real> class HubbardContinuum : public ::testing::Test {};
TYPED_TEST_SUITE(HubbardContinuum, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(HubbardContinuum, IndependentReferencesAndMomentumWitnesses)
{
  using R = TypeParam;
  R const pi = R{4} * std::atan(R{1}), tol = R{4096} * uni20::numeric_limits<R>::epsilon();
  // Independent Fourier-Bessel spinon values from test_hubbard_thermo.cpp.
  R const energy = uni20::parse_real<R>("0.180195419367165283255675569344877169275747983");
  R const p = uni20::parse_real<R>("0.147435530616030589362559903539606986354219786");
  for (R sign : {R{1}, R{-1}})
  {
    auto const lower = model::two_spinon_continuum(R{4}, sign * p);
    auto const upper = model::two_spinon_continuum(R{4}, sign * R{2} * p);
    ASSERT_TRUE(lower.converged && upper.converged);
    EXPECT_REAL_NEAR(*lower.lower_energy, energy, tol);
    EXPECT_REAL_NEAR(*upper.upper_energy, R{2} * energy, tol);
    for (auto const* state : {&lower, &upper})
      for (auto momenta : {state->lower_momenta, state->upper_momenta})
      {
        EXPECT_GE(momenta[0], R{0});
        EXPECT_LE(momenta[0], pi);
        EXPECT_GE(momenta[1], R{0});
        EXPECT_LE(momenta[1], pi);
        EXPECT_REAL_NEAR(model::detail::wrap(momenta[0] + momenta[1] - state->momentum), R{0}, tol);
      }
  }
  auto const edge = model::two_spinon_continuum(R{4}, pi);
  ASSERT_TRUE(edge.converged);
  EXPECT_EQ(*edge.lower_energy, R{0});
  EXPECT_REAL_NEAR(*edge.upper_energy, R{2} * uni20::parse_real<R>("1.24228148941956239984591553194945497655659639"),
                   tol);
}

TYPED_TEST(HubbardContinuum, ScatteringInteriorAndHeisenbergLimit)
{
  using R = TypeParam;
  R const pi = R{4} * std::atan(R{1}), tol = R{8192} * uni20::numeric_limits<R>::epsilon();
  for (R u : {R{1}, R{4}, R{20}})
  {
    R const q = pi * R{0.75};
    auto const edge = model::two_spinon_continuum(u, q);
    ASSERT_TRUE(edge.converged);
    for (int i : {1, 2, 3})
    {
      R const p = q * R(i) / R{4};
      auto const a = model::dispersion(model::Branch::spinon, u, p);
      auto const b = model::dispersion(model::Branch::spinon, u, q - p);
      ASSERT_TRUE(a.converged && b.converged);
      EXPECT_GE(*a.energy + *b.energy + tol, *edge.lower_energy);
      EXPECT_LE(*a.energy + *b.energy - tol, *edge.upper_energy);
    }
  }
  R const u = R{10000}, q = pi / R{3}, j = R{4} / u;
  auto const edge = model::two_spinon_continuum(u, q);
  ASSERT_TRUE(edge.converged);
  EXPECT_REAL_NEAR(*edge.lower_energy / j, pi / R{2} * std::sin(q), R{1e-6});
  EXPECT_REAL_NEAR(*edge.upper_energy / j, pi * std::sin(q / R{2}), R{1e-6});
}

TYPED_TEST(HubbardContinuum, SharedBudgetAndValidation)
{
  using R = TypeParam;
  model::Options<R> controls;
  controls.max_evaluations = 0;
  auto const zero = model::two_spinon_continuum(R{4}, R{0}, controls);
  ASSERT_TRUE(zero.converged);
  EXPECT_EQ(*zero.lower_energy, R{0});
  EXPECT_EQ(*zero.upper_energy, R{0});
  EXPECT_EQ(zero.evaluations, 0u);
  auto const lower = model::dispersion(model::Branch::spinon, R{4}, R{1});
  ASSERT_TRUE(lower.converged);
  controls.max_evaluations = lower.evaluations;
  auto const failed = model::two_spinon_continuum(R{4}, R{1}, controls);
  EXPECT_FALSE(failed.converged);
  EXPECT_FALSE(failed.lower_energy);
  EXPECT_FALSE(failed.upper_energy);
  EXPECT_EQ(failed.evaluations, controls.max_evaluations);
  EXPECT_EQ(failed.status, model::Status::quadrature_limit);
  EXPECT_THROW((void)model::two_spinon_continuum(R{0}, R{0}), std::invalid_argument);
  EXPECT_THROW((void)model::two_spinon_continuum(R{4}, R{4}), std::invalid_argument);
  EXPECT_THROW((void)model::two_spinon_continuum(R{4}, uni20::numeric_limits<R>::quiet_NaN()), std::invalid_argument);
}
} // namespace
