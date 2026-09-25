// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/hubbard_charge_continuum.hpp>

namespace
{
namespace model = bethe::hubbard::thermo;
template <typename R> class HubbardChargeContinuum : public ::testing::Test {};
TYPED_TEST_SUITE(HubbardChargeContinuum, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(HubbardChargeContinuum, IndependentOracle)
{
  using R = TypeParam;
  R const pi = R{4} * std::atan(R{1});
  // Development oracle is double precision; keep this broad channel test
  // affordable. Separate tests exercise native default tolerances.
  model::ChargeContinuumOptions<R> options;
  options.search.tolerance = R{1e-7};
  options.search.position_tolerance = R{1e-5};
  options.constituent.relative_tolerance = R{1e-10};
  struct Reference
  {
      model::ChargeChannel channel;
      R fraction;
      char const *lower, *upper;
  };
  for (auto ref : {Reference{model::ChargeChannel::spinon_holon, R{0}, "1.8856450004260146", "3.499612327587533"},
                   Reference{model::ChargeChannel::spinon_holon, R{0.25}, "1.4739946796582626", "4.47688442831685"},
                   Reference{model::ChargeChannel::holon_antiholon, R{0.5}, "4.009049756616836", "8.724870549966758"}})
  {
    auto const state =
        model::charge_continuum(ref.channel, R{4}, pi * ref.fraction, model::Convention::symmetric, options);
    ASSERT_TRUE(state.converged) << int(state.search_status) << "," << int(state.constituent_status) << " objectives "
                                 << state.objective_evaluations << " quadratures " << state.quadrature_evaluations;
    EXPECT_REAL_NEAR(*state.lower_energy, uni20::parse_real<R>(ref.lower), R{1e-6});
    EXPECT_REAL_NEAR(*state.upper_energy, uni20::parse_real<R>(ref.upper), R{1e-6});
    EXPECT_LE(state.quadrature_evaluations, options.max_quadrature_evaluations);
    // Independently differentiate the constituent sum at the interior minimum.
    // A coarse mesh winner is not sufficient: the stationary split must resolve.
    auto energy = [&](R p) {
      auto const first =
          ref.channel == model::ChargeChannel::holon_antiholon ? model::Branch::holon : model::Branch::spinon;
      auto const second =
          ref.channel == model::ChargeChannel::spinon_holon ? model::Branch::holon : model::Branch::antiholon;
      auto const a = model::dispersion(first, R{4}, p, model::Convention::symmetric, options.constituent);
      auto const b = model::dispersion(second, R{4}, model::detail::wrap(state.momentum - p),
                                       model::Convention::symmetric, options.constituent);
      EXPECT_TRUE(a.converged && b.converged);
      return a.energy.value_or(R{100}) + b.energy.value_or(R{100});
    };
    R const h = R{0.001}, p = state.lower_momenta[0];
    EXPECT_REAL_NEAR((energy(p + h) - energy(p - h)) / (R{2} * h), R{0}, R{1e-4});
  }
}
TEST(HubbardChargeContinuumCouplings, WeakAndStrongOracle)
{
  model::ChargeContinuumOptions<double> options;
  options.search.tolerance = 1e-7;
  options.search.position_tolerance = 1e-5;
  options.constituent.relative_tolerance = 1e-10;
  double const pi = 4 * std::atan(1.0);
  for (auto ref : {std::array<double, 3>{1, 1.3346630672738753, 4.202469304569685},
                   std::array<double, 3>{16, 6.4307959137581605, 9.728649273861945}})
  {
    auto const state = model::charge_continuum(model::ChargeChannel::spinon_holon, ref[0], pi / 4,
                                               model::Convention::symmetric, options);
    ASSERT_TRUE(state.converged) << int(state.search_status) << "," << int(state.constituent_status);
    EXPECT_NEAR(*state.lower_energy, ref[1], 1e-5);
    EXPECT_NEAR(*state.upper_energy, ref[2], 1e-5);
  }
}

TYPED_TEST(HubbardChargeContinuum, NativeDefaultTolerance)
{
  using R = TypeParam;
  auto const state = model::charge_continuum(model::ChargeChannel::holon_antiholon, R{4}, R{0});
  ASSERT_TRUE(state.converged) << int(state.search_status) << "," << int(state.constituent_status) << " objectives "
                               << state.objective_evaluations << " quadratures " << state.quadrature_evaluations;
  R const gap = uni20::parse_real<R>("0.643363511006452197366294884644872250891658422");
  R const tolerance = model::ChargeContinuumOptions<R>{}.search.tolerance * R{4};
  EXPECT_REAL_NEAR(*state.lower_energy, R{2} * gap, tolerance);
  EXPECT_REAL_NEAR(*state.upper_energy, R{2} * gap + R{8}, tolerance);
  EXPECT_EQ(state.delta_particles, 0);
  EXPECT_EQ(state.spin, uni20::half_int{0});
  EXPECT_LE(state.lower_error, tolerance);
  EXPECT_LE(state.upper_error, tolerance);
}

TYPED_TEST(HubbardChargeContinuum, ShiftAndConvention)
{
  using R = TypeParam;
  R const pi = R{4} * std::atan(R{1});
  model::ChargeContinuumOptions<R> o;
  o.search.tolerance = R{1e-7};
  o.search.position_tolerance = R{1e-5};
  o.constituent.relative_tolerance = R{1e-10};
  auto const holon =
      model::charge_continuum(model::ChargeChannel::spinon_holon, R{4}, pi / R{4}, model::Convention::unshifted, o);
  auto const anti = model::charge_continuum(model::ChargeChannel::spinon_antiholon, R{4}, -R{3} * pi / R{4},
                                            model::Convention::unshifted, o);
  ASSERT_TRUE(holon.converged && anti.converged);
  EXPECT_EQ(holon.delta_particles, -1);
  EXPECT_EQ(anti.delta_particles, 1);
  EXPECT_EQ(holon.spin, uni20::from_twice(1));
  EXPECT_REAL_NEAR(*holon.lower_symmetric_energy, *anti.lower_symmetric_energy, R{1e-6});
  EXPECT_REAL_NEAR(*holon.upper_symmetric_energy, *anti.upper_symmetric_energy, R{1e-6});
  EXPECT_REAL_NEAR(*anti.lower_energy - *holon.lower_energy, R{4}, R{1e-6});
  for (auto const* state : {&holon, &anti})
    for (auto const& momenta : {state->lower_momenta, state->upper_momenta})
      EXPECT_REAL_NEAR(model::detail::wrap(momenta[0] + momenta[1] - state->momentum), R{0},
                       R{64} * uni20::numeric_limits<R>::epsilon());
}

TYPED_TEST(HubbardChargeContinuum, FailureAndValidation)
{
  using R = TypeParam;
  model::ChargeContinuumOptions<R> o;
  o.max_quadrature_evaluations = 0;
  auto const failed =
      model::charge_continuum(model::ChargeChannel::spinon_holon, R{4}, R{1}, model::Convention::symmetric, o);
  EXPECT_FALSE(failed.converged);
  EXPECT_FALSE(failed.lower_energy);
  EXPECT_FALSE(failed.upper_symmetric_energy);
  EXPECT_EQ(failed.constituent_status, model::Status::quadrature_limit);
  EXPECT_EQ(failed.quadrature_evaluations, 0u);
  o.max_quadrature_evaluations = 50;
  auto const partial =
      model::charge_continuum(model::ChargeChannel::spinon_holon, R{4}, R{1}, model::Convention::symmetric, o);
  EXPECT_FALSE(partial.converged);
  EXPECT_FALSE(partial.lower_energy);
  EXPECT_FALSE(partial.upper_energy);
  EXPECT_EQ(partial.quadrature_evaluations, 50u);
  o = {};
  o.search.max_evaluations = 0;
  auto const search_failed =
      model::charge_continuum(model::ChargeChannel::holon_antiholon, R{4}, R{0}, model::Convention::symmetric, o);
  EXPECT_FALSE(search_failed.converged);
  EXPECT_FALSE(search_failed.lower_energy);
  EXPECT_EQ(search_failed.search_status, bethe::detail::ExtremaStatus::evaluation_limit);
  EXPECT_EQ(search_failed.quadrature_evaluations, 0u);
  EXPECT_THROW((void)model::charge_continuum(model::ChargeChannel::spinon_holon, R{0}, R{1}), std::invalid_argument);
  EXPECT_THROW((void)model::charge_continuum(model::ChargeChannel::spinon_holon, R{4}, R{4}), std::invalid_argument);
}
} // namespace
