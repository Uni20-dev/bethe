// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/lieb_liniger_thermal.hpp>
#include <bethe/observable_grid.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using uni20::Precision;
namespace model = bethe::lieb_liniger::thermal;
namespace math = uni20::math;

void check_precision(R const& value, Precision p)
{
  ASSERT_TRUE(value.initialized());
  EXPECT_FALSE(value.is_exact());
  EXPECT_EQ(uni20::precision_of(value), p);
  EXPECT_TRUE(uni20::isfinite(value));
}

void check_state(model::State<R> const& s, Precision p)
{
  for (auto const* x : {&s.interaction, &s.temperature, &s.chemical_potential, &s.cutoff, &s.nonlinear_residual,
                        &s.mesh_error, &s.cutoff_error})
    check_precision(*x, p);
  for (auto const* x : {&s.pressure, &s.density, &s.energy_per_length, &s.entropy_per_length})
  {
    ASSERT_EQ(x->has_value(), s.converged);
    if (*x) check_precision(**x, p);
  }
}
} // namespace

TEST(ThermalMpfr, ControlsAndStableFactorsHaveExplicitPrecision)
{
  auto const ambient = mpfr_get_default_prec();
  EXPECT_THROW((void)model::Options<R>(Precision::exact()), std::invalid_argument);
  EXPECT_THROW((void)model::DensityOptions<R>(Precision::exact()), std::invalid_argument);
  for (int bits : {80, 160, 333})
  {
    auto const p = Precision::bits(bits);
    model::Options<R> inner(p);
    model::DensityOptions<R> outer(p);
    EXPECT_EQ(inner.tolerance, 16384 * uni20::numeric_limits<R>::epsilon(p));
    EXPECT_EQ(outer.tolerance, 65536 * uni20::numeric_limits<R>::epsilon(p));
    check_precision(inner.tolerance, p);
    check_precision(outer.tolerance, p);
    for (R const energy : {R(0, p), R(1, p) / 3, R("1e1000", p)})
    {
      R const temperature = R(5, p) / 7;
      auto const a = bethe::detail::thermal_factors(energy, temperature);
      auto const b = bethe::detail::thermal_factors(-energy, temperature);
      for (auto const* x : {&a.log_weight, &a.filling, &a.entropy, &b.log_weight, &b.filling, &b.entropy})
        check_precision(*x, p);
      EXPECT_LT(math::abs(a.filling + b.filling - 1), 8 * uni20::numeric_limits<R>::epsilon(p));
      EXPECT_EQ(a.entropy, b.entropy);
    }
    auto const at_zero = bethe::detail::thermal_factors(R(0, p), R(2, p));
    EXPECT_EQ(at_zero.filling, R(1, p) / 2);
    EXPECT_EQ(at_zero.entropy, math::log(R(2, p)));
    EXPECT_THROW(bethe::detail::thermal_factors(R(1), R(2, p)), std::invalid_argument);
    EXPECT_THROW(bethe::detail::thermal_factors(R(1, Precision::bits(bits + 1)), R(2, p)), std::invalid_argument);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

#if UNI20_HAS_MPLAPACK_MPREAL
TEST(ThermalMpfr, RejectsMixedExactAndUnsetInputsBeforeSolving)
{
  auto const p = Precision::bits(160), other = Precision::bits(80);
  model::Options<R> options(p);
  EXPECT_THROW(model::equilibrium(R(4), R(1, p), R(-1, p), options), std::invalid_argument);
  EXPECT_THROW(model::equilibrium(R(4, p), R(1, other), R(-1, p), options), std::invalid_argument);
  EXPECT_THROW(model::equilibrium(R(4, p), R(1, p), R{}, options), std::invalid_argument);
  options.initial_cutoff = R(8, other);
  EXPECT_THROW(model::equilibrium(R(4, p), R(1, p), R(-1, p), options), std::invalid_argument);
  options.initial_cutoff.reset();
  options.tolerance = R("1e-20", other);
  EXPECT_THROW(model::equilibrium(R(4, p), R(1, p), R(-1, p), options), std::invalid_argument);
  model::DensityOptions<R> canonical(p);
  canonical.max_evaluations = 0;
  canonical.tolerance = R(1) / 100;
  EXPECT_THROW(model::at_density(R(4, p), R(1, p), R(1, p), canonical), std::invalid_argument);
}

TEST(ThermalMpfr, TonksFugacitySeriesShowsAccuracyBeyondBinary128)
{
  // Independent analytic Gaussian integrals of the fugacity series at T=1,
  // mu=-1. This uses no mesh, dressed operator, or finite-c solver.
  auto const reference_precision = Precision::bits(768);
  R const z = math::exp(R(-1, reference_precision));
  R power = z, pressure_sum(0, reference_precision), density_sum(0, reference_precision);
  for (int n = 1; n <= 800; ++n)
  {
    R const term = power / math::sqrt(R(n, reference_precision));
    density_sum += term;
    pressure_sum += term / n;
    power *= -z;
  }
  R const pi = 4 * math::atan(R(1, reference_precision));
  R const prefactor = 1 / (2 * math::sqrt(pi));
  std::array<R, 4> const reference{prefactor * pressure_sum, prefactor * density_sum, prefactor * pressure_sum / 2,
                                   prefactor * (3 * pressure_sum / 2 + density_sum)};
  R previous_error(1, reference_precision);
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto const p = Precision::bits(bits);
    model::Options<R> options(p);
    options.max_nodes = 512;
    // The 256/512 mesh difference is dominated by the coarser rule, not by
    // 256-bit arithmetic. Set its budget separately; the independent reference
    // below still requires an arithmetic-scale error, much smaller than this.
    if (bits == 256) options.tolerance = R("1e-60", p);
    // Finite-c corrections are far below the arithmetic errors being tested.
    auto const state = model::equilibrium(R("1e120", p), R(1, p), R(-1, p), options);
    ASSERT_TRUE(state.converged) << int(state.status) << " nodes=" << state.nodes << " cutoffs=" << state.cutoffs;
    check_state(state, p);
    std::array<R, 4> const actual{*state.pressure, *state.density, *state.energy_per_length, *state.entropy_per_length};
    R error(0, reference_precision);
    for (std::size_t i = 0; i < actual.size(); ++i)
      error = std::max(error, math::abs(uni20::scalar_at_precision<R>(actual[i], reference_precision) - reference[i]) /
                                  reference[i]);
    EXPECT_LT(error, uni20::scalar_at_precision<R>(8 * options.tolerance, reference_precision));
    EXPECT_LT(error, uni20::scalar_at_precision<R>(256 * uni20::numeric_limits<R>::epsilon(p), reference_precision));
    EXPECT_LT(error, previous_error * R("1e-15", reference_precision));
    previous_error = error;
    if (bits == 256) EXPECT_LT(error, R("1e-65", reference_precision));
  }
}

TEST(ThermalMpfr, InteractingScalingAndCanonicalRoundTrip)
{
  auto const p = Precision::bits(160);
  model::Options<R> options(p);
  options.tolerance = R("1e-27", p);
  auto const state = model::equilibrium(R(4, p), R(1, p), R(-1, p), options);
  ASSERT_TRUE(state.converged) << int(state.status);
  check_state(state, p);
  auto const scaled = model::equilibrium(R(8, p), R(4, p), R(-4, p), options);
  ASSERT_TRUE(scaled.converged);
  check_state(scaled, p);
  EXPECT_EQ(*scaled.pressure, 8 * *state.pressure);
  EXPECT_EQ(*scaled.density, 2 * *state.density);
  EXPECT_EQ(*scaled.energy_per_length, 8 * *state.energy_per_length);
  EXPECT_EQ(*scaled.entropy_per_length, 2 * *state.entropy_per_length);
  EXPECT_LT(math::abs(*state.pressure + *state.energy_per_length + *state.density - *state.entropy_per_length),
            8 * options.tolerance);
  model::DensityOptions<R> canonical(p);
  canonical.equilibrium = options;
  canonical.tolerance = R("1e-24", p);
  auto const recovered = model::at_density(R(4, p), R(1, p), *state.density, canonical);
  ASSERT_TRUE(recovered.converged) << int(recovered.status);
  ASSERT_TRUE(recovered.state);
  check_state(*recovered.state, p);
  check_precision(recovered.density_error, p);
  check_precision(recovered.requested_density, p);
  EXPECT_LT(math::abs(recovered.state->chemical_potential + 1), 32 * canonical.tolerance);
  EXPECT_LT(recovered.density_error, canonical.tolerance / 2);
}

TEST(ThermalMpfr, FailureResultsRetainPrecisionAndOmitObservables)
{
  auto const p = Precision::bits(160);
  for (int failure : {0, 1, 2})
  {
    model::Options<R> options(p);
    if (failure == 0) options.max_cutoffs = 0;
    if (failure == 1) options.max_iterations = 0;
    if (failure == 2) options.initial_nodes = options.max_nodes = 4;
    auto const state = model::equilibrium(R(4, p), R(1, p), R(-1, p), options);
    EXPECT_FALSE(state.converged);
    check_state(state, p);
    EXPECT_EQ(state.status, failure == 0   ? model::Status::cutoff_limit
                            : failure == 1 ? model::Status::iteration_limit
                                           : model::Status::mesh_limit);
  }
  model::DensityOptions<R> options(p);
  options.max_evaluations = 0;
  auto const state = model::at_density(R(4, p), R(1, p), R(1, p), options);
  EXPECT_FALSE(state.converged);
  EXPECT_FALSE(state.state);
  EXPECT_EQ(state.status, model::Status::density_limit);
  check_precision(state.requested_density, p);
  check_precision(state.density_error, p);
}

TEST(ThermalMpfr, PositiveObservableWithUnresolvableAbsoluteToleranceFails)
{
  auto const p = Precision::bits(128);
  model::Options<R> options(p);
  options.tolerance = R("1e-12", p);
  auto const reference = model::equilibrium(R(4, p), R(1, p), R(-1, p), options);
  ASSERT_TRUE(reference.converged);
  // Keep the physical pressure positive, but below the range at which its
  // requested absolute tolerance is representable. Do not alter MPFR's range.
  auto const exponent = 2 * ((static_cast<std::int64_t>(mpfr_get_emin()) + 32) / 3);
  R const temperature = math::scalbn(R(1, p), exponent);
  R const momentum_scale = math::sqrt(temperature);
  R const physical_pressure = *reference.pressure * momentum_scale * temperature;
  ASSERT_GT(physical_pressure, R(0, p));
  ASSERT_EQ(physical_pressure * options.tolerance, R(0, p));
  auto const state = model::equilibrium(4 * momentum_scale, temperature, -temperature, options);
  EXPECT_FALSE(state.converged);
  EXPECT_EQ(state.status, model::Status::precision_limit);
  EXPECT_FALSE(state.pressure);
}

TEST(ThermalMpfr, OwningMixedPrecisionGridMatchesCompleteStates)
{
  struct Request
  {
      R c, temperature, mu;
      model::Options<R> options;
  };
  auto const ambient = mpfr_get_default_prec();
  std::vector<Request> requests;
  for (int bits : {80, 160})
    for (bool fail : {false, true})
    {
      auto const p = Precision::bits(bits);
      model::Options<R> options(p);
      options.tolerance = R("1e-12", p);
      if (fail) options.max_iterations = 0;
      requests.push_back({R(4, p), R(2, p) / 3, R(-1, p), options});
    }
  std::vector<model::State<R>> expected;
  for (auto const& r : requests)
    expected.push_back(model::equilibrium(r.c, r.temperature, r.mu, r.options));
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<model::State<R>> actual;
    auto done = bethe::evaluate_observable_grid_async(
        requests.size(),
        [owned = requests](std::size_t i) {
          auto const& r = owned[i];
          return model::equilibrium(r.c, r.temperature, r.mu, r.options);
        },
        [&actual](std::size_t, model::State<R> state) { actual.push_back(std::move(state)); },
        {.scheduler = &scheduler, .max_in_flight = 2});
    (void)done.get_wait(scheduler);
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      auto const& a = actual[i];
      auto const& b = expected[i];
      check_state(a, requests[i].options.precision());
      EXPECT_EQ(a.interaction, b.interaction);
      EXPECT_EQ(a.temperature, b.temperature);
      EXPECT_EQ(a.chemical_potential, b.chemical_potential);
      EXPECT_EQ(a.pressure, b.pressure);
      EXPECT_EQ(a.density, b.density);
      EXPECT_EQ(a.energy_per_length, b.energy_per_length);
      EXPECT_EQ(a.entropy_per_length, b.entropy_per_length);
      EXPECT_EQ(a.cutoff, b.cutoff);
      EXPECT_EQ(a.nonlinear_residual, b.nonlinear_residual);
      EXPECT_EQ(a.mesh_error, b.mesh_error);
      EXPECT_EQ(a.cutoff_error, b.cutoff_error);
      EXPECT_EQ(a.nodes, b.nodes);
      EXPECT_EQ(a.iterations, b.iterations);
      EXPECT_EQ(a.cutoffs, b.cutoffs);
      EXPECT_EQ(a.converged, b.converged);
      EXPECT_EQ(a.status, b.status);
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
#endif
