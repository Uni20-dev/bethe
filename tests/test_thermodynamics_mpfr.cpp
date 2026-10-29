// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <array>
#include <bethe/lieb_liniger_thermo.hpp>
#include <bethe/prepared_observable_grid.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using uni20::Precision;
namespace model = bethe::lieb_liniger::thermo;
namespace math = uni20::math;
using namespace uni20::async;

void finite_precision(R const& value, Precision p)
{
  ASSERT_TRUE(value.initialized());
  EXPECT_FALSE(value.is_exact());
  EXPECT_EQ(value.precision(), p);
  EXPECT_TRUE(uni20::isfinite(value));
}

void check_background(model::Background<R> const& result, Precision p)
{
  for (auto const* value : {&result.interaction, &result.density, &result.mesh_error, &result.density_error})
    finite_precision(*value, p);
  if (result.converged)
  {
    ASSERT_TRUE(result.energy_per_length && result.chemical_potential && result.fermi_rapidity);
    finite_precision(*result.energy_per_length, p);
    finite_precision(*result.chemical_potential, p);
    finite_precision(*result.fermi_rapidity, p);
  }
  else
  {
    EXPECT_FALSE(result.energy_per_length);
    EXPECT_FALSE(result.chemical_potential);
    EXPECT_FALSE(result.fermi_rapidity);
  }
}
} // namespace

TEST(ThermodynamicsMpfr, QuadraturePolynomialMomentsAtIndependentPrecisions)
{
  auto const ambient = mpfr_get_default_prec();
  EXPECT_THROW(bethe::detail::gauss_legendre<R>(16, Precision::exact()), std::invalid_argument);
  for (int bits : {80, 160, 333})
  {
    auto const p = Precision::bits(bits);
    auto const rule = bethe::detail::gauss_legendre<R>(16, p);
    for (std::size_t j = 0; j < rule.x.size(); ++j)
    {
      finite_precision(rule.x[j], p);
      finite_precision(rule.w[j], p);
    }
    for (int power = 0; power < 32; ++power)
    {
      bethe::detail::CompensatedSum<R> sum(R(0, p));
      for (std::size_t j = 0; j < rule.x.size(); ++j)
      {
        R monomial(1, p);
        for (int k = 0; k < power; ++k)
          monomial *= rule.x[j];
        sum.add(rule.w[j] * monomial);
      }
      R const exact = power % 2 ? R(0, p) : R(2, p) / (power + 1);
      EXPECT_LT(math::abs(sum.value() - exact), 256 * uni20::numeric_limits<R>::epsilon(p));
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(ThermodynamicsMpfr, OptionsRequireFinitePrecision)
{
  auto const p = Precision::bits(160);
  EXPECT_THROW((void)model::Options<R>(Precision::exact()), std::invalid_argument);
  model::Options<R> options(p);
  finite_precision(options.tolerance, p);
  EXPECT_EQ(options.tolerance, 4096 * uni20::numeric_limits<R>::epsilon(p));
}

TEST(ThermodynamicsMpfr, RequiresMatchingPrecisionBeforeNumericalWork)
{
  auto const p = Precision::bits(160);
  model::Options<R> options(p);
  EXPECT_THROW(model::ground_state(R(4), R(1, p), options), std::invalid_argument);
  EXPECT_THROW(model::ground_state(R(4, p), R(1, Precision::bits(80)), options), std::invalid_argument);
  options.tolerance = R("1e-20", Precision::bits(80));
  EXPECT_THROW(model::ground_state(R(4, p), R(1, p), options), std::invalid_argument);
}

TEST(ThermodynamicsMpfr, SquareSolvePreservesSmallCorrections)
{
  auto const p = Precision::bits(333);
  R const one(1, p), gap("1e-50", p);
  std::vector<R> const a{one, one, one, one + gap};
  std::vector<R> rhs{2 * one, 2 * one + gap};
  ASSERT_TRUE(bethe::detail::newton_step(a, rhs, p));
  for (auto const& x : rhs)
  {
    finite_precision(x, p);
    EXPECT_EQ(x, one);
  }
  EXPECT_EQ(a.back(), one + gap);
}

TEST(ThermodynamicsMpfr, GroundEnergyGainsAccuracyAndPreservesDensityScaling)
{
  auto const p = Precision::bits(333);
  auto const reference = model::ground_state(R(4, p), R(1, p), model::Options<R>(p));
  ASSERT_TRUE(reference.converged) << int(reference.status) << " nodes=" << reference.nodes;
  check_background(reference, p);
  R previous_error(1, p);
  for (int bits : {128, 256})
  {
    auto const q = Precision::bits(bits);
    auto const state = model::ground_state(R(4, q), R(1, q), model::Options<R>(q));
    ASSERT_TRUE(state.converged) << int(state.status);
    check_background(state, q);
    R const error = math::abs(state.energy_per_length->at(p) - *reference.energy_per_length);
    EXPECT_GT(error, R(0, p));
    EXPECT_LT(error, previous_error * R("1e-20", p));
    EXPECT_LT(error, 8192 * uni20::numeric_limits<R>::epsilon(q).at(p));
    previous_error = error;
  }
  auto const scaled = model::ground_state(R(8, p), R(2, p), model::Options<R>(p));
  ASSERT_TRUE(scaled.converged);
  EXPECT_EQ(*scaled.energy_per_length, 8 * *reference.energy_per_length);
  EXPECT_EQ(*scaled.chemical_potential, 4 * *reference.chemical_potential);
  EXPECT_EQ(*scaled.fermi_rapidity, 2 * *reference.fermi_rapidity);
}

TEST(ThermodynamicsMpfr, ResolvesPhysicalCorrectionBeyondBinary128)
{
  auto const p = Precision::bits(333);
  R const gamma("1e50", p), pi = uni20::pi<R>.at(p);
  auto const state = model::ground_state(gamma, R(1, p), model::Options<R>(p));
  ASSERT_TRUE(state.converged) << int(state.status);
  check_background(state, p);
  // e/(pi^2/3) = 1 - 4/gamma + O(gamma^-2), and
  // mu/pi^2 = 1 - 16/(3*gamma) + O(gamma^-2).
  // Binary128 cannot resolve these departures from the Tonks limit.
  EXPECT_LT(math::abs((*state.energy_per_length / (pi * pi / 3) - 1) * gamma + 4), R("1e-43", p));
  EXPECT_LT(math::abs((*state.chemical_potential / (pi * pi) - 1) * gamma + R(16, p) / 3), R("1e-43", p));
}

TEST(ThermodynamicsMpfr, AsyncParameterSweepOwnsOptionsAndPreservesFailures)
{
  TbbScheduler scheduler(4);
  auto const ambient = mpfr_get_default_prec();
  std::vector<model::Background<R>> rows;
  auto pending = [&] {
    model::Options<R> options(Precision::bits(160));
    return bethe::evaluate_observable_grid_async(
        4,
        [options](std::size_t i) {
          auto local = options;
          if (i == 2) local.max_iterations = 0;
          return model::ground_state(R(8 + i, local.precision()), R(1, local.precision()), local);
        },
        [&](std::size_t i, auto row) {
          EXPECT_EQ(i, rows.size());
          rows.push_back(std::move(row));
        },
        {.scheduler = &scheduler, .max_in_flight = 2});
  }(); // options and the submitting callable are gone
  EXPECT_EQ(pending.get_wait(scheduler), 4u);
  ASSERT_EQ(rows.size(), 4u);
  for (std::size_t i = 0; i < rows.size(); ++i)
  {
    auto const p = Precision::bits(160);
    check_background(rows[i], p);
    EXPECT_EQ(rows[i].converged, i != 2);
    if (i == 2)
      EXPECT_EQ(rows[i].status, model::Status::density_limit);
    else
    {
      auto const serial = model::ground_state(R(8 + i, p), R(1, p), model::Options<R>(p));
      EXPECT_EQ(rows[i].energy_per_length, serial.energy_per_length);
      EXPECT_EQ(rows[i].chemical_potential, serial.chemical_potential);
      EXPECT_EQ(rows[i].fermi_rapidity, serial.fermi_rapidity);
      EXPECT_EQ(rows[i].mesh_error, serial.mesh_error);
      EXPECT_EQ(rows[i].density_error, serial.density_error);
      EXPECT_EQ(rows[i].iterations, serial.iterations);
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(ThermodynamicsMpfr, SharedBackgroundsRetainIndependentPrecisionsAndTinyGaps)
{
  TbbScheduler scheduler(4);
  auto const ambient = mpfr_get_default_prec();
  std::array<std::vector<model::Point<R>>, 2> curves;
  std::vector<Async<std::size_t>> pending;
  for (std::size_t k = 0; k < curves.size(); ++k)
  {
    auto const p = Precision::bits(k == 0 ? 128 : 256);
    auto background = bethe::prepare_background_async(
        [p] { return model::Solver<R>(R(8, p), R(1, p), model::Options<R>(p)); }, {&scheduler});
    pending.push_back(bethe::evaluate_observable_grid_async(
        background.read(), 6,
        [p](auto const& solver, std::size_t i) {
          auto const branch = i < 3 ? model::Branch::type_i : model::Branch::type_ii;
          R const momentum = i % 3 == 0 ? R(0, p) : i % 3 == 1 ? R("1e-50", p) : R(1, p);
          return solver.at_momentum(branch, momentum);
        },
        [&, k](auto, auto point) { curves[k].push_back(std::move(point)); }, {&scheduler, 128, 2}));
  } // both grids submitted, original background handles destroyed
  for (auto& result : pending)
    EXPECT_EQ(result.get_wait(scheduler), 6u);
  for (std::size_t k = 0; k < curves.size(); ++k)
  {
    auto const p = Precision::bits(k == 0 ? 128 : 256);
    model::Solver<R> serial(R(8, p), R(1, p), model::Options<R>(p));
    ASSERT_TRUE(serial.background().converged);
    for (auto const& point : curves[k])
    {
      ASSERT_TRUE(point.converged) << int(point.status);
      auto const expected = serial.at_momentum(point.branch, point.momentum);
      EXPECT_EQ(point.energy, expected.energy);
      EXPECT_EQ(point.rapidity, expected.rapidity);
      EXPECT_EQ(point.edge_distance, expected.edge_distance);
      EXPECT_EQ(point.energy_error, expected.energy_error);
      EXPECT_EQ(point.momentum_error, expected.momentum_error);
      EXPECT_EQ(point.iterations, expected.iterations);
      for (auto const* value : {&point.momentum, &point.energy_error, &point.momentum_error, &*point.energy,
                                &*point.rapidity, &*point.edge_distance})
        finite_precision(*value, p);
      if (point.momentum > R(0, p)) EXPECT_GT(*point.energy, R(0, p));
    }
    auto const end = serial.momentum_range(model::Branch::type_ii).second;
    auto const endpoint = serial.at_momentum(model::Branch::type_ii, end);
    ASSERT_TRUE(endpoint.converged);
    EXPECT_EQ(*endpoint.energy, R(0, p));
    finite_precision(*endpoint.energy, p);
    EXPECT_THROW(serial.at_momentum(model::Branch::type_i, R(1)), std::invalid_argument);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(ThermodynamicsMpfr, FailedBackgroundsAndInversionsRetainDiagnostics)
{
  auto const p = Precision::bits(160);
  for (bool mesh_limit : {false, true})
  {
    model::Options<R> options(p);
    if (mesh_limit)
      options.max_nodes = options.initial_nodes;
    else
      options.max_iterations = 0;
    model::Solver<R> solver(R(4, p), R(1, p), options);
    check_background(solver.background(), p);
    EXPECT_EQ(solver.background().status, mesh_limit ? model::Status::mesh_limit : model::Status::density_limit);
    auto const point = solver.at_momentum(model::Branch::type_i, R(0, p));
    EXPECT_FALSE(point.converged);
    EXPECT_FALSE(point.energy);
    finite_precision(point.energy_error, p);
    finite_precision(point.momentum_error, p);
  }
  model::Options<R> options(p);
  options.max_momentum_iterations = 0;
  model::Solver<R> solver(R(4, p), R(1, p), options);
  auto const point = solver.at_momentum(model::Branch::type_i, R(1, p));
  EXPECT_FALSE(point.energy);
  EXPECT_EQ(point.status, model::Status::momentum_limit);
}
