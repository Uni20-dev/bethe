// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <atomic>
#include <bethe/hubbard_charge_continuum.hpp>
#include <bethe/hubbard_continuum.hpp>
#include <bethe/hubbard_doped.hpp>
#include <bethe/observable_grid.hpp>
#include <bethe/prepared_observable_grid.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using uni20::Precision;
namespace model = bethe::hubbard::thermo;
namespace math = uni20::math;

void check(model::Point<R> const& s, Precision p)
{
  for (auto const* value :
       {&s.interaction, &s.momentum, &s.parameter, &s.energy_offset, &s.energy_error, &s.momentum_error})
  {
    ASSERT_TRUE(value->initialized());
    EXPECT_EQ(uni20::precision_of(*value), p);
  }
  for (auto const* value : {&s.energy, &s.symmetric_energy})
    if (*value)
    {
      EXPECT_EQ(uni20::precision_of(**value), p);
      EXPECT_TRUE(uni20::isfinite(**value));
    }
  EXPECT_EQ(s.converged, s.status == model::Status::converged);
  EXPECT_EQ(s.converged, s.energy.has_value());
}

auto fields(model::Point<R> const& s)
{
  return std::tie(s.branch, s.convention, s.interaction, s.momentum, s.parameter, s.energy_offset, s.energy_error,
                  s.momentum_error, s.energy, s.symmetric_energy, s.delta_particles, s.spin, s.evaluations,
                  s.iterations, s.converged, s.status);
}

// Independent original Fourier-Bessel integrals at U=4, lambda=2, k=1,
// mpmath 1.2.1 at 96 and 128 decimal digits. Both runs agree through all
// 90 digits shown. The weak-coupling minimum uses the positive odd K1 series.
// The production solver instead integrates nonoscillatory convolution forms.
char const* const references[] = {
    "0.180195419367165283255675569344877169275747982857939066888904265913419968700554935899809673",
    "0.147435530616030589362559903539606986354219786134238672485504679036370550127469923860951005",
    "3.62482643594807997383907104208212004914561353596139543345620414422224118406960936131405202",
    "0.0869244707432127794663054310826567415930967806182991896122678377545027768502062684148452409",
    "0.00000000000785641942083116472051349106943824118854939922451588081636261556995739008336698307366703898"};

TEST(HubbardThermoMpfr, IndependentLinesAndWeakGapGainAccuracyBeyondBinary128)
{
  auto const refp = Precision::bits(768);
  R previous(1, refp);
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto const p = Precision::bits(bits);
    model::Options<R> controls(p);
    auto const spin = model::spinon_at_rapidity(R(4, p), R(2, p), model::Convention::symmetric, controls);
    auto const charge =
        model::charge_at_bare_momentum(model::Branch::holon, R(4, p), R(1, p), model::Convention::symmetric, controls);
    R const pi = model::detail::pi(R(0, p));
    auto const weak =
        model::charge_at_bare_momentum(model::Branch::holon, R(1, p) / 4, pi, model::Convention::symmetric, controls);
    for (auto const* s : {&spin, &charge, &weak})
    {
      ASSERT_TRUE(s->converged) << int(s->status) << ", " << s->evaluations;
      check(*s, p);
    }
    std::array<R, 5> const actual{*spin.energy, spin.momentum, *charge.energy, charge.momentum, *weak.energy};
    R error(0, refp);
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      R const expected(references[i], refp);
      error = std::max(error, math::abs(actual[i].at(refp) - expected) / math::abs(expected));
    }
    EXPECT_LT(error, 2048 * uni20::numeric_limits<R>::epsilon(p).at(refp));
    EXPECT_LT(error, previous * R("1e-15", refp));
    if (bits == 256) EXPECT_LT(error, R("1e-70", refp));
    previous = error;
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(HubbardThermoMpfr, EndpointsInversionConventionsAndFailuresKeepPrecision)
{
  auto const p = Precision::bits(160);
  R const zero(0, p), u(4, p), pi = model::detail::pi(zero);
  model::Options<R> controls(p);
  for (R const& momentum : {zero, pi})
  {
    auto const s = model::dispersion(model::Branch::spinon, u, momentum, model::Convention::symmetric, controls);
    ASSERT_TRUE(s.converged);
    check(s, p);
    EXPECT_EQ(*s.energy, zero);
    EXPECT_EQ(s.evaluations, 0u);
  }
  auto const s =
      model::dispersion(model::Branch::spinon, u, R(references[1], p), model::Convention::symmetric, controls);
  ASSERT_TRUE(s.converged);
  check(s, p);
  EXPECT_LT(math::abs(*s.energy - R(references[0], p)), 4096 * uni20::numeric_limits<R>::epsilon(p));
  auto const h =
      model::charge_at_bare_momentum(model::Branch::holon, u, R(1, p), model::Convention::unshifted, controls);
  ASSERT_TRUE(h.converged);
  check(h, p);
  EXPECT_EQ(h.energy_offset, -u / 2);
  EXPECT_EQ(*h.energy, *h.symmetric_energy - u / 2);
  auto const edge = model::two_spinon_continuum(u, zero, controls);
  ASSERT_TRUE(edge.converged);
  EXPECT_EQ(*edge.lower_energy, zero);
  EXPECT_EQ(*edge.upper_energy, zero);
  for (auto const* v : {&edge.lower_energy_quad_error, &edge.upper_energy_quad_error, &edge.lower_momentum_error,
                        &edge.upper_momentum_error})
    EXPECT_EQ(uni20::precision_of(*v), p);
  controls.max_evaluations = 10;
  auto const limited = model::dispersion(model::Branch::holon, u, zero, model::Convention::symmetric, controls);
  check(limited, p);
  EXPECT_FALSE(limited.converged);
  EXPECT_EQ(limited.evaluations, 10u);
  EXPECT_EQ(limited.status, model::Status::quadrature_limit);
  for (R const& invalid : {R{}, R(1), R(1, Precision::bits(80)), R("nan", p)})
    EXPECT_THROW(model::spinon_at_rapidity(u, invalid, model::Convention::symmetric, controls), std::invalid_argument);
  controls.relative_tolerance = R(1, Precision::bits(80));
  EXPECT_THROW(model::dispersion(model::Branch::spinon, u, zero, model::Convention::symmetric, controls),
               std::invalid_argument);
  EXPECT_THROW((void)model::Options<R>(Precision::exact()), std::invalid_argument);
}

TEST(HubbardThermoMpfr, OwningMixedPrecisionQueriesPreserveCompleteResults)
{
  auto submit = [](uni20::async::IAsyncScheduler& scheduler, auto& output) {
    struct Query
    {
        R u, lambda;
        model::Options<R> controls;
    };
    std::vector<Query> requests;
    for (int bits : {80, 160})
    {
      auto const p = Precision::bits(bits);
      model::Options<R> c(p);
      requests.push_back({R(4, p), R(2, p), c});
      c.max_evaluations = 0;
      requests.push_back({R(4, p), R(1, p), c});
    }
    return bethe::evaluate_observable_grid_async(
        requests.size(),
        [owned = requests](std::size_t i) {
          auto const& q = owned[i];
          return model::spinon_at_rapidity(q.u, q.lambda, model::Convention::symmetric, q.controls);
        },
        [&output](std::size_t i, auto s) {
          EXPECT_EQ(i, output.size());
          output.push_back(std::move(s));
        },
        {.scheduler = &scheduler, .max_in_flight = 2});
  };
  auto const ambient = mpfr_get_default_prec();
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  std::vector<model::Point<R>> expected;
  auto done = submit(deferred, expected);
  EXPECT_EQ(done.get_wait(deferred), 4u);
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<model::Point<R>> actual;
    auto final = submit(scheduler, actual);
    EXPECT_EQ(final.get_wait(scheduler), expected.size());
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      EXPECT_TRUE(fields(actual[i]) == fields(expected[i]));
      check(actual[i], Precision::bits(i < 2 ? 80 : 160));
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(HubbardThermoMpfr, ExtremaControlsAndBudgetFailuresAreExplicit)
{
  auto const p = Precision::bits(160);
  R const zero(0, p), one(1, p);
  bethe::detail::ExtremaOptions<R> c(p);
  c.position_tolerance = R("1e-12", p);
  c.tolerance = R("1e-20", p);
  auto const result = bethe::detail::bounded_extrema(
      [zero](R x) { return std::optional{bethe::detail::ObjectiveSample<R>{(x - zero) * (x - zero), zero}}; }, -one,
      one, c);
  ASSERT_TRUE(result.converged);
  EXPECT_EQ(result.minimum->value, zero);
  EXPECT_EQ(result.maximum->value, one);
  EXPECT_EQ(uni20::precision_of(result.minimum->error), p);
  model::ChargeContinuumOptions<R> controls(p);
  controls.max_quadrature_evaluations = 0;
  auto const failed = model::charge_continuum(model::ChargeChannel::spinon_holon, R(4, p), zero,
                                              model::Convention::symmetric, controls);
  EXPECT_FALSE(failed.converged);
  EXPECT_FALSE(failed.lower_energy);
  EXPECT_EQ(failed.quadrature_evaluations, 0u);
  EXPECT_EQ(uni20::precision_of(failed.lower_error), p);
  EXPECT_EQ(uni20::precision_of(failed.energy_offset), p);
}

TEST(HubbardThermoMpfr, PhysicalChargeContinuumRetainsExplicitPrecision)
{
  auto const p = Precision::bits(160);
  model::ChargeContinuumOptions<R> controls(p);
  // Keep this global-search check bounded independently of arithmetic bits.
  // Accuracy growth of the constituent integrals is tested separately above.
  controls.search.tolerance = R("1e-8", p);
  controls.search.position_tolerance = R("1e-6", p);
  controls.constituent.relative_tolerance = R("1e-12", p);
  auto const result = model::charge_continuum(model::ChargeChannel::holon_antiholon, R(4, p), R(0, p),
                                              model::Convention::symmetric, controls);
  ASSERT_TRUE(result.converged) << int(result.search_status) << ", " << int(result.constituent_status);
  // At zero total momentum the charge-doublet edges are twice the gap and
  // twice the gap plus 8; the gap is an independent Fourier-Bessel reference.
  R const twice_gap("1.286727022012904394732589769289744501783316844", p);
  EXPECT_LT(math::abs(*result.lower_energy - twice_gap), R("1e-8", p));
  EXPECT_LT(math::abs(*result.upper_energy - twice_gap - 8), R("1e-8", p));
  for (auto const* v :
       {&result.interaction, &result.momentum, &result.energy_scale, &result.energy_offset, &*result.lower_energy,
        &*result.upper_energy, &*result.lower_symmetric_energy, &*result.upper_symmetric_energy, &result.lower_error,
        &result.upper_error, &result.lower_momentum_error, &result.upper_momentum_error, &result.lower_momenta[0],
        &result.lower_momenta[1], &result.upper_momenta[0], &result.upper_momenta[1]})
    EXPECT_EQ(uni20::precision_of(*v), p);
  EXPECT_LE(result.quadrature_evaluations, controls.max_quadrature_evaluations);
  EXPECT_EQ(result.delta_particles, 0);
}

#if UNI20_HAS_MPLAPACK_MPREAL
auto fields(model::DopedPoint<R> const& s)
{
  return std::tie(s.branch, s.convention, s.reference, s.momentum, s.parameter, s.energy_error, s.momentum_error,
                  s.energy, s.symmetric_energy, s.fermi_energy, s.delta_particles, s.spin, s.iterations, s.converged,
                  s.status);
}

TEST(HubbardThermoMpfr, OwningDopedPreparationsShareMixedPrecisionBackgrounds)
{
  using Outputs = std::array<std::vector<model::DopedPoint<R>>, 6>;
  auto submit = [](uni20::async::IAsyncScheduler& scheduler, Outputs& output, std::atomic<int>& preparations) {
    std::vector<uni20::async::Async<std::size_t>> finals;
    for (std::size_t family = 0; family < 3; ++family)
    {
      auto const p = Precision::bits(family == 1 ? 160 : 80);
      model::DopedOptions<R> controls(p);
      controls.tolerance = R("1e-12", p);
      if (family == 2) controls.max_background_iterations = 0;
      auto background = bethe::prepare_background_async(
          [p, controls, &preparations] {
            ++preparations;
            return model::DopedSolver<R>(R(40, p), R(1, p) / 3, controls);
          },
          {&scheduler});
      for (std::size_t channel = 0; channel < 2; ++channel)
      {
        auto const branch = channel == 0 ? model::DopedBranch::spinon : model::DopedBranch::holon;
        auto const index = 2 * family + channel;
        finals.push_back(bethe::evaluate_observable_grid_async(
            background.read(), 3,
            [branch](model::DopedSolver<R> const& solver, std::size_t i) {
              auto const [lo, hi] = solver.momentum_range(branch);
              R const momentum = i == 2 ? hi : lo + (hi - lo) * R(i) / 2;
              return solver.at_momentum(branch, momentum, model::Convention::unshifted, model::EnergyReference::fermi);
            },
            [&output, index](std::size_t i, auto state) {
              EXPECT_EQ(i, output[index].size());
              output[index].push_back(std::move(state));
            },
            {.scheduler = &scheduler, .max_in_flight = 2}));
      }
      // The preparation handle and local options expire before any final wait.
    }
    return finals;
  };
  auto const ambient = mpfr_get_default_prec();
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  Outputs expected;
  std::atomic<int> prepared{0};
  auto final = submit(deferred, expected, prepared);
  EXPECT_EQ(prepared, 0);
  for (auto& result : final)
    EXPECT_EQ(result.get_wait(deferred), 3u);
  EXPECT_EQ(prepared, 3);
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    Outputs actual;
    prepared = 0;
    auto finals = submit(scheduler, actual, prepared);
    for (auto& result : finals)
      EXPECT_EQ(result.get_wait(scheduler), 3u);
    EXPECT_EQ(prepared, 3);
    for (std::size_t family = 0; family < actual.size(); ++family)
    {
      ASSERT_EQ(actual[family].size(), 3u);
      ASSERT_EQ(expected[family].size(), 3u);
      auto const p = Precision::bits(family / 2 == 1 ? 160 : 80);
      for (std::size_t i = 0; i < 3; ++i)
      {
        auto const& state = actual[family][i];
        EXPECT_TRUE(fields(state) == fields(expected[family][i]));
        EXPECT_EQ(state.converged, family < 4);
        for (auto const* value : {&state.momentum, &state.parameter, &state.energy_error, &state.momentum_error})
          EXPECT_EQ(uni20::precision_of(*value), p);
        if (family >= 4)
        {
          EXPECT_EQ(state.status, model::DopedStatus::density_limit);
          EXPECT_FALSE(state.energy);
        }
        else
          for (auto const* energy : {&state.energy, &state.symmetric_energy, &state.fermi_energy})
          {
            ASSERT_TRUE(*energy);
            EXPECT_EQ(uni20::precision_of(**energy), p);
          }
      }
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(HubbardThermoMpfr, DopedFredholmReferencesAndPhysicalPrecisionGrowth)
{
  // Independent mpmath Gauss-Legendre/digamma equations at U=40, Q=1.
  // 64 nodes/96 digits and 96 nodes/112 digits agree beyond 85 decimal digits;
  // the requested density below is obtained from that background, not fitted.
  auto const refp = Precision::bits(768);
  R const density_ref("0.324321270322984098050102851457799211925925445763299821413665672162111298405245789656083022",
                      refp),
      mu_ref("-1.09260383304129321410274864209167693154960916564478503595799579631728287352309053032998345", refp),
      e0_ref("-0.54608999560441398159217752191704366471571732920728089904214447782014521020325842074560681", refp),
      es_ref("0.0272174980551578681104364045656371182104776147423962347936780945405613311251388154515130903", refp);
  R previous(1, refp);
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto const p = Precision::bits(bits);
    model::DopedOptions<R> controls(p);
    controls.max_nodes = 128;
    model::DopedSolver<R> solver(R(40, p), density_ref.at(p), controls);
    auto const& background = solver.background();
    ASSERT_TRUE(background.converged) << int(background.status) << ", " << background.nodes << ", "
                                      << background.iterations;
    for (auto const* value :
         {&background.interaction, &background.density, &background.mesh_error, &background.density_error})
      EXPECT_EQ(uni20::precision_of(*value), p);
    for (auto const* value : {&background.fermi_rapidity, &background.mu_unshifted, &background.mu_symmetric,
                              &background.energy_per_site_unshifted})
    {
      ASSERT_TRUE(*value);
      EXPECT_EQ(uni20::precision_of(**value), p);
    }
    auto const range = solver.momentum_range(model::DopedBranch::spinon);
    auto const spin = solver.at_momentum(model::DopedBranch::spinon, range.second / 2);
    ASSERT_TRUE(spin.converged) << int(spin.status);
    R const error = std::max({math::abs(background.fermi_rapidity->at(refp) - 1),
                              math::abs(background.mu_unshifted->at(refp) - mu_ref),
                              math::abs(background.energy_per_site_unshifted->at(refp) - e0_ref),
                              math::abs(spin.energy->at(refp) - es_ref)});
    EXPECT_LT(error, 4096 * uni20::numeric_limits<R>::epsilon(p).at(refp));
    EXPECT_LT(error, previous * R("1e-15", refp));
    if (bits == 256) EXPECT_LT(error, R("1e-70", refp));
    previous = error;
    for (auto const* value : {&spin.momentum, &spin.parameter, &spin.energy_error, &spin.momentum_error})
      EXPECT_EQ(uni20::precision_of(*value), p);
    auto const endpoint = solver.at_momentum(model::DopedBranch::spinon, range.first);
    ASSERT_TRUE(endpoint.converged);
    EXPECT_EQ(*endpoint.energy, R(0, p));
    EXPECT_EQ(uni20::precision_of(endpoint.parameter), p);
  }
}
#else
TEST(HubbardThermoMpfr, DopedMissingProviderIsAnExplicitConfigurationError)
{
  auto const p = Precision::bits(80);
  model::DopedOptions<R> controls(p);
  EXPECT_THROW((model::DopedSolver<R>(R(40, p), R(1, p) / 3, controls)), std::runtime_error);
}
#endif

TEST(HubbardThermoMpfr, DopedInputAndEarlyBudgetFailuresRetainPrecision)
{
  auto const p = Precision::bits(160);
  model::DopedOptions<R> controls(p);
  controls.max_background_iterations = 0;
  model::DopedSolver<R> solver(R(4, p), R(1, p) / 2, controls);
  auto const& background = solver.background();
  EXPECT_FALSE(background.converged);
  EXPECT_EQ(background.status, model::DopedStatus::density_limit);
  EXPECT_EQ(uni20::precision_of(background.mesh_error), p);
  auto const point = solver.at_momentum(model::DopedBranch::spinon, R(0, p));
  EXPECT_FALSE(point.energy);
  EXPECT_EQ(point.status, model::DopedStatus::density_limit);
  EXPECT_EQ(uni20::precision_of(point.parameter), p);
  EXPECT_EQ(uni20::precision_of(point.energy_error), p);
  EXPECT_THROW(solver.at_momentum(model::DopedBranch::spinon, R(0, Precision::bits(80))), std::invalid_argument);
  controls.tolerance = R(1, Precision::bits(80));
  EXPECT_THROW((model::DopedSolver<R>(R(4, p), R(1, p) / 2, controls)), std::invalid_argument);
}
} // namespace
