// SPDX-License-Identifier: GPL-3.0-or-later
#include "sine_gordon_vacuum_references.hpp"
#include "test_support.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/sine_gordon_vacuum.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using C = uni20::complex<R>;
using P = uni20::Precision;
namespace model = bethe::sine_gordon;
namespace math = uni20::math;
using State = model::VacuumState<R>;
void precision(State const& s, P p)
{
  for (auto const& x : {s.mass, s.length, s.coupling, s.scaled_length, s.contour_shift, s.verification_contour_shift,
                        s.cutoff, s.nonlinear_residual, s.kernel_error, s.mesh_error, s.cutoff_error, s.contour_error})
  {
    EXPECT_TRUE(x.initialized());
    EXPECT_EQ(uni20::precision_of(x), p);
  }
  for (auto const* x : {&s.scaling_function, &s.casimir_energy, &s.effective_central_charge})
  {
    EXPECT_EQ(x->has_value(), s.converged);
    if (*x) EXPECT_EQ(uni20::precision_of(**x), p);
  }
  EXPECT_EQ(s.converged, s.status == model::VacuumStatus::converged);
}
auto fields(State const& s)
{
  return std::tie(s.mass, s.length, s.coupling, s.scaled_length, s.contour_shift, s.verification_contour_shift,
                  s.cutoff, s.nonlinear_residual, s.kernel_error, s.mesh_error, s.cutoff_error, s.contour_error,
                  s.scaling_function, s.casimir_energy, s.effective_central_charge, s.intervals, s.iterations,
                  s.kernel_evaluations, s.cutoffs, s.converged, s.status);
}
TEST(SineGordonVacuumMpfr, IndependentDiracIntegralGainsAccuracyBeyondBinary128)
{
  auto hp = P::bits(768);
  R previous(1, hp);
  auto ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256, 333})
  {
    SCOPED_TRACE(bits);
    auto p = P::bits(bits);
    model::VacuumOptions<R> options(p);
    options.max_intervals = 8192;
    R worst(0, hp);
    for (auto const& reference : sg_vacuum_reference::free_dirac)
    {
      SCOPED_TRACE(reference.u);
      auto s = model::vacuum_energy(R(1, p), R(reference.u, p), R(1, p), options);
      ASSERT_TRUE(s.converged) << int(s.status) << " mesh " << s.intervals;
      precision(s, p);
      worst = std::max(worst, math::abs(s.scaling_function->at(hp) - R(reference.y, hp)));
      EXPECT_EQ(s.iterations, 0u);
      EXPECT_EQ(s.kernel_evaluations, 0u);
      EXPECT_LE(s.mesh_error, options.tolerance / 4);
      EXPECT_LE(s.cutoff_error, options.tolerance / 4);
      EXPECT_LE(s.contour_error, options.tolerance / 2);
      EXPECT_LT(math::abs(*s.scaling_function - s.length * *s.casimir_energy), options.tolerance);
    }
    EXPECT_LT(worst, 4 * options.tolerance.at(hp));
    EXPECT_LT(worst, previous * R("1e-15", hp));
    if (bits >= 160) EXPECT_LT(worst, R("1e-40", hp));
    if (bits >= 256) EXPECT_LT(worst, R("1e-68", hp));
    previous = worst;
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
TEST(SineGordonVacuumMpfr, FourierTableAgainstClosedTransforms)
{
  auto hp = P::bits(768);
  R previous(1, hp);
  for (int bits : {80, 160, 256, 333})
  {
    SCOPED_TRACE(bits);
    auto p = P::bits(bits);
    model::KernelOptions<R> options(p);
    options.max_evaluations = 1000000;
    R zero(0, p), eta("0.25", p), h("0.5", p), worst(0, hp);
    R pi = uni20::pi<R>.at(hp);
    for (char const* coupling : {"0.5", "2"})
    {
      std::size_t evaluations = 0;
      auto table = model::detail::kernel_table(R(coupling, p), eta, h, 4, options, evaluations);
      ASSERT_TRUE(table.converged) << coupling << " evals " << evaluations;
      EXPECT_LE(evaluations, options.max_evaluations);
      EXPECT_EQ(uni20::precision_of(table.error), p);
      for (std::size_t j = 0; j < 5; ++j)
      {
        for (bool shifted : {false, true})
        {
          C z(R(j, hp) / 2, R(shifted ? 1 : 0, hp) / 2);
          C expected;
          if (coupling == std::string("0.5"))
            expected = -1 / (2 * pi * math::cosh(z));
          else
            expected = j == 0 && !shifted ? C(1 / (2 * pi * pi), R(0, hp)) : z / (2 * pi * pi * math::sinh(z));
          C actual = shifted ? table.shifted[j] : C(table.real[j], zero);
          EXPECT_EQ(uni20::precision_of(actual), p);
          worst = std::max(worst, math::abs(actual.at(hp) - expected));
        }
      }
    }
    EXPECT_LT(worst, 8 * options.tolerance.at(hp));
    EXPECT_LT(worst, previous * R("1e-15", hp));
    if (bits >= 256) EXPECT_LT(worst, R("1e-70", hp));
    previous = worst;
  }
}
TEST(SineGordonVacuumMpfr, InteractingStatesAndScaling)
{
  auto p = P::bits(160);
  R one(1, p), coupling(2, p);
  model::VacuumOptions<R> c(p);
  c.tolerance = R("1e-9", p);
  for (char const* g : {"0.5", "2"})
  {
    auto state = model::vacuum_energy(one, one, R(g, p), c);
    ASSERT_TRUE(state.converged) << int(state.status);
    precision(state, p);
    auto expected = g == std::string("2") ? "-0.33346196574220655" : "-0.3860018985655546";
    EXPECT_LT(math::abs(*state.scaling_function - R(expected, p)), c.tolerance);
    auto scaled = model::vacuum_energy(R(4, p), R("0.25", p), R(g, p), c);
    ASSERT_TRUE(scaled.converged);
    EXPECT_EQ(*state.scaling_function, *scaled.scaling_function);
    EXPECT_EQ(4 * *state.casimir_energy, *scaled.casimir_energy);
  }
}
TEST(SineGordonVacuumMpfr, IndependentBudgetFailures)
{
  auto p = P::bits(160);
  R one(1, p), coupling(2, p);
  model::VacuumOptions<R> c(p);
  c.tolerance = R("1e-9", p);
  for (auto status : {model::VacuumStatus::kernel_limit, model::VacuumStatus::iteration_limit,
                      model::VacuumStatus::mesh_limit, model::VacuumStatus::cutoff_limit})
  {
    auto o = c;
    if (status == model::VacuumStatus::kernel_limit) o.kernel.max_evaluations = 0;
    if (status == model::VacuumStatus::iteration_limit) o.max_iterations = 0;
    if (status == model::VacuumStatus::mesh_limit) o.max_intervals = o.initial_intervals;
    if (status == model::VacuumStatus::cutoff_limit) o.max_cutoffs = 0;
    auto failed = model::vacuum_energy(one, one, coupling, o);
    precision(failed, p);
    EXPECT_EQ(failed.status, status);
  }
  // More bits do not waive the two successive mesh comparisons.
  model::VacuumOptions<R> fine(P::bits(333));
  auto limited = model::vacuum_energy(R(1, fine.precision()), R(1, fine.precision()), R(1, fine.precision()), fine);
  EXPECT_EQ(limited.status, model::VacuumStatus::mesh_limit);
  precision(limited, fine.precision());
}
TEST(SineGordonVacuumMpfr, RejectsUnsetExactMixedAndInvalidControls)
{
  auto p = P::bits(160);
  R one(1, p);
  model::VacuumOptions<R> c(p);
  c.max_cutoffs = 0;
  EXPECT_THROW((void)model::VacuumOptions<R>(P::exact()), std::invalid_argument);
  for (R const& bad : {R(uni20::uninitialized), R(1), R(1, P::bits(80)), R("nan", p)})
  {
    EXPECT_THROW(model::vacuum_energy(bad, one, one, c), std::invalid_argument);
    EXPECT_THROW(model::vacuum_energy(one, bad, one, c), std::invalid_argument);
    EXPECT_THROW(model::vacuum_energy(one, one, bad, c), std::invalid_argument);
    auto o = c;
    o.tolerance = bad;
    EXPECT_THROW(model::vacuum_energy(one, one, one, o), std::invalid_argument);
    o = c;
    o.contour_shift = bad;
    EXPECT_THROW(model::vacuum_energy(one, one, one, o), std::invalid_argument);
    o = c;
    o.initial_cutoff = bad;
    EXPECT_THROW(model::vacuum_energy(one, one, one, o), std::invalid_argument);
  }
  c.kernel = model::KernelOptions<R>(P::bits(80));
  EXPECT_THROW(model::vacuum_energy(one, one, one, c), std::invalid_argument);
}
TEST(SineGordonVacuumMpfr, OwningMixedPrecisionQueriesKeepCompleteStates)
{
  struct Query
  {
      R mass, length, coupling;
      model::VacuumOptions<R> controls;
  };
  auto evaluate = [](Query const& q) { return model::vacuum_energy(q.mass, q.length, q.coupling, q.controls); };
  std::vector<Query> queries;
  for (int bits : {80, 160})
  {
    auto p = P::bits(bits);
    model::VacuumOptions<R> c(p);
    c.tolerance = R("1e-6", p);
    queries.push_back({R(1, p), R(3, p), R("0.5", p), c});
    c.max_iterations = 0;
    queries.push_back({R(1, p), R(1, p), R(2, p), c});
  }
  std::vector<State> reference;
  for (auto const& q : queries)
    reference.push_back(evaluate(q));
  auto check = [&](uni20::async::IAsyncScheduler& scheduler) {
    std::vector<State> actual;
    auto final = [&] {
      auto owned = queries;
      auto count = owned.size();
      return bethe::evaluate_observable_grid_async(
          count, [owned = std::move(owned), evaluate](std::size_t i) { return evaluate(owned[i]); },
          [&actual](std::size_t i, State s) {
            EXPECT_EQ(i, actual.size());
            actual.push_back(std::move(s));
          },
          {.scheduler = &scheduler, .max_in_flight = 2});
    }();
    EXPECT_EQ(final.get_wait(scheduler), queries.size());
    ASSERT_EQ(actual.size(), reference.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      EXPECT_TRUE(fields(actual[i]) == fields(reference[i]));
      EXPECT_EQ(actual[i].converged, i % 2 == 0);
      if (actual[i].converged)
      {
        EXPECT_GT(actual[i].iterations, 0u);
        EXPECT_GT(actual[i].kernel_evaluations, 0u);
      }
      precision(actual[i], queries[i].controls.precision());
    }
  };
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  check(deferred);
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    check(scheduler);
  }
}
} // namespace

TEST(SineGordonVacuumMpfr, IndependentInteractingD3BeyondBinary128)
{
  using R = uni20::mpreal;
  using P = uni20::Precision;
  namespace model = bethe::sine_gordon;
  auto hp = P::bits(768);
  R reference(sg_vacuum_reference::d3_at_30, hp), previous(1, hp);
  for (int bits : {80, 160})
  {
    SCOPED_TRACE(bits);
    auto p = P::bits(bits);
    model::VacuumOptions<R> c(p);
    c.kernel.max_evaluations = 1000000;
    c.max_intervals = 8192;
    auto state = model::vacuum_energy(R(1, p), R(30, p), R(2, p), c);
    ASSERT_TRUE(state.converged) << int(state.status) << " mesh " << state.intervals;
    R error = uni20::math::abs(state.scaling_function->at(hp) - reference);
    EXPECT_LT(error, 4 * c.tolerance.at(hp));
    EXPECT_LT(error, previous * R("1e-15", hp));
    if (bits >= 160) EXPECT_LT(error, R("1e-40", hp));
    previous = error;
    EXPECT_GT(state.iterations, 0u);
  }
}
