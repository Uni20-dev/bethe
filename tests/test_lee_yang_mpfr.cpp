// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/lee_yang_excited.hpp>
#include <bethe/observable_grid.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using uni20::Precision;
namespace model = bethe::lee_yang;
namespace math = uni20::math;

template <class State> void check(State const& s, Precision p)
{
  for (R const* v : {&s.mass, &s.length, &s.scaled_length, &s.cutoff, &s.nonlinear_residual, &s.nonlinear_error,
                     &s.mesh_error, &s.cutoff_error, &s.direct_tail_bound})
  {
    EXPECT_TRUE(v->initialized());
    EXPECT_EQ(uni20::precision_of(*v), p);
  }
  EXPECT_EQ(s.converged, s.status == model::Status::converged);
  EXPECT_EQ(s.converged, s.scaling_function.has_value());
  EXPECT_EQ(s.converged, s.casimir_energy.has_value());
  for (auto const* v : {&s.scaling_function, &s.casimir_energy})
    if (*v) EXPECT_EQ(uni20::precision_of(**v), p);
}
template <class State> auto diagnostics(State const& s)
{
  return std::tie(s.mass, s.length, s.scaled_length, s.cutoff, s.nonlinear_residual, s.nonlinear_error, s.mesh_error,
                  s.cutoff_error, s.direct_tail_bound, s.intervals, s.iterations, s.cutoffs, s.kernel_products,
                  s.converged, s.status, s.scaling_function, s.casimir_energy);
}
void same(model::State<R> const& a, model::State<R> const& b)
{
  EXPECT_TRUE(diagnostics(a) == diagnostics(b));
  EXPECT_EQ(a.effective_central_charge, b.effective_central_charge);
}
void same(model::OneParticleState<R> const& a, model::OneParticleState<R> const& b)
{
  EXPECT_TRUE(diagnostics(a) == diagnostics(b));
  EXPECT_TRUE(std::tie(a.beta, a.pole_displacement, a.quantization_residual, a.source_error, a.root_iterations) ==
              std::tie(b.beta, b.pole_displacement, b.quantization_residual, b.source_error, b.root_iterations));
}

// Independent composite Gauss-Legendre oracle, direct complex S-matrix source,
// mpmath 1.2.1. 96/112 decimal digits, orders 64/80, and cutoffs 4/5 agree
// beyond 80 decimal digits. Production uses a trapezoid and reduced real sources.
char const* const reference[] = {
    "-0."
    "00000000000020700315803970804501898623814079162172476268090990402894985548605163866614558976260681450422940334135"
    "3",
    "30.0000000004661235303982110085257024618571758947677308341139611899836528067000078329391699015326242915",
    "0."
    "0000000000090210897973829821340583111561152188730833017778526208349825954359778154532008498577111446684736916651"
    "2"};

void check_independent_levels(std::initializer_list<int> precisions)
{
  auto const oracle = Precision::bits(768);
  R previous(1, oracle);
  auto const ambient = mpfr_get_default_prec();
  for (int bits : precisions)
  {
    SCOPED_TRACE(bits);
    auto p = Precision::bits(bits);
    model::OneParticleOptions<R> controls(p);
    auto const vacuum = model::ground_state(R(1, p), R(30, p), static_cast<model::Options<R> const&>(controls));
    auto const excited = model::one_particle(R(1, p), R(30, p), controls);
    ASSERT_TRUE(vacuum.converged) << int(vacuum.status);
    ASSERT_TRUE(excited.converged) << int(excited.status) << ", " << excited.kernel_products;
    check(vacuum, p);
    check(excited, p);
    R error(0, oracle);
    std::array<R, 3> const values{*vacuum.scaling_function, *excited.scaling_function, *excited.pole_displacement};
    for (std::size_t i = 0; i < values.size(); ++i)
      error = std::max(error, math::abs(values[i].at(oracle) - R(reference[i], oracle)));
    EXPECT_LT(error, 16 * controls.tolerance.at(oracle));
    EXPECT_LT(error, previous * R("1e-12", oracle));
    if (bits >= 160) EXPECT_LT(error, R("1e-40", oracle));
    if (bits == 256) EXPECT_LT(error, R("1e-65", oracle));
    previous = error;
    EXPECT_LE(vacuum.kernel_products, controls.max_kernel_products);
    EXPECT_LE(excited.kernel_products, controls.max_kernel_products);
    EXPECT_EQ(uni20::precision_of(*excited.beta), p);
    EXPECT_EQ(uni20::precision_of(excited.quantization_residual), p);
    EXPECT_EQ(uni20::precision_of(excited.source_error), p);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(LeeYangMpfr, IndependentLevelsGainAccuracyBeyondBinary128)
{
  check_independent_levels({80, 160});
}

TEST(LeeYangMpfr, IndependentLevelsAt256Bits)
{
  // The 1024-interval nonlinear refinement is retained in long-running CI.
  check_independent_levels({160, 256});
}

TEST(LeeYangMpfr, WorkLimitsAndRejectedInputsDoNotLosePrecision)
{
  auto const p = Precision::bits(160);
  R const mass(1, p), length(30, p);
  model::OneParticleOptions<R> controls(p);
  controls.max_kernel_products = 0;
  auto const vacuum = model::ground_state(mass, length, static_cast<model::Options<R> const&>(controls));
  auto const excited = model::one_particle(mass, length, controls);
  check(vacuum, p);
  check(excited, p);
  EXPECT_EQ(vacuum.status, model::Status::work_limit);
  EXPECT_EQ(excited.status, model::Status::work_limit);
  EXPECT_EQ(vacuum.kernel_products, 0u);
  EXPECT_EQ(excited.kernel_products, 0u);
  EXPECT_EQ(uni20::precision_of(excited.quantization_residual), p);
  EXPECT_EQ(uni20::precision_of(excited.source_error), p);
  controls.max_cutoffs = 0;
  auto const no_cutoffs = model::one_particle(mass, length, controls);
  check(no_cutoffs, p);
  EXPECT_EQ(no_cutoffs.status, model::Status::cutoff_limit);
  for (R const& invalid : {R{}, R(1), R(1, Precision::bits(80)), R("nan", p)})
    EXPECT_THROW(model::ground_state(invalid, length, static_cast<model::Options<R> const&>(controls)),
                 std::invalid_argument);
  controls.initial_cutoff = R(2, Precision::bits(80));
  EXPECT_THROW(model::one_particle(mass, length, controls), std::invalid_argument);
  controls.initial_cutoff.reset();
  controls.tolerance = R(1, Precision::bits(80));
  EXPECT_THROW(model::one_particle(mass, length, controls), std::invalid_argument);
  EXPECT_THROW((void)model::Options<R>(Precision::exact()), std::invalid_argument);
  EXPECT_THROW(model::kernel(R(1)), std::invalid_argument);
}

TEST(LeeYangMpfr, OwningMixedPrecisionSweepsPreserveCompleteLevels)
{
  using Output = std::pair<model::State<R>, model::OneParticleState<R>>;
  struct Query
  {
      R mass, length;
      model::OneParticleOptions<R> controls;
  };
  auto submit = [](uni20::async::IAsyncScheduler& scheduler, std::vector<Output>& output) {
    std::vector<Query> requests;
    for (int bits : {80, 160})
    {
      auto p = Precision::bits(bits);
      model::OneParticleOptions<R> c(p);
      c.tolerance = R("1e-10", p);
      requests.push_back({R(1, p), R(30, p), c});
      c.max_kernel_products = 0;
      requests.push_back({R(1, p), R(30, p), c});
    }
    auto const count = requests.size();
    return bethe::evaluate_observable_grid_async(
        count,
        [owned = std::move(requests)](std::size_t i) {
          auto const& q = owned[i];
          return Output{model::ground_state(q.mass, q.length, static_cast<model::Options<R> const&>(q.controls)),
                        model::one_particle(q.mass, q.length, q.controls)};
        },
        [&output](std::size_t i, Output value) {
          EXPECT_EQ(i, output.size());
          output.push_back(std::move(value));
        },
        {.scheduler = &scheduler, .max_in_flight = 2});
  };
  auto const ambient = mpfr_get_default_prec();
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  std::vector<Output> reference;
  auto done = submit(deferred, reference);
  EXPECT_EQ(done.get_wait(deferred), 4u);
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<Output> actual;
    auto final = submit(scheduler, actual);
    EXPECT_EQ(final.get_wait(scheduler), reference.size());
    ASSERT_EQ(actual.size(), reference.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      same(actual[i].first, reference[i].first);
      same(actual[i].second, reference[i].second);
      EXPECT_EQ(actual[i].first.converged, i % 2 == 0);
      EXPECT_EQ(actual[i].second.converged, i % 2 == 0);
      auto p = Precision::bits(i < 2 ? 80 : 160);
      check(actual[i].first, p);
      check(actual[i].second, p);
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
} // namespace
