// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <array>
#include <atomic>
#include <bethe/asep.hpp>
#include <bethe/prepared_observable_grid.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using C = uni20::complex<R>;
using uni20::Precision;
namespace model = bethe::asep;
namespace math = uni20::math;

auto fields(model::GapState<R> const& s)
{
  return std::tie(s.sites, s.particles, s.effective_particles, s.wave_index, s.iterations, s.continuation_steps,
                  s.seed_iterations, s.seed_newton_iterations, s.right_rate, s.left_rate, s.reached_ratio,
                  s.residual_norm, s.wave_base, s.scaled_roots, s.eigenvalue, s.gap, s.converged, s.analytic, s.status);
}

void check_precision(model::GapState<R> const& s, Precision p)
{
  for (auto const* value : {&s.right_rate, &s.left_rate, &s.reached_ratio, &s.residual_norm})
    EXPECT_EQ(uni20::precision_of(*value), p);
  EXPECT_EQ(uni20::precision_of(s.wave_base), p);
  for (auto const& root : s.scaled_roots)
    EXPECT_EQ(uni20::precision_of(root), p);
  if (s.eigenvalue) EXPECT_EQ(uni20::precision_of(*s.eigenvalue), p);
  if (s.gap) EXPECT_EQ(uni20::precision_of(*s.gap), p);
}

TEST(AsepMpfr, IndependentMarkovEigenvalueGainsAccuracy)
{
  // Independent ten-state L=5,N=2 Markov generator with right/left rates 1,1/2.
  auto high = Precision::bits(512);
  C const exact{
      R("-1.04039558904641505636142345412607413799326810221158347669691888089548565042267332323361234", high),
      R("0.159168821610180253202429215239399086838149271486788653791649970802963326547367324815530809", high)};
  R previous(1, high);
  auto ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto p = Precision::bits(bits);
    model::Options<R> options(p);
    auto s = model::relaxation_gap(5, 2, R(1, p), R(1, p) / 2, options);
    ASSERT_TRUE(s.converged) << int(s.status) << " steps " << s.continuation_steps;
    check_precision(s, p);
    auto error = math::abs(s.eigenvalue->at(high) - exact);
    EXPECT_LT(error, 8192 * uni20::numeric_limits<R>::epsilon(p).at(high));
    EXPECT_LT(error, previous * R("1e-15", high));
    if (bits == 256) EXPECT_LT(error, R("1e-70", high));
    previous = error;
    auto reflected = model::relaxation_gap(5, 3, R(3, p) / 2, R(3, p), options);
    ASSERT_TRUE(reflected.converged);
    EXPECT_EQ(*reflected.eigenvalue, 3 * *s.eigenvalue);
    model::PreparedGap<R> prepared(5, 2, options);
    EXPECT_TRUE(fields(prepared.at_rates(R(1, p), R(1, p) / 2)) == fields(s));
    auto totally = prepared.at_rates(R(1, p), R(0, p));
    ASSERT_TRUE(totally.converged);
    EXPECT_EQ(totally.eigenvalue, prepared.seed_state().eigenvalue);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(AsepMpfr, TinyBiasRetainsDriftAndAnalyticEndpoints)
{
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto p = Precision::bits(bits);
    model::Options<R> options(p);
    R const one(1, p), zero(0, p), eps = uni20::numeric_limits<R>::epsilon(p), pi = 4 * math::atan(one);
    R q = one - 8 * eps;
    auto near = model::relaxation_gap(5, 2, one, q, options);
    ASSERT_TRUE(near.converged) << int(near.status) << " steps " << near.continuation_steps;
    EXPECT_FALSE(near.analytic);
    check_precision(near, p);
    EXPECT_LT(math::abs(near.eigenvalue->imag() / (one - q) - math::sin(2 * pi / 5) / 3), 32768 * eps);
    C tiny(eps * eps, eps * eps);
    EXPECT_LT(math::abs(bethe::detail::complex_log1p(tiny) / tiny - one), 32 * eps);
    for (std::size_t n : {0u, 1u, 2u, 3u, 5u, 6u})
    {
      auto s = model::relaxation_gap(6, n, one, one, options);
      ASSERT_TRUE(s.converged);
      EXPECT_TRUE(s.analytic);
      check_precision(s, p);
      if (n == 0 || n == 6)
        EXPECT_FALSE(s.gap);
      else
      {
        ASSERT_TRUE(s.gap);
        EXPECT_LT(math::abs(*s.gap - one), 32 * eps);
        EXPECT_EQ(s.eigenvalue->imag(), zero);
      }
    }
  }
}

TEST(AsepMpfr, FailureDiagnosticsAndInputPrecision)
{
  auto p = Precision::bits(160);
  R one(1, p), half = one / 2;
  for (int mode = 0; mode < 4; ++mode)
  {
    model::Options<R> options(p);
    if (mode == 0) options.max_continuation_steps = 0;
    if (mode == 1) options.max_iterations = 0;
    if (mode == 2) options.seed_options.max_seed_iterations = 0;
    if (mode == 3) options.seed_options.max_iterations = 0;
    auto s = model::relaxation_gap(8, 3, one, half, options);
    EXPECT_FALSE(s.converged);
    EXPECT_FALSE(s.eigenvalue);
    EXPECT_FALSE(s.gap);
    EXPECT_EQ(s.status, mode == 0   ? model::Status::continuation_limit
                        : mode == 1 ? model::Status::iteration_limit
                                    : model::Status::seed_limit);
    check_precision(s, p);
    model::PreparedGap<R> prepared(8, 3, options);
    EXPECT_TRUE(fields(prepared.at_rates(one, half)) == fields(s));
    auto analytic = prepared.at_rates(one, one);
    EXPECT_TRUE(analytic.converged);
    EXPECT_TRUE(analytic.analytic);
    check_precision(analytic, p);
  }
  model::Options<R> options(p);
  EXPECT_THROW((model::Options<R>(Precision::exact())), std::invalid_argument);
  EXPECT_THROW(model::relaxation_gap(5, 2, R(1), half, options), std::invalid_argument);
  EXPECT_THROW(model::relaxation_gap(5, 2, one, R(0), options), std::invalid_argument);
  options.seed_options = bethe::tasep::Options<R>(Precision::bits(80));
  EXPECT_THROW((model::PreparedGap<R>(5, 2, options)), std::invalid_argument);
  options = model::Options<R>(p);
  options.seed_options.tolerance = R(-1, p);
  EXPECT_THROW(model::relaxation_gap(6, 2, one, one, options), std::invalid_argument);
}

TEST(AsepMpfr, RealBlockJacobianAtFiniteAndTinyBias)
{
  auto p = Precision::bits(160);
  R one(1, p), h = math::ldexp(one, -45);
  model::Options<R> options(p);
  for (R q : {one / 2, one - R("1e-30", p)})
  {
    auto state = model::relaxation_gap(5, 2, one, q, options);
    ASSERT_TRUE(state.converged);
    model::detail::System<R> system{5, 2, state.wave_index, q, one - q, state.wave_base};
    auto f = system.evaluate(state.scaled_roots, true);
    for (std::size_t column = 0; column < 4; ++column)
    {
      auto a = state.scaled_roots, b = a;
      C step(column % 2 ? R(0, p) : h, column % 2 ? h : R(0, p));
      a[column / 2] += step;
      b[column / 2] -= step;
      auto fa = system.evaluate(a), fb = system.evaluate(b);
      for (std::size_t row = 0; row < 4; ++row)
      {
        EXPECT_EQ(uni20::precision_of(f.jacobian[row * 4 + column]), p);
        EXPECT_LT(math::abs((fa.residual[row] - fb.residual[row]) / (2 * h) - f.jacobian[row * 4 + column]),
                  1000000 * h * h);
      }
    }
  }
}

TEST(AsepMpfr, OwningMixedPrecisionSharedSeedReaders)
{
  using State = model::GapState<R>;
  using Prepared = model::PreparedGap<R>;
  auto check = [&](uni20::async::IAsyncScheduler& scheduler) {
    std::array<std::vector<State>, 6> actual, expected;
    std::atomic<unsigned> preparations{0};
    std::vector<uni20::async::Async<std::size_t>> completions;
    for (std::size_t family = 0; family < 3; ++family)
    {
      auto p = Precision::bits(family == 0 ? 80 : family == 1 ? 160 : 256);
      model::Options<R> options(p);
      if (family == 2) options.seed_options.max_seed_iterations = 0;
      auto background = bethe::prepare_background_async(
          [options, &preparations] {
            ++preparations;
            return Prepared(5, 2, options);
          },
          {.scheduler = &scheduler});
      for (std::size_t b = 0; b < 2; ++b)
      {
        std::size_t index = 2 * family + b;
        R rate(b + 1, p);
        for (std::size_t i = 0; i < 3; ++i)
          expected[index].push_back(model::relaxation_gap(5, 2, rate, rate * R(i) / 2, options));
        completions.push_back(bethe::evaluate_observable_grid_async(
            background.read(), 3,
            [rate](Prepared const& prepared, std::size_t i) { return prepared.at_rates(rate, rate * R(i) / 2); },
            [&, index](std::size_t i, State s) {
              EXPECT_EQ(i, actual[index].size());
              actual[index].push_back(std::move(s));
            },
            {.scheduler = &scheduler, .max_in_flight = 2}));
      }
    }
    for (auto& done : completions)
      EXPECT_EQ(done.get_wait(scheduler), 3u);
    EXPECT_EQ(preparations.load(), 3u);
    for (std::size_t family = 0; family < 6; ++family)
    {
      ASSERT_EQ(actual[family].size(), 3u);
      for (std::size_t i = 0; i < 3; ++i)
      {
        EXPECT_TRUE(fields(actual[family][i]) == fields(expected[family][i]));
        check_precision(actual[family][i], uni20::precision_of(expected[family][i].right_rate));
      }
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
