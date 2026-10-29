// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/tasep.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using C = uni20::complex<R>;
using uni20::Precision;
namespace model = bethe::tasep;
namespace math = uni20::math;

void check_precision(model::GapState<R> const& s, Precision p)
{
  EXPECT_EQ(uni20::precision_of(s.rate), p);
  EXPECT_EQ(uni20::precision_of(s.residual_norm), p);
  EXPECT_EQ(uni20::precision_of(s.log_y), p);
  for (auto const& root : s.roots)
    EXPECT_EQ(uni20::precision_of(root), p);
  if (s.eigenvalue) EXPECT_EQ(uni20::precision_of(*s.eigenvalue), p);
  if (s.gap) EXPECT_EQ(uni20::precision_of(*s.gap), p);
}

TEST(TasepMpfr, IndependentMarkovEigenvalueGainsAccuracy)
{
  // Independent 90-decimal-digit eigensystem of the ten-state L=5,N=2
  // Markov generator; the native oracle checks every filling through L=10.
  auto high = Precision::bits(512);
  C const exact{
      R("-0.710290069123028672493762532365387684613821535715147285599991234743117333842877636225644555", high),
      R("0.326895513054695736451443863914497675039691838878110784103129825870014387978645791366909038", high)};
  R previous(1, high);
  auto ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto p = Precision::bits(bits);
    model::Options<R> options(p);
    auto state = model::relaxation_gap(5, 2, R(1, p), options);
    ASSERT_TRUE(state.converged) << int(state.status) << " seed " << state.seed_iterations << " updates "
                                 << state.iterations;
    check_precision(state, p);
    auto error = math::abs(state.eigenvalue->at(high) - exact);
    EXPECT_LT(error, R(8192, high) * uni20::numeric_limits<R>::epsilon(p).at(high));
    EXPECT_LT(error, previous * R("1e-15", high));
    if (bits == 256) EXPECT_LT(error, R("1e-70", high));
    previous = error;
    auto partner = model::relaxation_gap(5, 3, R(1, p), options);
    ASSERT_TRUE(partner.converged);
    EXPECT_EQ(partner.eigenvalue, state.eigenvalue);
    auto scaled = model::relaxation_gap(5, 2, R(3, p), options);
    ASSERT_TRUE(scaled.converged);
    EXPECT_EQ(*scaled.eigenvalue, 3 * *state.eigenvalue);
    EXPECT_LE(state.residual_norm, options.tolerance);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(TasepMpfr, AnalyticSectorsAndFailureDiagnosticsRetainPrecision)
{
  auto p = Precision::bits(160);
  model::Options<R> options(p);
  for (std::size_t n : {0u, 6u})
  {
    auto s = model::relaxation_gap(6, n, R(1, p), options);
    ASSERT_TRUE(s.converged);
    EXPECT_EQ(s.status, model::Status::stationary_only);
    EXPECT_FALSE(s.eigenvalue);
    EXPECT_FALSE(s.gap);
    check_precision(s, p);
  }
  for (std::size_t n : {1u, 5u})
  {
    auto s = model::relaxation_gap(6, n, R(1, p), options);
    ASSERT_TRUE(s.converged);
    check_precision(s, p);
    EXPECT_LT(math::abs(*s.gap - R(1, p) / 2), 32 * uni20::numeric_limits<R>::epsilon(p));
  }
  for (int failure : {0, 1})
  {
    options = model::Options<R>(p);
    if (failure == 0)
      options.max_seed_iterations = 0;
    else
      options.max_iterations = 0;
    auto s = model::relaxation_gap(8, 3, R(1, p), options);
    EXPECT_FALSE(s.converged);
    EXPECT_FALSE(s.eigenvalue);
    EXPECT_FALSE(s.gap);
    EXPECT_EQ(s.status, failure == 0 ? model::Status::seed_limit : model::Status::iteration_limit);
    check_precision(s, p);
  }
  options = model::Options<R>(p);
  EXPECT_THROW((model::Options<R>(Precision::exact())), std::invalid_argument);
  EXPECT_THROW(model::relaxation_gap(5, 2, R(1), options), std::invalid_argument);
  EXPECT_THROW(model::relaxation_gap(5, 2, R{}, options), std::invalid_argument);
  EXPECT_THROW(model::relaxation_gap(5, 2, R(1, Precision::bits(80)), options), std::invalid_argument);
  options.tolerance = R(1, Precision::bits(80));
  EXPECT_THROW(model::relaxation_gap(5, 2, R(1, p), options), std::invalid_argument);
}

TEST(TasepMpfr, RealBlockJacobianAndOwningMixedPrecisionRequests)
{
  auto p = Precision::bits(160);
  model::Options<R> options(p);
  auto state = model::relaxation_gap(5, 2, R(1, p), options);
  ASSERT_TRUE(state.converged);
  model::detail::System<R> system{5, 2};
  auto f = system.evaluate(state.roots, state.log_y, true);
  R h = math::ldexp(R(1, p), -45);
  for (std::size_t column = 0; column < 6; ++column)
  {
    auto a = state.roots, b = state.roots;
    auto ya = state.log_y, yb = state.log_y;
    C step(column % 2 ? R(0, p) : h, column % 2 ? h : R(0, p));
    if (column < 4)
    {
      a[column / 2] += step;
      b[column / 2] -= step;
    }
    else
    {
      ya += step;
      yb -= step;
    }
    auto fa = system.evaluate(a, ya, false), fb = system.evaluate(b, yb, false);
    for (std::size_t row = 0; row < 6; ++row)
    {
      EXPECT_EQ(uni20::precision_of(f.jacobian[row * 6 + column]), p);
      EXPECT_LT(math::abs((fa.residual[row] - fb.residual[row]) / (2 * h) - f.jacobian[row * 6 + column]),
                1000000 * h * h);
    }
  }
  struct Request
  {
      R rate;
      model::Options<R> options;
  };
  std::vector<Request> requests;
  for (int bits : {80, 160, 256})
  {
    auto q = Precision::bits(bits);
    model::Options<R> o(q);
    requests.push_back({R(1, q), o});
    o.max_iterations = 0;
    requests.push_back({R(1, q), o});
  }
  std::vector<model::GapState<R>> expected;
  for (auto const& r : requests)
    expected.push_back(model::relaxation_gap(5, 2, r.rate, r.options));
  auto fields = [](auto const& s) {
    return std::tie(s.sites, s.particles, s.effective_particles, s.iterations, s.seed_iterations, s.rate,
                    s.residual_norm, s.roots, s.log_y, s.eigenvalue, s.gap, s.converged, s.status);
  };
  auto check = [&](uni20::async::IAsyncScheduler& scheduler) {
    std::vector<model::GapState<R>> actual;
    auto done = [&] {
      auto owned = requests;
      return bethe::evaluate_observable_grid_async(
          owned.size(),
          [owned](std::size_t i) {
            auto const& r = owned[i];
            return model::relaxation_gap(5, 2, r.rate, r.options);
          },
          [&](std::size_t i, auto s) {
            EXPECT_EQ(i, actual.size());
            actual.push_back(std::move(s));
          },
          {.scheduler = &scheduler, .max_in_flight = 2});
    }();
    EXPECT_EQ(done.get_wait(scheduler), requests.size());
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      EXPECT_TRUE(fields(actual[i]) == fields(expected[i]));
      check_precision(actual[i], requests[i].options.precision());
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
