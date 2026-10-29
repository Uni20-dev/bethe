// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/biquadratic_qsystem.hpp>
#include <bethe/observable_grid.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using C = uni20::complex<R>;
using uni20::Precision;
namespace qs = bethe::xxz::quantum_group::qsystem;
namespace math = uni20::math;

void precision(qs::State<R> const& state, Precision p)
{
  for (auto const* v : {&state.delta, &state.residual_norm, &state.bethe_residual, &state.bethe_residual_bound,
                        &state.roots.residual_norm, &state.roots.reconstruction_error,
                        &state.roots.max_root_uncertainty, &state.roots.separation_ratio})
    EXPECT_EQ(uni20::precision_of(*v), p);
  if (state.energy) EXPECT_EQ(uni20::precision_of(*state.energy), p);
  for (auto const& values : {state.coefficients, state.partner})
    for (auto const& v : values)
      EXPECT_EQ(uni20::precision_of(v), p);
  for (auto const& values : {state.roots.roots, state.rapidities})
    for (auto const& v : values)
      EXPECT_EQ(uni20::precision_of(v), p);
}

auto fields(qs::State<R> const& s)
{
  return std::tie(s.sites, s.through_lines, s.iterations, s.delta, s.residual_norm, s.bethe_residual,
                  s.bethe_residual_bound, s.coefficients, s.partner, s.roots.roots, s.roots.status, s.roots.iterations,
                  s.roots.residual_norm, s.roots.reconstruction_error, s.roots.max_root_uncertainty,
                  s.roots.separation_ratio, s.rapidities, s.energy, s.status, s.converged);
}

TEST(OpenQSystemMpfr, IndependentComplexBranchGainsAccuracyBeyondBinary128)
{
  auto const oracle = Precision::bits(768);
  R previous(1, oracle);
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto p = Precision::bits(bits);
    bethe::SolverOptions<R> options(p);
    R error(0, oracle);
    std::vector<R> seed{R(3, p) / 2, -R(4, p) / 3};
    for (auto const& delta : {R(3, p) / 2, R(7, p) / 4, R(2, p), R(5, p) / 2})
    {
      auto const s = qs::solve<R>(4, delta, seed, options);
      ASSERT_TRUE(s.converged) << int(s.status) << " bits=" << bits;
      ASSERT_TRUE(s.energy);
      precision(s, p);
      seed = s.coefficients;
      R const d = delta.at(oracle);
      // Independent two-dimensional TL Hamiltonian eigenvalue.
      R const exact = -3 * d / 4 + math::sqrt(d * d + 2) / 2;
      error = std::max(error, math::abs(s.energy->at(oracle) - exact));
      ASSERT_EQ(s.roots.roots.size(), 2u);
      EXPECT_GT(math::abs(s.roots.roots[0].imag()), 1);
      EXPECT_LE(s.bethe_residual, s.bethe_residual_bound);
    }
    EXPECT_LT(error, 65536 * uni20::numeric_limits<R>::epsilon(p).at(oracle));
    EXPECT_LT(error, previous * R("1e-15", oracle));
    if (bits == 256) EXPECT_LT(error, R("1e-70", oracle));
    previous = error;
    auto const mapped = bethe::biquadratic::qsystem::solve<R>(4, std::vector<R>{R(3, p) / 2, -R(4, p) / 3}, options);
    ASSERT_TRUE(mapped.energy);
    EXPECT_EQ(uni20::precision_of(*mapped.energy), p);
    EXPECT_EQ(uni20::precision_of(*mapped.tl_energy), p);
    EXPECT_LT(math::abs(*mapped.energy - (-15 + math::sqrt(R(17, p))) / 2),
              65536 * uni20::numeric_limits<R>::epsilon(p));
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(OpenQSystemMpfr, IndependentWronskianAndJacobian)
{
  auto p = Precision::bits(160);
  R const zero(0, p), one(1, p), h = math::ldexp(one, -50);
  for (auto const& delta : {R(11, p) / 10, R(3, p) / 2, R(3, p)})
    for (std::size_t m : {1, 2, 3})
    {
      qs::System<R> system(6, m, delta);
      std::vector<R> c(m, one / 4), jac;
      auto const value = system.evaluate(c, &jac);
      EXPECT_EQ(uni20::precision_of(value.norm), p);
      for (std::size_t j = 0; j < m; ++j)
      {
        auto a = c, b = c;
        a[j] += h;
        b[j] -= h;
        auto plus = system.evaluate(a), minus = system.evaluate(b);
        for (std::size_t i = 0; i < m; ++i)
          EXPECT_LT(math::abs((plus.residual[i] - minus.residual[i]) / (2 * h) - jac[i * m + j]),
                    R(1000000, p) * h * h * std::max(one, math::abs(jac[i * m + j])));
      }
    }
  bethe::SolverOptions<R> options(p);
  R const delta = R(3, p) / 2;
  auto const s = qs::solve<R>(4, delta, std::vector<R>{delta, -R(4, p) / 3}, options);
  ASSERT_TRUE(s.converged);
  auto q = s.coefficients;
  q.push_back(one);
  auto at = [&](std::vector<R> const& c, C x) {
    C value(zero, zero);
    for (std::size_t k = c.size(); k-- > 0;)
      value = value * x + c[k];
    return value;
  };
  // Direct shifted hyperbolic Wronskian, independent of coefficient recurrences.
  C const u(one / 5, R(7, p) / 10);
  R const eta = math::acosh(delta);
  C const xp = math::cosh(2 * u + eta), xm = math::cosh(2 * u - eta), x = math::cosh(2 * u);
  C const lhs = at(s.partner, xp) * at(q, xm) - at(s.partner, xm) * at(q, xp);
  C rhs = 2 * math::sinh(eta) * math::sinh(2 * u);
  for (int k = 0; k < 4; ++k)
    rhs *= x - one;
  EXPECT_LT(math::abs(lhs - rhs), 8192 * uni20::numeric_limits<R>::epsilon(p) * std::max(one, math::abs(rhs)));
}

TEST(OpenQSystemMpfr, CompleteSmallModulesAndBoundedSearch)
{
  for (int bits : {80, 160})
  {
    SCOPED_TRACE(bits);
    auto p = Precision::bits(bits);
    bethe::SolverOptions<R> options(p);
    R const d = R(3, p) / 2, eps = uni20::numeric_limits<R>::epsilon(p);
    auto const pair = bethe::biquadratic::qsystem::spectrum<R>(4, 0, {.max_attempts = 100}, options);
    ASSERT_TRUE(pair.complete());
    ASSERT_EQ(pair.states.size(), 2u);
    for (int i = 0; i < 2; ++i)
    {
      R const exact = (-15 + (i ? 1 : -1) * math::sqrt(R(17, p))) / 2;
      ASSERT_TRUE(pair.states[i].energy);
      EXPECT_LT(math::abs(*pair.states[i].energy - exact), 65536 * eps);
      precision(pair.states[i].reference, p);
    }
    for (unsigned n : {4, 6})
    {
      auto const spectrum = qs::spectrum<R>(n, d, n - 2, {.max_attempts = 100}, options);
      ASSERT_TRUE(spectrum.complete());
      ASSERT_EQ(spectrum.states.size(), n - 1);
      for (unsigned i = 0; i < n - 1; ++i)
      {
        auto const& s = spectrum.states[i];
        R const exact = R(n - 1, p) * d / 4 - d + math::cos(uni20::pi<R>.at(p) * (n - 1 - i) / n);
        ASSERT_TRUE(s.energy);
        EXPECT_LT(math::abs(*s.energy - exact), 65536 * eps);
        precision(s, p);
      }
    }
    auto const bounded = qs::spectrum<R>(4, d, 0, {.max_attempts = 2}, options);
    EXPECT_LE(bounded.attempts, 2u);
    EXPECT_LE(bounded.states.size() + bounded.failed_attempts, bounded.attempts);
    auto const empty = qs::spectrum<R>(4, d, 0, {.max_attempts = 0}, options);
    EXPECT_EQ(empty.attempts, 0u);
    EXPECT_TRUE(empty.states.empty());
    EXPECT_FALSE(empty.complete());
  }
}

TEST(OpenQSystemMpfr, ExplicitPrecisionOnFailureAndInvalidInputs)
{
  auto p = Precision::bits(160);
  R const zero(0, p), d = R(3, p) / 2;
  bethe::SolverOptions<R> options(p);
  std::vector<R> const seed{d, -R(4, p) / 3};
  for (std::size_t budget : {0, 1})
  {
    options.max_iterations = budget;
    auto const s = qs::solve<R>(4, d, seed, options);
    EXPECT_EQ(s.status, qs::Status::iteration_limit);
    EXPECT_EQ(s.iterations, budget);
    EXPECT_FALSE(s.converged);
    precision(s, p);
  }
  auto const vacuum = qs::solve<R>(4, d, {}, options);
  EXPECT_TRUE(vacuum.converged);
  EXPECT_EQ(*vacuum.energy, R(9, p) / 8);
  precision(vacuum, p);
  options.max_iterations = 0;
  auto const pole = qs::solve<R>(4, d, std::vector<R>{-d}, options);
  EXPECT_FALSE(pole.energy);
  precision(pole, p);
  options = bethe::SolverOptions<R>(p);
  for (auto const& bad : {R(uni20::uninitialized), R{}, R(1), R(2, Precision::bits(80)), R("nan", p), R(1, p)})
    EXPECT_THROW(qs::solve<R>(4, bad, seed, options), std::invalid_argument);
  for (auto const& bad : {R(uni20::uninitialized), R{}, R(1), R(1, Precision::bits(80)), R("nan", p)})
  {
    EXPECT_THROW(qs::solve<R>(4, d, std::vector<R>{bad}, options), std::invalid_argument);
    options.residual_tolerance = bad;
    EXPECT_THROW(qs::solve<R>(4, d, seed, options), std::invalid_argument);
    options = bethe::SolverOptions<R>(p);
  }
  options.residual_tolerance = zero;
  EXPECT_THROW(qs::solve<R>(4, d, seed, options), std::invalid_argument);
  options.residual_tolerance = R(1, p);
  auto const loose = qs::solve<R>(4, d, std::vector<R>{zero}, options);
  EXPECT_FALSE(loose.converged);
  EXPECT_EQ(loose.status, qs::Status::inadmissible);
}

TEST(OpenQSystemMpfr, OwningMixedPrecisionComposedResults)
{
  struct Request
  {
      R delta;
      std::vector<R> seed;
      bethe::SolverOptions<R> options;
  };
  std::vector<Request> requests;
  for (int bits : {80, 160})
    for (bool fail : {false, true})
    {
      auto p = Precision::bits(bits);
      bethe::SolverOptions<R> options(p);
      if (fail) options.max_iterations = 0;
      requests.push_back({R(3, p) / 2, {R(3, p) / 2, -R(4, p) / 3}, options});
    }
  auto solve = [](Request const& r) { return qs::solve<R>(4, r.delta, r.seed, r.options); };
  std::vector<qs::State<R>> expected;
  for (auto const& r : requests)
    expected.push_back(solve(r));
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<qs::State<R>> actual;
    auto done = [&] {
      auto owned = requests;
      return bethe::evaluate_observable_grid_async(
          owned.size(), [owned, solve](std::size_t i) { return solve(owned[i]); },
          [&actual](std::size_t, qs::State<R> s) { actual.push_back(std::move(s)); },
          {.scheduler = &scheduler, .max_in_flight = 2});
    }();
    (void)done.get_wait(scheduler);
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      EXPECT_EQ(fields(actual[i]), fields(expected[i]));
      precision(actual[i], requests[i].options.precision());
    }
  }
}
} // namespace
