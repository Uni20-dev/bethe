// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/biquadratic_ferromagnetic.hpp>
#include <bethe/biquadratic_two_string.hpp>
#include <bethe/observable_grid.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using C = uni20::complex<R>;
using uni20::Precision;
namespace ts = bethe::xxz::quantum_group::two_string;
namespace math = uni20::math;

void precision(ts::State<R> const& s, Precision p)
{
  for (auto const* v : {&s.delta, &s.center, &s.log_deviation, &s.energy, &s.energy_shift, &s.residual_norm,
                        &s.phase_residual, &s.modulus_residual})
    EXPECT_EQ(uni20::precision_of(*v), p);
  for (auto const& v : s.rapidities)
    EXPECT_EQ(uni20::precision_of(v), p);
}

auto fields(ts::State<R> const& s)
{
  return std::tie(s.sites, s.delta, s.center, s.log_deviation, s.rapidities, s.energy, s.energy_shift, s.residual_norm,
                  s.phase_residual, s.modulus_residual, s.iterations, s.converged, s.status, s.string_label,
                  s.deviation_sign, s.real_labels);
}

// Independent real-component evaluation of the original complex Bethe equations.
// This deliberately does not require the separate MPC hyperbolic extension.
C sinh_components(C const& z)
{
  return {math::sinh(z.real()) * math::cos(z.imag()), math::cosh(z.real()) * math::sin(z.imag())};
}

void original_equations(ts::State<R> const& s)
{
  auto p = s.delta.precision();
  R const zero(0, p), one(1, p), eta = math::acosh(s.delta);
  R const d = s.deviation_sign * math::exp(-s.log_deviation);
  std::vector<C> roots;
  for (auto const& alpha : s.rapidities)
    roots.emplace_back(zero, alpha / 2);
  auto const pair = roots.size();
  roots.emplace_back((eta + d) / 2, s.center / 2);
  roots.emplace_back((eta + d) / 2, -s.center / 2);
  for (std::size_t i = 0; i < roots.size(); ++i)
  {
    C const drive = sinh_components(roots[i] + eta / 2) / sinh_components(roots[i] - eta / 2);
    C lhs(one, zero), rhs(one, zero);
    for (std::size_t k = 0; k < 2 * s.sites; ++k)
      lhs *= drive;
    for (std::size_t j = 0; j < roots.size(); ++j)
      if (i != j)
      {
        C const reflected =
            i >= pair && j >= pair ? C(math::sinh(d), zero) : sinh_components(roots[i] + roots[j] - eta);
        rhs *= sinh_components(roots[i] - roots[j] + eta) * sinh_components(roots[i] + roots[j] + eta) /
               (sinh_components(roots[i] - roots[j] - eta) * reflected);
      }
    // Reconstruction of eta+d still loses relative precision near a pole.
    EXPECT_LT(math::abs(lhs / rhs - C(one, zero)),
              8192 * R(s.sites, p) * uni20::numeric_limits<R>::epsilon(p) / math::abs(d));
  }
}

TEST(OpenTwoStringMpfr, IndependentFourSiteEnergyGainsAccuracy)
{
  auto const oracle = Precision::bits(768);
  auto const ambient = mpfr_get_default_prec();
  R previous(1, oracle);
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto p = Precision::bits(bits);
    bethe::SolverOptions<R> options(p);
    R error(0, oracle);
    for (auto const& delta : {R(3, p) / 2, R(2, p), R(5, p) / 2})
    {
      auto const s = ts::singlet<R>(4, delta, options);
      ASSERT_TRUE(s.converged) << int(s.status);
      precision(s, p);
      R const d = delta.at(oracle);
      R const exact = -3 * d / 4 + math::sqrt(d * d + 2) / 2;
      error = std::max(error, math::abs(s.energy.at(oracle) - exact));
      original_equations(s);
    }
    EXPECT_LT(error, 65536 * uni20::numeric_limits<R>::epsilon(p).at(oracle));
    EXPECT_LT(error, previous * R("1e-15", oracle));
    if (bits == 256) EXPECT_LT(error, R("1e-70", oracle));
    previous = error;
    auto const mapped = bethe::biquadratic::two_string::singlet<R>(4, options);
    ASSERT_TRUE(mapped.reference.converged);
    EXPECT_EQ(mapped.energy.precision(), p);
    EXPECT_EQ(mapped.tl_energy.precision(), p);
    EXPECT_LT(math::abs(mapped.energy - (-15 + math::sqrt(R(17, p))) / 2),
              65536 * uni20::numeric_limits<R>::epsilon(p));
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(OpenTwoStringMpfr, JacobianAndFiniteDeviationEquations)
{
  auto p = Precision::bits(160);
  R const one(1, p), h = math::ldexp(one, -50);
  for (auto const& delta : {R(11, p) / 10, R(3, p) / 2, R(3, p)})
    for (bool ideal : {false, true})
    {
      ts::detail::System<R> system(8, delta);
      auto x = system.seed();
      x.back() = R(4, p);
      std::vector<R> jac;
      auto value = system.evaluate(x, &jac, ideal);
      EXPECT_EQ(value.norm.precision(), p);
      EXPECT_EQ(value.phase_norm.precision(), p);
      EXPECT_EQ(value.modulus_norm.precision(), p);
      for (auto const& v : jac)
        EXPECT_EQ(v.precision(), p);
      for (std::size_t j = 0; j < x.size(); ++j)
      {
        auto a = x, b = x;
        a[j] += h;
        b[j] -= h;
        auto plus = system.evaluate(a, nullptr, ideal), minus = system.evaluate(b, nullptr, ideal);
        for (std::size_t i = 0; i < x.size(); ++i)
          EXPECT_LT(math::abs((plus.residual[i] - minus.residual[i]) / (2 * h) - jac[i * x.size() + j]),
                    R(1000000, p) * h * h * std::max(one, math::abs(jac[i * x.size() + j])));
      }
    }
  bethe::SolverOptions<R> options(p);
  for (std::size_t n : {6, 8, 16, 32})
  {
    auto const s = ts::singlet<R>(n, R(3, p) / 2, options);
    ASSERT_TRUE(s.converged) << n << ' ' << int(s.status);
    precision(s, p);
    original_equations(s);
  }
  for (std::size_t mode : {1, 2})
  {
    auto const s = bethe::biquadratic::ferromagnetic::bound_pair<R>(12, mode, options);
    ASSERT_TRUE(s.reference.converged) << mode << ' ' << int(s.reference.status);
    EXPECT_EQ(s.energy.precision(), p);
    EXPECT_EQ(s.tl_energy.precision(), p);
    EXPECT_EQ(s.tl_energy, -2 * s.reference.energy_shift);
    precision(s.reference, p);
    original_equations(s.reference);
  }
  auto const mixed = bethe::biquadratic::ferromagnetic::pair_defect<R>(10, 1, 1, options);
  ASSERT_TRUE(mixed.reference.converged);
  precision(mixed.reference, p);
  original_equations(mixed.reference);
}

TEST(OpenTwoStringMpfr, LogarithmicDeviationAndFailurePrecision)
{
  for (int bits : {80, 160, 256})
  {
    auto p = Precision::bits(bits);
    bethe::SolverOptions<R> options(p);
    R const delta = R(3, p) / 2;
    auto const bound = ts::bound_pair<R>(256, delta, 1, options);
    ASSERT_TRUE(bound.converged) << bits << ' ' << int(bound.status);
    precision(bound, p);
    EXPECT_GT(bound.log_deviation, R(100, p));
    EXPECT_LE(bound.residual_norm, options.residual_tolerance);
    for (std::size_t budget : {0, 1})
    {
      options.max_iterations = budget;
      auto const s = ts::singlet<R>(16, delta, options);
      EXPECT_FALSE(s.converged);
      EXPECT_EQ(s.status, bethe::xxz::quantum_group::SolveStatus::iteration_limit);
      EXPECT_EQ(s.iterations, budget);
      precision(s, p);
      auto x = s.rapidities;
      x.push_back(s.center);
      x.push_back(s.log_deviation);
      ts::detail::System<R> system(16, delta);
      EXPECT_EQ(s.residual_norm, system.evaluate(x).norm);
      EXPECT_EQ(s.energy, system.energy(x, delta));
    }
  }
}

TEST(OpenTwoStringMpfr, RejectsInvalidPrecisionBeforeSolving)
{
  auto p = Precision::bits(160);
  R const delta = R(3, p) / 2;
  bethe::SolverOptions<R> options(p);
  for (auto const& bad : {R(uni20::uninitialized), R{}, R(2), R(2, Precision::bits(80)), R("nan", p), R(1, p)})
    EXPECT_THROW(ts::singlet<R>(4, bad, options), std::invalid_argument);
  for (auto const& bad : {R(uni20::uninitialized), R{}, R(1), R(1, Precision::bits(80)), R("nan", p), R(0, p)})
  {
    options.residual_tolerance = bad;
    EXPECT_THROW(ts::singlet<R>(4, delta, options), std::invalid_argument);
    EXPECT_THROW(ts::bound_pair<R>(8, delta, 1, options), std::invalid_argument);
  }
  options = bethe::SolverOptions<R>(p);
  EXPECT_THROW(bethe::temperley_lieb::two_string::singlet<R>(4, R(uni20::uninitialized), options),
               std::invalid_argument);
  EXPECT_THROW(ts::pair_with_real_roots<R>(8, delta, std::vector<std::size_t>{2, 1}, 1, options),
               std::invalid_argument);
}

TEST(OpenTwoStringMpfr, OwningMixedPrecisionQueriesPreserveCompleteResults)
{
  struct Request
  {
      R delta;
      bethe::SolverOptions<R> options;
      bool pair;
  };
  std::vector<Request> requests;
  for (int bits : {80, 160})
    for (bool pair : {false, true})
      for (bool fail : {false, true})
      {
        auto p = Precision::bits(bits);
        bethe::SolverOptions<R> options(p);
        if (fail) options.max_iterations = 0;
        requests.push_back({R(3, p) / 2, options, pair});
      }
  auto solve = [](Request const& r) {
    return r.pair ? ts::bound_pair<R>(12, r.delta, 1, r.options) : ts::singlet<R>(8, r.delta, r.options);
  };
  std::vector<ts::State<R>> expected;
  for (auto const& r : requests)
    expected.push_back(solve(r));
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<ts::State<R>> actual;
    auto done = [&] {
      auto owned = requests;
      return bethe::evaluate_observable_grid_async(
          owned.size(), [owned, solve](std::size_t i) { return solve(owned[i]); },
          [&actual](std::size_t, ts::State<R> s) { actual.push_back(std::move(s)); },
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
