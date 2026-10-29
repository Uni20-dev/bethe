// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/polynomial_roots.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using C = uni20::complex<R>;
using uni20::Precision;
namespace roots = bethe::detail;
namespace math = uni20::math;
using Status = roots::PolynomialRootStatus;

void precision(roots::PolynomialRoots<R> const& s, Precision p)
{
  for (auto const* x : {&s.residual_norm, &s.reconstruction_error, &s.max_root_uncertainty, &s.separation_ratio})
    EXPECT_EQ(uni20::precision_of(*x), p);
  for (auto const& z : s.roots)
  {
    EXPECT_EQ(uni20::precision_of(z), p);
    EXPECT_EQ(uni20::precision_of(z.real()), p);
    EXPECT_EQ(uni20::precision_of(z.imag()), p);
  }
}

auto fields(roots::PolynomialRoots<R> const& s)
{
  return std::tie(s.roots, s.status, s.iterations, s.residual_norm, s.reconstruction_error, s.max_root_uncertainty,
                  s.separation_ratio);
}

TEST(PolynomialRootsMpfr, IndependentRootsGainAccuracyBeyondBinary128)
{
  auto const oracle = Precision::bits(768);
  R previous(1, oracle);
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto const p = Precision::bits(bits);
    roots::PolynomialRootOptions<R> options(p);
    R const zero(0, p);
    R error(0, oracle);
    for (int degree : {2, 3, 8, 17, 24})
    {
      SCOPED_TRACE(degree);
      // x^n-2 has independent polar roots, including irrational moduli.
      std::vector<R> coefficients(degree, zero);
      coefficients[0] = R(-2, p);
      auto const s = roots::recover_polynomial_roots<R>(coefficients, options);
      ASSERT_EQ(s.status, Status::resolved);
      ASSERT_EQ(s.roots.size(), std::size_t(degree));
      precision(s, p);
      R const radius = math::pow(R(2, oracle), R(1, oracle) / degree);
      std::vector<C> expected;
      for (int i = 0; i < degree; ++i)
      {
        R const a = 2 * uni20::pi<R>.at(oracle) * i / degree;
        expected.emplace_back(radius * math::cos(a), radius * math::sin(a));
      }
      for (auto const& z : s.roots)
      {
        auto best = expected.begin();
        for (auto i = expected.begin(); i != expected.end(); ++i)
          if (math::abs(z.at(oracle) - *i) < math::abs(z.at(oracle) - *best)) best = i;
        error = std::max(error, math::abs(z.at(oracle) - *best));
        expected.erase(best);
      }
      EXPECT_LE(s.residual_norm, options.tolerance);
      EXPECT_LE(s.reconstruction_error, options.tolerance);
      EXPECT_GT(s.separation_ratio, 8);
    }
    EXPECT_LT(error, 8192 * uni20::numeric_limits<R>::epsilon(p).at(oracle));
    EXPECT_LT(error, previous * R("1e-15", oracle));
    if (bits == 256) EXPECT_LT(error, R("1e-70", oracle));
    previous = error;
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(PolynomialRootsMpfr, ScalingVacuumAndBudgetDiagnostics)
{
  auto const p = Precision::bits(160);
  R const zero(0, p), one(1, p), eps = uni20::numeric_limits<R>::epsilon(p);
  roots::PolynomialRootOptions<R> options(p);
  auto const vacuum = roots::recover_polynomial_roots<R>({}, options);
  EXPECT_EQ(vacuum.status, Status::resolved);
  EXPECT_TRUE(vacuum.roots.empty());
  EXPECT_EQ(vacuum.iterations, 0u);
  precision(vacuum, p);
  for (int exponent : {-100, 100, 10000})
  {
    R const scale = math::ldexp(one, exponent);
    std::array<R, 3> const coefficients{-6 * scale * scale * scale, 11 * scale * scale, -6 * scale};
    auto const s = roots::recover_polynomial_roots<R>(coefficients, options);
    ASSERT_EQ(s.status, Status::resolved);
    precision(s, p);
    for (int i = 1; i <= 3; ++i)
    {
      R closest(100, p);
      for (auto const& z : s.roots)
        closest = std::min(closest, math::abs(z / scale - R(i, p)));
      EXPECT_LT(closest, 8192 * eps);
    }
  }
  std::array<R, 2> const coefficients{R(-2, p), zero};
  std::array<C, 2> const initial{C(-one, zero), C(one, zero)};
  options.max_iterations = 0;
  auto const initial_only = roots::recover_polynomial_roots<R>(coefficients, options, initial);
  EXPECT_EQ(initial_only.status, Status::iteration_limit);
  EXPECT_EQ(initial_only.iterations, 0u);
  EXPECT_EQ(initial_only.roots, (std::vector<C>(initial.begin(), initial.end())));
  EXPECT_GT(initial_only.residual_norm, R("0.1", p));
  precision(initial_only, p);
  options.max_iterations = 1;
  auto const single = roots::recover_polynomial_roots<R>(coefficients, options, initial);
  EXPECT_EQ(single.status, Status::iteration_limit);
  EXPECT_EQ(single.iterations, 1u);
  precision(single, p);
}

TEST(PolynomialRootsMpfr, ClustersAndDuplicatedRootsAreNotReportedAsResolved)
{
  auto const p = Precision::bits(160);
  R const zero(0, p), one(1, p);
  roots::PolynomialRootOptions<R> options(p);
  for (auto const& c : {std::vector<R>{zero, zero}, std::vector<R>{one, R(-2, p)}})
  {
    auto const s = roots::recover_polynomial_roots<R>(c, options);
    EXPECT_NE(s.status, Status::resolved);
    precision(s, p);
  }
  std::array<R, 2> const c{-one, zero};
  std::array<C, 2> const duplicated{C(one, zero), C(one, zero)};
  options.max_iterations = 0;
  auto const failed = roots::recover_polynomial_roots<R>(c, options, duplicated);
  EXPECT_EQ(failed.residual_norm, 0);
  EXPECT_GT(failed.reconstruction_error, 1);
  EXPECT_EQ(failed.status, Status::iteration_limit);
  options.max_iterations = 200;
  EXPECT_EQ(roots::recover_polynomial_roots<R>(c, options, duplicated).status, Status::unresolved_cluster);
  // Coefficients distinguish the roots here; recovery must not collapse them.
  R const gap = math::ldexp(one, -30);
  auto const nearby = roots::recover_polynomial_roots<R>(std::array<R, 2>{one - gap * gap, R(-2, p)}, options);
  ASSERT_EQ(nearby.status, Status::resolved);
  EXPECT_GT(math::abs(nearby.roots[0] - nearby.roots[1]), gap);
}

TEST(PolynomialRootsMpfr, InvalidAndOutOfRangeInputsKeepExplicitDiagnostics)
{
  auto const p = Precision::bits(160);
  R const zero(0, p), one(1, p);
  roots::PolynomialRootOptions<R> options(p);
  std::array<R, 2> const c{-one, zero};
  for (R const& bad : {R{}, R(1), R(1, Precision::bits(80)), R("nan", p)})
    EXPECT_THROW(roots::recover_polynomial_roots<R>(std::array<R, 1>{bad}, options), std::invalid_argument);
  for (C const& bad : {C{}, C(1), C(R(1, Precision::bits(80)), R(0, Precision::bits(80))), C(R("nan", p), zero)})
    EXPECT_THROW(roots::recover_polynomial_roots<R>(c, options, std::array<C, 2>{bad, C(one, zero)}),
                 std::invalid_argument);
  for (R const& bad : {R{}, R(1), R(1, Precision::bits(80)), R("nan", p), zero})
  {
    options.tolerance = bad;
    EXPECT_THROW(roots::recover_polynomial_roots<R>(c, options), std::invalid_argument);
  }
  options = roots::PolynomialRootOptions<R>(p);
  options.max_degree = 1;
  EXPECT_THROW(roots::recover_polynomial_roots<R>(c, options), std::length_error);
  EXPECT_THROW((void)roots::PolynomialRootOptions<R>(Precision::exact()), std::invalid_argument);
  options.max_degree = 64;
  std::array<R, 2> const extreme{math::next_up(zero), math::next_down(R("inf", p))};
  auto const s = roots::recover_polynomial_roots<R>(extreme, options);
  EXPECT_EQ(s.status, Status::nonfinite_or_unrepresentable);
  EXPECT_TRUE(s.roots.empty());
  precision(s, p);
}

TEST(PolynomialRootsMpfr, OwningMixedPrecisionResultsAgreeAcrossSchedulers)
{
  struct Request
  {
      roots::PolynomialRootOptions<R> options;
      std::vector<R> coefficients;
  };
  std::vector<Request> requests;
  for (int bits : {80, 160})
    for (bool failed : {false, true})
    {
      auto const p = Precision::bits(bits);
      roots::PolynomialRootOptions<R> options(p);
      if (failed) options.max_iterations = 0;
      requests.push_back({options, {R(-2, p), R(0, p), R(0, p)}});
    }
  auto solve = [](Request const& r) { return roots::recover_polynomial_roots<R>(r.coefficients, r.options); };
  std::vector<roots::PolynomialRoots<R>> expected;
  for (auto const& r : requests)
    expected.push_back(solve(r));
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<roots::PolynomialRoots<R>> actual;
    auto done = [&] {
      auto owned = requests;
      return bethe::evaluate_observable_grid_async(
          owned.size(), [owned, solve](std::size_t i) { return solve(owned[i]); },
          [&actual](std::size_t, roots::PolynomialRoots<R> s) { actual.push_back(std::move(s)); },
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
