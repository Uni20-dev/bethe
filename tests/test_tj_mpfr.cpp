// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/tj.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using C = uni20::complex<R>;
using uni20::Precision;
namespace model = bethe::tj;
namespace math = uni20::math;

void check_precision(model::State<R> const& s, Precision p)
{
  auto check = [&](R const& v) {
    ASSERT_TRUE(v.initialized());
    EXPECT_EQ(uni20::precision_of(v), p);
    EXPECT_TRUE(uni20::isfinite(v));
  };
  for (auto const* v : {&s.energy, &s.momentum, &s.residual_norm})
    check(*v);
  for (auto const& v : s.level_residuals)
    check(v);
  for (auto const& level : s.rapidities)
    for (auto const& v : level)
      check(v);
}

TEST(TJMpfr, SeedsAnalyticBranchesAndZeroBudgetKeepPrecision)
{
  auto const p = Precision::bits(160);
  R const zero(0, p);
  model::SolverOptions<R> options(p);
  model::detail::GroundSystem<R> system(8, 3, 3, zero);
  auto const seed = system.seed();
  EXPECT_TRUE(system.physical(seed));
  auto const roots = system.expand(seed);
  EXPECT_EQ(roots[0][2], zero);
  for (auto const& level : roots)
    for (auto const& r : level)
      EXPECT_EQ(uni20::precision_of(r), p);
  options.max_iterations = 0;
  auto const partial = model::ground_state<R>(8, 3, 3, options);
  EXPECT_FALSE(partial.converged);
  EXPECT_EQ(partial.status, model::SolveStatus::iteration_limit);
  EXPECT_EQ(partial.iterations, 0u);
  check_precision(partial, p);
  for (auto counts : {std::array<std::size_t, 3>{7, 0, 0}, {7, 7, 0}, {7, 0, 4}})
  {
    auto const s = model::ground_state<R>(counts[0], counts[1], counts[2], options);
    EXPECT_TRUE(s.converged);
    EXPECT_EQ(s.iterations, 0u);
    EXPECT_EQ(s.branch, model::Branch::polarized_free);
    check_precision(s, p);
  }
  for (R tolerance : {R{}, R(1), R("1e-20", Precision::bits(80)), R(0, p), R(-1, p)})
  {
    options.residual_tolerance = tolerance;
    EXPECT_THROW((void)model::ground_state<R>(7, 0, 0, options), std::invalid_argument);
  }
}

TEST(TJMpfr, FreeMomentumAndNoHoleReduction)
{
  auto const p = Precision::bits(256), reference = Precision::bits(768);
  model::SolverOptions<R> options(p);
  auto const free = model::ground_state<R>(7, 4, 0, options);
  R energy(0, reference), momentum(0, reference);
  R const pi = uni20::pi<R>.at(reference);
  for (auto mode : free.free_modes)
  {
    R const k = 2 * pi * R(mode, reference) / 7;
    energy -= 2 * math::cos(k);
    momentum += k;
  }
  EXPECT_LT(math::abs(free.energy.at(reference) - energy), R("1e-70", reference));
  EXPECT_LT(math::abs(free.momentum.at(reference) - momentum), R("1e-70", reference));
  check_precision(free, p);
  auto const full = model::ground_state<R>(6, 3, 3, options);
  ASSERT_TRUE(full.converged);
  EXPECT_EQ(full.branch, model::Branch::no_holes_xxx);
  EXPECT_LT(math::abs(full.energy.at(reference) + 5 + math::sqrt(R(13, reference))), R("1e-65", reference));
  check_precision(full, p);
  options.max_iterations = 0;
  auto const incomplete = model::ground_state<R>(6, 3, 3, options);
  EXPECT_FALSE(incomplete.converged);
  EXPECT_EQ(incomplete.iterations, 0u);
  check_precision(incomplete, p);
  // A finite O(1) energy must survive the extensive XXX offset cancellation.
  std::size_t const n = 1000000000000000000ULL;
  auto const near_polarized = model::ground_state<R>(n, n - 1, 1, options);
  EXPECT_TRUE(near_polarized.converged);
  EXPECT_EQ(near_polarized.energy, R(-4, p));
  check_precision(near_polarized, p);
}

TEST(TJMpfr, DopedAnalyticEnergyGainsAccuracyBeyondBinary128)
{
  auto const reference = Precision::bits(768);
  R const exact = -3 - math::sqrt(R(5, reference));
  R previous(1, reference);
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto const p = Precision::bits(bits);
    model::SolverOptions<R> options(p);
    auto const s = model::ground_state<R>(5, 3, 1, options);
    ASSERT_TRUE(s.converged);
    check_precision(s, p);
    R const error = math::abs(s.energy.at(reference) - exact);
    EXPECT_LT(error, (2048 * uni20::numeric_limits<R>::epsilon(p)).at(reference));
    EXPECT_LT(error, previous * R("1e-15", reference));
    if (bits == 256) EXPECT_LT(error, R("1e-65", reference));
    previous = error;
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(TJMpfr, OriginalMultiplicativeEquationsAndSpinReversal)
{
  auto const p = Precision::bits(160);
  R const zero(0, p), one(1, p), eps = uni20::numeric_limits<R>::epsilon(p);
  auto e = [](R const& x, R const& width) { return C(x, width / 2) / C(x, -width / 2); };
  for (auto counts : {std::array<std::size_t, 3>{3, 1, 1}, {5, 3, 1}, {8, 3, 3}, {12, 5, 3}})
  {
    auto const [n, up, down] = counts;
    SCOPED_TRACE(n);
    auto const s = model::ground_state<R>(n, up, down, model::SolverOptions<R>(p));
    ASSERT_TRUE(s.converged);
    check_precision(s, p);
    R energy(2 * s.holes, p);
    for (std::size_t a = 0; a < 2; ++a)
      for (std::size_t j = 0; j < s.rapidities[a].size(); ++j)
      {
        R const& x = s.rapidities[a][j];
        C lhs(one, zero), rhs(one, zero);
        if (a == 0)
        {
          for (std::size_t k = 0; k < n; ++k)
            lhs *= e(x, one);
          for (std::size_t k = 0; k < s.rapidities[0].size(); ++k)
            if (j != k) rhs *= e(x - s.rapidities[0][k], R(2, p));
          for (auto const& y : s.rapidities[1])
            rhs /= e(x - y, one);
          energy -= one / (x * x + one / 4);
        }
        else
          for (auto const& y : s.rapidities[0])
            rhs *= e(x - y, one);
        EXPECT_LT(math::abs(lhs - rhs), 512 * R(n, p) * eps);
      }
    EXPECT_LT(math::abs(energy - s.energy), 128 * R(n, p) * eps);
    auto const reversed = model::ground_state<R>(n, down, up, model::SolverOptions<R>(p));
    EXPECT_EQ(reversed.energy, s.energy);
    EXPECT_EQ(reversed.rapidities, s.rapidities);
  }
}

TEST(TJMpfr, OwningMixedPrecisionAsyncGridsPreserveCompleteStates)
{
  struct Request
  {
      std::size_t n, up, down;
      model::SolverOptions<R> options;
  };
  std::vector<Request> requests;
  for (int bits : {80, 160})
    for (int branch : {0, 1, 2, 3})
    {
      model::SolverOptions<R> options(Precision::bits(bits));
      if (branch == 3) options.max_iterations = 0;
      requests.push_back({branch == 1 ? 6u : 5u, 3u, branch == 2 ? 0u : branch == 1 ? 3u : 1u, options});
    }
  std::vector<model::State<R>> expected;
  for (auto const& r : requests)
    expected.push_back(model::ground_state<R>(r.n, r.up, r.down, r.options));
  auto const ambient = mpfr_get_default_prec();
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<model::State<R>> actual;
    auto done = bethe::evaluate_observable_grid_async(
        requests.size(),
        [owned = requests](std::size_t index) {
          auto const& r = owned[index];
          return model::ground_state<R>(r.n, r.up, r.down, r.options);
        },
        [&actual](std::size_t, model::State<R> state) { actual.push_back(std::move(state)); },
        {.scheduler = &scheduler, .max_in_flight = 2});
    (void)done.get_wait(scheduler);
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t j = 0; j < actual.size(); ++j)
    {
      auto const& a = actual[j];
      auto const& b = expected[j];
      check_precision(a, requests[j].options.precision());
      EXPECT_EQ(a.sites, b.sites);
      EXPECT_EQ(a.num_up, b.num_up);
      EXPECT_EQ(a.num_down, b.num_down);
      EXPECT_EQ(a.holes, b.holes);
      EXPECT_EQ(a.branch, b.branch);
      EXPECT_EQ(a.quantum_numbers, b.quantum_numbers);
      EXPECT_EQ(a.rapidities, b.rapidities);
      EXPECT_EQ(a.free_modes, b.free_modes);
      EXPECT_EQ(a.energy, b.energy);
      EXPECT_EQ(a.momentum, b.momentum);
      EXPECT_EQ(a.residual_norm, b.residual_norm);
      EXPECT_EQ(a.level_residuals, b.level_residuals);
      EXPECT_EQ(a.momentum_index, b.momentum_index);
      EXPECT_EQ(a.iterations, b.iterations);
      EXPECT_EQ(a.converged, b.converged);
      EXPECT_EQ(a.spin_reversed, b.spin_reversed);
      EXPECT_EQ(a.status, b.status);
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
} // namespace
