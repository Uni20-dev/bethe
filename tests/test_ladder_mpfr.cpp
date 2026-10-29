// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/ladder.hpp>
#include <bethe/prepared_observable_grid.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using C = uni20::complex<R>;
using uni20::Precision;
namespace model = bethe::ladder;
namespace math = uni20::math;

auto fields(model::State<R> const& s)
{
  auto const& b = s.highest_weight;
  return std::tie(s.rungs, s.singlets, s.iterations, s.branches, s.tableaux, s.populations, s.rung_coupling,
                  s.magnetic_field, s.magnetization, s.energy, s.converged, s.descendant, s.analytic, s.status, b.shape,
                  b.labels, b.rapidities, b.energy, b.residual, b.momentum_index, b.iterations, b.converged, b.status);
}
void check_precision(model::State<R> const& s, Precision p)
{
  auto check = [&](R const& x) {
    ASSERT_TRUE(x.initialized());
    EXPECT_EQ(uni20::precision_of(x), p);
    EXPECT_TRUE(uni20::isfinite(x));
  };
  for (auto const* x : {&s.rung_coupling, &s.magnetic_field, &s.highest_weight.energy, &s.highest_weight.residual})
    check(*x);
  if (s.energy) check(*s.energy);
  for (auto const& level : s.highest_weight.rapidities)
    for (auto const& x : level)
      check(x);
}
void compare(model::SectorScan<R> const& actual, model::SectorScan<R> const& expected, Precision p)
{
  EXPECT_EQ(std::tie(actual.iterations, actual.branches, actual.tableaux, actual.complete, actual.status),
            std::tie(expected.iterations, expected.branches, expected.tableaux, expected.complete, expected.status));
  ASSERT_EQ(actual.sectors.size(), expected.sectors.size());
  for (std::size_t i = 0; i < actual.sectors.size(); ++i)
  {
    check_precision(actual.sectors[i], p);
    EXPECT_EQ(fields(actual.sectors[i]), fields(expected.sectors[i]));
  }
}

TEST(LadderMpfr, TypedProductBranchesZeroBudgetsAndValidation)
{
  auto const p = Precision::bits(160);
  model::SolverOptions<R> options(p);
  options.max_branches = options.max_iterations = 0;
  auto const product = model::ground_state<R>(7, R(6, p), options);
  ASSERT_TRUE(product.converged && product.analytic);
  EXPECT_EQ(*product.energy, R(-105, p) / 4);
  check_precision(product, p);
  auto const polarized = model::magnetization_ground_state<R>(7, -7, R(-2, p), R(-9, p), options);
  ASSERT_TRUE(polarized.converged && polarized.analytic);
  check_precision(polarized, p);
  model::PreparedSectors<R> prepared(6, options);
  EXPECT_FALSE(prepared.complete());
  EXPECT_EQ(prepared.branches(), 0u);
  auto const empty = prepared.sector_ground_states(R(1, p));
  compare(empty, model::sector_ground_states<R>(6, R(1, p), options), p);
  EXPECT_FALSE(empty.complete);
  for (auto const& s : empty.sectors)
    EXPECT_FALSE(s.energy.has_value());
  auto const bypass = prepared.ground_state(R(6, p));
  EXPECT_TRUE(bypass.converged && bypass.analytic);
  check_precision(bypass, p);
  for (R invalid : {R{}, R(1), R(1, Precision::bits(80))})
  {
    EXPECT_THROW((void)prepared.ground_state(invalid), std::invalid_argument);
    EXPECT_THROW((void)prepared.ground_state(R(1, p), invalid), std::invalid_argument);
  }
  options.residual_tolerance = R(1);
  EXPECT_THROW((void)model::ground_state<R>(7, R(6, p), options), std::invalid_argument);
}

TEST(LadderMpfr, DescendantAlgebraicEnergyGainsAccuracyBeyondBinary128)
{
  auto const reference = Precision::bits(768);
  R const exact = -R(1, reference) / 2 - math::sqrt(R(5, reference));
  R previous(1, reference);
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto const p = Precision::bits(bits);
    model::SolverOptions<R> options(p);
    auto const s = model::sector_ground_state<R>(6, 4, R(0, p), options);
    ASSERT_TRUE(s.converged);
    ASSERT_TRUE(s.energy);
    EXPECT_TRUE(s.descendant);
    EXPECT_EQ(s.highest_weight.shape, (model::detail::Shape{4, 2, 0, 0}));
    check_precision(s, p);
    R const error = math::abs(s.energy->at(reference) - exact);
    EXPECT_LT(error, (4096 * uni20::numeric_limits<R>::epsilon(p)).at(reference));
    EXPECT_LT(error, previous * R("1e-15", reference));
    if (bits == 256) EXPECT_LT(error, R("1e-65", reference));
    previous = error;
    R const field(1, p);
    auto const shifted = model::magnetization_ground_state<R>(6, 1, R(0, p), field, options, 4);
    ASSERT_TRUE(shifted.converged);
    EXPECT_LT(math::abs(*shifted.energy + field - *s.energy), 4096 * uni20::numeric_limits<R>::epsilon(p));
    check_precision(shifted, p);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(LadderMpfr, OriginalComplexEquationsAndShiftedSeas)
{
  auto const p = Precision::bits(160);
  R const zero(0, p), one(1, p), eps = uni20::numeric_limits<R>::epsilon(p);
  auto e = [](R const& x, R const& width) { return C(x, width / 2) / C(x, -width / 2); };
  for (auto const& shape : {model::detail::Shape{2, 1, 1, 0}, {3, 2, 1, 0}})
    for (auto const& shift : model::detail::shifts(shape))
    {
      auto const b = model::detail::solve<R>(shape, shift, 32 * eps, 10000);
      ASSERT_TRUE(b.converged);
      std::size_t const n = shape[0] + shape[1] + shape[2] + shape[3];
      R energy(n, p);
      for (std::size_t a = 0; a < b.rapidities.size(); ++a)
        for (auto const& root : b.rapidities[a])
        {
          C product(-one, zero);
          if (a == 0)
          {
            for (std::size_t j = 0; j < n; ++j)
              product *= e(root, -one);
            energy -= one / (root * root + one / 4);
          }
          for (auto const& x : b.rapidities[a])
            product *= e(root - x, 2 * one);
          if (a)
            for (auto const& x : b.rapidities[a - 1])
              product *= e(root - x, -one);
          if (a + 1 < b.rapidities.size())
            for (auto const& x : b.rapidities[a + 1])
              product *= e(root - x, -one);
          EXPECT_LT(math::abs(product - C(one, zero)), 1024 * R(n, p) * eps);
        }
      EXPECT_LT(math::abs(energy - b.energy), 1024 * R(n, p) * eps);
    }
}

TEST(LadderMpfr, PreparedQueriesPreserveCompleteAndIncompleteScans)
{
  auto const p = Precision::bits(160);
  for (int kind : {0, 1, 2, 3})
  {
    model::SolverOptions<R> options(p);
    if (kind == 1) options.max_iterations = 0;
    if (kind == 2) options.max_branches = 1;
    if (kind == 3) options.max_branches = 3;
    model::PreparedSectors<R> prepared(6, options);
    EXPECT_EQ(prepared.options().precision(), p);
    for (auto const& values : {std::pair{0, 0}, {1, 2}, {-1, -2}, {6, 1}})
    {
      R const rung(values.first, p), field(values.second, p);
      compare(prepared.sector_ground_states(rung, field), model::sector_ground_states(6, rung, field, options), p);
      auto const s = prepared.ground_state(rung, field);
      EXPECT_EQ(fields(s), fields(model::ground_state(6, rung, field, options)));
      check_precision(s, p);
    }
  }
}

TEST(LadderMpfr, AsyncPreparedGridsOwnMixedPrecisionContexts)
{
  auto const ambient = mpfr_get_default_prec();
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::array<std::vector<model::State<R>>, 2> actual;
    auto submit = [&](std::size_t slot, int bits) {
      auto const p = Precision::bits(bits);
      model::SolverOptions<R> options(p);
      auto background =
          bethe::prepare_background_async([options] { return model::PreparedSectors<R>(4, options); }, {&scheduler});
      return bethe::evaluate_observable_grid_async(
          background.read(), 3,
          [p](auto const& prepared, std::size_t i) { return prepared.ground_state(R(1, p), R(i, p) / 7); },
          [&, slot](std::size_t, auto state) { actual[slot].push_back(std::move(state)); },
          {.scheduler = &scheduler, .max_in_flight = 2});
    };
    auto first = submit(0, 80), second = submit(1, 160);
    EXPECT_EQ(first.get_wait(scheduler), 3u);
    EXPECT_EQ(second.get_wait(scheduler), 3u);
    for (std::size_t slot = 0; slot < 2; ++slot)
    {
      auto const p = Precision::bits(slot ? 160 : 80);
      ASSERT_EQ(actual[slot].size(), 3u);
      for (std::size_t i = 0; i < 3; ++i)
      {
        check_precision(actual[slot][i], p);
        auto const expected = model::ground_state<R>(4, R(1, p), R(i, p) / 7, model::SolverOptions<R>(p));
        EXPECT_EQ(fields(actual[slot][i]), fields(expected));
      }
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
} // namespace
