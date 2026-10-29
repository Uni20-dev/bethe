// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/lieb_liniger.hpp>
#include <bethe/observable_grid.hpp>
#include <bethe/su_fermions.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using C = uni20::complex<R>;
using uni20::Precision;
namespace model = bethe::su_fermions;
namespace math = uni20::math;

void check_precision(model::State<R> const& s, Precision p)
{
  auto check = [&](R const& value) {
    ASSERT_TRUE(value.initialized());
    EXPECT_EQ(uni20::precision_of(value), p);
    EXPECT_TRUE(uni20::isfinite(value));
  };
  for (auto const* value : {&s.length, &s.interaction, &s.residual_norm, &s.target_residual_norm, &s.momentum})
    check(*value);
  for (auto const& level : s.rapidities)
    for (auto const& value : level)
      check(value);
  for (auto const& value : s.level_residuals)
    check(value);
  if (s.energy) check(*s.energy);
  if (s.reached_interaction) check(*s.reached_interaction);
}

void check_equations(model::State<R> const& s, Precision p)
{
  R const zero(0, p), one(1, p), c = *s.reached_interaction;
  R const tolerance = 4096 * R(s.particles, p) * uni20::numeric_limits<R>::epsilon(p);
  auto ratio = [](R const& d, R const& w) { return C(d, w) / C(d, -w); };
  auto const& roots = s.rapidities;
  for (std::size_t a = 0; a < roots.size(); ++a)
    for (std::size_t j = 0; j < roots[a].size(); ++j)
    {
      C lhs(one, zero), rhs(one, zero);
      if (a == 0)
      {
        lhs = math::exp(C(zero, roots[0][j] * s.length));
        for (auto const& r : roots[1])
          rhs *= ratio(roots[0][j] - r, c / 2);
      }
      else
      {
        for (auto const& r : roots[a - 1])
          lhs *= ratio(roots[a][j] - r, c / 2);
        if (a + 1 < roots.size())
          for (auto const& r : roots[a + 1])
            lhs *= ratio(roots[a][j] - r, c / 2);
        for (std::size_t k = 0; k < roots[a].size(); ++k)
          if (k != j) rhs *= ratio(roots[a][j] - roots[a][k], c);
      }
      EXPECT_LT(math::abs(lhs - rhs), tolerance);
      EXPECT_EQ(roots[a][j], -roots[a][roots[a].size() - 1 - j]);
    }
}
} // namespace

TEST(SUFermionsMpfr, ExplicitOptionsAndNestedCentralZeros)
{
  auto const p = Precision::bits(160);
  EXPECT_THROW((void)model::SolverOptions<R>(Precision::exact()), std::invalid_argument);
  model::SolverOptions<R> options(p);
  EXPECT_EQ(options.precision(), p);
  EXPECT_EQ(uni20::precision_of(options.residual_tolerance), p);
  model::detail::GroundSystem<R> system({5, 2, 1}, R(0, p));
  auto const x = system.seed();
  ASSERT_TRUE(system.physical(x, R(320, p)));
  auto const roots = system.expand(x, R(320, p));
  EXPECT_EQ(roots[0][2], R(0, p));
  EXPECT_EQ(roots[2][0], R(0, p));
  for (auto const& level : roots)
    for (auto const& value : level)
      EXPECT_EQ(uni20::precision_of(value), p);
  auto const evaluation = system.evaluate(x, R(320, p));
  EXPECT_EQ(uni20::precision_of(evaluation.norm), p);
  for (auto const& value : evaluation.levels)
    EXPECT_EQ(uni20::precision_of(value), p);
}

TEST(SUFermionsMpfr, TwoBodyJumpConditionGainsAccuracyBeyondBinary128)
{
  auto const reference_precision = Precision::bits(768);
  R const pi = 4 * math::atan(R(1, reference_precision)), exact = pi * pi / 2;
  R previous(1, reference_precision);
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto const p = Precision::bits(bits);
    auto const s = model::ground_state<R>(std::vector<std::size_t>{1, 1, 0}, R(1, p),
                                          uni20::scalar_at_precision<R>(pi, p), model::SolverOptions<R>(p));
    ASSERT_TRUE(s.converged) << int(s.status);
    check_precision(s, p);
    // Independent contact jump condition: at c=pi,L=1, k=+/-pi/2.
    R const error = math::abs(uni20::scalar_at_precision<R>(*s.energy, reference_precision) - exact);
    EXPECT_LT(error, uni20::scalar_at_precision<R>(4096 * uni20::numeric_limits<R>::epsilon(p), reference_precision));
    EXPECT_LT(error, previous * R("1e-15", reference_precision));
    if (bits == 256) EXPECT_LT(error, R("1e-65", reference_precision));
    previous = error;
  }
}

TEST(SUFermionsMpfr, NestedEquationsAndIndependentBosonReduction)
{
  auto const p = Precision::bits(160);
  for (auto const& pop :
       {std::vector<std::size_t>{1, 1, 1}, std::vector<std::size_t>{3, 1, 1}, std::vector<std::size_t>{1, 1, 1, 1}})
    for (char const* coupling : {"0.01", "1", "100"})
    {
      R const length(1, p), c(coupling, p);
      auto const s = model::ground_state<R>(pop, length, c, model::SolverOptions<R>(p));
      ASSERT_TRUE(s.converged) << int(s.status) << " c=" << coupling;
      check_precision(s, p);
      check_equations(s, p);
      if (pop.front() == 1)
      {
        auto const b = bethe::lieb_liniger::ground_state(pop.size(), length, c, bethe::SolverOptions<R>(p));
        ASSERT_TRUE(b.converged);
        EXPECT_LT(math::abs(*s.energy - b.energy), 8192 * uni20::numeric_limits<R>::epsilon(p) * (1 + b.energy));
      }
    }
}

TEST(SUFermionsMpfr, FreeVacuumPermutationAndIncompleteStages)
{
  auto const p = Precision::bits(160);
  model::SolverOptions<R> options(p);
  for (auto const& pop :
       {std::vector<std::size_t>{0, 0, 0}, std::vector<std::size_t>{2, 1, 4}, std::vector<std::size_t>{0, 6, 0}})
  {
    auto const s = model::ground_state<R>(pop, R(2, p), R(0, p), options);
    ASSERT_TRUE(s.converged && s.free);
    check_precision(s, p);
    R const pi = 4 * math::atan(R(1, p));
    R exact(0, p);
    for (auto n : pop)
      exact += pi * pi * R(n, p) * (R(n, p) * R(n, p) + (n % 2 ? -1 : 2)) / 12;
    EXPECT_LT(math::abs(*s.energy - exact), 256 * uni20::numeric_limits<R>::epsilon(p) * (1 + exact));
  }
  auto const a = model::ground_state<R>(std::vector<std::size_t>{3, 1, 1}, R(1, p), R(1, p), options);
  auto const b = model::ground_state<R>(std::vector<std::size_t>{1, 0, 3, 1}, R(2, p), R(1, p) / 2, options);
  ASSERT_TRUE(a.converged && b.converged);
  EXPECT_EQ(*a.energy, 4 * *b.energy);
  EXPECT_EQ(b.component_order, (std::vector<std::size_t>{2, 0, 3}));
  options.max_stages = 1;
  auto const partial = model::ground_state<R>(std::vector<std::size_t>{3, 1, 1}, R(1, p), R(1, p), options);
  ASSERT_FALSE(partial.converged);
  ASSERT_TRUE(partial.energy && partial.reached_interaction);
  EXPECT_GT(*partial.reached_interaction, partial.interaction);
  EXPECT_EQ(partial.status, model::SolveStatus::stage_limit);
  check_precision(partial, p);
  check_equations(partial, p);
  options.max_stages = 0;
  auto const empty = model::ground_state<R>(std::vector<std::size_t>{3, 1, 1}, R(1, p), R(1, p), options);
  EXPECT_FALSE(empty.energy);
  EXPECT_TRUE(empty.rapidities.empty());
  check_precision(empty, p);
  EXPECT_THROW(model::ground_state<R>(std::vector<std::size_t>{1, 1}, R(1), R(1, p), options), std::invalid_argument);
  EXPECT_THROW(model::ground_state<R>(std::vector<std::size_t>{1, 1}, R(1, p), R(1, Precision::bits(80)), options),
               std::invalid_argument);
  options.residual_tolerance = R{};
  EXPECT_THROW(model::ground_state<R>(std::vector<std::size_t>{0}, R(1, p), R(0, p), options), std::invalid_argument);
}

TEST(SUFermionsMpfr, OwningMixedPrecisionGridsPreserveCompleteStates)
{
  struct Request
  {
      std::vector<std::size_t> populations;
      R length, c;
      model::SolverOptions<R> options;
  };
  std::vector<Request> requests;
  for (int bits : {80, 160})
    for (int kind : {0, 1, 2, 3})
    {
      auto const p = Precision::bits(bits);
      model::SolverOptions<R> options(p);
      if (kind == 1) options.max_stages = 1;
      if (kind == 2) options.max_iterations = 0;
      requests.push_back(
          {kind == 3 ? std::vector<std::size_t>{0, 0} : std::vector<std::size_t>{1, 1, 1}, R(1, p), R(1, p), options});
    }
  std::vector<model::State<R>> expected;
  for (auto const& r : requests)
    expected.push_back(model::ground_state<R>(r.populations, r.length, r.c, r.options));
  auto const ambient = mpfr_get_default_prec();
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<model::State<R>> actual;
    auto done = bethe::evaluate_observable_grid_async(
        requests.size(),
        [owned = requests](std::size_t i) {
          auto const& r = owned[i];
          return model::ground_state<R>(r.populations, r.length, r.c, r.options);
        },
        [&actual](std::size_t, model::State<R> s) { actual.push_back(std::move(s)); },
        {.scheduler = &scheduler, .max_in_flight = 2});
    (void)done.get_wait(scheduler);
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t j = 0; j < actual.size(); ++j)
    {
      auto const& a = actual[j];
      auto const& b = expected[j];
      check_precision(a, requests[j].options.precision());
      EXPECT_EQ(a.populations, b.populations);
      EXPECT_EQ(a.component_order, b.component_order);
      EXPECT_EQ(a.particles, b.particles);
      EXPECT_EQ(a.iterations, b.iterations);
      EXPECT_EQ(a.stages, b.stages);
      EXPECT_EQ(a.length, b.length);
      EXPECT_EQ(a.interaction, b.interaction);
      EXPECT_EQ(a.residual_norm, b.residual_norm);
      EXPECT_EQ(a.target_residual_norm, b.target_residual_norm);
      EXPECT_EQ(a.momentum, b.momentum);
      EXPECT_EQ(a.momentum_index, b.momentum_index);
      EXPECT_EQ(a.quantum_numbers, b.quantum_numbers);
      EXPECT_EQ(a.rapidities, b.rapidities);
      EXPECT_EQ(a.level_residuals, b.level_residuals);
      EXPECT_EQ(a.free_modes, b.free_modes);
      EXPECT_EQ(a.energy, b.energy);
      EXPECT_EQ(a.reached_interaction, b.reached_interaction);
      EXPECT_EQ(a.free, b.free);
      EXPECT_EQ(a.converged, b.converged);
      EXPECT_EQ(a.status, b.status);
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
