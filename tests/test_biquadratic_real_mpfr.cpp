// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/biquadratic_ferromagnetic.hpp>
#include <bethe/observable_grid.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using C = uni20::complex<R>;
using uni20::Precision;
namespace math = uni20::math;
namespace model = bethe::biquadratic;
namespace qg = bethe::xxz::quantum_group;

void precision(model::State<R> const& s, Precision p)
{
  auto const& r = s.reference;
  for (auto const* v : {&s.energy, &s.tl_energy, &r.delta, &r.energy, &r.energy_shift, &r.residual_norm})
  {
    EXPECT_EQ(uni20::precision_of(*v), p);
    EXPECT_TRUE(uni20::isfinite(*v));
  }
  for (auto const& values : {r.angles, r.rapidities})
    for (auto const& v : values)
    {
      EXPECT_EQ(uni20::precision_of(v), p);
      EXPECT_TRUE(uni20::isfinite(v));
    }
}

auto fields(model::State<R> const& s)
{
  auto const& r = s.reference;
  return std::tie(s.energy, s.tl_energy, s.through_lines, s.multiplicity, s.exchange, r.sites, r.delta, r.angles,
                  r.rapidities, r.quantum_numbers, r.energy, r.residual_norm, r.iterations, r.converged, r.status,
                  r.energy_shift);
}

TEST(BiquadraticRealMpfr, AnalyticModesAndZeroBudgetKeepPrecisionWithoutAProvider)
{
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256})
  {
    auto const p = Precision::bits(bits);
    bethe::SolverOptions<R> options(p);
    options.max_iterations = 0;
    auto const vacuum = model::ferromagnetic::solve_real<R>(8, {}, options);
    ASSERT_TRUE(vacuum.reference.converged);
    EXPECT_EQ(vacuum.energy, 7);
    EXPECT_EQ(vacuum.tl_energy, 0);
    precision(vacuum, p);
    auto const failed = model::ground_state<R>(8, options);
    EXPECT_FALSE(failed.reference.converged);
    EXPECT_EQ(failed.reference.status, qg::SolveStatus::iteration_limit);
    EXPECT_EQ(failed.reference.iterations, 0u);
    precision(failed, p);
    auto const ferro = model::ferromagnetic::ground_space<R>(8, p);
    EXPECT_EQ(uni20::precision_of(ferro.energy), p);
    auto const level = model::ferromagnetic::one_defect_level<R>(4, 3, p);
    for (auto const* v : {&level.energy, &level.gap, &level.wave_number})
      EXPECT_EQ(uni20::precision_of(*v), p);
    EXPECT_LT(math::abs(level.gap - (3 - math::sqrt(R(2, p)))), 32 * uni20::numeric_limits<R>::epsilon(p));
    EXPECT_EQ(model::ferromagnetic::spectral_gap<R>(4, p), level.gap);
    EXPECT_EQ(level.multiplicity, 8u);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(BiquadraticRealMpfr, InvalidInputsFailBeforeNumericalWork)
{
  auto const p = Precision::bits(160);
  bethe::SolverOptions<R> options(p);
  options.max_iterations = 0;
  for (R const& delta : {R{}, R(3), R("1.5", Precision::bits(80)), R("nan", p), R(1, p)})
    EXPECT_THROW(qg::solve_real<R>(4, delta, {}, options), std::invalid_argument);
  for (R const& tolerance : {R{}, R(1), R("1e-20", Precision::bits(80)), R(0, p), R("nan", p)})
  {
    options.residual_tolerance = tolerance;
    EXPECT_THROW(model::ground_state<R>(4, options), std::invalid_argument);
  }
  EXPECT_THROW(model::ferromagnetic::ground_space<R>(4, Precision::exact()), std::invalid_argument);
  EXPECT_THROW(model::ferromagnetic::spectral_gap<R>(4, Precision::exact()), std::invalid_argument);
  EXPECT_THROW(model::ferromagnetic::one_defect_level<R>(4, 3, Precision::exact()), std::invalid_argument);
#if !UNI20_HAS_MPLAPACK_MPREAL
  options = bethe::SolverOptions<R>(p);
  EXPECT_THROW(model::ground_state<R>(4, options), std::runtime_error);
#endif
}

#if UNI20_HAS_MPLAPACK_MPREAL
TEST(BiquadraticRealMpfr, SpinHamiltonianEnergiesGainAccuracyBeyondBinary128)
{
  auto const oracle = Precision::bits(768);
  R const reference = -(15 + math::sqrt(R(17, oracle))) / 2;
  R previous(1, oracle);
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto const p = Precision::bits(bits);
    bethe::SolverOptions<R> options(p);
    auto const two = model::ground_state<R>(2, options);
    auto const four = model::ground_state<R>(4, options);
    ASSERT_TRUE(two.reference.converged);
    ASSERT_TRUE(four.reference.converged);
    precision(two, p);
    precision(four, p);
    EXPECT_LT(math::abs(two.energy + 4), 256 * uni20::numeric_limits<R>::epsilon(p));
    R const error = math::abs(four.energy.at(oracle) - reference);
    EXPECT_LT(error, 4096 * uni20::numeric_limits<R>::epsilon(p).at(oracle));
    EXPECT_LT(error, previous * R("1e-15", oracle));
    if (bits == 256) EXPECT_LT(error, R("1e-70", oracle));
    previous = error;
    auto const opposite = model::ferromagnetic::solve_real<R>(4, four.reference.quantum_numbers, options);
    EXPECT_EQ(opposite.energy, -four.energy);
    EXPECT_EQ(opposite.reference.rapidities, four.reference.rapidities);
    precision(opposite, p);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(BiquadraticRealMpfr, OriginalMultiplicativeEquationsAndAnalyticBand)
{
  auto const p = Precision::bits(160);
  R const zero(0, p), one(1, p), delta = R(3, p) / 2, eta = math::acosh(delta);
  R const eps = uni20::numeric_limits<R>::epsilon(p);
  bethe::SolverOptions<R> options(p);
  auto ratio = [&](R alpha, R width) {
    C const numerator(math::sinh(width) * math::cos(alpha / 2), math::cosh(width) * math::sin(alpha / 2));
    return numerator / math::conj(numerator);
  };
  for (std::size_t n : {4, 8, 16})
  {
    auto const s = model::ground_state<R>(n, options);
    ASSERT_TRUE(s.reference.converged);
    precision(s, p);
    auto const& roots = s.reference.rapidities;
    R energy = R(n - 1, p) * delta / 4;
    for (std::size_t i = 0; i < roots.size(); ++i)
    {
      C lhs(one, zero), rhs(one, zero);
      auto const drive = ratio(roots[i], eta / 2);
      for (std::size_t k = 0; k < 2 * n; ++k)
        lhs *= drive;
      for (std::size_t j = 0; j < roots.size(); ++j)
        if (i != j) rhs *= ratio(roots[i] - roots[j], eta) * ratio(roots[i] + roots[j], eta);
      EXPECT_LT(math::abs(lhs - rhs), 1024 * R(n, p) * eps);
      energy -= (delta * delta - 1) / (delta - math::cos(roots[i]));
    }
    EXPECT_LT(math::abs(energy - s.reference.energy), 4096 * R(n, p) * eps);
  }
  auto const all = model::ferromagnetic::real_excitations<R>(7, 5, {.count = 100}, options);
  ASSERT_TRUE(all.converged());
  ASSERT_EQ(all.levels.size(), 6u);
  for (std::size_t i = 0; i < all.levels.size(); ++i)
  {
    auto const& level = all.levels[i];
    precision(level.state, p);
    ASSERT_TRUE(level.gap);
    EXPECT_EQ(uni20::precision_of(*level.gap), p);
    auto const exact = model::ferromagnetic::one_defect_level<R>(7, 6 - i, p);
    EXPECT_LT(math::abs(*level.gap - exact.gap), 1024 * eps);
    EXPECT_EQ(level.state.reference.quantum_numbers.front(), uni20::half_int(std::int64_t(6 - i)));
  }
}

TEST(BiquadraticRealMpfr, OwningMixedPrecisionSweepsAndConcurrentScans)
{
  struct Request
  {
      bethe::SolverOptions<R> options;
      bool ferro;
  };
  std::vector<Request> requests;
  std::vector<model::State<R>> expected;
  auto solve = [](Request const& r) {
    auto s = model::ground_state<R>(8, r.options);
    return r.ferro ? model::ferromagnetic::detail::reverse(std::move(s)) : s;
  };
  for (int bits : {80, 160})
    for (bool failed : {false, true})
      for (bool ferro : {false, true})
      {
        bethe::SolverOptions<R> options(Precision::bits(bits));
        if (failed) options.max_iterations = 0;
        requests.push_back({options, ferro});
        expected.push_back(solve(requests.back()));
      }
  auto const ambient = mpfr_get_default_prec();
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<model::State<R>> actual;
    auto done = [&] {
      auto owned = requests;
      return bethe::evaluate_observable_grid_async(
          owned.size(), [owned, solve](std::size_t i) { return solve(owned[i]); },
          [&actual](std::size_t, model::State<R> s) { actual.push_back(std::move(s)); },
          {.scheduler = &scheduler, .max_in_flight = 2});
    }();
    (void)done.get_wait(scheduler);
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      EXPECT_EQ(fields(actual[i]), fields(expected[i]));
      precision(actual[i], requests[i].options.precision());
    }
    for (bool failed : {false, true})
    {
      bethe::SolverOptions<R> options(Precision::bits(160));
      if (failed) options.max_iterations = 0;
      auto const serial = model::real_excitations<R>(8, 4, {.count = 6}, options);
      auto const parallel = model::real_excitations<R>(
          8, 4, {.count = 6, .execution = {.scheduler = &scheduler, .max_in_flight = 2}}, options);
      EXPECT_EQ(parallel.candidate_count, serial.candidate_count);
      EXPECT_EQ(parallel.converged_count, serial.converged_count);
      EXPECT_EQ(fields(parallel.ground_state), fields(serial.ground_state));
      ASSERT_EQ(parallel.first_unconverged.has_value(), serial.first_unconverged.has_value());
      if (parallel.first_unconverged) EXPECT_EQ(fields(*parallel.first_unconverged), fields(*serial.first_unconverged));
      ASSERT_EQ(parallel.levels.size(), serial.levels.size());
      for (std::size_t i = 0; i < parallel.levels.size(); ++i)
      {
        EXPECT_EQ(fields(parallel.levels[i].state), fields(serial.levels[i].state));
        EXPECT_EQ(parallel.levels[i].gap, serial.levels[i].gap);
        precision(parallel.levels[i].state, options.precision());
      }
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
#endif
} // namespace
