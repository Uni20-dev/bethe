// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/takhtajan_babujian.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using C = uni20::complex<R>;
namespace model = bethe::takhtajan_babujian;
namespace math = uni20::math;

void check_precision(model::State<R> const& s, uni20::Precision p)
{
  auto check = [&](R const& value) {
    EXPECT_TRUE(value.initialized());
    EXPECT_FALSE(value.is_exact());
    EXPECT_EQ(uni20::precision_of(value), p);
    EXPECT_TRUE(uni20::isfinite(value));
  };
  for (auto const* value : {&s.energy, &s.momentum, &s.residual_norm, &s.phase_residual, &s.modulus_residual})
    check(*value);
  for (auto const& value : s.centers)
    check(value);
  for (auto const& value : s.deviations)
    check(value);
  for (auto const& value : s.rapidities)
  {
    check(value.real());
    check(value.imag());
  }
}
} // namespace

TEST(TBMpfr, ExplicitPrecisionSeedsIncludeCentralZero)
{
  auto const p = uni20::Precision::bits(160);
  model::detail::GroundSystem<R> system(6, R(0, p));
  auto const seed = system.seed();
  EXPECT_TRUE(system.physical(seed));
  auto const [centers, deviations] = system.expand(seed);
  ASSERT_EQ(centers.size(), 3u);
  EXPECT_EQ(centers[1], R(0, p));
  for (auto const& values : {seed, centers, deviations})
    for (auto const& value : values)
      EXPECT_EQ(uni20::precision_of(value), p);
  auto const evaluation = system.evaluate(seed);
  EXPECT_EQ(uni20::precision_of(evaluation.phase), p);
  EXPECT_EQ(uni20::precision_of(evaluation.modulus), p);
}

#if UNI20_HAS_MPLAPACK_MPREAL
TEST(TBMpfr, FourSiteAnalyticEnergyGainsAccuracyBeyondBinary128)
{
  auto const reference_precision = uni20::Precision::bits(768);
  R const exact = -11 - math::sqrt(R(41, reference_precision));
  R previous(1, reference_precision);
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto const p = uni20::Precision::bits(bits);
    model::SolverOptions<R> options(p);
    auto const state = model::ground_state<R>(4, options);
    ASSERT_TRUE(state.converged) << int(state.status) << " iterations=" << state.iterations;
    check_precision(state, p);
    EXPECT_LE(state.residual_norm, options.residual_tolerance);
    R const error = math::abs(uni20::scalar_at_precision<R>(state.energy, reference_precision) - exact);
    EXPECT_LT(error, uni20::scalar_at_precision<R>(2048 * uni20::numeric_limits<R>::epsilon(p), reference_precision));
    EXPECT_LT(error, previous * R("1e-15", reference_precision));
    if (bits == 256) EXPECT_LT(error, R("1e-65", reference_precision));
    previous = error;
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(TBMpfr, OriginalComplexEquationsAndReflection)
{
  auto const p = uni20::Precision::bits(160);
  R const zero(0, p), one(1, p), eps = uni20::numeric_limits<R>::epsilon(p);
  C const unit(one, zero), i(zero, one), four(R(4, p), zero);
  for (std::size_t n : {4, 6, 12})
  {
    SCOPED_TRACE(n);
    auto const state = model::ground_state<R>(n, model::SolverOptions<R>(p));
    ASSERT_TRUE(state.converged);
    check_precision(state, p);
    C energy(zero, zero), momentum = unit;
    for (std::size_t j = 0; j < n; ++j)
    {
      C const z = state.rapidities[j], q = (z + i) / (z - i);
      C product = unit;
      for (std::size_t site = 0; site < n; ++site)
        product *= q;
      for (std::size_t k = 0; k < n; ++k)
        if (j != k) product *= (z - state.rapidities[k] - i) / (z - state.rapidities[k] + i);
      EXPECT_LT(math::abs(product.real() - 1), 512 * R(n, p) * eps);
      EXPECT_LT(math::abs(product.imag()), 512 * R(n, p) * eps);
      energy -= four / (unit + z * z);
      momentum *= q;
    }
    EXPECT_LT(math::abs(energy.real() - state.energy), 128 * R(n, p) * eps);
    EXPECT_LT(math::abs(energy.imag()), 128 * R(n, p) * eps);
    EXPECT_LT(math::abs(momentum.real() - 1), 128 * R(n, p) * eps);
    EXPECT_LT(math::abs(momentum.imag()), 128 * R(n, p) * eps);
    EXPECT_EQ(state.momentum, zero);
    for (std::size_t j = 0; j < n / 2; ++j)
    {
      EXPECT_EQ(state.centers[j], -state.centers[n / 2 - 1 - j]);
      EXPECT_EQ(state.deviations[j], state.deviations[n / 2 - 1 - j]);
      EXPECT_EQ(state.string_quantum_numbers[j].twice(), 2 * std::int64_t(j) - std::int64_t(n / 2 - 1));
    }
  }
}

TEST(TBMpfr, InvalidControlsAndIncompleteStates)
{
  auto const p = uni20::Precision::bits(160);
  model::SolverOptions<R> options(p);
  options.max_iterations = 0;
  auto const initial = model::ground_state<R>(6, options);
  EXPECT_FALSE(initial.converged);
  EXPECT_EQ(initial.status, model::SolveStatus::iteration_limit);
  EXPECT_EQ(initial.iterations, 0u);
  check_precision(initial, p);
  for (R tolerance : {R{}, R(1), R("1e-20", uni20::Precision::bits(80)), R(0, p), R(-1, p)})
  {
    options.residual_tolerance = tolerance;
    EXPECT_THROW((void)model::ground_state<R>(4, options), std::invalid_argument);
  }
}

TEST(TBMpfr, OwningMixedPrecisionGridsMatchCompleteStates)
{
  struct Request
  {
      std::size_t sites;
      model::SolverOptions<R> options;
  };
  std::vector<Request> requests;
  for (int bits : {80, 160})
    for (bool fail : {false, true})
    {
      model::SolverOptions<R> options(uni20::Precision::bits(bits));
      if (fail) options.max_iterations = 0;
      requests.push_back({fail ? 6u : 4u, options});
    }
  std::vector<model::State<R>> expected;
  for (auto const& request : requests)
    expected.push_back(model::ground_state<R>(request.sites, request.options));
  auto const ambient = mpfr_get_default_prec();
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<model::State<R>> actual;
    auto done = bethe::evaluate_observable_grid_async(
        requests.size(),
        [owned = requests](std::size_t index) {
          return model::ground_state<R>(owned[index].sites, owned[index].options);
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
      EXPECT_EQ(a.iterations, b.iterations);
      EXPECT_EQ(a.string_quantum_numbers, b.string_quantum_numbers);
      EXPECT_EQ(a.centers, b.centers);
      EXPECT_EQ(a.deviations, b.deviations);
      EXPECT_EQ(a.rapidities, b.rapidities);
      EXPECT_EQ(a.energy, b.energy);
      EXPECT_EQ(a.momentum, b.momentum);
      EXPECT_EQ(a.momentum_index, b.momentum_index);
      EXPECT_EQ(a.residual_norm, b.residual_norm);
      EXPECT_EQ(a.phase_residual, b.phase_residual);
      EXPECT_EQ(a.modulus_residual, b.modulus_residual);
      EXPECT_EQ(a.converged, b.converged);
      EXPECT_EQ(a.status, b.status);
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
#endif
