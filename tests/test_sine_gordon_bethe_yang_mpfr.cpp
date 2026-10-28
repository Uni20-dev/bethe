// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/sine_gordon_bethe_yang.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using C = uni20::complex<R>;
using uni20::Precision;
namespace sg = bethe::sine_gordon;
namespace math = uni20::math;
// Independent rapidity-space p=2 phase and closed complex kernel. Generated
// with reference_sine_gordon_excitations.py at 128 and 160 decimal digits;
// all six output fields agree within 1e-120. Production uses Fourier integrals.
char const* const phase_p2 = "0."
                             "21695200769630950421273978546518952993433638101065038565253108557308211615399945348037891"
                             "940363253371032766020689986186262";
char const* const pair_rapidity = "0."
                                  "291770417218905301263082140230901300783511514593893326056881773976976334577794130100"
                                  "29003765353674236942549006555545113324";
char const* const pair_energy = "2."
                                "08573561878514048803097793223750795709003826540309156866332797938201260358145525632036"
                                "16916129861034106042386768200734652";
char const* const kernel_p2_real = "0."
                                   "04698435098352462286189247477730282383025855612644623064181508907938153276288066552"
                                   "4046851816073805700771802695655888094614";
char const* const kernel_p2_imag = "-0."
                                   "00213066271844231787665042412916624296574750708185896903780303670298725471289818338"
                                   "95728047339489055385343777449173513941922";

void check(sg::BetheYangResult<R> const& state, Precision p)
{
  for (auto const* value : {&state.residual, &state.phase_error, &state.quadrature_error, &state.tail_bound})
  {
    ASSERT_TRUE(value->initialized());
    EXPECT_EQ(uni20::precision_of(*value), p);
  }
  EXPECT_EQ(state.converged, state.status == sg::BetheYangStatus::converged);
  EXPECT_EQ(state.converged, state.energy.has_value());
  EXPECT_EQ(state.converged, state.rapidity.has_value());
  for (auto const* value : {&state.energy, &state.rapidity})
    if (*value)
    {
      EXPECT_EQ(uni20::precision_of(**value), p);
      EXPECT_TRUE(uni20::isfinite(**value));
    }
}
auto fields(sg::BetheYangResult<R> const& s)
{
  return std::tie(s.rapidity, s.energy, s.residual, s.phase_error, s.quadrature_error, s.tail_bound, s.iterations,
                  s.evaluations, s.phase_status, s.status, s.converged);
}

TEST(SineGordonBetheYangMpfr, IndependentPhaseKernelAndLevelsGainAccuracy)
{
  auto const reference = Precision::bits(768);
  R previous(1, reference);
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto const p = Precision::bits(bits);
    sg::KernelOptions<R> kernel_options(p);
    auto const phase = sg::soliton_phase(R("0.7", p), R(2, p), kernel_options);
    ASSERT_TRUE(phase.converged) << int(phase.status);
    auto const kernel = sg::scattering_kernel(C(R("0.7", p), R("0.2", p)), R(2, p), kernel_options);
    ASSERT_TRUE(kernel.converged) << int(kernel.status);
    sg::BetheYangOptions<R> options(p);
    auto const pair = sg::same_charge_pair(R(1, p), R(10, p), R(2, p), uni20::from_twice(1), options);
    ASSERT_TRUE(pair.converged) << int(pair.status) << ", " << pair.iterations;
    check(pair, p);
    std::array<R, 5> const actual{*phase.phase, kernel.value->real(), kernel.value->imag(), *pair.rapidity,
                                  *pair.energy};
    std::array<char const*, 5> const expected{phase_p2, kernel_p2_real, kernel_p2_imag, pair_rapidity, pair_energy};
    R error(0, reference);
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      EXPECT_EQ(uni20::precision_of(actual[i]), p);
      error = std::max(error, math::abs(actual[i].at(reference) - R(expected[i], reference)));
    }
    EXPECT_LT(error, 32768 * uni20::numeric_limits<R>::epsilon(p).at(reference));
    EXPECT_LT(error, previous * R("1e-15", reference));
    if (bits == 256) EXPECT_LT(error, R("1e-70", reference));
    previous = error;
    // Reconstruct the original exponentiated equation, including S_ss's minus sign.
    auto const phase_at_root = sg::soliton_phase(2 * *pair.rapidity, R(2, p), kernel_options);
    ASSERT_TRUE(phase_at_root.converged);
    R const counting = 10 * math::sinh(*pair.rapidity) + *phase_at_root.phase;
    EXPECT_LT(math::abs(math::sin(counting)), 2 * options.tolerance);
    EXPECT_LT(math::abs(math::cos(counting) + 1), 2 * options.tolerance);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(SineGordonBetheYangMpfr, FreeAttractiveAndNearFreeIdentitiesPreservePrecision)
{
  auto const p = Precision::bits(160);
  R const zero(0, p), one(1, p), theta("0.7", p), pi = bethe::detail::pi_like(zero);
  sg::KernelOptions<R> options(p);
  auto const attractive = sg::soliton_phase(theta, one / 2, options);
  ASSERT_TRUE(attractive.converged);
  EXPECT_LT(math::abs(*attractive.phase + math::atan(math::sinh(theta))), options.tolerance);
  auto const negative = sg::soliton_phase(-theta, one / 2, options);
  ASSERT_TRUE(negative.converged);
  EXPECT_EQ(*negative.phase, -*attractive.phase);
  auto const free = sg::scattering_kernel(C(theta, zero), one, options);
  ASSERT_TRUE(free.converged);
  EXPECT_EQ(free.evaluations, 0u);
  EXPECT_EQ(free.value->real(), zero);
  EXPECT_EQ(uni20::precision_of(free.value->imag()), p);
  R const delta = 8 * uni20::numeric_limits<R>::epsilon(p);
  options.tolerance = 4096 * uni20::numeric_limits<R>::epsilon(p) * delta;
  auto const near = sg::soliton_phase(theta, one + delta, options);
  ASSERT_TRUE(near.converged);
  EXPECT_GT(*near.phase, zero);
  EXPECT_LT(math::abs(*near.phase / delta - pi / 2 * math::tanh(theta / 2)),
            8192 * uni20::numeric_limits<R>::epsilon(p));
  auto const p333 = Precision::bits(333);
  sg::BetheYangOptions<R> pair_options(p333);
  auto const pair = sg::same_charge_pair(R(1, p333), R(10, p333), R(1, p333), uni20::from_twice(1), pair_options);
  ASSERT_TRUE(pair.converged);
  check(pair, p333);
  EXPECT_EQ(pair.evaluations, 0u);
  R const pi333 = bethe::detail::pi_like(R(0, p333));
  EXPECT_LT(math::abs(*pair.energy - 2 * math::hypot(R(1, p333), pi333 / 10)),
            128 * uni20::numeric_limits<R>::epsilon(p333));
}

TEST(SineGordonBetheYangMpfr, FailuresAndRejectedControlsRetainTheirMeaning)
{
  auto const p = Precision::bits(160);
  R const zero(0, p), one(1, p);
  sg::BetheYangOptions<R> controls(p);
  for (int failure : {0, 1, 2})
  {
    auto c = controls;
    if (failure == 0) c.max_iterations = 0;
    if (failure == 1) c.phase.max_evaluations = 0;
    if (failure == 2) c.phase.max_cutoffs = 0;
    auto const result = sg::same_charge_pair(one, R(10, p), R(2, p), uni20::from_twice(1), c);
    check(result, p);
    EXPECT_FALSE(result.converged);
    EXPECT_EQ(result.status, failure == 0 ? sg::BetheYangStatus::iteration_limit : sg::BetheYangStatus::phase_limit);
    if (failure > 0)
      EXPECT_EQ(result.phase_status,
                failure == 1 ? sg::KernelStatus::quadrature_limit : sg::KernelStatus::cutoff_limit);
  }
  for (R const& invalid : {R{}, R(1), R(1, Precision::bits(80)), R("nan", p)})
    EXPECT_THROW(sg::same_charge_pair(one, invalid, R(2, p), uni20::from_twice(1), controls), std::invalid_argument);
  EXPECT_THROW(sg::BetheYangOptions<R>(Precision::exact()), std::invalid_argument);
  EXPECT_THROW(sg::KernelOptions<R>(Precision::exact()), std::invalid_argument);
  controls.phase = sg::KernelOptions<R>(Precision::bits(80));
  EXPECT_THROW(sg::same_charge_pair(one, R(10, p), R(2, p), uni20::from_twice(1), controls), std::invalid_argument);
  controls = sg::BetheYangOptions<R>(p);
  controls.phase.max_levels = 25;
  EXPECT_THROW(sg::same_charge_pair(one, R(10, p), R(2, p), uni20::from_twice(1), controls), std::invalid_argument);
}

TEST(SineGordonBetheYangMpfr, OwningMixedPrecisionLevelsAgreeWithDeferredExecution)
{
  auto submit = [](uni20::async::IAsyncScheduler& scheduler, auto& output) {
    struct Request
    {
        R mass, length, coupling;
        sg::BetheYangOptions<R> controls;
    };
    std::vector<Request> requests;
    for (int bits : {80, 160})
    {
      auto const p = Precision::bits(bits);
      sg::BetheYangOptions<R> controls(p);
      controls.tolerance = R("1e-12", p);
      requests.push_back({R(1, p), R(10, p), R(2, p), controls});
      controls.max_iterations = 0;
      requests.push_back({R(1, p), R(10, p), R(2, p), controls});
    }
    return bethe::evaluate_observable_grid_async(
        requests.size(),
        [owned = requests](std::size_t i) {
          auto const& r = owned[i];
          return sg::same_charge_pair(r.mass, r.length, r.coupling, uni20::from_twice(1), r.controls);
        },
        [&output](std::size_t i, auto value) {
          EXPECT_EQ(i, output.size());
          output.push_back(std::move(value));
        },
        {.scheduler = &scheduler, .max_in_flight = 2});
  };
  auto const ambient = mpfr_get_default_prec();
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  std::vector<sg::BetheYangResult<R>> expected;
  auto ready = submit(deferred, expected);
  EXPECT_EQ(ready.get_wait(deferred), 4u);
  for (int threads : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(threads);
    std::vector<sg::BetheYangResult<R>> actual;
    auto done = submit(scheduler, actual);
    EXPECT_EQ(done.get_wait(scheduler), 4u);
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      EXPECT_TRUE(fields(actual[i]) == fields(expected[i]));
      EXPECT_EQ(actual[i].converged, i % 2 == 0);
      check(actual[i], Precision::bits(i < 2 ? 80 : 160));
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
} // namespace
