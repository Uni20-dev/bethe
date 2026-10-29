// SPDX-License-Identifier: GPL-3.0-or-later
#include "sine_gordon_vacuum_references.hpp"
#include "test_support.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/sine_gordon_excited.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>
namespace
{
using R = uni20::mpreal;
using C = uni20::complex<R>;
using P = uni20::Precision;
namespace sg = bethe::sine_gordon;
namespace math = uni20::math;
using State = sg::TwoSolitonState<R>;
void precision(State const& s, P p)
{
  for (auto const& x : {s.mass, s.length, s.coupling, s.scaled_length, s.contour_shift, s.verification_contour_shift,
                        s.cutoff, s.nonlinear_residual, s.kernel_error, s.mesh_error, s.cutoff_error, s.contour_error,
                        s.quantization_residual, s.hole_error, s.source_quadrature_error, s.source_tail_bound})
  {
    EXPECT_TRUE(x.initialized());
    EXPECT_EQ(uni20::precision_of(x), p);
  }
  for (auto const* x : {&s.scaling_function, &s.casimir_energy, &s.rapidity})
  {
    EXPECT_EQ(x->has_value(), s.converged);
    if (*x) EXPECT_EQ(uni20::precision_of(**x), p);
  }
  EXPECT_FALSE(s.effective_central_charge);
  EXPECT_EQ(s.converged, s.status == sg::VacuumStatus::converged);
}
auto fields(State const& s)
{
  return std::tie(s.mass, s.length, s.coupling, s.scaled_length, s.contour_shift, s.verification_contour_shift,
                  s.cutoff, s.nonlinear_residual, s.kernel_error, s.mesh_error, s.cutoff_error, s.contour_error,
                  s.scaling_function, s.casimir_energy, s.effective_central_charge, s.intervals, s.iterations,
                  s.kernel_evaluations, s.cutoffs, s.converged, s.status, s.rapidity, s.number, s.quantization_residual,
                  s.hole_error, s.source_quadrature_error, s.source_tail_bound, s.root_iterations,
                  s.source_evaluations);
}
TEST(SineGordonExcitedMpfr, FreeDiracLevelsGainAccuracyBeyondBinary128)
{
  auto hp = P::bits(768);
  R previous(1, hp), pi = uni20::pi<R>.at(hp);
  auto ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256, 333})
  {
    SCOPED_TRACE(bits);
    auto p = P::bits(bits);
    sg::TwoSolitonOptions<R> options(p);
    options.max_intervals = 8192;
    options.max_root_iterations = 8192;
    R worst(0, hp);
    for (auto const& reference : {sg_vacuum_reference::free_dirac[1]})
      for (int twice : {1, 3})
      {
        SCOPED_TRACE(reference.u);
        SCOPED_TRACE(twice);
        auto s = sg::two_soliton_level(R(1, p), R(reference.u, p), R(1, p), uni20::from_twice(twice), options);
        ASSERT_TRUE(s.converged) << int(s.status) << " mesh " << s.intervals;
        precision(s, p);
        R u(reference.u, hp), t = pi * twice;
        R expected = R(reference.y, hp) + 2 * math::sqrt(u * u + t * t);
        worst = std::max(worst, math::abs(s.scaling_function->at(hp) - expected));
        EXPECT_LT(math::abs(s.rapidity->at(hp) - math::asinh(t / u)), options.tolerance.at(hp));
        EXPECT_LE(s.quantization_residual, options.tolerance / 8);
        EXPECT_LE(s.hole_error, options.tolerance / 8);
        EXPECT_EQ(s.iterations, 0u);
        EXPECT_GT(s.root_iterations, 0u);
        EXPECT_EQ(s.source_evaluations, 0u);
      }
    EXPECT_LT(worst, 4 * options.tolerance.at(hp));
    EXPECT_LT(worst, previous * R("1e-15", hp));
    if (bits >= 160) EXPECT_LT(worst, R("1e-38", hp));
    if (bits >= 256) EXPECT_LT(worst, R("1e-65", hp));
    previous = worst;
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
TEST(SineGordonExcitedMpfr, FreeDiracExtremeScalesAt333Bits)
{
  auto p = P::bits(333), hp = P::bits(768);
  R pi = uni20::pi<R>.at(hp);
  sg::TwoSolitonOptions<R> options(p);
  options.max_intervals = 8192;
  options.max_root_iterations = 8192;
  for (std::size_t index : {0u, 2u})
    for (int twice : {1, 3})
    {
      auto const& reference = sg_vacuum_reference::free_dirac[index];
      SCOPED_TRACE(reference.u);
      SCOPED_TRACE(twice);
      auto s = sg::two_soliton_level(R(1, p), R(reference.u, p), R(1, p), uni20::from_twice(twice), options);
      ASSERT_TRUE(s.converged) << int(s.status) << " trials " << s.root_iterations;
      precision(s, p);
      R u(reference.u, hp), target = pi * twice;
      R expected = R(reference.y, hp) + 2 * math::sqrt(u * u + target * target);
      EXPECT_LT(math::abs(s.scaling_function->at(hp) - expected), 4 * options.tolerance.at(hp));
      EXPECT_LT(math::abs(s.rapidity->at(hp) - math::asinh(target / u)), options.tolerance.at(hp));
    }
}
TEST(SineGordonExcitedMpfr, SourcesAgainstIndependentClosedKernel)
{
  auto hp = P::bits(768);
  R zero(0, hp), one(1, hp), pi = uni20::pi<R>.at(hp), previous(1, hp);
  C izero(zero, zero), ii(zero, one);
  auto kernel = [&](C const& z) { return z == izero ? C(1 / (2 * pi * pi), zero) : z / (2 * pi * pi * math::sinh(z)); };
  auto a = bethe::detail::gauss_legendre<R>(192, hp), b = bethe::detail::gauss_legendre<R>(256, hp);
  auto phase = [&](C const& z, auto const& rule) {
    bethe::detail::CompensatedSum<C> sum(izero);
    for (std::size_t i = 0; i < rule.x.size(); ++i)
      sum.add(pi * z * rule.w[i] * kernel(z * (1 + rule.x[i]) / 2));
    return sum.value();
  };
  std::vector<C> source, counting;
  for (int j = 0; j <= 8; ++j)
  {
    R x = R(3, hp) * (j - 4) / 4;
    C z(x, R("0.5", hp)), shift(R("0.3", hp), zero);
    auto pa = phase(z - shift, a) + phase(z + shift, a), pb = phase(z - shift, b) + phase(z + shift, b);
    ASSERT_LT(math::abs(pa - pb), R("1e-110", hp));
    source.push_back(-ii * pb);
    counting.push_back(kernel(C(R("0.3", hp) - x, R("-0.5", hp))));
  }
  C z(R("0.6", hp), zero);
  R expected_phase = phase(z, b).real(), expected_derivative = 2 * pi * kernel(z).real();
  for (int bits : {80, 160, 256, 333})
  {
    SCOPED_TRACE(bits);
    auto p = P::bits(bits);
    sg::KernelOptions<R> o(p);
    o.max_evaluations = 1000000;
    std::size_t count = 0;
    auto t = sg::detail::hole_table(R(2, p), R("0.5", p), R(3, p), 8, R("0.3", p), o, count);
    ASSERT_TRUE(t.converged) << count;
    EXPECT_LE(count, o.max_evaluations);
    R worst = std::max(math::abs(t.phase.at(hp) - expected_phase),
                       math::abs(t.phase_derivative.at(hp) - expected_derivative));
    for (auto const& x : {t.phase, t.phase_derivative, t.error, t.quadrature_error, t.tail_bound})
      EXPECT_EQ(uni20::precision_of(x), p);
    for (int j = 0; j <= 8; ++j)
    {
      EXPECT_EQ(uni20::precision_of(t.source[j]), p);
      EXPECT_EQ(uni20::precision_of(t.counting[j]), p);
      worst =
          std::max({worst, math::abs(t.source[j].at(hp) - source[j]), math::abs(t.counting[j].at(hp) - counting[j])});
    }
    EXPECT_LT(worst, 8 * o.tolerance.at(hp));
    EXPECT_LT(worst, previous * R("1e-15", hp));
    if (bits >= 256) EXPECT_LT(worst, R("1e-68", hp));
    previous = worst;
  }
}
TEST(SineGordonExcitedMpfr, FailuresAndInvalidPrecision)
{
  auto p = P::bits(160);
  R one(1, p);
  sg::TwoSolitonOptions<R> c(p);
  c.tolerance = R("1e-7", p);
  for (int which = 0; which < 6; ++which)
  {
    SCOPED_TRACE(which);
    auto o = c;
    auto expected = sg::VacuumStatus::iteration_limit;
    if (which == 0) o.max_root_iterations = 0;
    if (which == 1) o.max_iterations = 0;
    if (which == 2)
    {
      o.kernel.max_evaluations = 0;
      expected = sg::VacuumStatus::kernel_limit;
    }
    if (which == 3)
    {
      o.max_cutoffs = 0;
      expected = sg::VacuumStatus::cutoff_limit;
    }
    if (which == 4)
    {
      o.max_intervals = o.initial_intervals;
      expected = sg::VacuumStatus::mesh_limit;
    }
    if (which == 5) o.max_root_iterations = 1;
    auto s = sg::two_soliton_level(one, one, R(which >= 4 ? 1 : 2, p), uni20::from_twice(1), o);
    EXPECT_EQ(s.status, expected);
    precision(s, p);
  }
  EXPECT_THROW((void)sg::TwoSolitonOptions<R>(P::exact()), std::invalid_argument);
  c.max_cutoffs = 0;
  for (auto const& bad : {R(uni20::uninitialized), R(1), R(1, P::bits(80)), R("nan", p)})
  {
    EXPECT_THROW(sg::two_soliton_level(bad, one, one, uni20::from_twice(1), c), std::invalid_argument);
    EXPECT_THROW(sg::two_soliton_level(one, bad, one, uni20::from_twice(1), c), std::invalid_argument);
    EXPECT_THROW(sg::two_soliton_level(one, one, bad, uni20::from_twice(1), c), std::invalid_argument);
    auto o = c;
    o.tolerance = bad;
    EXPECT_THROW(sg::two_soliton_level(one, one, one, uni20::from_twice(1), o), std::invalid_argument);
    o = c;
    o.contour_shift = bad;
    EXPECT_THROW(sg::two_soliton_level(one, one, one, uni20::from_twice(1), o), std::invalid_argument);
    o = c;
    o.initial_cutoff = bad;
    EXPECT_THROW(sg::two_soliton_level(one, one, one, uni20::from_twice(1), o), std::invalid_argument);
  }
  c.kernel = sg::KernelOptions<R>(P::bits(80));
  EXPECT_THROW(sg::two_soliton_level(one, one, one, uni20::from_twice(1), c), std::invalid_argument);
}
TEST(SineGordonExcitedMpfr, OwningMixedPrecisionResults)
{
  struct Query
  {
      R length, coupling;
      sg::TwoSolitonOptions<R> controls;
  };
  std::vector<Query> qs;
  for (int bits : {80, 160})
  {
    auto p = P::bits(bits);
    sg::TwoSolitonOptions<R> c(p);
    qs.push_back({R(1, p), R(1, p), c});
    c.kernel.max_evaluations = 0;
    qs.push_back({R(1, p), R(2, p), c});
  }
  auto evaluate = [](Query const& q) {
    return sg::two_soliton_level(R(1, q.controls.precision()), q.length, q.coupling, uni20::from_twice(1), q.controls);
  };
  std::vector<State> expected;
  for (auto const& q : qs)
    expected.push_back(evaluate(q));
  auto check = [&](uni20::async::IAsyncScheduler& scheduler) {
    std::vector<State> actual;
    auto final = [&] {
      auto owned = qs;
      auto count = owned.size();
      return bethe::evaluate_observable_grid_async(
          count, [owned = std::move(owned), evaluate](std::size_t i) { return evaluate(owned[i]); },
          [&](std::size_t i, State s) {
            EXPECT_EQ(i, actual.size());
            actual.push_back(std::move(s));
          },
          {.scheduler = &scheduler, .max_in_flight = 2});
    }();
    EXPECT_EQ(final.get_wait(scheduler), qs.size());
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      EXPECT_TRUE(fields(actual[i]) == fields(expected[i]));
      EXPECT_EQ(actual[i].converged, i % 2 == 0);
      precision(actual[i], qs[i].controls.precision());
    }
  };
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  check(deferred);
  for (int n : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(n);
    check(scheduler);
  }
}
TEST(SineGordonExcitedMpfr, InteractingLevelsAndScaling)
{
  auto p = P::bits(160);
  sg::TwoSolitonOptions<R> c(p);
  c.tolerance = R("1e-7", p);
  for (int twice : {1, 3})
  {
    auto level = sg::two_soliton_level(R(1, p), R(1, p), R(2, p), uni20::from_twice(twice), c);
    ASSERT_TRUE(level.converged) << int(level.status);
    precision(level, p);
    R expected(twice == 1 ? "5.11055662981005" : "17.14609634255262", p);
    EXPECT_LT(math::abs(*level.scaling_function - expected), c.tolerance);
    auto scaled = sg::two_soliton_level(R(4, p), R("0.25", p), R(2, p), uni20::from_twice(twice), c);
    ASSERT_TRUE(scaled.converged);
    EXPECT_EQ(level.scaling_function, scaled.scaling_function);
    EXPECT_EQ(4 * *level.casimir_energy, *scaled.casimir_energy);
  }
}
} // namespace
