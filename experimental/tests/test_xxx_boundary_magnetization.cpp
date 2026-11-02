// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "exact_spectrum.hpp"
#include "test_support.hpp"
#include "xxx_open_wave_oracle.hpp"
#include <bethe/experimental/xxx_boundary_magnetization.hpp>
#include <condition_variable>
#include <thread>
#include <uni20/async/tbb_scheduler.hpp>

namespace boundary = bethe::experimental::xxx::open;
namespace model = bethe::heisenberg::open;
using namespace uni20::async;
namespace
{
template <uni20::Real R> void check_waves(bethe::SolverOptions<R> const& solver)
{
  using C = uni20::complex<R>;
  R const zero = bethe::detail::solver_zero(solver), one = uni20::scalar_like(zero, 1);
  R const tol = 4194304 * uni20::numeric_limits<R>::epsilon(zero);
  DebugScheduler serial;
  TbbScheduler parallel(2);
  for (unsigned n : {2, 4, 6, 8})
  {
    SCOPED_TRACE(n);
    auto result = boundary::boundary_real_triplet_magnetization<R>(n, solver, 250000, {&serial});
    auto pending = boundary::boundary_real_triplet_magnetization_async<R>(n, solver, 250000, {&parallel, 1, 2});
    auto const& concurrent = pending.get_wait(parallel);
    ASSERT_TRUE(result.converged());
    ASSERT_TRUE(concurrent.converged());
    ASSERT_TRUE(result.projected_norm);
    ASSERT_EQ(result.pairs.size(), result.states.size() * (result.states.size() + 1) / 2);
    EXPECT_EQ(result.projected_norm, concurrent.projected_norm);
    boundary::BoundarySpinFormFactors<R> calculator(n, result.spectrum.ground_state);
    auto ground = test_support::open_aba_wave<R>(n, result.spectrum.ground_state.rapidities, zero);
    std::vector<C> raised(1U << n, C(zero, zero));
    for (unsigned bits = 0; bits < (1U << n); ++bits)
      if (bits & 1U) raised[bits ^ 1U] = ground[bits];
    std::vector<std::vector<C>> waves;
    std::vector<C> couplings;
    std::vector<boundary::BoundarySpinState<R>> states;
    for (auto const& row : result.spectrum.records)
    {
      ASSERT_TRUE(row.state);
      states.push_back(calculator.prepare(*row.state));
      ASSERT_TRUE(states.back().ready());
      EXPECT_EQ(uni20::precision_of(*states.back().diagonal()), uni20::precision_of(zero));
      EXPECT_EQ(uni20::precision_of(*states.back().raising_amplitude()), uni20::precision_of(zero));
      ASSERT_TRUE(states.back().solve_backward_error());
      EXPECT_LE(*states.back().solve_backward_error(), 4096 * n * uni20::numeric_limits<R>::epsilon(zero));
      waves.push_back(test_support::open_aba_wave<R>(n, row.state->rapidities, zero));
      C coupling(zero, zero);
      for (unsigned bits = 0; bits < (1U << n); ++bits)
        coupling += uni20::conj(waves.back()[bits]) * raised[bits];
      couplings.push_back(coupling);
    }
    bool negative = false;
    for (std::size_t k = 0; k < result.pairs.size(); ++k)
    {
      auto const& pair = result.pairs[k];
      auto const i = pair.left_state - 1, j = pair.right_state - 1;
      ASSERT_TRUE(pair.term && pair.element && pair.element->value);
      EXPECT_TRUE(pair.diagnostic.empty());
      EXPECT_EQ(pair.left_state, concurrent.pairs[k].left_state);
      EXPECT_EQ(pair.right_state, concurrent.pairs[k].right_state);
      EXPECT_EQ(pair.term->frequency, concurrent.pairs[k].term->frequency);
      EXPECT_EQ(pair.term->coefficient, concurrent.pairs[k].term->coefficient);
      auto element = calculator(states[i], states[j]), reverse = calculator(states[j], states[i]);
      ASSERT_TRUE(element.value && reverse.value);
      EXPECT_REAL_NEAR(*element.value, *reverse.value, tol);
      EXPECT_EQ(*pair.element->value, *element.value);
      C spin(zero, zero);
      for (unsigned bits = 0; bits < (1U << n); ++bits)
        spin += uni20::conj(waves[i][bits]) * waves[j][bits] * ((bits & 1U) ? -one / 2 : one / 2);
      C const expected = uni20::conj(couplings[i]) * couplings[j] * spin;
      R const actual = *states[i].raising_amplitude() * *states[j].raising_amplitude() * *element.value;
      EXPECT_REAL_NEAR(uni20::math::abs(expected - actual), zero, tol);
      negative = negative || pair.term->coefficient < zero;
    }
    if (n >= 4) EXPECT_TRUE(negative); // A squared-weight-only implementation cannot pass.
    for (int t : {0, 1, 5, 10, 50, 100})
    {
      R const time = uni20::scalar_like(zero, t);
      std::vector<C> projected(1U << n, C(zero, zero));
      for (std::size_t i = 0; i < waves.size(); ++i)
      {
        R const phase = -result.spectrum.records[i].line->gap * time;
        C const rotation(uni20::math::cos(phase), uni20::math::sin(phase));
        for (unsigned bits = 0; bits < (1U << n); ++bits)
          projected[bits] += uni20::math::sqrt(2 * one) * rotation * couplings[i] * waves[i][bits];
      }
      bethe::detail::CompensatedSum<R> norm(zero), first(zero), last(zero), total(zero);
      for (unsigned bits = 0; bits < (1U << n); ++bits)
      {
        R const probability =
            projected[bits].real() * projected[bits].real() + projected[bits].imag() * projected[bits].imag();
        norm.add(probability);
        first.add(((bits & 1U) ? -one / 2 : one / 2) * probability);
        last.add(((bits & (1U << (n - 1))) ? -one / 2 : one / 2) * probability);
        total.add((uni20::scalar_like(zero, n) / 2 - uni20::scalar_like(zero, std::popcount(bits))) * probability);
      }
      EXPECT_REAL_NEAR(norm.value(), *result.projected_norm, tol);
      EXPECT_REAL_NEAR(total.value(), *result.projected_norm, tol);
      EXPECT_REAL_NEAR(result.at(time), first.value(), tol);
      EXPECT_REAL_NEAR(result.at(time, true), last.value(), tol);
      EXPECT_EQ(result.at(time), concurrent.at(time));
      EXPECT_EQ(result.at(time), result.at(-time));
      EXPECT_LE(uni20::math::abs(result.at(time)), *result.projected_norm / 2 + tol);
      if (n == 2)
      {
        EXPECT_REAL_NEAR(result.at(time), one / 2, tol);
        EXPECT_REAL_NEAR(result.at(time, true), one / 2, tol);
      }
    }
  }
}
AsyncTask final_value(ReadBuffer<boundary::BoundaryMagnetization<double>> input, WriteBuffer<double> output)
{
  auto const& result = co_await input;
  co_await output = result.at(0);
}
} // namespace

template <typename R> class ExperimentalBoundaryMagnetization : public ::testing::Test {};
TYPED_TEST_SUITE(ExperimentalBoundaryMagnetization, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(ExperimentalBoundaryMagnetization, SignedCyclicProductsWavesAndConcurrency) { check_waves<TypeParam>({}); }

TEST(ExperimentalBoundaryMagnetizationED, CompleteEvolutionAndProjectionBound)
{
  using C = std::complex<double>;
  for (unsigned n : {2, 4, 6, 8})
  {
    auto result = boundary::boundary_real_triplet_magnetization<double>(n);
    ASSERT_TRUE(result.converged());
    auto gs = test_support::exact_eigensystem(n, n / 2, 0, false);
    auto ex = test_support::exact_eigensystem(n, n / 2 - 1, 0, false);
    auto g = std::min_element(gs.energies.begin(), gs.energies.end()) - gs.energies.begin();
    auto gd = gs.basis.size(), ed = ex.basis.size();
    std::vector<double> raised(1U << n), couplings(ed);
    double static_end = 0;
    for (std::size_t i = 0; i < gd; ++i)
    {
      if (gs.basis[i] & 1U) raised[gs.basis[i] ^ 1U] = gs.vectors[i * gd + g] * std::sqrt(2.);
      static_end += gs.vectors[i * gd + g] * gs.vectors[i * gd + g] * ((gs.basis[i] & 1U) ? -.5 : .5) *
                    ((gs.basis[i] & (1U << (n - 1))) ? -.5 : .5);
    }
    for (std::size_t k = 0; k < ed; ++k)
      for (std::size_t i = 0; i < ed; ++i)
        couplings[k] += ex.vectors[i * ed + k] * raised[ex.basis[i]];
    for (double time : {0., 1., 10., 100.})
    {
      std::vector<C> wave(ed);
      for (std::size_t k = 0; k < ed; ++k)
      {
        C const phase = std::exp(C(0, -(ex.energies[k] - gs.energies[g]) * time));
        for (std::size_t i = 0; i < ed; ++i)
          wave[i] += phase * couplings[k] * ex.vectors[i * ed + k];
      }
      double first = 0, last = 0;
      for (std::size_t i = 0; i < ed; ++i)
      {
        first += std::norm(wave[i]) * ((ex.basis[i] & 1U) ? -.5 : .5);
        last += std::norm(wave[i]) * ((ex.basis[i] & (1U << (n - 1))) ? -.5 : .5);
      }
      if (n <= 4)
      {
        EXPECT_NEAR(result.at(time), first, 3e-11);
        EXPECT_NEAR(result.at(time, true), last, 3e-11);
      }
      else
      {
        ASSERT_TRUE(result.truncation_bound());
        EXPECT_LE(std::abs(result.at(time) - first), *result.truncation_bound() + 3e-11);
        EXPECT_LE(std::abs(result.at(time, true) - last), *result.truncation_bound() + 3e-11);
      }
    }
    if (n <= 4)
    {
      EXPECT_NEAR(result.at(0), .5, 3e-11);
      EXPECT_NEAR(result.at(0, true), -2 * static_end, 3e-11);
    }
  }
}

TEST(ExperimentalBoundaryMagnetizationExecution, OwnedPreparationPreflightFailuresAndComposition)
{
  for (auto order : {DebugSchedulerOrder::fifo, DebugSchedulerOrder::reverse, DebugSchedulerOrder::random})
  {
    DebugScheduler scheduler({.order = order, .random_seed = 91});
    Async<double> final;
    {
      auto pending = boundary::boundary_real_triplet_magnetization_async<double>(4, {}, 6, {&scheduler});
      scheduler.schedule(final_value(pending.read(), final.write()));
    }
    EXPECT_NEAR(final.get_wait(scheduler), .5, 3e-12);
  }
  EXPECT_THROW((void)boundary::boundary_real_triplet_magnetization_async<double>(5), std::invalid_argument);
  EXPECT_THROW((void)boundary::boundary_real_triplet_magnetization_async<double>(64, {}, 139655), std::length_error);
  EXPECT_THROW((void)boundary::boundary_real_triplet_magnetization_async<double>(256), std::length_error);
  EXPECT_THROW((void)boundary::boundary_real_triplet_magnetization_async<double>(4, {}, 0), std::length_error);
  EXPECT_THROW((void)boundary::boundary_real_triplet_magnetization_async<double>(4, {.residual_tolerance = 0}),
               std::invalid_argument);
  auto failed = boundary::boundary_real_triplet_magnetization<double>(8, {.max_iterations = 0});
  EXPECT_FALSE(failed.converged());
  EXPECT_FALSE(failed.matrix_ready);
  EXPECT_FALSE(failed.omitted_norm());
  EXPECT_FALSE(failed.projected_norm);
  ASSERT_EQ(failed.pairs.size(), 55u);
  for (auto const& p : failed.pairs)
  {
    EXPECT_FALSE(p.element);
    EXPECT_FALSE(p.term);
    EXPECT_FALSE(p.diagnostic.empty());
  }
  EXPECT_THROW((void)failed.at(0), std::logic_error);
  auto ground = model::ground_state<double>(4);
  boundary::BoundarySpinFormFactors<double> calculator(4, ground), other(4, ground);
  auto excited = model::sector_ground_state<double>(4, uni20::half_int{1});
  auto prepared = calculator.prepare(excited);
  ASSERT_TRUE(prepared.ready());
  auto copy = calculator;
  EXPECT_TRUE(copy(prepared, prepared).converged());
  EXPECT_THROW((void)other(prepared, prepared), std::invalid_argument);
  auto moved = std::move(prepared);
  EXPECT_FALSE(prepared.ready());
  EXPECT_TRUE(moved.ready());
  EXPECT_FALSE(calculator(prepared, moved).value);
  excited.rapidities[0] += .001;
  EXPECT_FALSE(calculator.prepare(excited).ready());
  auto dimer = boundary::boundary_real_triplet_magnetization<double>(2);
  EXPECT_THROW((void)dimer.at(std::numeric_limits<double>::infinity()), std::invalid_argument);
  dimer.projected_norm = 1.01;
  EXPECT_FALSE(dimer.omitted_norm());
  EXPECT_FALSE(dimer.truncation_bound());
  auto unresolved = boundary::boundary_detail::submit_magnetization<double>(
      4, {}, 6, {},
      [](auto const&, auto const&, auto const&, std::stop_token) { return boundary::BoundarySpinElement<double>{}; });
  auto const& bad = unresolved.get_wait();
  EXPECT_FALSE(bad.converged());
  EXPECT_FALSE(bad.matrix_ready);
  ASSERT_EQ(bad.pairs.size(), 6u);
  for (auto const& p : bad.pairs)
  {
    ASSERT_TRUE(p.element);
    EXPECT_FALSE(p.term);
    EXPECT_FALSE(p.diagnostic.empty());
  }
  EXPECT_THROW((void)bad.at(0), std::logic_error);
  std::stop_source stop;
  stop.request_stop();
  auto cancelled =
      boundary::boundary_real_triplet_magnetization_async<double>(4, {}, 6, {nullptr, 1, 2, stop.get_token()});
  EXPECT_THROW((void)cancelled.get_wait(), bethe::ExecutionCancelled);
}

TEST(ExperimentalBoundaryMagnetizationExecution, SpectralFamilySurvivesIndependentDensityFailure)
{
  auto reference = boundary::boundary_real_triplet_magnetization<double>(4);
  ASSERT_TRUE(reference.converged());
  for (bool rejected_density : {false, true})
  {
    SCOPED_TRACE(rejected_density);
    auto spectrum = reference.spectrum;
    if (rejected_density)
      spectrum.density = {}; // Failed independent ground-density calculation.
    else
      spectrum.density.bond_correlation = 0.; // Density-based first-moment diagnostic fails.
    ASSERT_TRUE(spectrum.family_converged());
    ASSERT_FALSE(spectrum.sum_rules_consistent());
    ASSERT_FALSE(spectrum.converged());
    DebugScheduler scheduler;
    Async<boundary::BoundarySpectrum<double>> input(std::move(spectrum));
    Async<boundary::BoundaryMagnetization<double>> output;
    std::size_t evaluated = 0;
    auto task = boundary::boundary_detail::magnetization_task<double>(
        scheduler, input.read(), reference.pairs.size(), {&scheduler},
        [&](auto const& calculator, auto const& x, auto const& y, std::stop_token) {
          ++evaluated;
          return calculator(x, y);
        },
        output.write());
    bethe::detail::schedule_owned(scheduler, std::move(task));
    auto const& result = output.get_wait(scheduler);
    ASSERT_TRUE(result.converged());
    EXPECT_TRUE(result.matrix_ready);
    EXPECT_FALSE(result.spectrum.converged()); // The independent diagnostic remains visible.
    EXPECT_EQ(result.spectrum.density.converged(), !rejected_density);
    EXPECT_EQ(evaluated, reference.pairs.size());
    EXPECT_EQ(result.projected_norm, reference.projected_norm);
    for (double time : {0., 1., 10., 100.})
    {
      EXPECT_EQ(result.at(time), reference.at(time));
      EXPECT_EQ(result.at(time, true), reference.at(time, true));
    }
  }
}

TEST(ExperimentalBoundaryMagnetizationExecution, PairFailureAndCancellationDrainAdmittedWork)
{
  TbbScheduler scheduler(2);
  for (bool cancel : {false, true})
  {
    std::mutex mutex;
    std::condition_variable cv;
    bool held = false, failed = false, release = false, finished = false;
    std::stop_source stop;
    std::size_t index = 0;
    auto pending = boundary::boundary_detail::submit_magnetization<double>(
        4, {}, 6, {&scheduler, 1, 2, stop.get_token()}, [&](auto const&, auto const&, auto const&, std::stop_token) {
          std::unique_lock lock(mutex);
          if (index++ == 0)
          {
            held = true;
            cv.notify_all();
            EXPECT_TRUE(cv.wait_for(lock, std::chrono::seconds(5), [&] { return release; }));
            finished = true;
          }
          else
          {
            EXPECT_TRUE(cv.wait_for(lock, std::chrono::seconds(5), [&] { return held; }));
            failed = true;
            cv.notify_all();
            if (cancel)
              stop.request_stop();
            else
              throw std::bad_alloc{};
          }
          return boundary::BoundarySpinElement<double>{};
        });
    std::jthread observer([&] {
      std::unique_lock lock(mutex);
      EXPECT_TRUE(cv.wait_for(lock, std::chrono::seconds(5), [&] { return held && failed; }));
      EXPECT_FALSE(pending.read().await_ready());
      release = true;
      cv.notify_all();
    });
    if (cancel)
      EXPECT_THROW((void)pending.get_wait(scheduler), bethe::ExecutionCancelled);
    else
      EXPECT_THROW((void)pending.get_wait(scheduler), std::bad_alloc);
    observer.join();
    EXPECT_TRUE(finished);
  }
}
#if UNI20_HAS_MPREAL
TEST(ExperimentalBoundaryMagnetizationMpfr, ExplicitPrecisionAndMismatch)
{
  auto ambient = mpfr_get_default_prec();
  for (int bits : {128, 256})
    check_waves(bethe::SolverOptions<uni20::mpreal>(uni20::Precision::bits(bits)));
  auto result = boundary::boundary_real_triplet_magnetization<uni20::mpreal>(
      4, bethe::SolverOptions<uni20::mpreal>(uni20::Precision::bits(160)));
  EXPECT_THROW((void)result.at(uni20::mpreal(0, uni20::Precision::bits(80))), std::invalid_argument);
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
#endif
