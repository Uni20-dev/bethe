// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "exact_spectrum.hpp"
#include "test_support.hpp"
#include <bethe/experimental/xxx_boundary_dynamics.hpp>
#include <condition_variable>
#include <thread>
#include <uni20/async/tbb_scheduler.hpp>

namespace boundary = bethe::experimental::xxx::open;
namespace model = bethe::heisenberg::open;
using namespace uni20::async;

namespace
{
template <uni20::Real R> void check_scan(bethe::SolverOptions<R> const& options)
{
  R const zero = bethe::detail::solver_zero(options), one = uni20::scalar_like(zero, 1);
  R const eps = uni20::numeric_limits<R>::epsilon(zero), tol = 524288 * eps;
  DebugScheduler serial;
  TbbScheduler parallel(2);
  for (unsigned n : {2, 4, 6, 8})
  {
    SCOPED_TRACE(n);
    auto reference = boundary::boundary_real_triplet_spectrum<R>(n, options, 1000, {&serial});
    auto pending = boundary::boundary_real_triplet_spectrum_async<R>(n, options, 1000, {&parallel, 1, 2});
    auto const& result = pending.get_wait(parallel);
    ASSERT_TRUE(result.converged());
    ASSERT_EQ(result.records.size(), bethe::heisenberg::real_excitation_count(n, uni20::half_int{1}));
    ASSERT_EQ(result.lines.size(), result.records.size());
    ASSERT_EQ(result.end_to_end_lines.size(), result.lines.size());
    EXPECT_EQ(result.ground_reflection_parity, n / 2 % 2 ? -1 : 1);
    for (std::size_t i = 0; i < result.records.size(); ++i)
    {
      auto const& row = result.records[i];
      ASSERT_TRUE(row.state && row.factor && row.line);
      EXPECT_EQ(row.candidate_id, reference.records[i].candidate_id);
      EXPECT_EQ(row.quantum_numbers, row.state->quantum_numbers);
      if (i) EXPECT_LT(result.records[i - 1].quantum_numbers, row.quantum_numbers);
      EXPECT_TRUE(row.diagnostic.empty());
      EXPECT_EQ(row.line->state_id, i + 1);
      EXPECT_EQ(row.line->gap, reference.lines[i].gap);
      EXPECT_EQ(row.line->weight, reference.lines[i].weight);
      EXPECT_EQ(row.line->weight, *row.factor->weight / 2);
      auto const& cross = result.end_to_end_lines[i];
      EXPECT_EQ(cross.coefficient.real(), result.ground_reflection_parity * row.reflection_parity * row.line->weight);
      EXPECT_EQ(cross.coefficient.imag(), zero);
      EXPECT_EQ(uni20::precision_of(cross.coefficient), uni20::precision_of(zero));
    }
    EXPECT_EQ(result.moments.weight, reference.moments.weight);
    EXPECT_EQ(result.end_to_end_moments.value.real(), reference.end_to_end_moments.value.real());
    EXPECT_EQ(result.end_to_end_moments.value.imag(), reference.end_to_end_moments.value.imag());
    EXPECT_GE(-*result.density.zz() - result.moments.first_moment, -tol);
    if (n <= 4)
    {
      EXPECT_REAL_NEAR(result.moments.weight, one / 4, tol);
      EXPECT_REAL_NEAR(result.end_to_end_moments.value.real(), n == 2 ? -one / 4 : -one / 12, tol);
    }
    else
    {
      ASSERT_TRUE(result.omitted_weight());
      EXPECT_GT(*result.omitted_weight(), zero);
      ASSERT_TRUE(result.omitted_first_moment());
      EXPECT_GT(*result.omitted_first_moment(), zero);
    }
    for (int t : {0, 1, 10, 100})
    {
      R const time = uni20::scalar_like(zero, t);
      auto const c = bethe::local_time_correlation<R>(result.lines, time);
      EXPECT_LE(uni20::math::abs(c), result.moments.weight + tol);
      auto const reverse = bethe::local_time_correlation<R>(result.lines, -time);
      EXPECT_REAL_NEAR(uni20::math::abs(reverse - uni20::conj(c)), zero, tol);
      if (n == 2)
        EXPECT_REAL_NEAR(
            uni20::math::abs(c - uni20::complex<R>(uni20::math::cos(time), -uni20::math::sin(time)) / (4 * one)), zero,
            tol);
    }
  }
}

AsyncTask count_lines(ReadBuffer<boundary::BoundarySpectrum<double>> input, WriteBuffer<std::size_t> output)
{
  auto const& spectrum = co_await input;
  co_await output = spectrum.lines.size();
}
} // namespace

template <typename R> class ExperimentalBoundaryDynamics : public ::testing::Test {};
TYPED_TEST_SUITE(ExperimentalBoundaryDynamics, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(ExperimentalBoundaryDynamics, SerialConcurrentNormAndEndpointRelations) { check_scan<TypeParam>({}); }

TEST(ExperimentalBoundaryDynamicsED, LocalAndEndToEndTimeAndFullMoments)
{
  for (unsigned n : {2, 4, 6, 8})
  {
    SCOPED_TRACE(n);
    auto result = boundary::boundary_real_triplet_spectrum<double>(n);
    ASSERT_TRUE(result.converged());
    auto gs = test_support::exact_eigensystem(n, n / 2, 0, false);
    auto ex = test_support::exact_eigensystem(n, n / 2 - 1, 0, false);
    auto g = std::min_element(gs.energies.begin(), gs.energies.end()) - gs.energies.begin();
    std::vector<double> left(1U << n, 0), right(1U << n, 0);
    double exact_static = 0;
    for (std::size_t i = 0; i < gs.basis.size(); ++i)
    {
      auto bits = gs.basis[i];
      double a = gs.vectors[i * gs.basis.size() + g];
      if (bits & 1U) left[bits ^ 1U] = a;
      if (bits & (1U << (n - 1))) right[bits ^ (1U << (n - 1))] = a;
      exact_static += a * a * ((bits & 1U) ? -.5 : .5) * ((bits & (1U << (n - 1))) ? -.5 : .5);
    }
    std::vector<bethe::LocalSpectralLine<double>> full, retained;
    std::vector<bethe::CrossSpectralLine<double>> full_cross, retained_cross;
    for (std::size_t k = 0; k < ex.energies.size(); ++k)
    {
      double a = 0, b = 0;
      for (std::size_t i = 0; i < ex.basis.size(); ++i)
      {
        a += ex.vectors[i * ex.basis.size() + k] * left[ex.basis[i]];
        b += ex.vectors[i * ex.basis.size() + k] * right[ex.basis[i]];
      }
      double gap = ex.energies[k] - gs.energies[g];
      full.push_back({k + 1, gap, a * a / 2});
      full_cross.push_back({k + 1, gap, {a * b / 2, 0}});
      if (std::any_of(result.lines.begin(), result.lines.end(),
                      [&](auto const& line) { return std::abs(line.gap - gap) < 1e-9; }))
      {
        retained.push_back(full.back());
        retained_cross.push_back(full_cross.back());
      }
    }
    auto m = bethe::local_spectral_moments<double>(full);
    EXPECT_NEAR(m.weight, .25, 3e-12);
    EXPECT_NEAR(m.first_moment, -*result.density.zz(), 3e-12);
    auto cm = bethe::cross_spectral_moments<double>(full_cross);
    EXPECT_NEAR(cm.value.real(), exact_static, 3e-12);
    if (n > 2) EXPECT_NEAR(cm.first_moment.real(), 0, 3e-12);
    for (double time : {0., .25, 1., 10., 100.})
    {
      auto local = bethe::local_time_correlation<double>(result.lines, time);
      auto cross = bethe::cross_time_correlation<double>(result.end_to_end_lines, time);
      EXPECT_NEAR(std::abs(local - bethe::local_time_correlation<double>(retained, time)), 0, 3e-10);
      EXPECT_NEAR(std::abs(cross - bethe::cross_time_correlation<double>(retained_cross, time)), 0, 3e-10);
      if (auto omitted = result.omitted_weight())
      {
        EXPECT_LE(std::abs(local - bethe::local_time_correlation<double>(full, time)), *omitted + 3e-10);
        EXPECT_LE(std::abs(cross - bethe::cross_time_correlation<double>(full_cross, time)), *omitted + 3e-10);
      }
    }
  }
}

TEST(ExperimentalBoundaryDynamicsExecution, OwnershipFailuresPreflightAndComposition)
{
  for (auto order : {DebugSchedulerOrder::fifo, DebugSchedulerOrder::reverse, DebugSchedulerOrder::random})
  {
    DebugScheduler scheduler({.order = order, .random_seed = 41});
    Async<std::size_t> final;
    {
      bethe::SolverOptions<double> options;
      auto pending = boundary::boundary_real_triplet_spectrum_async<double>(8, options, 100, {&scheduler});
      options.max_iterations = 0;
      scheduler.schedule(count_lines(pending.read(), final.write()));
    }
    EXPECT_EQ(final.get_wait(scheduler), 10u);
  }
  DebugScheduler scheduler;
  EXPECT_THROW((void)boundary::boundary_real_triplet_spectrum_async<double>(5), std::invalid_argument);
  EXPECT_THROW((void)boundary::boundary_real_triplet_spectrum_async<double>(256, {}, 8255), std::length_error);
  EXPECT_THROW((void)boundary::boundary_real_triplet_spectrum_async<double>(8, {}, 100, {&scheduler, 1, 0}),
               std::invalid_argument);
  EXPECT_THROW((void)boundary::boundary_real_triplet_spectrum_async<double>(8, {.residual_tolerance = 0}),
               std::invalid_argument);
  auto failed = boundary::boundary_real_triplet_spectrum<double>(8, {.max_iterations = 0}, 100, {&scheduler});
  EXPECT_FALSE(failed.converged());
  EXPECT_FALSE(failed.omitted_weight());
  EXPECT_FALSE(failed.omitted_first_moment());
  EXPECT_TRUE(failed.lines.empty());
  ASSERT_EQ(failed.records.size(), 10u);
  for (auto const& row : failed.records)
  {
    EXPECT_FALSE(row.state);
    EXPECT_FALSE(row.factor);
    EXPECT_FALSE(row.line);
    EXPECT_FALSE(row.candidate_id.empty());
    EXPECT_FALSE(row.diagnostic.empty());
  }
  auto mixed = boundary::boundary_detail::submit_spectrum<double>(
      4, {}, 100, {&scheduler}, [](auto n, auto labels, auto const& prepared, auto options, std::stop_token stop) {
        if (labels.front() == uni20::half_int{1}) options.max_iterations = 0;
        return boundary::boundary_detail::evaluate(n, std::move(labels), prepared, options, stop);
      });
  auto const& m = mixed.get_wait(scheduler);
  EXPECT_FALSE(m.converged());
  ASSERT_EQ(m.records.size(), 3u);
  ASSERT_TRUE(m.records[0].state);
  EXPECT_FALSE(m.records[0].state->converged);
  EXPECT_FALSE(m.records[0].line);
  EXPECT_EQ(m.lines.size(), 2u);
  auto exceeded = boundary::boundary_real_triplet_spectrum<double>(4);
  exceeded.moments.weight = .251;
  EXPECT_FALSE(exceeded.omitted_weight());
  EXPECT_FALSE(exceeded.converged());
  exceeded.moments.first_moment = 1;
  EXPECT_FALSE(exceeded.omitted_first_moment());
  EXPECT_FALSE(exceeded.sum_rules_consistent());
  std::stop_source stop;
  stop.request_stop();
  auto cancelled =
      boundary::boundary_real_triplet_spectrum_async<double>(4, {}, 100, {&scheduler, 1, 2, stop.get_token()});
  EXPECT_THROW((void)cancelled.get_wait(scheduler), bethe::ExecutionCancelled);
}

TEST(ExperimentalBoundaryDynamicsExecution, FailureAndCancellationDrainAdmittedWork)
{
  TbbScheduler scheduler(2);
  for (bool cancel : {false, true})
  {
    std::mutex mutex;
    std::condition_variable cv;
    bool held = false, failed = false, release = false, finished = false;
    std::stop_source stop;
    auto pending = boundary::boundary_detail::submit_spectrum<double>(
        4, {}, 100, {&scheduler, 1, 2, stop.get_token()},
        [&](auto, auto labels, auto const&, auto const&, std::stop_token) {
          std::unique_lock lock(mutex);
          if (labels.front() == uni20::half_int{1})
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
          return boundary::BoundaryRecord<double>{};
        });
    std::jthread observer([&] {
      std::unique_lock lock(mutex);
      EXPECT_TRUE(cv.wait_for(lock, std::chrono::seconds(5), [&] { return held && failed; }));
      auto reader = pending.read();
      EXPECT_FALSE(reader.await_ready());
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
TEST(ExperimentalBoundaryDynamicsMpfr, OwnedExplicitPrecisionWorkers)
{
  auto ambient = mpfr_get_default_prec();
  for (int bits : {128, 256})
    check_scan(bethe::SolverOptions<uni20::mpreal>(uni20::Precision::bits(bits)));
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
#endif
