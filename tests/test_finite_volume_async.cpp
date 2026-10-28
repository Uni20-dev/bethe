// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "../apps/independent-solves.hpp"
#include "test_support.hpp"
#include <bethe/lee_yang_excited.hpp>
#include <bethe/sine_gordon_excited.hpp>
#include <future>
#include <tuple>

namespace
{
namespace ly = bethe::lee_yang;
namespace sg = bethe::sine_gordon;
using namespace uni20::async;

template <class R> auto diagnostics(ly::Diagnostics<R> const& s)
{
  return std::tie(s.mass, s.length, s.scaled_length, s.cutoff, s.nonlinear_residual, s.nonlinear_error, s.mesh_error,
                  s.cutoff_error, s.direct_tail_bound, s.intervals, s.iterations, s.cutoffs, s.kernel_products,
                  s.converged, s.status);
}
template <class R> void same(ly::State<R> const& a, ly::State<R> const& b)
{
  EXPECT_TRUE(diagnostics(a) == diagnostics(b));
  EXPECT_TRUE(std::tie(a.scaling_function, a.casimir_energy, a.effective_central_charge) ==
              std::tie(b.scaling_function, b.casimir_energy, b.effective_central_charge));
}
template <class R> void same(ly::OneParticleState<R> const& a, ly::OneParticleState<R> const& b)
{
  EXPECT_TRUE(diagnostics(a) == diagnostics(b));
  EXPECT_TRUE(std::tie(a.scaling_function, a.casimir_energy, a.beta, a.pole_displacement, a.quantization_residual,
                       a.source_error, a.root_iterations) == std::tie(b.scaling_function, b.casimir_energy, b.beta,
                                                                      b.pole_displacement, b.quantization_residual,
                                                                      b.source_error, b.root_iterations));
}
template <class R> void same(sg::VacuumState<R> const& a, sg::VacuumState<R> const& b)
{
  EXPECT_TRUE(std::tie(a.mass, a.length, a.coupling, a.scaled_length, a.contour_shift, a.verification_contour_shift,
                       a.cutoff, a.scaling_function, a.casimir_energy, a.effective_central_charge, a.nonlinear_residual,
                       a.kernel_error, a.mesh_error, a.cutoff_error, a.contour_error, a.intervals, a.iterations,
                       a.kernel_evaluations, a.cutoffs, a.converged, a.status) ==
              std::tie(b.mass, b.length, b.coupling, b.scaled_length, b.contour_shift, b.verification_contour_shift,
                       b.cutoff, b.scaling_function, b.casimir_energy, b.effective_central_charge, b.nonlinear_residual,
                       b.kernel_error, b.mesh_error, b.cutoff_error, b.contour_error, b.intervals, b.iterations,
                       b.kernel_evaluations, b.cutoffs, b.converged, b.status));
}
template <class R> void same(sg::TwoSolitonState<R> const& a, sg::TwoSolitonState<R> const& b)
{
  same(static_cast<sg::VacuumState<R> const&>(a), static_cast<sg::VacuumState<R> const&>(b));
  EXPECT_TRUE(std::tie(a.rapidity, a.number, a.quantization_residual, a.hole_error, a.source_quadrature_error,
                       a.source_tail_bound, a.root_iterations, a.source_evaluations) ==
              std::tie(b.rapidity, b.number, b.quantization_residual, b.hole_error, b.source_quadrature_error,
                       b.source_tail_bound, b.root_iterations, b.source_evaluations));
}

template <class Evaluate> void check_grid(Evaluate evaluate)
{
  using Value = std::invoke_result_t<Evaluate, std::size_t>;
  std::vector<Value> reference;
  for (std::size_t i = 0; i != 4; ++i)
    reference.push_back(evaluate(i));
  EXPECT_TRUE(reference[0].converged);
  EXPECT_TRUE(reference[1].converged);
  EXPECT_FALSE(reference[2].converged);
  EXPECT_FALSE(reference[3].converged);
  auto check = [&](IAsyncScheduler& scheduler) {
    std::vector<Value> actual;
    auto complete = [&] {
      auto owned = evaluate;
      return bethe::evaluate_observable_grid_async(4, std::move(owned),
                                                   [&](std::size_t i, Value value) {
                                                     EXPECT_EQ(i, actual.size());
                                                     actual.push_back(std::move(value));
                                                   },
                                                   {.scheduler = &scheduler, .max_in_flight = 2});
    }(); // the submitted callable owns its parameter/options copies
    EXPECT_EQ(complete.get_wait(scheduler), 4u);
    ASSERT_EQ(actual.size(), reference.size());
    for (std::size_t i = 0; i != actual.size(); ++i)
      same(actual[i], reference[i]);
  };
  DebugScheduler deferred;
  check(deferred);
  for (int threads : {1, 4})
  {
    TbbScheduler scheduler(threads);
    check(scheduler);
  }
}

template <typename R> class FiniteVolumeAsync : public ::testing::Test {};
TYPED_TEST_SUITE(FiniteVolumeAsync, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(FiniteVolumeAsync, LeeYangCompleteOwningStatesAndFailures)
{
  using R = TypeParam;
  ly::OneParticleOptions<R> options;
  options.tolerance = R{1e-9};
  check_grid([options](std::size_t i) {
    auto o = options;
    if (i == 2) o.max_iterations = 0;
    if (i == 3) o.max_cutoffs = 0;
    return ly::ground_state(R{1}, R(5 + i), static_cast<ly::Options<R> const&>(o));
  });
  check_grid([options](std::size_t i) {
    auto o = options;
    if (i == 2) o.max_root_iterations = 0;
    if (i == 3) o.max_kernel_products = 0;
    return ly::one_particle(R{1}, R(5 + i), o);
  });
}

TYPED_TEST(FiniteVolumeAsync, SineGordonCompleteOwningStatesAndFailures)
{
  using R = TypeParam;
  sg::TwoSolitonOptions<R> options;
  options.tolerance = R{1e-6};
  check_grid([options](std::size_t i) {
    auto o = options;
    if (i == 2) o.max_iterations = 0;
    if (i == 3) o.max_cutoffs = 0;
    return sg::vacuum_energy(R{1}, R(1 + i), R{2}, static_cast<sg::VacuumOptions<R> const&>(o));
  });
  check_grid([options](std::size_t i) {
    auto o = options;
    if (i == 2) o.max_root_iterations = 0;
    if (i == 3) o.max_cutoffs = 0;
    return sg::two_soliton_level(R{1}, R(1 + i), R{2}, uni20::from_twice(1 + 2 * int(i % 2)), o);
  });
}

TEST(IndependentSolves, SerialAndConcurrentMoveOnlyAndRepeatedTypes)
{
  for (int threads : {1, 2, 8})
  {
    auto const values =
        bethe::cli::independent_solves([x = std::make_unique<int>(7)] { return std::make_unique<int>(*x); },
                                       [] { return std::make_unique<int>(9); }, threads);
    EXPECT_EQ(*values.first, 7);
    EXPECT_EQ(*values.second, 9);
  }
  EXPECT_THROW(bethe::cli::independent_solves([] { return 1; }, [] { return 2; }, 0), std::invalid_argument);
}

TEST(IndependentSolves, FailedFirstWaitsForHeldSecond)
{
  using namespace std::chrono_literals;
  std::promise<void> entered, release, failed;
  auto entered_future = entered.get_future().share(), release_future = release.get_future().share();
  auto failed_future = failed.get_future();
  auto work = std::async(std::launch::async, [&] {
    return bethe::cli::independent_solves(
        [&]() -> int {
          if (entered_future.wait_for(5s) != std::future_status::ready) throw std::logic_error("second never entered");
          failed.set_value();
          throw std::runtime_error("first solve failed");
        },
        [&] {
          entered.set_value();
          if (release_future.wait_for(5s) != std::future_status::ready) throw std::logic_error("release timed out");
          return std::string("second finished");
        },
        2);
  });
  EXPECT_EQ(failed_future.wait_for(5s), std::future_status::ready);
  EXPECT_EQ(work.wait_for(20ms), std::future_status::timeout);
  release.set_value();
  try
  {
    (void)work.get();
    FAIL() << "missing first failure";
  }
  catch (std::runtime_error const& e)
  {
    EXPECT_STREQ(e.what(), "first solve failed");
  }
}
} // namespace
