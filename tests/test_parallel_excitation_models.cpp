// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/biquadratic_ferromagnetic.hpp>
#include <bethe/lieb_liniger_open.hpp>
#include <bethe/q_boson.hpp>
#include <bethe/xxz_quantum_group_critical.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
class RecordingScheduler : public uni20::async::DebugScheduler {
  public:
    std::size_t tasks = 0, largest_batch = 0;

  private:
    void execute_batch_impl(uni20::async::LightweightTaskBatch const& batch) override
    {
      tasks += batch.size();
      largest_batch = std::max(largest_batch, batch.size());
      for (std::size_t i = batch.size(); i > 0; --i)
        batch(i - 1);
    }
};

template <typename State> void equal_state(State const& a, State const& b)
{
  EXPECT_EQ(a.energy, b.energy);
  if constexpr (requires { a.reference; })
  {
    EXPECT_EQ(a.tl_energy, b.tl_energy);
    EXPECT_EQ(a.through_lines, b.through_lines);
    EXPECT_EQ(a.multiplicity, b.multiplicity);
    equal_state(a.reference, b.reference);
  }
  else
  {
    EXPECT_EQ(a.quantum_numbers, b.quantum_numbers);
    EXPECT_EQ(a.converged, b.converged);
    EXPECT_EQ(a.status, b.status);
    EXPECT_EQ(a.residual_norm, b.residual_norm);
    EXPECT_EQ(a.iterations, b.iterations);
    if constexpr (requires { a.momenta; }) EXPECT_EQ(a.momenta, b.momenta);
    if constexpr (requires { a.rapidities; }) EXPECT_EQ(a.rapidities, b.rapidities);
    if constexpr (requires { a.angles; }) EXPECT_EQ(a.angles, b.angles);
    if constexpr (requires { a.scaled_roots; }) EXPECT_EQ(a.scaled_roots, b.scaled_roots);
    if constexpr (requires { a.energy_shift; }) EXPECT_EQ(a.energy_shift, b.energy_shift);
    if constexpr (requires { a.modes; }) EXPECT_EQ(a.modes, b.modes);
    if constexpr (requires { a.deviations; }) EXPECT_EQ(a.deviations, b.deviations);
    if constexpr (requires { a.momentum; }) EXPECT_EQ(a.momentum, b.momentum);
    if constexpr (requires { a.momentum_index; }) EXPECT_EQ(a.momentum_index, b.momentum_index);
  }
}

template <typename Scan> void equal_scan(Scan const& a, Scan const& b)
{
  EXPECT_EQ(a.candidate_count, b.candidate_count);
  EXPECT_EQ(a.converged_count, b.converged_count);
  EXPECT_EQ(a.converged(), b.converged());
  equal_state(a.ground_state, b.ground_state);
  ASSERT_EQ(a.first_unconverged.has_value(), b.first_unconverged.has_value());
  if (a.first_unconverged) equal_state(*a.first_unconverged, *b.first_unconverged);
  ASSERT_EQ(a.levels.size(), b.levels.size());
  for (std::size_t i = 0; i < a.levels.size(); ++i)
  {
    SCOPED_TRACE(i);
    equal_state(a.levels[i].state, b.levels[i].state);
    EXPECT_EQ(a.levels[i].gap, b.levels[i].gap);
  }
}

template <typename Solve> void check_schedulers(Solve solve)
{
  using namespace uni20::async;
  DebugScheduler serial({.order = DebugSchedulerOrder::fifo});
  TbbScheduler parallel(4);
  auto* original = get_global_scheduler();
  for (std::size_t iterations : {10000, 0})
  {
    SCOPED_TRACE(iterations);
    auto const reference = solve({.count = 3, .max_candidates = 100, .execution = {&serial, 1}}, iterations);
    EXPECT_EQ(reference.converged(), iterations != 0);
    if (iterations == 0) EXPECT_TRUE(reference.first_unconverged);
    RecordingScheduler recording;
    auto const reversed = solve({.count = 3, .max_candidates = 100, .execution = {&recording, 3}}, iterations);
    EXPECT_EQ(recording.tasks, reference.candidate_count); // scheduler selection was not ignored
    EXPECT_LE(recording.largest_batch, 3U);
    equal_scan(reference, reversed);
    auto const concurrent = solve({.count = 3, .max_candidates = 100, .execution = {&parallel, 7}}, iterations);
    equal_scan(reference, concurrent);
  }
  EXPECT_EQ(get_global_scheduler(), original);
}
} // namespace

template <typename T> class ParallelExcitationModels : public ::testing::Test {};
TYPED_TEST_SUITE(ParallelExcitationModels, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(ParallelExcitationModels, LiebLinigerPeriodicAndOpen)
{
  using Real = TypeParam;
  check_schedulers([](bethe::RealExcitationOptions scan, std::size_t iterations) {
    return bethe::lieb_liniger::real_excitations<Real>(3, Real{4}, Real{2}, 1, scan, {.max_iterations = iterations});
  });
  check_schedulers([](bethe::RealExcitationOptions scan, std::size_t iterations) {
    return bethe::lieb_liniger::open::real_excitations<Real>(3, Real{4}, Real{2}, 3, scan,
                                                             {.max_iterations = iterations});
  });
}

TYPED_TEST(ParallelExcitationModels, QBoson)
{
  using Real = TypeParam;
  check_schedulers([](bethe::RealExcitationOptions scan, std::size_t iterations) {
    return bethe::q_boson::real_excitations<Real>(4, 3, Real{1}, scan, {.max_iterations = iterations});
  });
}

TYPED_TEST(ParallelExcitationModels, QuantumGroupCritical)
{
  using Real = TypeParam;
  check_schedulers([](bethe::RealExcitationOptions scan, std::size_t iterations) {
    return bethe::xxz::quantum_group::critical::real_excitations<Real>(8, Real{0.6}, 4, scan,
                                                                       {.max_iterations = iterations});
  });
}

TYPED_TEST(ParallelExcitationModels, BiquadraticSignsAndWindow)
{
  using Real = TypeParam;
  check_schedulers([](bethe::RealExcitationOptions scan, std::size_t iterations) {
    return bethe::biquadratic::real_excitations<Real>(8, 2, scan, {.max_iterations = iterations});
  });
  check_schedulers([](bethe::RealExcitationOptions scan, std::size_t iterations) {
    return bethe::biquadratic::ferromagnetic::real_excitations<Real>(9, 5, scan, {.max_iterations = iterations});
  });
  check_schedulers([](bethe::RealExcitationOptions scan, std::size_t iterations) {
    return bethe::biquadratic::ferromagnetic::real_excitations_window<Real>(9, 5, 3, scan,
                                                                            {.max_iterations = iterations});
  });
}
