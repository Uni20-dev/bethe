// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <atomic>
#include <bethe/xxx_structure_factor.hpp>
#include <bethe/xxz_excitations.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
struct State
{
    std::vector<uni20::half_int> quantum_numbers;
    double energy = 0;
    bool converged = true;
};
using Scan = bethe::RealExcitationScan<State>;
class RecordingScheduler : public uni20::async::DebugScheduler {
  public:
    std::vector<std::size_t> sizes;

  private:
    void execute_batch_impl(uni20::async::LightweightTaskBatch const& batch) override
    {
      sizes.push_back(batch.size());
      for (std::size_t i = batch.size(); i > 0; --i)
        batch(i - 1);
    }
};
auto scan(bethe::ExecutionOptions execution, bethe::detail::EnergyOrder order)
{
  return bethe::detail::scan_real_combinations<Scan>(
      8, 2, 0, {5, 100, execution},
      [](auto const& labels) {
        // Tied energies and multiple failed candidates test deterministic merges.
        return State{labels, double(labels[0].twice() % 3), labels[1].twice() % 4 != 0};
      },
      [] { return State{{}, -1}; }, order);
}
void equal(Scan const& a, Scan const& b)
{
  EXPECT_EQ(a.candidate_count, b.candidate_count);
  EXPECT_EQ(a.converged_count, b.converged_count);
  ASSERT_TRUE(a.first_unconverged && b.first_unconverged);
  EXPECT_EQ(a.first_unconverged->quantum_numbers, b.first_unconverged->quantum_numbers);
  ASSERT_EQ(a.levels.size(), b.levels.size());
  for (std::size_t i = 0; i < a.levels.size(); ++i)
  {
    EXPECT_EQ(a.levels[i].state.quantum_numbers, b.levels[i].state.quantum_numbers);
    EXPECT_EQ(a.levels[i].state.energy, b.levels[i].state.energy);
    EXPECT_EQ(a.levels[i].gap, b.levels[i].gap);
  }
}
} // namespace

TEST(ParallelScan, BoundedHeapFailureAndTieOrdering)
{
  using namespace uni20::async;
  DebugScheduler fifo({.order = DebugSchedulerOrder::fifo});
  DebugScheduler reverse;
  DebugScheduler random({.order = DebugSchedulerOrder::random, .random_seed = 42});
  TbbScheduler parallel(4);
  for (auto order : {bethe::detail::EnergyOrder::ascending, bethe::detail::EnergyOrder::descending})
  {
    auto const reference = scan({&fifo, 1}, order);
    for (auto* scheduler : std::vector<IAsyncScheduler*>{&reverse, &random, &parallel})
      for (std::size_t batch : {1, 3, 128})
        equal(reference, scan({scheduler, batch}, order));
  }
}

TEST(ParallelScan, PreflightVacuumAndException)
{
  uni20::async::TbbScheduler scheduler(4);
  std::atomic<int> calls = 0;
  auto solve = [&](auto const& labels) {
    ++calls;
    return State{labels, 0};
  };
  auto ground = [] { return State{{}, 0}; };
  EXPECT_THROW((bethe::detail::scan_real_combinations<Scan>(8, 2, 0, {5, 1, {&scheduler}}, solve, ground)),
               std::length_error);
  EXPECT_THROW((bethe::detail::scan_real_combinations<Scan>(8, 2, 0, {5, 100, {&scheduler, 0}}, solve, ground)),
               std::invalid_argument);
  EXPECT_EQ(calls, 0);
  auto const vacuum = bethe::detail::scan_real_combinations<Scan>(0, 0, 0, {1, 1, {&scheduler}}, solve, ground);
  EXPECT_EQ(vacuum.levels.size(), 1);
  EXPECT_EQ(calls, 0); // reference reused, including the empty quantum-number set
  std::atomic<int> active = 0;
  auto fail = [&](auto const&) -> State {
    struct Guard
    {
        std::atomic<int>& active;
        Guard(std::atomic<int>& a) : active(a) { ++active; }
        ~Guard() { --active; }
    } guard(active);
    throw std::runtime_error("injected solve failure");
  };
  EXPECT_THROW((bethe::detail::scan_real_combinations<Scan>(8, 2, 0, {5, 100, {&scheduler}}, fail, ground)),
               std::runtime_error);
  EXPECT_EQ(active, 0); // no work remains after synchronous exception propagation
}

TEST(ParallelScan, BoundedBatchesAndActiveSchedulerFallback)
{
  RecordingScheduler active, explicit_scheduler;
  uni20::async::ScopedScheduler scope(&active);
  auto const reference = scan({nullptr, 6}, bethe::detail::EnergyOrder::ascending);
  EXPECT_EQ(active.sizes, (std::vector<std::size_t>{6, 6, 6, 6, 4}));
  auto const result = scan({&explicit_scheduler, 8}, bethe::detail::EnergyOrder::ascending);
  EXPECT_EQ(explicit_scheduler.sizes, (std::vector<std::size_t>{8, 8, 8, 4}));
  EXPECT_EQ(active.sizes.size(), 5); // explicit selection did not touch the active scheduler
  EXPECT_EQ(uni20::async::get_global_scheduler(), &active);
  equal(reference, result);
}

template <typename T> class ParallelStructureFactor : public ::testing::Test {};
TYPED_TEST_SUITE(ParallelStructureFactor, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(ParallelStructureFactor, ExactAgreementAcrossSchedulers)
{
  using namespace uni20::async;
  DebugScheduler fifo({.order = DebugSchedulerOrder::fifo});
  DebugScheduler random({.order = DebugSchedulerOrder::random, .random_seed = 123});
  TbbScheduler parallel(4);
  auto* original = get_global_scheduler();
  auto const reference = bethe::heisenberg::two_spinon_structure_factor<TypeParam>(12, {}, 100, {&fifo, 1});
  ASSERT_TRUE(reference.converged());
  for (auto* scheduler : std::vector<IAsyncScheduler*>{&random, &parallel})
  {
    auto const result = bethe::heisenberg::two_spinon_structure_factor<TypeParam>(12, {}, 100, {scheduler, 7});
    ASSERT_TRUE(result.converged());
    ASSERT_EQ(result.lines.size(), reference.lines.size());
    for (std::size_t i = 0; i < result.lines.size(); ++i)
    {
      EXPECT_EQ(result.lines[i].state_id, reference.lines[i].state_id);
      EXPECT_EQ(result.lines[i].momentum_index, reference.lines[i].momentum_index);
      EXPECT_EQ(result.lines[i].gap, reference.lines[i].gap);
      EXPECT_EQ(result.lines[i].weight, reference.lines[i].weight);
      EXPECT_EQ(result.scan.levels[i].state.rapidities, reference.scan.levels[i].state.rapidities);
      EXPECT_EQ(result.form_factors[i].pivot_margin, reference.form_factors[i].pivot_margin);
    }
    EXPECT_EQ(result.weight_fraction, reference.weight_fraction);
    EXPECT_EQ(result.first_moment_fraction, reference.first_moment_fraction);
    for (std::size_t q = 0; q < 12; ++q)
    {
      EXPECT_EQ(result.moments[q].weight, reference.moments[q].weight);
      EXPECT_EQ(result.moments[q].first_moment, reference.moments[q].first_moment);
    }
  }
  EXPECT_EQ(get_global_scheduler(), original);
}

TYPED_TEST(ParallelStructureFactor, XXZOpenAndPeriodicSolvers)
{
  using Real = TypeParam;
  uni20::async::DebugScheduler serial;
  uni20::async::TbbScheduler parallel(4);
  auto check = [&](auto solve) {
    auto const a = solve(bethe::ExecutionOptions{&serial, 1});
    auto const b = solve(bethe::ExecutionOptions{&parallel, 7});
    ASSERT_TRUE(a.converged());
    ASSERT_TRUE(b.converged());
    ASSERT_EQ(a.levels.size(), b.levels.size());
    EXPECT_EQ(a.ground_state.rapidities, b.ground_state.rapidities);
    for (std::size_t i = 0; i < a.levels.size(); ++i)
    {
      EXPECT_EQ(a.levels[i].state.quantum_numbers, b.levels[i].state.quantum_numbers);
      EXPECT_EQ(a.levels[i].state.rapidities, b.levels[i].state.rapidities);
      EXPECT_EQ(a.levels[i].state.energy, b.levels[i].state.energy);
      EXPECT_EQ(a.levels[i].gap, b.levels[i].gap);
    }
  };
  // The existing finite-real excitation scanner supports 0 <= Delta <= 1;
  // negative/gapped ground-state and thermodynamic APIs are separate paths.
  for (Real delta : {Real{0}, Real{0.5}, Real{1}})
  {
    check([&](auto execution) {
      return bethe::xxz::real_excitations<Real>(8, delta, uni20::half_int{1}, {100, 100, execution});
    });
    check([&](auto execution) {
      return bethe::xxz::open::real_excitations<Real>(8, delta, uni20::half_int{1}, {100, 100, execution});
    });
  }
}
