// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <algorithm>
#include <bethe/detail/continuous_work.hpp>
#include <bethe/execution.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mutex>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <uni20/async/async.hpp>
#include <uni20/common/half_int.hpp>
#include <utility>
#include <vector>

namespace bethe
{
/// Scan the entire supported real-root family in a model-defined window, retaining
/// at most count lowest converged states. The supported family's minimum is
/// included; this need not be the true sector minimum if complex roots are missing.
struct RealExcitationOptions
{
    std::size_t count = 10;
    /// Reject larger families before allocating roots or solving any state.
    std::size_t max_candidates = 10000;
    ExecutionOptions execution = {};
};

namespace detail
{
template <typename T> struct EnergyValue
{
    using type = T;
};
template <typename T> struct EnergyValue<std::optional<T>>
{
    using type = T;
};
} // namespace detail

template <typename State> struct RealExcitation
{
    State state;
    /// E-E0 at the selected precision; absent if the ground reference failed.
    std::optional<typename detail::EnergyValue<decltype(State{}.energy)>::type> gap;
};

template <typename State> struct RealExcitationScan
{
    State ground_state;
    /// Energy ordered; exact energy ties use lexicographic quantum numbers.
    /// Degenerate, distinct states are not merged. Failed solves are excluded.
    std::vector<RealExcitation<State>> levels;
    std::size_t candidate_count = 0;
    std::size_t converged_count = 0;
    /// A diagnostic example; converged_count accounts for ALL candidates.
    std::optional<State> first_unconverged;

    /// Completeness only within the explicitly supported quantum-number window,
    /// never within the full symmetry sector or the full physical spectrum.
    bool family_converged() const { return candidate_count == converged_count; }
    bool converged() const { return family_converged() && ground_state.converged; }
};

namespace detail
{
enum class EnergyOrder
{
  ascending,
  descending
};

struct StateEnergy
{
    template <typename State> auto operator()(State const& state) const { return state.energy; }
};

inline std::size_t bounded_binomial(std::size_t slots, std::size_t m, std::size_t limit)
{
  if (m > slots) return 0;
  auto const k = std::min(m, slots - m);
  std::size_t count = 1;
  if (limit == 0) throw std::length_error("real-root family exceeds max_candidates");
  for (std::size_t i = 1; i <= k; ++i)
  {
    auto numerator = slots - k + i;
    auto const divisor = std::gcd(numerator, i);
    numerator /= divisor;
    // The remaining denominator divides the previous binomial coefficient.
    count /= i / divisor;
    if (count > limit / numerator) throw std::length_error("real-root family exceeds max_candidates");
    count *= numerator;
  }
  return count;
}

/// Advance ordered unit-spaced labels lexicographically within an inclusive
/// upper bound. The caller supplies a valid combination. Empty sets have one
/// combination, so advancing an empty span returns false.
template <typename Integer> bool advance_combination(std::span<Integer> labels, Integer last)
{
  std::size_t i = labels.size();
  while (i > 0 && labels[i - 1] == last)
  {
    --i;
    if (i > 0) --last;
  }
  if (i == 0) return false;
  ++labels[i - 1];
  for (std::size_t j = i; j < labels.size(); ++j)
  {
    labels[j] = labels[j - 1];
    ++labels[j];
  }
  return true;
}

/// Synchronous core, called within one owning Async operation. The caller has
/// checked the budget and initialized result.candidate_count. transform enriches
/// a converged level outside the collection lock, before releasing its runner;
/// the returned record must retain .state and .gap for ordering and diagnostics.
/// This fuses per-state follow-up work without a second queue or per-state task.
template <typename Result, typename Solve, typename Ground, typename Energy, typename Transform>
void scan_real_combinations_into(uni20::async::IAsyncScheduler& scheduler, std::size_t slots, std::size_t m,
                                 std::int64_t first, RealExcitationOptions const& scan, Solve& solve, Ground& ground,
                                 EnergyOrder order, Energy& energy, Result& result, Transform transform)
{
  if (scan.execution.stop.stop_requested()) throw ScanCancelled{};
  result.ground_state = ground();
  if (scan.execution.stop.stop_requested()) throw ScanCancelled{};
  using State = decltype(result.ground_state);
  std::vector<uni20::half_int> numbers(m);
  for (std::size_t i = 0; i < m; ++i)
    numbers[i] = uni20::from_twice(first + 2 * static_cast<std::int64_t>(i));
  auto less = [order, &energy](auto const& left, auto const& right) {
    auto const a = energy(left.state), b = energy(right.state);
    if (a != b) return order == EnergyOrder::ascending ? a < b : a > b;
    return left.state.quantum_numbers < right.state.quantum_numbers;
  };
  auto const keep = std::min(scan.count, result.candidate_count);
  result.levels.reserve(keep);
  std::mutex mutex;
  std::size_t next = 0, first_failure = result.candidate_count;
  auto completed = continuous_work(
      scheduler, scan.execution.max_in_flight, result.candidate_count,
      [&](std::size_t, std::stop_token stop) {
        std::size_t id;
        std::vector<uni20::half_int> labels;
        {
          // Label generation and ID assignment are atomic together. Arrival at
          // the collector never defines a candidate's enumeration identity.
          std::lock_guard lock(mutex);
          id = next++;
          labels = numbers;
          if (next < result.candidate_count)
            advance_combination<uni20::half_int>(numbers,
                                                 uni20::from_twice(first + 2 * static_cast<std::int64_t>(slots - 1)));
        }
        auto state = labels == result.ground_state.quantum_numbers ? result.ground_state : solve(labels);
        // A cancelled scan never returns partial data. Do not start a follow-up
        // kernel after cancellation (or another runner's execution failure).
        if (stop.stop_requested()) return;
        if (state.converged)
        {
          RealExcitation<State> solved{.state = std::move(state), .gap = std::nullopt};
          if (result.ground_state.converged) solved.gap = energy(solved.state) - energy(result.ground_state);
          auto level = transform(std::move(solved), std::as_const(result.ground_state));
          std::lock_guard lock(mutex);
          ++result.converged_count;
          if (result.levels.size() < keep)
          {
            result.levels.push_back(std::move(level));
            std::push_heap(result.levels.begin(), result.levels.end(), less);
          }
          else if (less(level, result.levels.front()))
          {
            std::pop_heap(result.levels.begin(), result.levels.end(), less);
            result.levels.back() = std::move(level);
            std::push_heap(result.levels.begin(), result.levels.end(), less);
          }
        }
        else
        {
          std::lock_guard lock(mutex);
          if (id < first_failure)
          {
            first_failure = id;
            result.first_unconverged = std::move(state);
          }
        }
      },
      scan.execution.stop);
  // Cancellation is an execution outcome, not numerical nonconvergence. Active
  // kernels finish before it is observed; no partly scanned family is returned.
  if (completed != result.candidate_count || scan.execution.stop.stop_requested()) throw ScanCancelled{};
  std::sort_heap(result.levels.begin(), result.levels.end(), less);
}

/// One owning coroutine, one joined batch of continuous runners. Neither the
/// coroutine inventory nor completed intermediate queues grow with family size.
template <typename Result, typename Solve, typename Ground, typename Energy>
uni20::async::AsyncTask scan_real_combinations_task(uni20::async::IAsyncScheduler& scheduler, std::size_t slots,
                                                    std::size_t m, std::int64_t first, RealExcitationOptions scan,
                                                    Solve solve, Ground ground, EnergyOrder order, Energy energy,
                                                    Result result, uni20::async::WriteBuffer<Result> output)
{
  // Enroll the exception sink before any preparation or allocation can fail.
  auto destination = co_await output;
  scan_real_combinations_into(scheduler, slots, m, first, scan, solve, ground, order, energy, result,
                              [](auto level, auto const&) { return level; });
  destination = std::move(result);
}

/// Submit a complete finite-family scan. Preflight errors throw before any
/// scientific work. Callbacks/configuration are decay-copied into the task:
/// references inside a callback remain the caller's lifetime responsibility.
/// The scheduler is borrowed until completion, even if the result is discarded.
template <typename Result, typename Solve, typename Ground, typename Energy = StateEnergy>
uni20::async::Async<Result> scan_real_combinations_async(std::size_t slots, std::size_t m, std::int64_t first,
                                                         RealExcitationOptions scan, Solve solve, Ground ground,
                                                         EnergyOrder order = EnergyOrder::ascending, Energy energy = {},
                                                         Result seed = {})
{
  if (scan.count == 0 || scan.max_candidates == 0)
    throw std::invalid_argument("excitation count and max_candidates must be positive");
  auto& scheduler = scan.execution.continuous_scheduler();
  seed.candidate_count = bounded_binomial(slots, m, scan.max_candidates);
  if (seed.candidate_count == 0) throw std::invalid_argument("no supported real-root configurations at this precision");
  uni20::async::Async<Result> output;
  // Only ONE task is admitted. All later work is within execute_batch's joining
  // scope, including allocation/worker failure; no partially submitted graph.
  auto task = scan_real_combinations_task(scheduler, slots, m, first, std::move(scan), std::move(solve),
                                          std::move(ground), order, std::move(energy), std::move(seed), output.write());
  schedule_owned(scheduler, std::move(task));
  return output;
}

/// Move a privately owned result only after waiting on the explicitly resolved
/// scheduler. Uni20's move_from_wait uses the active scheduler if not yet ready.
template <typename T> T take_scan_result(uni20::async::Async<T>& result, uni20::async::IAsyncScheduler& scheduler)
{
  return take_result(result, scheduler);
}

template <typename Result, typename Solve, typename Ground, typename Energy = StateEnergy>
Result scan_real_combinations(std::size_t slots, std::size_t m, std::int64_t first, RealExcitationOptions const& scan,
                              Solve&& solve, Ground&& ground, EnergyOrder order = EnergyOrder::ascending,
                              Energy energy = {})
{
  auto& scheduler = scan.execution.resolved_scheduler();
  auto result = scan_real_combinations_async<Result>(slots, m, first, scan, std::forward<Solve>(solve),
                                                     std::forward<Ground>(ground), order, std::move(energy));
  return take_scan_result(result, scheduler);
}

} // namespace detail
} // namespace bethe
