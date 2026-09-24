// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
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
};

template <typename State> struct RealExcitation
{
    State state;
    /// E-E0 at the selected precision; absent if the ground reference failed.
    std::optional<decltype(State{}.energy)> gap;
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

template <typename Result, typename Solve, typename Ground, typename Energy = StateEnergy>
Result scan_real_combinations(std::size_t slots, std::size_t m, std::int64_t first, RealExcitationOptions const& scan,
                              Solve&& solve, Ground&& ground, EnergyOrder order = EnergyOrder::ascending,
                              Energy energy = {})
{
  if (scan.count == 0 || scan.max_candidates == 0)
    throw std::invalid_argument("excitation count and max_candidates must be positive");
  Result result;
  using State = decltype(result.ground_state);
  result.candidate_count = bounded_binomial(slots, m, scan.max_candidates);
  if (result.candidate_count == 0)
    throw std::invalid_argument("no supported real-root configurations at this precision");
  result.ground_state = ground();
  std::vector<uni20::half_int> numbers(m);
  for (std::size_t i = 0; i < m; ++i)
    numbers[i] = uni20::from_twice(first + 2 * static_cast<std::int64_t>(i));

  // A bounded max heap avoids retaining roots for every solved candidate.
  // An optional energy projection keeps small gaps and level ordering from
  // being lost to a shared extensive offset. Other models keep their existing
  // total-energy convention through the default projection.
  auto less = [order, &energy](auto const& left, auto const& right) {
    auto const a = energy(left.state), b = energy(right.state);
    if (a != b) return order == EnergyOrder::ascending ? a < b : a > b;
    return left.state.quantum_numbers < right.state.quantum_numbers;
  };
  auto const keep = std::min(scan.count, result.candidate_count);
  result.levels.reserve(keep);
  for (;;)
  {
    auto state = numbers == result.ground_state.quantum_numbers ? result.ground_state : solve(numbers);
    if (state.converged)
    {
      ++result.converged_count;
      RealExcitation<State> level{.state = std::move(state), .gap = std::nullopt};
      if (result.ground_state.converged) level.gap = energy(level.state) - energy(result.ground_state);
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
    else if (!result.first_unconverged)
      result.first_unconverged = std::move(state);

    // Lexicographic combinations of M slots, including the unique M=0 set.
    std::size_t i = m;
    while (i > 0 && numbers[i - 1].twice() == first + 2 * static_cast<std::int64_t>(slots - m + i - 1))
      --i;
    if (i == 0) break;
    numbers[i - 1] = uni20::from_twice(numbers[i - 1].twice() + 2);
    for (std::size_t j = i; j < m; ++j)
      numbers[j] = uni20::from_twice(numbers[j - 1].twice() + 2);
  }
  std::sort_heap(result.levels.begin(), result.levels.end(), less);
  return result;
}

} // namespace detail
} // namespace bethe
