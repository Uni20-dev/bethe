// SPDX-License-Identifier: GPL-3.0-or-later
// Frozen batched scanner from public 5909c21, for regression tests only.
#pragma once
#include <bethe/real_excitations.hpp>
namespace scan_reference
{
using bethe::RealExcitation;
using bethe::RealExcitationOptions;
using bethe::detail::advance_combination;
using bethe::detail::bounded_binomial;
using bethe::detail::EnergyOrder;
using bethe::detail::StateEnergy;
template <typename Result, typename Solve, typename Ground, typename Energy = StateEnergy>
Result scan_real_combinations(std::size_t slots, std::size_t m, std::int64_t first, RealExcitationOptions const& scan,
                              Solve&& solve, Ground&& ground, EnergyOrder order = EnergyOrder::ascending,
                              Energy energy = {})
{
  if (scan.count == 0 || scan.max_candidates == 0)
    throw std::invalid_argument("excitation count and max_candidates must be positive");
  auto& scheduler = scan.execution.resolved_scheduler();
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
  for (std::size_t offset = 0; offset < result.candidate_count;)
  {
    auto const size = std::min(scan.execution.batch_size, result.candidate_count - offset);
    std::vector<std::vector<uni20::half_int>> labels;
    labels.reserve(size);
    for (std::size_t i = 0; i < size; ++i)
    {
      labels.push_back(numbers);
      // Empty sets have one combination; no slots-1 arithmetic is needed.
      if (offset + i + 1 < result.candidate_count)
        advance_combination<uni20::half_int>(numbers,
                                             uni20::from_twice(first + 2 * static_cast<std::int64_t>(slots - 1)));
    }
    std::vector<std::optional<State>> states(size);
    scheduler.execute_batch(size, [&](std::size_t i) {
      states[i] = labels[i] == result.ground_state.quantum_numbers ? result.ground_state : solve(labels[i]);
    });
    // Merge in enumeration order, independent of completion order. In particular
    // the first failure and equal-energy tie breaking retain their serial meaning.
    for (auto& slot : states)
    {
      auto state = std::move(*slot);
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
    }
    offset += size;
  }
  std::sort_heap(result.levels.begin(), result.levels.end(), less);
  return result;
}
} // namespace scan_reference
