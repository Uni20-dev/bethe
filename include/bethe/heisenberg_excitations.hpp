// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include <bethe/heisenberg_open.hpp>

#include <numeric>
#include <optional>
#include <utility>

namespace bethe::heisenberg
{
/// Scan the entire supported real-root family at a fixed total spin, retaining
/// at most count lowest converged multiplets. The sector minimum is included.
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
    uni20::half_int spin;
    State ground_state;
    /// Energy ordered; exact energy ties use lexicographic quantum numbers.
    /// Degenerate, distinct multiplets are not merged. Failed solves are excluded.
    std::vector<RealExcitation<State>> levels;
    std::size_t candidate_count = 0;
    std::size_t converged_count = 0;
    /// A diagnostic example; converged_count accounts for ALL candidates.
    std::optional<State> first_unconverged;

    /// Completeness only within the explicitly supported quantum-number window,
    /// never within the full spin sector or the full physical spectrum.
    bool family_converged() const { return candidate_count == converged_count; }
    bool converged() const { return family_converged() && ground_state.converged; }
};

/// Both boundaries have C(N-M,M) supported configurations, M=N/2-S.
/// Count using integer arithmetic, rejecting limit/size_t overflow without
/// allocating or iterating through the configurations. Negative S is invalid.
[[nodiscard]] inline std::size_t real_excitation_count(std::size_t sites, uni20::half_int spin,
                                                       std::size_t limit = std::numeric_limits<std::size_t>::max())
{
  auto const n = detail::checked_sites(sites);
  if (spin.twice() < 0 || spin.twice() > n || (n - spin.twice()) % 2 != 0)
    throw std::invalid_argument("S must lie in [0,N/2] with the same half-integer parity as N/2");
  auto const m = static_cast<std::size_t>((n - spin.twice()) / 2);
  auto const slots = sites - m;
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

namespace detail
{
template <typename State, typename Solve, typename Ground>
RealExcitationScan<State> scan_real_excitations(std::size_t sites, uni20::half_int spin,
                                                RealExcitationOptions const& scan, bool periodic, Solve&& solve,
                                                Ground&& ground)
{
  if (scan.count == 0 || scan.max_candidates == 0)
    throw std::invalid_argument("excitation count and max_candidates must be positive");
  RealExcitationScan<State> result;
  result.spin = spin;
  result.candidate_count = real_excitation_count(sites, spin, scan.max_candidates);
  result.ground_state = ground();
  auto const m = static_cast<std::size_t>((static_cast<std::int64_t>(sites) - spin.twice()) / 2);
  auto const slots = sites - m;
  auto const first = periodic ? -static_cast<std::int64_t>(slots - 1) : std::int64_t{2};
  QuantumNumbers numbers(m);
  for (std::size_t i = 0; i < m; ++i)
    numbers[i] = uni20::from_twice(first + 2 * static_cast<std::int64_t>(i));

  // A bounded max heap avoids retaining roots for every solved candidate.
  auto less = [](auto const& left, auto const& right) {
    if (left.state.energy != right.state.energy) return left.state.energy < right.state.energy;
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
      if (result.ground_state.converged) level.gap = level.state.energy - result.ground_state.energy;
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

/// Periodic real-root highest-weight multiplets at total spin S (Sz=S).
/// Exhaustive within the conventional all-1-string window, not a full spectrum.
/// O(count*M+N) retained storage; every candidate still requires a solve.
template <uni20::Real Real = double>
[[nodiscard]] RealExcitationScan<RealState<Real>> real_excitations(std::size_t sites, uni20::half_int spin,
                                                                   RealExcitationOptions const& scan = {},
                                                                   SolverOptions<Real> const& solver = {})
{
  return detail::scan_real_excitations<RealState<Real>>(
      sites, spin, scan, true, [&](QuantumNumbers const& numbers) { return solve_real<Real>(sites, numbers, solver); },
      [&] { return ground_state<Real>(sites, solver); });
}

namespace open
{
/// Free-end analogue: positive real roots with integer labels in [1,N-M].
/// No lattice momentum. String states and SU(2) descendants are not enumerated.
template <uni20::Real Real = double>
[[nodiscard]] RealExcitationScan<RealState<Real>> real_excitations(std::size_t sites, uni20::half_int spin,
                                                                   RealExcitationOptions const& scan = {},
                                                                   SolverOptions<Real> const& solver = {})
{
  return heisenberg::detail::scan_real_excitations<RealState<Real>>(
      sites, spin, scan, false,
      [&](QuantumNumbers const& numbers) { return open::solve_real<Real>(sites, numbers, solver); },
      [&] { return open::ground_state<Real>(sites, solver); });
}
} // namespace open
} // namespace bethe::heisenberg
