// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include <bethe/heisenberg_open.hpp>
#include <bethe/real_excitations.hpp>

#include <numeric>
#include <optional>
#include <utility>

namespace bethe::heisenberg
{
using RealExcitationOptions = bethe::RealExcitationOptions;
template <typename State> using RealExcitation = bethe::RealExcitation<State>;
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
  return bethe::detail::bounded_binomial(slots, m, limit);
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
  (void)real_excitation_count(sites, spin, scan.max_candidates);
  auto const m = static_cast<std::size_t>((static_cast<std::int64_t>(sites) - spin.twice()) / 2);
  auto const slots = sites - m;
  auto const first = periodic ? -static_cast<std::int64_t>(slots - 1) : std::int64_t{2};
  auto result = bethe::detail::scan_real_combinations<RealExcitationScan<State>>(
      slots, m, first, scan, std::forward<Solve>(solve), std::forward<Ground>(ground));
  result.spin = spin;
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
