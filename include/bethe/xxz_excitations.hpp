// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/real_excitations.hpp>
#include <bethe/xxz.hpp>
#include <bethe/xxz_open.hpp>

namespace bethe::xxz
{
using RealExcitationOptions = bethe::RealExcitationOptions;
template <typename State> using RealExcitation = bethe::RealExcitation<State>;
template <typename State> struct RealExcitationScan : bethe::RealExcitationScan<State>
{
    uni20::half_int sz;
    decltype(State{}.delta) delta = 0;
    RealQuantumNumberWindow window;
};

/// Count combinations in the supported anisotropy-dependent window, not all
/// states of the Sz sector. Preflight is allocation-free and overflow-checked.
template <uni20::Real Real = double>
[[nodiscard]] std::size_t real_excitation_count(std::size_t sites, Real delta, uni20::half_int sz,
                                                std::size_t limit = std::numeric_limits<std::size_t>::max())
{
  auto const window = xxz::real_quantum_number_window(sites, delta, sz);
  return bethe::detail::bounded_binomial(window.slots, detail::sector_roots(sites, sz), limit);
}

/// Scan every configuration in the supported finite-real window, returning
/// the lowest count converged states. Includes the sector minimum; NOT the
/// complete Sz spectrum. At Delta=1 these are the XXX highest-weight states
/// (or spin-reversed partners), not their infinite-root descendants.
/// Gaps use the GLOBAL ground reference, not the selected sector minimum.
template <uni20::Real Real = double>
[[nodiscard]] RealExcitationScan<RealState<Real>> real_excitations(std::size_t sites, Real delta, uni20::half_int sz,
                                                                   RealExcitationOptions const& scan = {},
                                                                   SolverOptions<Real> const& solver = {})
{
  auto const window = xxz::real_quantum_number_window(sites, delta, sz);
  auto const m = detail::sector_roots(sites, sz);
  auto result = bethe::detail::scan_real_combinations<RealExcitationScan<RealState<Real>>>(
      window.slots, m, window.first.twice(), scan,
      [&](QuantumNumbers const& numbers) { return xxz::solve_real<Real>(sites, delta, numbers, solver); },
      [&] { return xxz::ground_state<Real>(sites, delta, solver); });
  result.sz = sz;
  result.delta = delta;
  result.window = window;
  for (auto& level : result.levels)
  {
    level.state.sz = sz;
    level.state.spin_reversed = sz.twice() < 0;
  }
  if (result.first_unconverged)
  {
    result.first_unconverged->sz = sz;
    result.first_unconverged->spin_reversed = sz.twice() < 0;
  }
  return result;
}
namespace open
{
using xxz::RealExcitation;
using xxz::RealExcitationOptions;
using xxz::RealExcitationScan;

template <uni20::Real Real = double>
[[nodiscard]] std::size_t real_excitation_count(std::size_t sites, Real delta, uni20::half_int sz,
                                                std::size_t limit = std::numeric_limits<std::size_t>::max())
{
  auto const window = open::real_quantum_number_window(sites, delta, sz);
  return bethe::detail::bounded_binomial(window.slots, xxz::detail::sector_roots(sites, sz), limit);
}

/// Positive finite-root, free-end XXZ states at fixed Sz; NOT the complete
/// sector spectrum. Includes the sector minimum. Gaps use the global ground.
template <uni20::Real Real = double>
[[nodiscard]] RealExcitationScan<RealState<Real>> real_excitations(std::size_t sites, Real delta, uni20::half_int sz,
                                                                   RealExcitationOptions const& scan = {},
                                                                   SolverOptions<Real> const& solver = {})
{
  auto const window = open::real_quantum_number_window(sites, delta, sz);
  auto const m = xxz::detail::sector_roots(sites, sz);
  auto result = bethe::detail::scan_real_combinations<RealExcitationScan<RealState<Real>>>(
      window.slots, m, window.first.twice(), scan,
      [&](QuantumNumbers const& numbers) { return open::solve_real<Real>(sites, delta, numbers, solver); },
      [&] {
        auto const numbers =
            open::sector_ground_quantum_numbers(sites, uni20::from_twice(static_cast<std::int64_t>(sites % 2)));
        return open::solve_real<Real>(sites, delta, numbers, solver);
      });
  result.sz = sz;
  result.delta = delta;
  result.window = window;
  for (auto& level : result.levels)
  {
    level.state.sz = sz;
    level.state.spin_reversed = sz.twice() < 0;
  }
  if (result.first_unconverged)
  {
    result.first_unconverged->sz = sz;
    result.first_unconverged->spin_reversed = sz.twice() < 0;
  }
  return result;
}
} // namespace open
} // namespace bethe::xxz
