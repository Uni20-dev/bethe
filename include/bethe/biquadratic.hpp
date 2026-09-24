// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/temperley_lieb.hpp>

namespace bethe::biquadratic
{
enum class Exchange
{
  antiferromagnetic,
  ferromagnetic
};

/// Free-end spin-1 chain. Default API: H=-sum (S_i.S_(i+1))^2.
/// Selected real-root states allow odd/even N; ground-state helpers require even N.
/// The ferromagnetic API reverses physical and TL energies, not the XXZ reference.
/// TL labels are NOT physical spin. The AF ell=0 ground state is a unique singlet.
/// Check reference.converged before interpreting energy as an eigenvalue.
template <uni20::Real Real> struct State
{
    Real energy{}, tl_energy{};
    std::size_t through_lines = 0;
    /// Copies per TL eigenvector, not a count of distinct energies or SU(2) multiplets.
    /// nullopt means uint64 overflow, never zero or an approximate multiplicity.
    std::optional<std::uint64_t> multiplicity;
    xxz::quantum_group::RealState<Real> reference;
    Exchange exchange = Exchange::antiferromagnetic;
};
template <uni20::Real Real> using GroundState = State<Real>;

namespace detail
{
template <uni20::Real Real> State<Real> from_tl(temperley_lieb::OpenState<Real> tl)
{
  return {.energy = tl.energy - Real(tl.reference.sites - 1),
          .tl_energy = tl.energy,
          .through_lines = tl.through_lines,
          .multiplicity = temperley_lieb::spin_chain_multiplicity(3, tl.through_lines),
          .reference = std::move(tl.reference)};
}
} // namespace detail

template <uni20::Real Real = double>
[[nodiscard]] State<Real> solve_real(std::size_t sites, std::span<uni20::half_int const> numbers,
                                     SolverOptions<Real> const& options = {})
{
  return detail::from_tl(temperley_lieb::solve_real<Real>(sites, Real{3}, numbers, options));
}

/// Lowest energy in one TL module, even ell in [0,N]. Not a physical-spin sector minimum.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> sector_ground_state(std::size_t sites, std::size_t through_lines,
                                              SolverOptions<Real> const& options = {})
{
  return detail::from_tl(temperley_lieb::sector_ground_state(sites, Real{3}, through_lines, options));
}

/// e_i=(S_i.S_(i+1))^2-1=3*P_singlet: H=H_TL-(N-1).
/// Even N>=2, free ends, singlet ground state.
/// There is no mapping to the zero-boundary-field XXZ solver or its physical Sz labels.
template <uni20::Real Real = double>
[[nodiscard]] GroundState<Real> ground_state(std::size_t sites, SolverOptions<Real> const& options = {})
{
  return sector_ground_state<Real>(sites, 0, options);
}

using temperley_lieb::real_excitation_count;
template <uni20::Real Real> using RealExcitationScan = temperley_lieb::RealExcitationScan<State<Real>>;

/// Positive finite-real family in one TL module, with exact representation weights.
/// Includes the module minimum; gaps use the global singlet ground state.
template <uni20::Real Real = double>
[[nodiscard]] RealExcitationScan<Real> real_excitations(std::size_t sites, std::size_t through_lines,
                                                        RealExcitationOptions const& scan = {},
                                                        SolverOptions<Real> const& solver = {})
{
  return temperley_lieb::detail::map_real_scan(
      temperley_lieb::real_excitations(sites, Real{3}, through_lines, scan, solver),
      [](auto state) { return detail::from_tl(std::move(state)); }, Real{1});
}
} // namespace bethe::biquadratic
