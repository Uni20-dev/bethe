// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/temperley_lieb.hpp>

namespace bethe::biquadratic
{
/// Even free-end spin-1 chain, H=-sum (S_i.S_(i+1))^2, coefficient -1.
/// Unique singlet ground state (TL through-lines=0, multiplicity=1).
/// Check reference.converged before interpreting energy as a ground energy.
template <uni20::Real Real> struct GroundState
{
    Real energy{}, tl_energy{};
    xxz::quantum_group::GroundState<Real> reference;
};

/// e_i=(S_i.S_(i+1))^2-1=3*P_singlet: H=H_TL-(N-1).
/// Only even N>=2, free ends, and the singlet ground state are supported.
/// There is no mapping to the zero-boundary-field XXZ solver or its physical Sz labels.
template <uni20::Real Real = double>
[[nodiscard]] GroundState<Real> ground_state(std::size_t sites, SolverOptions<Real> const& options = {})
{
  auto tl = temperley_lieb::open_ground_state(sites, Real{3}, options);
  return {.energy = tl.energy - Real(sites - 1), .tl_energy = tl.energy, .reference = std::move(tl.reference)};
}
} // namespace bethe::biquadratic
