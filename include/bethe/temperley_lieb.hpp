// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz_quantum_group.hpp>
#include <optional>

namespace bethe::temperley_lieb
{
/// Multiplicity of an open TL through-line module in the generic spin-chain
/// singlet-projector representation with local dimension d>=2.
/// m_0=1, m_1=d, m_(ell+1)=d*m_ell-m_(ell-1). Not a physical SU(2) spin label
/// or a universal rule for other TL representations (e.g. RSOS quotients).
/// Exact uint64, or nullopt on overflow. Coincident levels in different modules
/// acquire the sum of their multiplicities, not just one module's value.
[[nodiscard]] inline std::optional<std::uint64_t> spin_chain_multiplicity(std::uint64_t local_dimension,
                                                                          std::size_t through_lines)
{
  if (local_dimension < 2) throw std::invalid_argument("TL spin representation requires local dimension >= 2");
  if (through_lines == 0) return 1;
  if (local_dimension == 2)
  {
    if (through_lines >= std::numeric_limits<std::uint64_t>::max()) return std::nullopt;
    return std::uint64_t(through_lines) + 1;
  }
  std::uint64_t previous = 1, current = local_dimension;
  for (std::size_t i = 1; i < through_lines; ++i)
  {
    // (d-1)*current + (current-previous) avoids an overflowing d*current
    // intermediate when the final subtraction would still fit.
    auto const remainder = current - previous;
    if (current > (std::numeric_limits<std::uint64_t>::max() - remainder) / (local_dimension - 1)) return std::nullopt;
    auto const next = (local_dimension - 1) * current + remainder;
    previous = current;
    current = next;
  }
  return current;
}

/// The even-chain zero-through-line ground level of H_TL=-sum e_i,
/// e_i^2=loop_weight*e_i, loop_weight>2. Representation multiplicities are
/// separate from this spectral problem; no physical spin or momentum is inferred.
template <uni20::Real Real> struct OpenGroundState
{
    Real loop_weight{}, energy{};
    xxz::quantum_group::GroundState<Real> reference;
};

template <uni20::Real Real = double>
[[nodiscard]] OpenGroundState<Real> open_ground_state(std::size_t sites, Real loop_weight,
                                                      SolverOptions<Real> const& options = {})
{
  if (!uni20::isfinite(loop_weight) || loop_weight <= Real{2})
    throw std::invalid_argument("massive TL ground state requires finite loop weight > 2");
  auto reference = xxz::quantum_group::ground_state(sites, loop_weight / Real{2}, options);
  Real const energy = Real{2} * reference.energy - Real(sites - 1) * (loop_weight / Real{4});
  if (!uni20::isfinite(energy)) throw std::overflow_error("TL energy overflow");
  return {.loop_weight = loop_weight, .energy = energy, .reference = std::move(reference)};
}
} // namespace bethe::temperley_lieb
