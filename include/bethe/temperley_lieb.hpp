// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/real_excitations.hpp>
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

/// A finite-real-root level of H_TL=-sum e_i,
/// e_i^2=loop_weight*e_i, loop_weight>2. Representation multiplicities are
/// separate from this spectral problem; no physical spin or momentum is inferred.
template <uni20::Real Real> struct OpenState
{
    Real loop_weight{}, energy{};
    std::size_t through_lines = 0;
    xxz::quantum_group::RealState<Real> reference;
};
template <uni20::Real Real> using OpenGroundState = OpenState<Real>;

namespace detail
{
template <uni20::Real Real> void validate_loop_weight(Real loop_weight)
{
  if (!uni20::isfinite(loop_weight) || loop_weight <= Real{2})
    throw std::invalid_argument("massive TL solver requires finite loop weight > 2");
}

template <uni20::Real Real>
OpenState<Real> from_reference(xxz::quantum_group::RealState<Real> reference, Real loop_weight)
{
  Real const energy = Real{2} * reference.energy_shift;
  if (!uni20::isfinite(energy)) throw std::overflow_error("TL energy overflow");
  auto const ell = reference.sites - 2 * reference.quantum_numbers.size();
  return {.loop_weight = loop_weight, .energy = energy, .through_lines = ell, .reference = std::move(reference)};
}
} // namespace detail

template <uni20::Real Real = double>
[[nodiscard]] OpenState<Real> solve_real(std::size_t sites, Real loop_weight, std::span<uni20::half_int const> numbers,
                                         SolverOptions<Real> const& options = {})
{
  detail::validate_loop_weight(loop_weight);
  return detail::from_reference(xxz::quantum_group::solve_real<Real>(sites, loop_weight / Real{2}, numbers, options),
                                loop_weight);
}

template <uni20::Real Real = double>
[[nodiscard]] OpenState<Real> sector_ground_state(std::size_t sites, Real loop_weight, std::size_t through_lines,
                                                  SolverOptions<Real> const& options = {})
{
  detail::validate_loop_weight(loop_weight);
  return detail::from_reference(
      xxz::quantum_group::sector_ground_state(sites, loop_weight / Real{2}, through_lines, options), loop_weight);
}

template <uni20::Real Real = double>
[[nodiscard]] OpenGroundState<Real> open_ground_state(std::size_t sites, Real loop_weight,
                                                      SolverOptions<Real> const& options = {})
{
  return sector_ground_state(sites, loop_weight, 0, options);
}

/// Completeness refers only to the positive finite-real family in one TL module.
/// Distinct Bethe solutions are kept separate, including numerical energy ties.
template <typename State> struct RealExcitationScan
{
    State ground_state;
    std::vector<bethe::RealExcitation<State>> levels;
    std::size_t candidate_count = 0, converged_count = 0;
    std::optional<State> first_unconverged;
    bool family_converged() const { return candidate_count == converged_count; }
    bool converged() const { return family_converged() && ground_state.reference.converged; }
};

namespace detail
{
template <typename Scan, typename Map, typename Real> auto map_real_scan(Scan input, Map map, Real energy_scale)
{
  RealExcitationScan<decltype(map(std::move(input.ground_state)))> result;
  result.ground_state = map(std::move(input.ground_state));
  result.candidate_count = input.candidate_count;
  result.converged_count = input.converged_count;
  result.levels.reserve(input.levels.size());
  for (auto& level : input.levels)
    result.levels.push_back({.state = map(std::move(level.state)),
                             .gap = level.gap ? std::optional<Real>{energy_scale * *level.gap} : std::nullopt});
  if (input.first_unconverged) result.first_unconverged = map(std::move(*input.first_unconverged));
  return result;
}
} // namespace detail

/// Allocation-free preflight: choose(N-M,M), M=(N-ell)/2; not the TL module dimension.
[[nodiscard]] inline std::size_t real_excitation_count(std::size_t sites, std::size_t through_lines,
                                                       std::size_t limit = std::numeric_limits<std::size_t>::max())
{
  auto const m = xxz::quantum_group::detail::sector_roots(sites, through_lines);
  return bethe::detail::bounded_binomial(sites - m, m, limit);
}

/// Scan increasing integer labels in [1,N-M], retaining the lowest count converged
/// levels, including the module minimum. Gaps are to the GLOBAL ground state.
/// No complex roots or quantum-group descendants: weights are attached separately.
template <uni20::Real Real = double>
[[nodiscard]] RealExcitationScan<OpenState<Real>>
real_excitations(std::size_t sites, Real loop_weight, std::size_t through_lines, RealExcitationOptions const& scan = {},
                 SolverOptions<Real> const& solver = {})
{
  detail::validate_loop_weight(loop_weight);
  auto const m = xxz::quantum_group::detail::sector_roots(sites, through_lines);
  Real const delta = loop_weight / Real{2};
  auto result = bethe::detail::scan_real_combinations<bethe::RealExcitationScan<xxz::quantum_group::RealState<Real>>>(
      sites - m, m, 2, scan,
      [&](xxz::QuantumNumbers const& numbers) {
        return xxz::quantum_group::solve_real<Real>(sites, delta, numbers, solver);
      },
      [&] { return xxz::quantum_group::ground_state(sites, delta, solver); });
  return detail::map_real_scan(
      std::move(result), [loop_weight](auto state) { return detail::from_reference(std::move(state), loop_weight); },
      Real{2});
}
} // namespace bethe::temperley_lieb
