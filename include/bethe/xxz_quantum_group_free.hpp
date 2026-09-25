// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <array>
#include <bethe/real_excitations.hpp>
#include <bethe/xxz_common.hpp>

namespace bethe::xxz::quantum_group::free
{
/// One Jordan block of H=sum(SxSx+SySy)+i/2*(Sz_1-Sz_N).
/// Distinct records can have equal energies; block_size is not degeneracy.
template <uni20::Real Real> struct Block
{
    Real energy{};
    std::vector<std::size_t> modes; // Occupied dispersive k, epsilon_k=cos(pi*k/N).
    std::size_t zero_occupation = 0, down = 0, block_size = 1;
};

/// Occupation description of a single block, not its spin-basis Jordan vectors.
/// Odd N: one ordinary zero mode. Even N: one two-dimensional defective
/// zero-mode space; occupation 1 produces J_2(0), occupations 0 and 2 J_1(0).
template <uni20::Real Real = double>
Block<Real> block(std::size_t sites, std::span<std::size_t const> modes, std::size_t zero_occupation)
{
  xxz::detail::checked_sites(sites);
  auto const zeros = sites % 2 ? 1u : 2u;
  if (zero_occupation > zeros) throw std::invalid_argument("invalid zero-mode occupation");
  bethe::detail::CompensatedSum<Real> energy;
  Real const pi = Real{4} * std::atan(Real{1});
  for (std::size_t i = 0; i < modes.size(); ++i)
  {
    auto const k = modes[i];
    if (k == 0 || k >= sites || (sites % 2 == 0 && k == sites / 2) || (i && k <= modes[i - 1]))
      throw std::invalid_argument("dispersive modes must increase in 1..N-1, excluding even-N k=N/2");
    // Enforce particle-hole paired energies using the same native cosine.
    auto const reflected = k > sites / 2;
    Real const value = std::cos(pi * Real(reflected ? sites - k : k) / Real(sites));
    energy.add(reflected ? -value : value);
  }
  return {energy.value(),
          {modes.begin(), modes.end()},
          zero_occupation,
          modes.size() + zero_occupation,
          zeros == 2 && zero_occupation == 1 ? 2u : 1u};
}

struct SectorOptions
{
    std::size_t max_blocks = 100000;
    std::size_t max_mode_entries = 1000000;
};

/// Complete fixed-down-spin sector, ordered by zero occupation then by mode
/// combination, NOT by energy. Coincident energies are deliberately not merged.
/// Reject oversized output before allocating: no partial spectrum is returned.
/// Reuses the common bounded counting and lexicographic combination helpers.
template <uni20::Real Real = double>
std::vector<Block<Real>> sector(std::size_t sites, std::size_t down, SectorOptions const& options = {})
{
  xxz::detail::checked_sites(sites);
  if (down > sites) throw std::invalid_argument("down-spin count must lie in 0..N");
  std::size_t const zeros = sites % 2 ? 1 : 2, slots = sites - zeros;
  std::array<std::size_t, 3> counts{};
  std::size_t total = 0, entries = 0;
  for (std::size_t z = 0; z <= zeros; ++z)
    if (down >= z && down - z <= slots)
    {
      auto const m = down - z;
      std::size_t count;
      try
      {
        count = bethe::detail::bounded_binomial(slots, m, options.max_blocks - total);
      }
      catch (std::length_error const&)
      {
        throw std::length_error("free quantum-group sector exceeds max_blocks");
      }
      if (m && count > (options.max_mode_entries - entries) / m)
        throw std::length_error("free quantum-group sector exceeds max_mode_entries");
      entries += count * m;
      total += count;
      counts[z] = count;
    }
  std::vector<Block<Real>> result;
  result.reserve(total);
  for (std::size_t z = 0; z <= zeros; ++z)
    if (counts[z])
    {
      std::vector<std::size_t> indices(down - z), modes(down - z);
      std::iota(indices.begin(), indices.end(), std::size_t{1});
      do
      {
        for (std::size_t i = 0; i < indices.size(); ++i)
          modes[i] = indices[i] + (zeros == 2 && indices[i] >= sites / 2 ? 1 : 0);
        result.push_back(block<Real>(sites, modes, z));
      }
      while (bethe::detail::advance_combination<std::size_t>(indices, slots));
    }
  return result;
}
} // namespace bethe::xxz::quantum_group::free
