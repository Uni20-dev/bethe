// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/biquadratic.hpp>
#include <bethe/temperley_lieb_qsystem.hpp>

namespace bethe::biquadratic::qsystem
{
template <uni20::Real Real> struct State
{
    std::optional<Real> energy, tl_energy;
    std::size_t through_lines = 0;
    std::optional<std::uint64_t> multiplicity;
    xxz::quantum_group::qsystem::State<Real> reference;
    Exchange exchange = Exchange::antiferromagnetic;
};

template <uni20::Real Real> [[nodiscard]] State<Real> from_reference(xxz::quantum_group::qsystem::State<Real> reference)
{
  if (reference.delta != Real{3} / Real{2}) throw std::invalid_argument("biquadratic Q-system requires Delta=3/2");
  auto tl = temperley_lieb::qsystem::from_reference(std::move(reference));
  auto const energy = tl.energy ? std::optional<Real>{*tl.energy - Real(tl.reference.sites - 1)} : std::nullopt;
  return {.energy = energy,
          .tl_energy = tl.energy,
          .through_lines = tl.through_lines,
          .multiplicity = temperley_lieb::spin_chain_multiplicity(3, tl.through_lines),
          .reference = std::move(tl.reference)};
}

template <uni20::Real Real = double>
[[nodiscard]] State<Real> solve(std::size_t sites, std::span<Real const> seed, SolverOptions<Real> const& options = {})
{
  return from_reference(xxz::quantum_group::qsystem::solve<Real>(sites, Real{3} / Real{2}, seed, options));
}

template <uni20::Real Real> struct Spectrum
{
    std::vector<State<Real>> states;
    std::optional<std::size_t> expected_count;
    std::size_t attempts = 0, failed_attempts = 0;
    bool complete() const { return expected_count && states.size() == *expected_count; }
};

template <uni20::Real Real = double>
[[nodiscard]] Spectrum<Real> spectrum(std::size_t sites, std::size_t through_lines,
                                      xxz::quantum_group::qsystem::SearchOptions const& search = {},
                                      SolverOptions<Real> const& options = {})
{
  auto input = xxz::quantum_group::qsystem::spectrum<Real>(sites, Real{3} / Real{2}, through_lines, search, options);
  Spectrum<Real> result{.states = {},
                        .expected_count = input.expected_count,
                        .attempts = input.attempts,
                        .failed_attempts = input.failed_attempts};
  for (auto& state : input.states)
    result.states.push_back(from_reference(std::move(state)));
  return result;
}
} // namespace bethe::biquadratic::qsystem
