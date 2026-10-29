// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/temperley_lieb_two_string.hpp>

namespace bethe::biquadratic::two_string
{
template <uni20::Real Real> struct State
{
    Real energy{}, tl_energy{};
    std::size_t through_lines = 0;
    std::uint64_t multiplicity = 1;
    xxz::quantum_group::two_string::State<Real> reference;
};

/// Selected complex-root singlet branch, not an enumeration of singlets.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> singlet(std::size_t sites, SolverOptions<Real> const& options = {})
{
  auto tl = temperley_lieb::two_string::singlet<Real>(sites, uni20::scalar_at_precision<Real>(3, options.precision()),
                                                      options);
  return {.energy = tl.energy - Real(sites - 1), .tl_energy = tl.energy, .reference = std::move(tl.reference)};
}
} // namespace bethe::biquadratic::two_string
