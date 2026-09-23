// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/temperley_lieb.hpp>
#include <bethe/xxz_open_two_string.hpp>

namespace bethe::temperley_lieb::two_string
{
template <uni20::Real Real> struct State
{
    Real loop_weight{}, energy{};
    xxz::quantum_group::two_string::State<Real> reference;
};

template <uni20::Real Real = double>
[[nodiscard]] State<Real> singlet(std::size_t sites, Real loop_weight, SolverOptions<Real> const& options = {})
{
  temperley_lieb::detail::validate_loop_weight(loop_weight);
  auto reference = xxz::quantum_group::two_string::singlet<Real>(sites, loop_weight / Real{2}, options);
  Real const energy = Real{2} * reference.energy - Real(sites - 1) * loop_weight / Real{4};
  return {.loop_weight = loop_weight, .energy = energy, .reference = std::move(reference)};
}
} // namespace bethe::temperley_lieb::two_string
