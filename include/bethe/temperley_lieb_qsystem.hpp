// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/temperley_lieb.hpp>
#include <bethe/xxz_open_qsystem.hpp>

namespace bethe::temperley_lieb::qsystem
{
template <uni20::Real Real> struct State
{
    Real loop_weight{};
    std::optional<Real> energy;
    std::size_t through_lines = 0;
    xxz::quantum_group::qsystem::State<Real> reference;
};

template <uni20::Real Real> [[nodiscard]] State<Real> from_reference(xxz::quantum_group::qsystem::State<Real> reference)
{
  Real const loop = Real{2} * reference.delta;
  auto const ell = reference.through_lines;
  auto const energy =
      reference.energy ? std::optional<Real>{Real{2} * *reference.energy - Real(reference.sites - 1) * loop / Real{4}}
                       : std::nullopt;
  return {.loop_weight = loop, .energy = energy, .through_lines = ell, .reference = std::move(reference)};
}

template <uni20::Real Real = double>
[[nodiscard]] State<Real> solve(std::size_t sites, Real loop_weight, std::span<Real const> seed,
                                SolverOptions<Real> const& options = {})
{
  temperley_lieb::detail::validate_loop_weight(loop_weight);
  return from_reference(xxz::quantum_group::qsystem::solve<Real>(sites, loop_weight / Real{2}, seed, options));
}
} // namespace bethe::temperley_lieb::qsystem
