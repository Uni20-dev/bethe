// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "result-output.hpp"

namespace bethe::cli
{
// Shared decay-mode schema. Model-specific continuation/root diagnostics live
// in metadata and the model's own root table, not in energy-labelled columns.
template <typename State> void relaxation_table(ResultOutput& output, State const& state, std::string status)
{
  using Real = typename decltype(state.gap)::value_type;
  using Optional = std::optional<Real>;
  Optional const real = state.eigenvalue ? Optional(state.eigenvalue->real()) : std::nullopt;
  Optional const imag = state.eigenvalue ? Optional(state.eigenvalue->imag()) : std::nullopt;
  output.table(
      "relaxation", "Leading relaxation mode",
      [&](auto& table) {
        table.append(state.effective_particles != 0, real, imag, state.gap, imag, state.residual_norm, state.iterations,
                     state.seed_iterations, state.converged, status);
      },
      column<bool>("has_mode"), column<Optional>("lambda_real"), column<Optional>("lambda_imag"),
      column<Optional>("gap"), column<Optional>("frequency"), column<Real>("residual"),
      column<std::size_t>("iterations"), column<std::size_t>("seed_iterations"), column<bool>("converged"),
      column<std::string>("status"));
}
} // namespace bethe::cli
