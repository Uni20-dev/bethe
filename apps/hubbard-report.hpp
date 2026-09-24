// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "result-output.hpp"
#include <bethe/hubbard_common.hpp>

namespace bethe::cli
{
template <uni20::Real Real, typename State>
int print_hubbard_state(State const& state, uni20::half_int sz, std::string_view precision,
                        DataOutputOptions const& output_options, Real tolerance, uni20::run_context& context,
                        bool roots)
{
  constexpr bool periodic = requires { state.momentum_index; };
  auto const model_description =
      periodic ? "periodic Hubbard, t=1, U*n_up*n_down" : "free-end Hubbard, t=1, U*n_up*n_down";
  auto const status = state.converged                                 ? "converged"
                      : state.status == hubbard::SolveStatus::stalled ? "line search stalled; unconverged estimate"
                                                                      : "iteration limit reached; unconverged estimate";
  auto const method = state.free_fermion ? "exact free fermions"
                                         : (periodic ? "Lieb-Wu (damped Newton + continuation)"
                                                     : "open-chain Bethe ansatz (damped Newton + continuation)");
  std::string mapping;
  if (state.shiba_transformed) mapping += "Shiba (down-spin particle-hole); ";
  if (state.particle_hole_transformed) mapping += "full particle-hole; ";
  if (state.spin_reversed) mapping += "spin reversal; ";
  if (mapping.empty())
    mapping = "none";
  else
    mapping.resize(mapping.size() - 2);
  bethe::cli::RunReport report(context, periodic ? "Hubbard (periodic) - sector ground state"
                                                 : "Hubbard (free ends) - sector ground state");
  report.status(state.converged ? semantic_glyph::success : semantic_glyph::warning, status)
      .result(state.converged, status)
      .field("model", "Model", model_description)
      .field("sites", "Sites", state.sites)
      .field("particles", "Particles", state.particles)
      .field("down_spins", "Down spins", state.down_spins)
      .field("sz", "Sz", sz, {.fractions = true})
      .field("u", "U", state.interaction)
      .field("precision", "Precision", precision)
      .field("method", "Method", method)
      .field("symmetry_mapping", "Symmetry mapping", mapping)
      .field("residual_tolerance", "Residual tolerance", tolerance)
      .field("iterations", "Iterations", state.iterations)
      .field("completed_continuation_stages", "Completed continuation stages", state.continuation_steps)
      .field("charge_residual", "Charge residual", state.charge_residual)
      .field("spin_residual", "Spin residual", state.spin_residual)
      .field("residual_norm", "Residual norm", state.residual_norm);
  if constexpr (periodic)
    report.field("momentum_index", "Momentum index", state.momentum_index)
        .field("momentum_p", "Momentum P", state.momentum);
  report.field("total_energy", "Total energy", state.energy)
      .field("energy_per_site", "Energy per site", state.energy / Real(state.sites));
  if (state.auxiliary_roots())
  {
    report.field("roots_and_residuals", "Roots and residuals", "auxiliary sector (not physical-sector Bethe roots)")
        .field("root_particles", "Root particles", state.root_particles)
        .field("root_down_spins", "Root down spins", state.root_down_spins)
        .field("root_u", "Root U", state.root_interaction)
        .field("energy_offset", "Energy offset", state.energy_offset);
    if constexpr (periodic) report.field("momentum_index_offset", "Momentum index offset", state.momentum_offset);
  }
  std::vector<std::string> names{"states"};
  if (roots)
  {
    names.push_back(state.free_fermion ? "free_modes" : "charge_roots");
    if (!state.free_fermion) names.push_back("spin_roots");
  }
  ResultOutput output(report, output_options, names);
  auto const state_columns = std::tuple{
      column<std::size_t>("state_id"),   column<uni20::half_int>("sz", "Sz"), column<Real>("energy", "Energy"),
      column<Real>("charge_residual"),   column<Real>("spin_residual"),       column<Real>("residual"),
      column<std::size_t>("iterations"), column<std::size_t>("stages"),       column<bool>("auxiliary_roots"),
      column<bool>("converged"),         column<std::string>("status")};
  std::apply(
      [&](auto... columns) {
        if constexpr (periodic)
          output.table(
              "states", "State",
              [&](auto& table) {
                table.append(0, sz, state.energy, state.charge_residual, state.spin_residual, state.residual_norm,
                             state.iterations, state.continuation_steps, state.auxiliary_roots(), state.converged,
                             std::string(status), state.momentum_index, state.momentum);
              },
              columns..., column<std::size_t>("momentum_index"), column<Real>("p", "P"));
        else
          output.table(
              "states", "State",
              [&](auto& table) {
                table.append(0, sz, state.energy, state.charge_residual, state.spin_residual, state.residual_norm,
                             state.iterations, state.continuation_steps, state.auxiliary_roots(), state.converged,
                             std::string(status));
              },
              columns...);
      },
      state_columns);
  if (roots)
  {
    std::string const prefix = state.auxiliary_roots() ? "Auxiliary: " : "";
    if (state.free_fermion)
      output.table(
          "free_modes", prefix + "Free-fermion occupied momenta (no Bethe labels)",
          [&](auto& table) {
            for (std::size_t j = 0; j < state.charge_momenta.size(); ++j)
              table.append(0, j, state.charge_momenta[j]);
          },
          column<std::size_t>("state_id"), column<std::size_t>("index", "Index"), column<Real>("k"));
    else
    {
      output.table(
          "charge_roots", prefix + (periodic ? "Charge momenta" : "Charge wave numbers (standing waves)"),
          [&](auto& table) {
            for (std::size_t j = 0; j < state.charge_momenta.size(); ++j)
              table.append(0, j, state.charge_momenta[j], state.quantum_numbers.charge[j]);
          },
          column<std::size_t>("state_id"), column<std::size_t>("index", "Index"), column<Real>("k"),
          column<uni20::half_int>("quantum_number", "I"));
      output.table(
          "spin_roots", prefix + "Spin rapidities (conventional Lambda)",
          [&](auto& table) {
            for (std::size_t j = 0; j < state.spin_rapidities.size(); ++j)
              table.append(0, j, state.spin_rapidities[j], state.quantum_numbers.spin[j]);
          },
          column<std::size_t>("state_id"), column<std::size_t>("index", "Index"), column<Real>("rapidity", "Lambda"),
          column<uni20::half_int>("quantum_number", "J"));
    }
  }
  output.finish();
  if (!state.converged) std::cerr << "Hubbard solve: " << status << "; consider a larger budget or higher precision.\n";
  return state.converged ? 0 : 2;
}
} // namespace bethe::cli
