// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "report-common.hpp"
#include <bethe/hubbard_common.hpp>

namespace bethe::cli
{
template <uni20::Real Real, typename State>
int print_hubbard_state(State const& state, uni20::half_int sz, std::string_view precision, std::string_view format,
                        Real tolerance, std::string_view cpu_time, bool roots)
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
  bool const pretty = format == "pretty" || (format == "auto" && terminal::is_a_terminal(stdout));
  if (pretty)
  {
    using bethe::cli::semantic_glyph;
    using bethe::cli::table_alignment;
    bethe::cli::report_builder report(periodic ? "Hubbard (periodic) - sector ground state"
                                               : "Hubbard (free ends) - sector ground state");
    report.status(state.converged ? semantic_glyph::success : semantic_glyph::warning, status)
        .field("Model", model_description)
        .field("Sites", state.sites)
        .field("Particles", state.particles)
        .field("Down spins", state.down_spins)
        .field("Sz", uni20::to_string_fraction(sz))
        .field("U", uni20::format_real(state.interaction))
        .field("Precision", precision)
        .field("Method", method)
        .field("Symmetry mapping", mapping)
        .field("Residual tolerance", uni20::format_real(tolerance))
        .field("CPU time", cpu_time)
        .field("Iterations", state.iterations)
        .field("Completed continuation stages", state.continuation_steps)
        .field("Charge residual", uni20::format_real(state.charge_residual))
        .field("Spin residual", uni20::format_real(state.spin_residual))
        .field("Residual norm", uni20::format_real(state.residual_norm));
    if constexpr (periodic)
      report.field("Momentum index", state.momentum_index).field("Momentum P", uni20::format_real(state.momentum));
    report.field("Total energy", uni20::format_real(state.energy))
        .field("Energy per site", uni20::format_real(state.energy / Real(state.sites)));
    if (state.auxiliary_roots())
    {
      report.field("Roots and residuals", "auxiliary sector (not physical-sector Bethe roots)")
          .field("Root particles", state.root_particles)
          .field("Root down spins", state.root_down_spins)
          .field("Root U", uni20::format_real(state.root_interaction))
          .field("Energy offset", uni20::format_real(state.energy_offset));
      if constexpr (periodic) report.field("Momentum index offset", state.momentum_offset);
    }
    if (roots)
    {
      auto& charge =
          report.table(std::string(state.auxiliary_roots() ? "Auxiliary: " : "") +
                       (state.free_fermion ? "Free-fermion occupied momenta (no Bethe labels)"
                                           : (periodic ? "Charge momenta" : "Charge wave numbers (standing waves)")));
      charge.header_separator().column("Index").column("k", table_alignment::decimal);
      if (!state.free_fermion) charge.column("I");
      for (std::size_t j = 0; j < state.charge_momenta.size(); ++j)
        if (state.free_fermion)
          charge.row(j, uni20::format_real(state.charge_momenta[j]));
        else
          charge.row(j, uni20::format_real(state.charge_momenta[j]),
                     uni20::to_string_fraction(state.quantum_numbers.charge[j]));
      if (!state.free_fermion)
      {
        auto& spin = report.table(std::string(state.auxiliary_roots() ? "Auxiliary: " : "") +
                                  "Spin rapidities (conventional Lambda)");
        spin.header_separator().column("Index").column("Lambda", table_alignment::decimal).column("J");
        for (std::size_t a = 0; a < state.spin_rapidities.size(); ++a)
          spin.row(a, uni20::format_real(state.spin_rapidities[a]),
                   uni20::to_string_fraction(state.quantum_numbers.spin[a]));
      }
    }
    bethe::cli::print_report(report);
  }
  else
  {
    std::cout << "Sites: " << state.sites << "\nModel: " << model_description << '\n'
              << "Particles: " << state.particles << "\nDown spins: " << state.down_spins
              << "\nSz: " << uni20::to_string_fraction(sz) << '\n'
              << "U: " << uni20::format_real(state.interaction) << "\nPrecision: " << precision << '\n'
              << "Method: " << method << "\nSymmetry mapping: " << mapping
              << "\nResidual tolerance: " << uni20::format_real(tolerance) << '\n'
              << "Status: " << status << "\nIterations: " << state.iterations << "\nCPU time: " << cpu_time << '\n'
              << "Completed continuation stages: " << state.continuation_steps << '\n'
              << "Charge residual: " << uni20::format_real(state.charge_residual) << '\n'
              << "Spin residual: " << uni20::format_real(state.spin_residual) << '\n'
              << "Residual norm: " << uni20::format_real(state.residual_norm) << '\n';
    if constexpr (periodic)
      std::cout << "Momentum index: " << state.momentum_index << "\nMomentum: " << uni20::format_real(state.momentum)
                << '\n';
    std::cout << "Total energy: " << uni20::format_real(state.energy) << '\n'
              << "Energy per site: " << uni20::format_real(state.energy / Real(state.sites)) << '\n';
    if (state.auxiliary_roots())
    {
      std::cout << "Roots and residuals: auxiliary sector (not physical-sector Bethe roots)\n"
                << "Root particles: " << state.root_particles << "\nRoot down spins: " << state.root_down_spins
                << "\nRoot U: " << uni20::format_real(state.root_interaction)
                << "\nEnergy offset: " << uni20::format_real(state.energy_offset) << '\n';
      if constexpr (periodic) std::cout << "Momentum index offset: " << state.momentum_offset << '\n';
    }
    if (roots)
    {
      if (state.auxiliary_roots()) std::cout << "# Auxiliary-sector roots/occupations follow\n";
      std::cout << (state.free_fermion ? "# index k (free-fermion occupations; no Bethe labels)\n" : "# index k I\n");
      for (std::size_t j = 0; j < state.charge_momenta.size(); ++j)
      {
        std::cout << j << ' ' << uni20::format_real(state.charge_momenta[j]);
        if (!state.free_fermion) std::cout << ' ' << state.quantum_numbers.charge[j];
        std::cout << '\n';
      }
      if (!state.free_fermion)
      {
        std::cout << "# index Lambda J\n";
        for (std::size_t a = 0; a < state.spin_rapidities.size(); ++a)
          std::cout << a << ' ' << uni20::format_real(state.spin_rapidities[a]) << ' ' << state.quantum_numbers.spin[a]
                    << '\n';
      }
    }
  }
  if (!state.converged) std::cerr << "Hubbard solve: " << status << "; consider a larger budget or higher precision.\n";
  return state.converged ? 0 : 2;
}
} // namespace bethe::cli
