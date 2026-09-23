// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include "excitation-report.hpp"
#include "heisenberg-cli.hpp"

namespace bethe::cli
{
template <uni20::Real Real>
report_builder report_header(std::size_t sites, std::string_view precision,
                             heisenberg::SolverOptions<Real> const& options, std::string_view mode,
                             std::string_view cpu_time, bool periodic = true)
{
  report_builder report("Heisenberg XXX - " + std::string(mode));
  report.field("Model", periodic ? "periodic spin-1/2, J=1, h=0" : "open spin-1/2, free ends, J=1, h=0")
      .field("Sites", sites)
      .field("Precision", precision)
      .field("Residual tolerance", uni20::format_real(options.residual_tolerance))
      .field("CPU time", cpu_time);
  return report;
}

template <uni20::Real Real>
bool print_spinons(std::size_t sites, std::string_view precision, heisenberg::SolverOptions<Real> const& solver,
                   std::vector<heisenberg::SpinonState<Real>> const& branch, bool roots, std::string_view cpu_time,
                   DataOutputOptions const& options, int argc, char** argv)
{
  using State = heisenberg::RealState<Real>;
  auto report = report_header(sites, precision, solver, "one-spinon branch", cpu_time);
  report.field("Sz", "0.5")
      .field("Spinon momentum", "k = pi/2 - 2*pi*hole/N")
      .field("Lattice momentum", "P = 2*pi*momentum_index/N (mod 2*pi)")
      .field("Bulk reference", "e_inf = 1/4 - log(2)")
      .field("Thermodynamic curve", "epsilon_inf = (pi/2)*sin(k)")
      .field("Finite-size energy", "E - N*e_inf retains finite-size corrections");
  std::vector<State const*> states;
  std::size_t converged = 0;
  for (auto const& point : branch)
  {
    states.push_back(&point.state);
    converged += point.state.converged;
  }
  add_scan_status(report, converged, branch.size());
  report.field("Status", converged == branch.size() ? "converged" : "incomplete; unconverged estimates");
  auto names = spin_tables<State>(roots, false);
  names.push_back("spinons");
  ResultOutput output(report, options, "bethe-xxx-pbc", argc, argv, names);
  spin_rows<Real>(output, "states", "Lattice energies and momenta", states, {});
  output.table(
      "spinons", "Spinon dispersion",
      [&](auto& table) {
        for (std::size_t i = 0; i < branch.size(); ++i)
        {
          auto const& point = branch[i];
          table.append(i, point.hole, point.spinon_momentum, point.bulk_subtracted_energy,
                       heisenberg::spinon_energy(point.spinon_momentum));
        }
      },
      column<std::size_t>("state_id"), column<uni20::half_int>("hole", "Hole"), column<Real>("k"),
      column<Real>("bulk_subtracted_energy", "E - N*e_inf"), column<Real>("epsilon_inf", "epsilon_inf"));
  spin_roots<Real>(output, states, roots, false);
  output.finish();
  return converged == branch.size();
}
} // namespace bethe::cli
