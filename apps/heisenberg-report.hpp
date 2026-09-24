// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include "excitation-report.hpp"
#include "heisenberg-cli.hpp"

namespace bethe::cli
{
template <uni20::Real Real>
RunReport report_header(std::size_t sites, std::string_view precision, heisenberg::SolverOptions<Real> const& options,
                        std::string_view mode, uni20::run_context& context, bool periodic = true)
{
  RunReport report(context, "Heisenberg XXX - " + std::string(mode));
  report.field("model", "Model", periodic ? "periodic spin-1/2, J=1, h=0" : "open spin-1/2, free ends, J=1, h=0")
      .field("sites", "Sites", sites)
      .field("precision", "Precision", precision)
      .field("residual_tolerance", "Residual tolerance", options.residual_tolerance);
  return report;
}

template <uni20::Real Real>
bool print_spinons(std::size_t sites, std::string_view precision, heisenberg::SolverOptions<Real> const& solver,
                   std::vector<heisenberg::SpinonState<Real>> const& branch, bool roots, uni20::run_context& context,
                   DataOutputOptions const& options)
{
  using State = heisenberg::RealState<Real>;
  auto report = report_header(sites, precision, solver, "one-spinon branch", context);
  report.field("sz", "Sz", uni20::from_twice(1))
      .field("spinon_momentum", "Spinon momentum", "k = pi/2 - 2*pi*hole/N")
      .field("lattice_momentum", "Lattice momentum", "P = 2*pi*momentum_index/N (mod 2*pi)")
      .field("bulk_reference", "Bulk reference", "e_inf = 1/4 - log(2)")
      .field("thermodynamic_curve", "Thermodynamic curve", "epsilon_inf = (pi/2)*sin(k)")
      .field("finite_size_energy", "Finite-size energy", "E - N*e_inf retains finite-size corrections");
  std::vector<State const*> states;
  std::size_t converged = 0;
  for (auto const& point : branch)
  {
    states.push_back(&point.state);
    converged += point.state.converged;
  }
  add_scan_status(report, converged, branch.size());
  report.result(converged == branch.size(),
                converged == branch.size() ? "converged" : "incomplete; unconverged estimates");
  auto names = spin_tables<State>(roots, false);
  names.push_back("spinons");
  ResultOutput output(report, options, names);
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
