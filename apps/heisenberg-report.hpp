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

template <uni20::Real Real, typename State>
bool print_state(std::size_t sites, std::string_view precision, heisenberg::SolverOptions<Real> const& options,
                 State const& state, std::string_view mode, bool roots, std::string_view cpu_time)
{
  constexpr bool periodic = requires { state.momentum; };
  return print_state_report<Real>(report_header(sites, precision, options, mode, cpu_time, periodic), sites, state,
                                  roots);
}

template <uni20::Real Real, typename State>
bool print_sectors(std::size_t sites, std::string_view precision, heisenberg::SolverOptions<Real> const& options,
                   std::vector<State> const& states, bool roots, std::string_view cpu_time)
{
  constexpr bool periodic = requires(State state) { state.momentum; };
  return print_sector_report(report_header(sites, precision, options, "sector minima", cpu_time, periodic), states,
                             roots);
}

template <uni20::Real Real>
bool print_spinons(std::size_t sites, std::string_view precision, heisenberg::SolverOptions<Real> const& options,
                   std::vector<heisenberg::SpinonState<Real>> const& branch, bool roots, std::string_view cpu_time)
{
  auto report = report_header(sites, precision, options, "one-spinon branch", cpu_time);
  report.field("Sz", "1/2")
      .field("Spinon momentum", "k = pi/2 - 2*pi*hole/N")
      .field("Lattice momentum", "P = 2*pi*momentum_index/N (mod 2*pi)")
      .field("Bulk reference", "e_inf = 1/4 - log(2)")
      .field("Thermodynamic curve", "epsilon_inf = (pi/2)*sin(k)")
      .field("Finite-size energy", "E - N*e_inf retains finite-size corrections");
  auto& dispersion = report.table("Spinon dispersion");
  dispersion.header_separator()
      .column("Hole")
      .column("k", table_alignment::decimal)
      .column("E - N*e_inf", table_alignment::decimal)
      .column("epsilon_inf", table_alignment::decimal);
  auto& energies = report.table("Lattice energies and momenta");
  energies.header_separator()
      .column("Hole")
      .column("Momentum index")
      .column("P", table_alignment::decimal)
      .column("Energy", table_alignment::decimal);
  auto& diagnostics = report.table("Convergence");
  diagnostics.header_separator()
      .column("Hole")
      .column("Residual")
      .column("Iterations")
      .column("Status", table_alignment::left);
  std::size_t converged = 0;
  for (auto const& point : branch)
  {
    auto const& state = point.state;
    auto const hole = uni20::to_string_fraction(point.hole);
    dispersion.row(hole, uni20::format_real(point.spinon_momentum), uni20::format_real(point.bulk_subtracted_energy),
                   uni20::format_real(heisenberg::spinon_energy(point.spinon_momentum)));
    energies.row(hole, state.momentum_index, uni20::format_real(state.momentum), uni20::format_real(state.energy));
    diagnostics.row(hole, uni20::format_real(state.residual_norm), state.iterations,
                    convergence_status(state.converged));
    converged += state.converged;
    if (roots) add_roots(report, state, "Rapidities: hole=" + hole);
  }
  add_scan_status(report, converged, branch.size());
  print_report(report);
  return converged == branch.size();
}

template <uni20::Real Real, typename State>
bool print_excitations(std::size_t sites, std::string_view precision, heisenberg::SolverOptions<Real> const& options,
                       heisenberg::RealExcitationScan<State> const& scan, bool roots, std::string_view cpu_time,
                       bool pretty)
{
  constexpr bool periodic = requires(State state) { state.momentum; };
  return print_excitation_report(report_header(sites, precision, options, "real-root excitations", cpu_time, periodic),
                                 scan,
                                 {.family = "restricted real-root highest-weight multiplets; NOT a complete spectrum",
                                  .sector_label = "S",
                                  .sector = scan.spin,
                                  .multiplet_size = scan.spin.twice() + 1},
                                 roots, cpu_time, pretty);
}
} // namespace bethe::cli
