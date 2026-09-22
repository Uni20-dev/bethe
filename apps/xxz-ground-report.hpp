// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include "report-common.hpp"
#include <bethe/xxz_ground_state.hpp>

namespace bethe::cli::xxz_ground_report
{
inline std::string_view status_code(bethe::xxz::GroundSolveStatus status)
{
  using Status = bethe::xxz::GroundSolveStatus;
  switch (status)
  {
    case Status::converged:
      return "converged";
    case Status::iteration_limit:
      return "iteration_limit";
    case Status::stalled:
      return "stalled";
  }
  throw std::logic_error("unknown periodic XXZ solve status");
}

template <typename State> std::string_view status_text(State const& state)
{
  if (state.converged) return "converged";
  if constexpr (requires { state.status; })
    if (state.status == bethe::xxz::GroundSolveStatus::stalled) return "line search stalled; unconverged estimate";
  return "iteration limit reached; unconverged estimate";
}

template <uni20::Real Real> void print_roots(bethe::xxz::GroundState<Real> const& state)
{
  if (state.residual_convention == bethe::xxz::GroundResidualConvention::logarithmic_phase)
    return cli::print_roots(state);
  std::cout << "# index rapidity I lambda\n";
  for (std::size_t i = 0; i < state.rapidities.size(); ++i)
    std::cout << i << ' ' << uni20::format_real(state.rapidities[i]) << ' ' << state.quantum_numbers[i] << ' '
              << uni20::format_real(state.log_rapidities[i]) << '\n';
}

template <uni20::Real Real>
void add_roots(report_builder& report, bethe::xxz::GroundState<Real> const& state, std::string const& title)
{
  if (state.log_rapidities.empty()) return cli::add_roots(report, state, title);
  auto& table = report.table(title);
  table.header_separator()
      .column("Index")
      .column("I")
      .column("Rapidity z", table_alignment::decimal)
      .column("Lambda", table_alignment::decimal);
  for (std::size_t i = 0; i < state.rapidities.size(); ++i)
    table.row(i, uni20::to_string_fraction(state.quantum_numbers[i]), uni20::format_real(state.rapidities[i]),
              uni20::format_real(state.log_rapidities[i]));
}

template <uni20::Real Real, typename State>
bool print_state(report_builder report, std::size_t sites, State const& state, bool roots)
{
  return cli::print_state_report<Real>(std::move(report), sites, state, roots);
}

template <uni20::Real Real>
bool print_state(report_builder report, std::size_t sites, bethe::xxz::GroundState<Real> const& state, bool roots)
{
  if (state.residual_convention == bethe::xxz::GroundResidualConvention::logarithmic_phase)
    return cli::print_state_report<Real>(std::move(report), sites, state, roots);
  report.status(state.converged ? semantic_glyph::success : semantic_glyph::warning, std::string(status_text(state)))
      .field("Sz", uni20::to_string_fraction(state.sz))
      .field("Reference vacuum", state.spin_reversed ? "all down (spin reversed)" : "all up")
      .field("Momentum index", state.momentum_index)
      .field("Momentum P", uni20::format_real(state.momentum))
      .field("Iterations", state.iterations)
      .field("Residual norm", uni20::format_real(state.residual_norm))
      .field("Total energy", uni20::format_real(state.energy))
      .field("Energy per site", uni20::format_real(state.energy / Real(sites)));
  if (roots) xxz_ground_report::add_roots(report, state, "Rapidities");
  cli::print_report(report);
  return state.converged;
}

template <uni20::Real Real>
bool print_sectors(report_builder report, std::vector<bethe::xxz::GroundState<Real>> const& states, bool roots)
{
  if (states.front().residual_convention == bethe::xxz::GroundResidualConvention::logarithmic_phase)
    return cli::print_sector_report(std::move(report), states, roots);
  report.field("Momentum convention", "P = 2*pi*momentum_index/N (mod 2*pi)");
  auto& energies = report.table("Sector energies");
  energies.header_separator()
      .column("Sz")
      .column("Momentum index")
      .column("Momentum P", table_alignment::decimal)
      .column("Energy", table_alignment::decimal);
  auto& diagnostics = report.table("Convergence");
  diagnostics.header_separator()
      .column("Sz")
      .column("Residual")
      .column("Iterations")
      .column("Status", table_alignment::left);
  std::size_t converged = 0;
  for (auto const& state : states)
  {
    auto const sz = uni20::to_string_fraction(state.sz);
    energies.row(sz, state.momentum_index, uni20::format_real(state.momentum), uni20::format_real(state.energy));
    diagnostics.row(sz, uni20::format_real(state.residual_norm), state.iterations, status_code(state.status));
    if (roots)
      xxz_ground_report::add_roots(report, state,
                                   "Rapidities: Sz=" + sz + (state.spin_reversed ? " (spin reversed)" : ""));
    converged += state.converged;
  }
  cli::add_scan_status(report, converged, states.size());
  cli::print_report(report);
  return converged == states.size();
}

inline int finish(bool converged)
{
  if (!converged)
    std::cerr << "The calculation did not converge; each nonconverged energy is an unconverged estimate.\n";
  return converged ? 0 : 2;
}
} // namespace bethe::cli::xxz_ground_report
