// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include "report-common.hpp"
#include <bethe/xxz_open.hpp>

namespace bethe::cli::xxz_open_report
{
inline std::string_view residual_description(bethe::xxz::open::GroundResidualConvention convention)
{
  using Convention = bethe::xxz::open::GroundResidualConvention;
  switch (convention)
  {
    case Convention::logarithmic_phase:
      return "max|F|/(2*N), logarithmic phase";
    case Convention::massive_regularized:
      return "normalized bulk and regularized boundary equations";
    case Convention::negative_rank_scaled:
      return "rank-subtracted equations divided by N*s; s=sqrt((1+Delta)/(1-Delta))";
  }
  throw std::logic_error("unknown open XXZ residual convention");
}

inline std::string_view status_code(bethe::xxz::open::GroundSolveStatus status)
{
  using Status = bethe::xxz::open::GroundSolveStatus;
  switch (status)
  {
    case Status::converged:
      return "converged";
    case Status::iteration_limit:
      return "iteration_limit";
    case Status::stalled:
      return "stalled";
  }
  throw std::logic_error("unknown open XXZ solve status");
}

template <typename State> std::string_view status_text(State const& state)
{
  if (state.converged) return "converged";
  if constexpr (requires { state.status; })
    if (state.status == bethe::xxz::open::GroundSolveStatus::stalled)
      return "line search stalled; unconverged estimate";
  return "iteration limit reached; unconverged estimate";
}

template <uni20::Real Real> std::string_view boundary_kind(bethe::xxz::open::BoundaryRoot<Real> const& root)
{
  if (root.inverse_square > Real{0}) return "real";
  if (root.inverse_square < Real{0}) return "imaginary";
  return "infinity";
}

template <uni20::Real Real> void print_roots(bethe::xxz::open::GroundState<Real> const& state)
{
  if (state.residual_convention == bethe::xxz::open::GroundResidualConvention::negative_rank_scaled)
  {
    std::cout << "# index rapidity I lambda\n";
    for (std::size_t i = 0; i < state.rapidities.size(); ++i)
      std::cout << i << ' ' << uni20::format_real(state.rapidities[i]) << ' ' << state.quantum_numbers[i] << ' '
                << uni20::format_real(state.log_rapidities[i]) << '\n';
  }
  else
    cli::print_roots(state);
  if (state.boundary_root)
  {
    auto const& root = *state.boundary_root;
    std::cout << "# boundary I inverse_square log_distance kind\n"
              << "# boundary " << root.quantum_number << ' ' << uni20::format_real(root.inverse_square) << ' '
              << uni20::format_real(root.log_distance) << ' ' << boundary_kind(root) << '\n';
  }
}

template <uni20::Real Real>
void add_roots(report_builder& report, bethe::xxz::open::GroundState<Real> const& state, std::string const& title)
{
  if (!state.log_rapidities.empty())
  {
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
  else if (state.rapidities.empty() && state.boundary_root)
    report.table(title + " (no bulk roots)");
  else
    cli::add_roots(report, state, title);
  if (state.boundary_root)
  {
    auto const& root = *state.boundary_root;
    report.table(title + " - boundary root")
        .header_separator()
        .column("I")
        .column("Inverse square y", table_alignment::decimal)
        .column("Log distance w", table_alignment::decimal)
        .column("Kind of z", table_alignment::left)
        .row(uni20::to_string_fraction(root.quantum_number), uni20::format_real(root.inverse_square),
             uni20::format_real(root.log_distance), boundary_kind(root));
  }
}

template <uni20::Real Real, typename State>
bool print_state(report_builder report, std::size_t sites, State const& state, bool roots)
{
  return cli::print_state_report<Real>(std::move(report), sites, state, roots);
}

template <uni20::Real Real>
bool print_state(report_builder report, std::size_t sites, bethe::xxz::open::GroundState<Real> const& state, bool roots)
{
  if (state.residual_convention == bethe::xxz::open::GroundResidualConvention::logarithmic_phase)
    return cli::print_state_report<Real>(std::move(report), sites, state, roots);
  report.status(state.converged ? semantic_glyph::success : semantic_glyph::warning, std::string(status_text(state)))
      .field("Sz", uni20::to_string_fraction(state.sz))
      .field("Reference vacuum", state.spin_reversed ? "all down (spin reversed)" : "all up")
      .field("Iterations", state.iterations)
      .field("Root Delta", uni20::format_real(state.root_delta))
      .field("Boundary root", state.boundary_root ? boundary_kind(*state.boundary_root) : "none")
      .field("Residual norm", uni20::format_real(state.residual_norm))
      .field("Total energy", uni20::format_real(state.energy))
      .field("Energy per site", uni20::format_real(state.energy / Real(sites)));
  if (roots) xxz_open_report::add_roots(report, state, "Rapidities");
  cli::print_report(report);
  return state.converged;
}

template <uni20::Real Real>
bool print_sectors(report_builder report, std::vector<bethe::xxz::open::GroundState<Real>> const& states, bool roots)
{
  if (states.front().residual_convention == bethe::xxz::open::GroundResidualConvention::logarithmic_phase)
    return cli::print_sector_report(std::move(report), states, roots);
  auto& energies = report.table("Sector energies");
  energies.header_separator().column("Sz").column("Energy", table_alignment::decimal);
  auto& diagnostics = report.table("Convergence");
  diagnostics.header_separator()
      .column("Sz")
      .column("Residual")
      .column("Iterations")
      .column("Root Delta")
      .column("Status", table_alignment::left)
      .column("Boundary root", table_alignment::left);
  std::size_t converged = 0;
  for (auto const& state : states)
  {
    auto const sz = uni20::to_string_fraction(state.sz);
    energies.row(sz, uni20::format_real(state.energy));
    diagnostics.row(sz, uni20::format_real(state.residual_norm), state.iterations, uni20::format_real(state.root_delta),
                    status_code(state.status), state.boundary_root ? boundary_kind(*state.boundary_root) : "none");
    if (roots)
      xxz_open_report::add_roots(report, state,
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
} // namespace bethe::cli::xxz_open_report
