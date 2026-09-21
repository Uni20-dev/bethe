// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include <bethe/heisenberg.hpp>
#include <uni20/common/presentation.hpp>

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

// Presentation belongs to the front end, not the numerical library. Pass all
// real-valued cells through format_real before handing them to the renderer:
// this retains max_digits10 for the selected type, including binary128.
namespace bethe::cli
{
namespace presentation = uni20::presentation;
using presentation::report_builder;
using presentation::semantic_glyph;
using presentation::table_alignment;

inline auto convergence_status(bool converged)
{
  return converged ? presentation::style("Green")(semantic_glyph::success, "converged")
                   : presentation::style("Yellow;Bold")(semantic_glyph::warning, "unconverged estimate");
}

// Our tables have one labeled column per cell, with no spans or separator rows.
// Prefer vertical records when a table would exceed the terminal width. Never
// ask the table renderer to hard-wrap a high-precision numeric token.
inline void print_report(report_builder const& report)
{
  auto policy = presentation::terminal_policy(stdout);
  policy.wrap_width = std::nullopt;
  auto const columns = static_cast<std::size_t>(terminal::columns(stdout));
  report_builder header(report.title());
  for (auto const& [glyph, label] : report.statuses())
    header.status(glyph, label);
  for (auto const& [key, value] : report.fields())
    header.field(key, value);
  std::cout << presentation::render_terminal(header, policy);

  for (auto const& table : report.tables())
  {
    report_builder block;
    block.table("") = table;
    auto const plain = presentation::render_plain(block, policy);
    bool fits = true;
    for (std::string_view remaining = plain; !remaining.empty();)
    {
      auto const newline = remaining.find('\n');
      fits = fits && presentation::display_width(remaining.substr(0, newline), policy) <= columns;
      if (newline == std::string_view::npos)
        break;
      remaining.remove_prefix(newline + 1);
    }
    std::cout << '\n';
    if (fits)
    {
      std::cout << presentation::render_terminal(block, policy);
      continue;
    }

    presentation::styled_text records;
    records.append(presentation::style("Cyan")(table.title())).append("\n");
    std::size_t label_width = 0;
    for (auto const& column : table.columns())
      label_width = std::max(label_width, presentation::display_width(column.heading, policy));
    bool first = true;
    for (auto const& entry : table.entries())
    {
      auto const& row = std::get<std::vector<presentation::table_cell>>(entry);
      if (!first)
        records.append("\n");
      first = false;
      for (std::size_t i = 0; i < row.size(); ++i)
        records.append("  ")
            .append(presentation::style("LightGray")(
                presentation::pad_right(table.columns()[i].heading, label_width + 2, policy)))
            .append(row[i].content)
            .append("\n");
    }
    std::cout << presentation::render_terminal(records, policy);
  }
}

template <uni20::Real Real>
report_builder report_header(std::size_t sites, std::string_view precision,
                             heisenberg::SolverOptions<Real> const& options, std::string_view mode)
{
  report_builder report("Heisenberg XXX - " + std::string(mode));
  report.field("Model", "periodic spin-1/2, J=1, h=0")
      .field("Sites", sites)
      .field("Precision", precision)
      .field("Residual tolerance", uni20::format_real(options.residual_tolerance));
  return report;
}

template <uni20::Real Real>
void add_roots(report_builder& report, heisenberg::RealState<Real> const& state, std::string title)
{
  if (state.rapidities.empty())
  {
    report.table(std::move(title) + " (none; fully polarized)");
    return;
  }
  auto& table = report.table(std::move(title));
  table.header_separator().column("Index").column("I").column("Rapidity", table_alignment::decimal);
  for (std::size_t i = 0; i < state.rapidities.size(); ++i)
    table.row(i, uni20::to_string_fraction(state.quantum_numbers[i]), uni20::format_real(state.rapidities[i]));
}

template <uni20::Real Real>
bool print_state(std::size_t sites, std::string_view precision, heisenberg::SolverOptions<Real> const& options,
                 heisenberg::RealState<Real> const& state, std::string_view mode, bool roots)
{
  auto report = report_header(sites, precision, options, mode);
  report
      .status(state.converged ? semantic_glyph::success : semantic_glyph::warning,
              state.converged ? "converged" : "iteration limit reached; unconverged estimate")
      .field("Sz", uni20::to_string_fraction(state.sz))
      .field("Reference vacuum", state.spin_reversed ? "all down (spin reversed)" : "all up")
      .field("Momentum index", state.momentum_index)
      .field("Momentum P", uni20::format_real(state.momentum))
      .field("Iterations", state.iterations)
      .field("Residual norm", uni20::format_real(state.residual_norm))
      .field("Total energy", uni20::format_real(state.energy))
      .field("Energy per site", uni20::format_real(state.energy / static_cast<Real>(sites)));
  if (roots)
    add_roots(report, state, "Rapidities");
  print_report(report);
  return state.converged;
}

inline void add_scan_status(report_builder& report, std::size_t converged, std::size_t total)
{
  bool const all = converged == total;
  report.status(all ? semantic_glyph::success : semantic_glyph::warning,
                std::to_string(converged) + "/" + std::to_string(total) + " states converged" +
                    (all ? "" : "; remaining energies are unconverged estimates"));
}

template <uni20::Real Real>
bool print_sectors(std::size_t sites, std::string_view precision, heisenberg::SolverOptions<Real> const& options,
                   std::vector<heisenberg::RealState<Real>> const& states, bool roots)
{
  auto report = report_header(sites, precision, options, "sector minima");
  report.field("Momentum convention", "P = 2*pi*momentum_index/N (mod 2*pi)");
  auto& energies = report.table("Sector energies and momenta");
  energies.header_separator()
      .column("Sz")
      .column("Momentum index")
      .column("P", table_alignment::decimal)
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
    diagnostics.row(sz, uni20::format_real(state.residual_norm), state.iterations, convergence_status(state.converged));
    converged += state.converged;
    if (roots)
      add_roots(report, state, "Rapidities: Sz=" + sz + (state.spin_reversed ? " (spin reversed)" : ""));
  }
  add_scan_status(report, converged, states.size());
  print_report(report);
  return converged == states.size();
}

template <uni20::Real Real>
bool print_spinons(std::size_t sites, std::string_view precision, heisenberg::SolverOptions<Real> const& options,
                   std::vector<heisenberg::SpinonState<Real>> const& branch, bool roots)
{
  auto report = report_header(sites, precision, options, "one-spinon branch");
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
    if (roots)
      add_roots(report, state, "Rapidities: hole=" + hole);
  }
  add_scan_status(report, converged, branch.size());
  print_report(report);
  return converged == branch.size();
}
} // namespace bethe::cli
