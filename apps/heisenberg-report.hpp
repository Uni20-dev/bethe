// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include "heisenberg-cli.hpp"
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
      if (newline == std::string_view::npos) break;
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
      if (!first) records.append("\n");
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

template <typename State> void add_roots(report_builder& report, State const& state, std::string title)
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

template <uni20::Real Real, typename State>
bool print_state(std::size_t sites, std::string_view precision, heisenberg::SolverOptions<Real> const& options,
                 State const& state, std::string_view mode, bool roots, std::string_view cpu_time)
{
  constexpr bool periodic = requires { state.momentum; };
  auto report = report_header(sites, precision, options, mode, cpu_time, periodic);
  report
      .status(state.converged ? semantic_glyph::success : semantic_glyph::warning,
              state.converged ? "converged" : "iteration limit reached; unconverged estimate")
      .field("Sz", uni20::to_string_fraction(state.sz))
      .field("Reference vacuum", state.spin_reversed ? "all down (spin reversed)" : "all up");
  if constexpr (periodic)
    report.field("Momentum index", state.momentum_index).field("Momentum P", uni20::format_real(state.momentum));
  report.field("Iterations", state.iterations)
      .field("Residual norm", uni20::format_real(state.residual_norm))
      .field("Total energy", uni20::format_real(state.energy))
      .field("Energy per site", uni20::format_real(state.energy / static_cast<Real>(sites)));
  if (roots) add_roots(report, state, "Rapidities");
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

template <uni20::Real Real, typename State>
bool print_sectors(std::size_t sites, std::string_view precision, heisenberg::SolverOptions<Real> const& options,
                   std::vector<State> const& states, bool roots, std::string_view cpu_time)
{
  constexpr bool periodic = requires(State state) { state.momentum; };
  auto report = report_header(sites, precision, options, "sector minima", cpu_time, periodic);
  if constexpr (periodic) report.field("Momentum convention", "P = 2*pi*momentum_index/N (mod 2*pi)");
  auto& energies = report.table(periodic ? "Sector energies and momenta" : "Sector energies");
  energies.header_separator().column("Sz");
  if constexpr (periodic) energies.column("Momentum index").column("P", table_alignment::decimal);
  energies.column("Energy", table_alignment::decimal);
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
    if constexpr (periodic)
      energies.row(sz, state.momentum_index, uni20::format_real(state.momentum), uni20::format_real(state.energy));
    else
      energies.row(sz, uni20::format_real(state.energy));
    diagnostics.row(sz, uni20::format_real(state.residual_norm), state.iterations, convergence_status(state.converged));
    converged += state.converged;
    if (roots) add_roots(report, state, "Rapidities: Sz=" + sz + (state.spin_reversed ? " (spin reversed)" : ""));
  }
  add_scan_status(report, converged, states.size());
  print_report(report);
  return converged == states.size();
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

inline std::string quantum_number_text(heisenberg::QuantumNumbers const& numbers)
{
  std::string result;
  for (auto number : numbers)
  {
    if (!result.empty()) result += ',';
    result += uni20::to_string_fraction(number);
  }
  return result.empty() ? "-" : result;
}

template <uni20::Real Real, typename State>
bool print_excitations(std::size_t sites, std::string_view precision, heisenberg::SolverOptions<Real> const& options,
                       heisenberg::RealExcitationScan<State> const& scan, bool roots, std::string_view cpu_time,
                       bool pretty)
{
  constexpr bool periodic = requires(State state) { state.momentum; };
  constexpr std::string_view family = "restricted real-root highest-weight multiplets; NOT a complete spectrum";
  auto const ordering =
      scan.family_converged() ? "complete within supported family" : "incomplete; failed candidates excluded";
  auto const spin = uni20::to_string_fraction(scan.spin);
  auto const& ground = scan.ground_state;
  auto gap_text = [](auto const& level) { return level.gap ? uni20::format_real(*level.gap) : "unavailable"; };
  if (!pretty)
  {
    std::cout << "# CPU time: " << cpu_time << '\n'
              << "# Family: " << family << '\n'
              << "# S: " << spin << '\n'
              << "# Multiplet size: " << scan.spin.twice() + 1 << '\n'
              << "# Candidates: " << scan.candidate_count << '\n'
              << "# Converged candidates: " << scan.converged_count << '\n'
              << "# Returned multiplets: " << scan.levels.size() << '\n'
              << "# Ordering: " << ordering << '\n'
              << "# Ground energy: " << uni20::format_real(ground.energy) << '\n'
              << "# Ground converged: " << ground.converged << '\n'
              << "# Ground residual: " << uni20::format_real(ground.residual_norm) << '\n'
              << "# Ground iterations: " << ground.iterations << '\n';
    if (!ground.converged) std::cout << "# Ground reference is an unconverged estimate; gaps unavailable\n";
    if (scan.first_unconverged)
    {
      auto const& failed = *scan.first_unconverged;
      std::cout << "# First failed I: " << quantum_number_text(failed.quantum_numbers) << '\n'
                << "# First failed residual: " << uni20::format_real(failed.residual_norm) << '\n'
                << "# First failed iterations: " << failed.iterations << '\n';
    }
    std::cout << "# level S";
    if constexpr (periodic) std::cout << " momentum_index P";
    std::cout << " energy gap residual iterations converged I\n";
    for (std::size_t i = 0; i < scan.levels.size(); ++i)
    {
      auto const& level = scan.levels[i];
      auto const& state = level.state;
      std::cout << i + 1 << ' ' << spin;
      if constexpr (periodic) std::cout << ' ' << state.momentum_index << ' ' << uni20::format_real(state.momentum);
      std::cout << ' ' << uni20::format_real(state.energy) << ' ' << gap_text(level) << ' '
                << uni20::format_real(state.residual_norm) << ' ' << state.iterations << ' ' << state.converged << ' '
                << quantum_number_text(state.quantum_numbers) << '\n';
      if (roots)
      {
        std::cout << "# Roots: level=" << i + 1 << "; S=" << spin << '\n';
        print_roots(state);
      }
    }
    return scan.converged();
  }

  auto report = report_header(sites, precision, options, "real-root excitations", cpu_time, periodic);
  report.field("Family", family)
      .field("S", spin)
      .field("Multiplet size", scan.spin.twice() + 1)
      .field("Candidates", scan.candidate_count)
      .field("Returned multiplets", scan.levels.size())
      .field("Ordering", ordering)
      .field("Ground energy", uni20::format_real(ground.energy))
      .field("Ground status", ground.converged ? "converged" : "unconverged estimate")
      .field("Ground residual", uni20::format_real(ground.residual_norm))
      .field("Ground iterations", ground.iterations)
      .field("Gap reference", ground.converged ? "E-E0; global ground state" : "unavailable; ground solve failed");
  report.status(scan.family_converged() ? semantic_glyph::success : semantic_glyph::warning,
                std::to_string(scan.converged_count) + "/" + std::to_string(scan.candidate_count) +
                    " candidates converged" + (scan.family_converged() ? "" : "; failed candidates excluded"));
  if (!ground.converged) report.status(semantic_glyph::warning, "ground reference failed; gaps unavailable");
  if (scan.first_unconverged)
  {
    auto const& failed = *scan.first_unconverged;
    report.field("First failed I", quantum_number_text(failed.quantum_numbers))
        .field("First failed residual", uni20::format_real(failed.residual_norm))
        .field("First failed iterations", failed.iterations);
  }
  auto& energies = report.table("Real-root multiplet energies");
  energies.header_separator().column("Level").column("S");
  if constexpr (periodic) energies.column("Momentum index").column("P", table_alignment::decimal);
  energies.column("Energy", table_alignment::decimal).column("E-E0", table_alignment::decimal);
  auto& diagnostics = report.table("Convergence and quantum numbers");
  diagnostics.header_separator().column("Level").column("Residual").column("Iterations").column("Status").column("I");
  for (std::size_t i = 0; i < scan.levels.size(); ++i)
  {
    auto const& level = scan.levels[i];
    auto const& state = level.state;
    if constexpr (periodic)
      energies.row(i + 1, spin, state.momentum_index, uni20::format_real(state.momentum),
                   uni20::format_real(state.energy), gap_text(level));
    else
      energies.row(i + 1, spin, uni20::format_real(state.energy), gap_text(level));
    diagnostics.row(i + 1, uni20::format_real(state.residual_norm), state.iterations,
                    convergence_status(state.converged), quantum_number_text(state.quantum_numbers));
    if (roots) add_roots(report, state, "Rapidities: level=" + std::to_string(i + 1) + "; S=" + spin);
  }
  print_report(report);
  return scan.converged();
}
} // namespace bethe::cli
