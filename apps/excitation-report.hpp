// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "report-common.hpp"
#include <span>

namespace bethe::cli
{
inline std::string quantum_number_text(std::span<uni20::half_int const> numbers)
{
  std::string result;
  for (auto number : numbers)
  {
    if (!result.empty()) result += ',';
    result += uni20::to_string_fraction(number);
  }
  return result.empty() ? "-" : result;
}

struct ExcitationReportInfo
{
    std::string_view family;
    std::string_view sector_label;
    uni20::half_int sector;
    std::optional<std::int64_t> multiplet_size;
};

template <typename Scan>
bool print_excitation_report(report_builder report, Scan const& scan, ExcitationReportInfo const& info, bool roots,
                             std::string_view cpu_time, bool pretty)
{
  using State = decltype(scan.ground_state);
  constexpr bool periodic = requires(State state) { state.momentum; };
  auto const family = info.family;
  auto const returned_label = info.multiplet_size ? "Returned multiplets" : "Returned states";
  auto const ordering =
      scan.family_converged() ? "complete within supported family" : "incomplete; failed candidates excluded";
  auto const spin = uni20::to_string_fraction(info.sector);
  auto const& ground = scan.ground_state;
  auto gap_text = [](auto const& level) { return level.gap ? uni20::format_real(*level.gap) : "unavailable"; };
  if (!pretty)
  {
    std::cout << "# CPU time: " << cpu_time << '\n'
              << "# Family: " << family << '\n'
              << "# " << info.sector_label << ": " << spin << '\n';
    if (info.multiplet_size) std::cout << "# Multiplet size: " << *info.multiplet_size << '\n';
    std::cout << "# Candidates: " << scan.candidate_count << '\n'
              << "# Converged candidates: " << scan.converged_count << '\n'
              << "# " << returned_label << ": " << scan.levels.size() << '\n'
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
    std::cout << "# level " << info.sector_label;
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
        std::cout << "# Roots: level=" << i + 1 << "; " << info.sector_label << "=" << spin << '\n';
        print_roots(state);
      }
    }
    return scan.converged();
  }

  report.field("Family", family).field(std::string(info.sector_label), spin);
  if (info.multiplet_size) report.field("Multiplet size", *info.multiplet_size);
  report.field("Candidates", scan.candidate_count)
      .field(returned_label, scan.levels.size())
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
  auto& energies = report.table(info.multiplet_size ? "Real-root multiplet energies" : "Real-root state energies");
  energies.header_separator().column("Level").column(std::string(info.sector_label));
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
    if (roots)
      add_roots(report, state,
                "Rapidities: level=" + std::to_string(i + 1) + "; " + std::string(info.sector_label) + "=" + spin);
  }
  print_report(report);
  return scan.converged();
}
} // namespace bethe::cli
