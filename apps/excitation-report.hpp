// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "spin-output.hpp"
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
                             DataOutputOptions const& options, std::string program, int argc, char** argv)
{
  using State = std::remove_cvref_t<decltype(scan.ground_state)>;
  using Real = std::remove_cvref_t<decltype(scan.ground_state.energy)>;
  auto const family = info.family;
  auto const returned_label = info.multiplet_size ? "Returned multiplets" : "Returned states";
  auto const ordering =
      scan.family_converged() ? "complete within supported family" : "incomplete; failed candidates excluded";
  auto const spin = uni20::to_string(info.sector);
  auto const& ground = scan.ground_state;
  report.field("Family", family).field(std::string(info.sector_label), spin);
  if (info.multiplet_size) report.field("Multiplet size", *info.multiplet_size);
  report.field("Candidates", scan.candidate_count)
      .field("Converged candidates", scan.converged_count)
      .field(returned_label, scan.levels.size())
      .field("Ordering", ordering)
      .field("Ground energy", uni20::format_real(ground.energy))
      .field("Ground converged", ground.converged ? 1 : 0)
      .field("Ground status", ground.converged ? "converged" : "unconverged estimate")
      .field("Ground residual", uni20::format_real(ground.residual_norm))
      .field("Ground iterations", ground.iterations)
      .field("Gap reference",
             ground.converged ? "E-E0; global ground state" : "unavailable; ground solve failed; gaps unavailable");
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
  report.field("Status", scan.converged() ? "converged" : "incomplete scan or ground reference");
  ResultOutput output(report, options, std::move(program), argc, argv,
                      spin_tables<State>(roots, true, true, scan.first_unconverged.has_value()));
  std::vector<State const*> states;
  std::vector<std::optional<Real>> gaps;
  for (auto const& level : scan.levels)
  {
    states.push_back(&level.state);
    gaps.push_back(level.gap);
  }
  spin_rows<Real>(output, "states", info.multiplet_size ? "Real-root multiplet energies" : "Real-root state energies",
                  states, gaps);
  spin_rows<Real>(output, "reference", "Ground reference", std::vector{&ground}, {}, states.size());
  states.push_back(&ground);
  if (scan.first_unconverged)
  {
    spin_rows<Real>(output, "failed", "First failed state (unranked estimate)", std::vector{&*scan.first_unconverged},
                    {}, states.size());
    states.push_back(&*scan.first_unconverged);
  }
  spin_roots<Real>(output, states, roots, true);
  output.finish();
  return scan.converged();
}
} // namespace bethe::cli
