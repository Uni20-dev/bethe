// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "report.hpp"

namespace bethe::apps::biquadratic
{
template <uni20::Real Real> int run_ferro_analytic(Arguments const& args, cli::RunReport report)
{
  auto computation = report.context().computation();
  auto const ground = model::ferromagnetic::ground_space<Real>(args.sites);
  std::vector<model::ferromagnetic::OneDefectLevel<Real>> levels;
  if (args.one_defect)
  {
    if (args.sites - 1 > args.max_candidates.value_or(10000))
      throw std::length_error("one-defect band exceeds max_candidates; raise --max-candidates");
    levels.reserve(args.sites - 1);
    for (std::size_t j = args.sites - 1; j > 0; --j)
      levels.push_back(model::ferromagnetic::one_defect_level<Real>(args.sites, j));
  }
  else if (args.through_lines && *args.through_lines == args.sites - 2)
    levels.push_back(model::ferromagnetic::one_defect_level<Real>(args.sites, args.sites - 1));
  computation.finish();
  report
      .field("calculation", "Calculation",
             args.one_defect  ? "complete one-defect TL module"
             : levels.empty() ? "exact ferro ground space"
                              : "exact one-defect module minimum")
      .field("coverage", "Coverage",
             args.one_defect  ? "all N-1 eigenvalues of ell=N-2; NOT the full excited spectrum"
             : levels.empty() ? "ground space only; zero modes are not individually enumerated"
                              : "first positive level above the entire ground space")
      .field("wave_number", "Wave number", "k=pi*j/N is an OBC standing-wave coordinate, not lattice momentum")
      .field("ground_energy", "Ground energy", ground.energy)
      .field("ground_multiplicity", "Ground-space dimension", ground.multiplicity, {.missing = "overflow (>uint64)"})
      .field("gap_reference", "Gap reference", "E-(N-1); exact degenerate ferro ground space")
      .field("multiplicity_meaning", "Multiplicity meaning",
             "physical states per TL eigenvector; not SU(2) multiplets or accidental-degeneracy sums")
      .result(true, "exact spectral rules");
  cli::ResultOutput output(report, args.output, table_names(args, {"states", "reference"}));
  output.table(
      "states", "Analytic ferromagnetic levels",
      [&](auto& table) {
        if (levels.empty())
          table.append(std::size_t{0}, ground.through_lines, std::size_t{0}, ground.multiplicity,
                       std::optional<std::size_t>{}, std::optional<Real>{}, ground.energy, Real{0});
        for (std::size_t i = 0; i < levels.size(); ++i)
        {
          auto const& level = levels[i];
          table.append(i, level.through_lines, std::size_t{1}, level.multiplicity,
                       std::optional<std::size_t>{level.mode}, std::optional<Real>{level.wave_number}, level.energy,
                       level.gap);
        }
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("through_lines"),
      cli::column<std::size_t>("defects"), cli::column<std::optional<std::uint64_t>>("multiplicity"),
      cli::column<std::optional<std::size_t>>("mode", "j"),
      cli::column<std::optional<Real>>("wave_number", "Standing-wave k"), cli::column<Real>("energy", "Energy"),
      cli::column<Real>("gap", "E-E0"));
  output.table(
      "reference", "Exact ferromagnetic ground space",
      [&](auto& table) {
        table.append(std::max(std::size_t{1}, levels.size()), ground.through_lines, ground.energy, ground.multiplicity);
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("through_lines"),
      cli::column<Real>("energy", "Energy"), cli::column<std::optional<std::uint64_t>>("multiplicity"));
  write_spin_content(output, args, std::max(std::size_t{1}, levels.size()) + 1,
                     [&](std::size_t i) { return i < levels.size() ? levels[i].through_lines : ground.through_lines; });
  output.finish();
  return 0;
}

} // namespace bethe::apps::biquadratic
