// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "report.hpp"

namespace bethe::apps::biquadratic
{
template <uni20::Real Real>
int run_real(Arguments const& args, bethe::SolverOptions<Real> const& options, cli::RunReport report)
{
  auto& context = report.context();
  auto computation = context.computation();
  if (args.excitations)
  {
    auto const ell = args.through_lines.value_or(args.ferromagnetic ? args.sites - 2 : 2);
    bethe::RealExcitationOptions const enumeration{.count = *args.excitations,
                                                   .max_candidates = args.max_candidates.value_or(10000)};
    auto const scan =
        args.real_window     ? model::ferromagnetic::real_excitations_window<Real>(args.sites, ell, *args.real_window,
                                                                                   enumeration, options)
        : args.ferromagnetic ? model::ferromagnetic::real_excitations<Real>(args.sites, ell, enumeration, options)
                             : model::real_excitations<Real>(args.sites, ell, enumeration, options);
    computation.finish();
    auto const& ground = scan.ground_state;
    report.field("calculation", "Calculation", "restricted real-root excitations")
        .field("family", "Family", "positive finite roots; complex-root levels excluded; NOT the complete TL spectrum")
        .field("tl_through_lines", "TL through-lines", ell)
        .field("tl_defects", "TL defects M", (args.sites - ell) / 2)
        .field("multiplicity_per_tl_eigenvector", "Multiplicity per TL eigenvector",
               bethe::temperley_lieb::spin_chain_multiplicity(3, ell), {.missing = "overflow (>uint64)"})
        .field("multiplicity_meaning", "Multiplicity meaning",
               "physical states, not SU(2) multiplets; not an accidental-degeneracy sum")
        .field("candidates", "Candidates", scan.candidate_count)
        .field("converged_candidates", "Converged candidates", scan.converged_count)
        .field("returned_levels", "Returned levels", scan.levels.size())
        .field("ordering", "Ordering",
               !scan.family_converged() ? "incomplete; failed candidates excluded"
               : args.real_window       ? "complete within selected label window; NOT a global excitation ranking"
                                        : "complete within supported family")
        .field("ground_energy", "Ground energy", ground.energy)
        .field("ground_status", "Ground status", status(ground.reference.status))
        .field("ground_residual", "Ground residual", ground.reference.residual_norm)
        .field("ground_iterations", "Ground iterations", ground.reference.iterations)
        .field("gap_reference", "Gap reference",
               args.ferromagnetic           ? "E-(N-1); exact degenerate ferro ground space"
               : ground.reference.converged ? "E-E0; global singlet ground state"
                                            : "unavailable; ground solve failed")
        .result(scan.converged(), scan.converged() ? "converged" : "incomplete scan or ground reference");
    if (args.real_window)
    {
      auto const m = (args.sites - ell) / 2, last = args.sites - m;
      report.field("real_window_width", "Real-label window width", *args.real_window)
          .field("real_window_first", "First allowed I", last - *args.real_window + 1)
          .field("real_window_last", "Last allowed I", last)
          .field("coverage", "Coverage",
                 "selected real-root scattering window only; bound and mixed-string branches excluded")
          .field("energy_ordering", "Energy ordering",
                 "direct excitation energy; independent of the extensive total-energy offset");
    }
    if (scan.first_unconverged)
      report
          .field("first_failed_i", "First failed I",
                 cli::quantum_number_text(scan.first_unconverged->reference.quantum_numbers))
          .field("first_failed_status", "First failed status", status(scan.first_unconverged->reference.status))
          .field("first_failed_residual", "First failed residual", scan.first_unconverged->reference.residual_norm);
    std::vector<model::State<Real> const*> states;
    std::vector<std::optional<Real>> gaps;
    for (auto const& level : scan.levels)
    {
      states.push_back(&level.state);
      gaps.push_back(level.gap);
    }
    write_output(report, args, states, gaps, &ground, scan.first_unconverged ? &*scan.first_unconverged : nullptr);
    return finish(scan.converged());
  }
  if (args.sectors)
  {
    std::vector<model::State<Real>> states;
    states.push_back(model::ground_state<Real>(args.sites, options));
    bool converged = states.front().reference.converged;
    for (std::size_t ell = 2; ell <= args.sites; ell += 2)
    {
      states.push_back(model::sector_ground_state<Real>(args.sites, ell, options));
      converged = converged && states.back().reference.converged;
    }
    computation.finish();
    report.field("calculation", "Calculation", "TL module minima (not physical-spin sectors)")
        .field("multiplicity_meaning", "Multiplicity meaning",
               "physical states per TL eigenvector, not SU(2) multiplets")
        .result(converged, converged ? "converged" : "incomplete; unconverged estimates");
    std::vector<model::State<Real> const*> rows;
    std::vector<std::optional<Real>> gaps;
    for (auto const& state : states)
    {
      rows.push_back(&state);
      gaps.push_back(state.reference.converged && states.front().reference.converged
                         ? std::optional<Real>{Real{2} * (state.reference.energy - states.front().reference.energy)}
                         : std::nullopt);
    }
    write_output(report, args, rows, gaps);
    return finish(converged);
  }
  auto const state =
      args.numbers ? (args.ferromagnetic ? model::ferromagnetic::solve_real<Real>(args.sites, *args.numbers, options)
                                         : model::solve_real<Real>(args.sites, *args.numbers, options))
                   : model::sector_ground_state<Real>(args.sites, args.through_lines.value_or(0), options);
  computation.finish();
  auto const& reference = state.reference;
  report
      .status(reference.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning,
              status(reference.status))
      .field("calculation", "Calculation",
             args.numbers               ? "specified real-root TL level"
             : state.through_lines == 0 ? "even-chain singlet ground state"
                                        : "TL module minimum")
      .field("tl_through_lines", "TL through-lines", state.through_lines)
      .field("tl_defects", "TL defects M", (args.sites - state.through_lines) / 2)
      .field("multiplicity",
             !args.ferromagnetic && state.through_lines == 0 ? "Ground-state multiplicity"
                                                             : "Multiplicity per TL eigenvector",
             state.multiplicity, {.missing = "overflow (>uint64)"})
      .result(reference.converged, status(reference.status))
      .field("total_energy", "Total energy", state.energy)
      .field("energy_per_site", "Energy per site", state.energy / Real(args.sites))
      .field("tl_energy_sum_e_i", args.ferromagnetic ? "TL energy (+sum e_i)" : "TL energy (-sum e_i)", state.tl_energy)
      .field("xxz_reference_energy", "XXZ reference energy", reference.energy)
      .field("residual_norm", "Residual norm", reference.residual_norm)
      .field("iterations", "Iterations", reference.iterations);
  if (state.through_lines == 0) report.field("total_spin", "Total spin", 0);
  if (state.through_lines != 0)
    report.field("multiplicity_meaning", "Multiplicity meaning", "physical states, not a physical-spin label");
  if (args.ferromagnetic)
    report.field("gap_reference", "Gap reference", "E-(N-1); exact degenerate ferro ground space");
  write_output<Real>(report, args, {&state},
                     {args.ferromagnetic && reference.converged ? std::optional<Real>{state.tl_energy} : std::nullopt});
  return finish(reference.converged);
}
} // namespace bethe::apps::biquadratic
