// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "report.hpp"
#include <bethe/biquadratic_two_string.hpp>

namespace bethe::apps::biquadratic
{
template <uni20::Real Real>
int run_singlet(Arguments const& args, bethe::SolverOptions<Real> const& options, cli::RunReport report)
{
  auto& context = report.context();
  auto computation = context.computation();
  auto const state = model::two_string::singlet<Real>(args.sites, options);
  auto const ground = model::ground_state<Real>(args.sites, options);
  computation.finish();
  auto const& r = state.reference;
  bool const complete = r.converged && ground.reference.converged;
  Real const eta = std::acosh(r.delta), d = std::exp(-r.log_deviation);
  std::optional<Real> gap;
  if (complete) gap = Real{2} * (r.energy - ground.reference.energy);
  report.field("calculation", "Calculation", "selected complex-root singlet excitation")
      .field("family", "Family", "one positive-deviation two-string; real I=1,...,N/2-2; string label 1")
      .field("ordering", "Ordering", "targeted branch, not an exhaustive search or a global first-excitation guarantee")
      .field("tl_through_lines", "TL through-lines", 0)
      .field("multiplicity_meaning", "Multiplicity meaning", "one physical singlet per TL eigenvector")
      .field("residual_convention", "Residual convention",
             "max phase/log-modulus equation residual divided by 2N; not an energy-error bound")
      .field("string_coordinate", "String coordinate", "u=(eta+d)/2 +/- i*a/2; d=exp(-L)>0; L is authoritative")
      .field("string_deviation", "String deviation",
             eta + d == eta ? "unresolved in rounded u; retained by L=-log(d)"
                            : "resolved in rounded u; L=-log(d) also retained")
      .field("gap_reference", "Gap reference",
             ground.reference.converged ? "E-E0; global singlet ground state" : "unavailable; ground solve failed")
      .result(complete, complete ? "converged" : "incomplete; unconverged estimate");
  std::vector<std::string> names{"states", "reference", "string"};
  if (args.roots) names.push_back("roots");
  cli::ResultOutput output(report, args.output, table_names(args, names));
  output.table(
      "states", "Selected two-string singlet",
      [&](auto& table) {
        table.append(std::size_t{0}, state.through_lines, state.multiplicity, state.energy, gap, state.tl_energy,
                     r.energy, r.residual_norm, r.phase_residual, r.modulus_residual, r.iterations, r.converged,
                     std::string(status(r.status)));
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("through_lines"),
      cli::column<std::uint64_t>("multiplicity"), cli::column<Real>("energy", "Energy"),
      cli::column<std::optional<Real>>("gap", "E-E0"), cli::column<Real>("tl_energy"),
      cli::column<Real>("reference_energy"), cli::column<Real>("residual"), cli::column<Real>("phase_residual"),
      cli::column<Real>("modulus_residual"), cli::column<std::size_t>("iterations"), cli::column<bool>("converged"),
      cli::column<std::string>("status", "Status"));
  output.table(
      "reference", "Global ground reference",
      [&](auto& table) {
        table.append(std::size_t{1}, ground.energy, ground.reference.residual_norm, ground.reference.converged,
                     std::string(status(ground.reference.status)));
      },
      cli::column<std::size_t>("state_id"), cli::column<Real>("energy", "Energy"), cli::column<Real>("residual"),
      cli::column<bool>("converged"), cli::column<std::string>("status", "Status"));
  output.table(
      "string", "Two-string parameters (L remains valid below root-coordinate resolution)",
      [&](auto& table) {
        table.append(std::size_t{0}, r.center, r.log_deviation, d > Real{0} ? std::optional<Real>{d} : std::nullopt);
      },
      cli::column<std::size_t>("state_id"), cli::column<Real>("center", "a"),
      cli::column<Real>("log_deviation", "L=-log(d)"),
      cli::column<std::optional<Real>>("deviation", "d (null if underflow)"));
  if (args.roots)
    output.table(
        "roots", "Reference roots (rounded u; use L for the string deviation)",
        [&](auto& table) {
          for (std::size_t j = 0; j < r.rapidities.size(); ++j)
            table.append(std::size_t{0}, j, std::string("real sea"), Real{0}, r.rapidities[j] / Real{2});
          for (std::size_t j = 0; j < 2; ++j)
            table.append(std::size_t{0}, r.rapidities.size() + j, std::string("two-string"), (eta + d) / Real{2},
                         (j ? -r.center : r.center) / Real{2});
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("index"), cli::column<std::string>("kind"),
        cli::column<Real>("u_real"), cli::column<Real>("u_imag"));
  write_spin_content(output, args, 2,
                     [&](std::size_t i) { return i == 0 ? state.through_lines : ground.through_lines; });
  output.finish();
  if (!complete) std::cerr << "Two-string singlet or ground reference unconverged; no verified gap.\n";
  return complete ? 0 : 2;
}

} // namespace bethe::apps::biquadratic
