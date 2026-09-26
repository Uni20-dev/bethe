// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "excitation-report.hpp"
#include "options.hpp"
#include <bethe/biquadratic.hpp>
#include <bethe/biquadratic_ferromagnetic.hpp>
#include <map>

namespace bethe::apps::biquadratic
{
namespace model = bethe::biquadratic;
inline char const* status(bethe::xxz::quantum_group::SolveStatus value)
{
  using Status = bethe::xxz::quantum_group::SolveStatus;
  switch (value)
  {
    case Status::converged:
      return "converged";
    case Status::iteration_limit:
      return "iteration limit reached; unconverged estimate";
    case Status::singular_jacobian:
      return "singular or ill-conditioned Jacobian; unconverged estimate";
    case Status::stalled:
      return "line search or representable precision stalled; unconverged estimate";
  }
  return "unknown";
}
inline std::vector<std::string> table_names(Arguments const& args, std::vector<std::string> names)
{
  if (args.spin_content) names.push_back("spin_content");
  return names;
}

// State IDs also cover references and failed estimates. Spin content belongs to
// the module independently of convergence; it never certifies a numerical level.
template <typename ThroughLines>
void write_spin_content(cli::ResultOutput& output, Arguments const& args, std::size_t count, ThroughLines ell_of)
{
  if (!args.spin_content) return;
  output.table(
      "spin_content", "Physical SU(2) content per TL eigenvector (not spectral weights)",
      [&](auto& table) {
        std::map<std::size_t, std::optional<std::vector<std::uint64_t>>> cache;
        for (std::size_t i = 0; i < count; ++i)
        {
          auto const ell = ell_of(i);
          auto [found, inserted] = cache.try_emplace(ell);
          if (inserted) found->second = bethe::temperley_lieb::spin_one_multiplets(ell);
          auto const& counts = found->second;
          if (!counts)
            table.append(i, ell, std::optional<uni20::half_int>{}, std::optional<std::uint64_t>{},
                         std::optional<std::uint64_t>{}, std::string("total dimension exceeds uint64"));
          else
            for (std::size_t spin = 0; spin < counts->size(); ++spin)
              if ((*counts)[spin])
                table.append(i, ell, std::optional<uni20::half_int>{uni20::half_int(std::int64_t(spin))},
                             std::optional<std::uint64_t>{(*counts)[spin]},
                             std::optional<std::uint64_t>{(2 * spin + 1) * (*counts)[spin]}, std::string("exact"));
        }
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("through_lines"),
      cli::column<std::optional<uni20::half_int>>("spin", "S"), cli::column<std::optional<std::uint64_t>>("multiplets"),
      cli::column<std::optional<std::uint64_t>>("magnetic_states", "(2S+1)*multiplets"),
      cli::column<std::string>("status"));
}

template <uni20::Real Real>
void write_output(cli::RunReport& report, Arguments const& args, std::vector<model::State<Real> const*> const& states,
                  std::vector<std::optional<Real>> const& gaps, model::State<Real> const* reference = nullptr,
                  model::State<Real> const* failed = nullptr)
{
  std::vector<std::string> names{"states", "quantum_numbers"};
  if (reference) names.push_back("reference");
  if (failed) names.push_back("failed");
  if (args.roots) names.push_back("roots");
  cli::ResultOutput output(report, args.output, table_names(args, names));
  auto write_states = [&](std::string name, std::string title, auto const& rows, std::size_t offset, bool ranked) {
    output.table(
        name, title,
        [&](auto& table) {
          for (std::size_t i = 0; i < rows.size(); ++i)
          {
            auto const& s = *rows[i];
            auto const& r = s.reference;
            table.append(offset + i, s.through_lines, s.multiplicity, s.energy, ranked ? gaps[i] : std::nullopt,
                         s.tl_energy, r.energy, r.residual_norm, r.iterations, r.converged,
                         std::string(status(r.status)));
          }
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("through_lines", "Through-lines"),
        cli::column<std::optional<std::uint64_t>>("multiplicity", "Multiplicity"),
        cli::column<Real>("energy", "Energy"), cli::column<std::optional<Real>>("gap", "E-E0"),
        cli::column<Real>("tl_energy"), cli::column<Real>("reference_energy"),
        cli::column<Real>("residual", "Residual"), cli::column<std::size_t>("iterations", "Iterations"),
        cli::column<bool>("converged"), cli::column<std::string>("status", "Status"));
  };
  write_states("states",
               args.excitations ? (args.ferromagnetic ? "Real-root TL levels (complex levels excluded)"
                                                      : "Real-root TL levels (module minimum included)")
               : args.sectors   ? "TL module minima"
                                : "State",
               states, 0, true);
  auto all = states;
  if (reference)
  {
    write_states("reference", "Ground reference", std::vector{reference}, all.size(), false);
    all.push_back(reference);
  }
  if (failed)
  {
    write_states("failed", "First failed state (unranked estimate)", std::vector{failed}, all.size(), false);
    all.push_back(failed);
  }
  output.table(
      "quantum_numbers", "Reference quantum numbers (not physical spin-1 labels)",
      [&](auto& table) {
        for (std::size_t i = 0; i < all.size(); ++i)
          for (std::size_t j = 0; j < all[i]->reference.quantum_numbers.size(); ++j)
            table.append(i, j, all[i]->reference.quantum_numbers[j]);
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("index", "Index"),
      cli::column<uni20::half_int>("quantum_number", "I"));
  if (args.roots)
    output.table(
        "roots", "Reference XXZ roots (not physical spin-1 quantum numbers)",
        [&](auto& table) {
          for (std::size_t i = 0; i < all.size(); ++i)
            for (std::size_t j = 0; j < all[i]->reference.rapidities.size(); ++j)
              table.append(i, j, all[i]->reference.quantum_numbers[j], all[i]->reference.rapidities[j],
                           all[i]->reference.angles[j]);
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("index", "Index"),
        cli::column<uni20::half_int>("quantum_number", "I"), cli::column<Real>("rapidity", "alpha"),
        cli::column<Real>("angle", "x=Theta_1/2"));
  write_spin_content(output, args, all.size(), [&](std::size_t i) { return all[i]->through_lines; });
  output.finish();
}

template <uni20::Real Real>
auto preamble(uni20::run_context& context, Arguments const& args, bethe::SolverOptions<Real> const& options)
{
  cli::RunReport report(context, "Spin-1 pure biquadratic chain (free ends)");
  report
      .field("hamiltonian", "Hamiltonian", args.ferromagnetic ? "H=+sum_i (S_i.S_(i+1))^2" : "H=-sum_i (S_i.S_(i+1))^2")
      .field("sites", "Sites", args.sites)
      .field("spin", "Spin", 1)
      .field("tl_loop_weight", "TL loop weight", 3)
      .field("xxz_delta", "XXZ Delta", Real{1.5})
      .field("xxz_reference", "XXZ reference", "spin-half exchange 1; +sqrt(5)/4*(sz_1-sz_N)")
      .field("precision", "Precision", args.precision)
      .field("residual_tolerance", "Residual tolerance", options.residual_tolerance);
  if (args.spin_content)
    report.field("spin_content_convention", "Spin content",
                 "physical SU(2) multiplets per TL eigenvector, independent of solver convergence; "
                 "not spectral weights; null when total dimension exceeds uint64");
  if (args.ferromagnetic)
    report.field("tl_convention", "TL convention", "tl_energy=+sum e_i=E-(N-1); auxiliary XXZ reference is unchanged")
        .field("defect_meaning", "Defect meaning", "M=(N-ell)/2; TL singlet defects, not physical spin flips")
        .field("spectral_gap", "Exact gap above ground space", model::ferromagnetic::spectral_gap<Real>(args.sites));
  return report;
}

inline int finish(bool converged)
{
  if (!converged) std::cerr << "Biquadratic solve incomplete; consider a larger budget or higher precision.\n";
  return converged ? 0 : 2;
}

} // namespace bethe::apps::biquadratic
