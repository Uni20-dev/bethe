// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "report.hpp"
#include <bethe/biquadratic_qsystem.hpp>

namespace bethe::apps::biquadratic
{
inline char const* status(bethe::xxz::quantum_group::qsystem::Status value)
{
  using S = bethe::xxz::quantum_group::qsystem::Status;
  switch (value)
  {
    case S::converged:
      return "converged; numerical admissibility checks passed";
    case S::iteration_limit:
      return "iteration limit reached; unconverged estimate";
    case S::singular_jacobian:
      return "singular Jacobian; unconverged estimate";
    case S::stalled:
      return "line search stalled; unconverged estimate";
    case S::unresolved_roots:
      return "unresolved polynomial roots; unverified estimate";
    case S::inadmissible:
      return "root admissibility or original Bethe-equation check failed";
  }
  return "unknown";
}

template <uni20::Real Real> std::vector<Real> parse_q_seed(std::string_view text)
{
  if (text == "none") return {};
  std::vector<Real> result;
  for (;;)
  {
    auto const comma = text.find(',');
    result.push_back(uni20::parse_real<Real>(std::string(text.substr(0, comma))));
    if (comma == std::string_view::npos) return result;
    text.remove_prefix(comma + 1);
  }
}

template <uni20::Real Real>
int run_qsystem(Arguments const& args, bethe::SolverOptions<Real> const& options, cli::RunReport report)
{
  auto& context = report.context();
  if (args.sites > 8)
    std::cerr << "Warning: Q-system spectrum validation covers N<=8. Larger chains are experimental; "
                 "memory, search cost and string conditioning may prevent convergence or completeness.\n";
  std::vector<model::qsystem::State<Real>> states;
  auto const ell = args.through_lines.value_or(args.ferromagnetic ? args.sites - 2 : 0);
  bool complete = false;
  if (args.q_seed)
  {
    states.push_back(context.measure([&] {
      auto const seed = parse_q_seed<Real>(*args.q_seed);
      return args.ferromagnetic ? model::ferromagnetic::qsystem::solve<Real>(args.sites, seed, options)
                                : model::qsystem::solve<Real>(args.sites, seed, options);
    }));
    complete = states.front().reference.converged;
    report.field("calculation", "Calculation", "selected Q-system branch; not necessarily a lowest level");
  }
  else
  {
    auto scan = context.measure([&] {
      bethe::xxz::quantum_group::qsystem::SearchOptions search{.max_attempts = args.max_attempts.value_or(4000)};
      return args.ferromagnetic ? model::ferromagnetic::qsystem::spectrum<Real>(args.sites, ell, search, options)
                                : model::qsystem::spectrum<Real>(args.sites, ell, search, options);
    });
    complete = scan.complete();
    report.field("calculation", "Calculation", "Q-system spectrum search (real and complex roots)")
        .field("expected_module_dimension", "Expected module dimension", scan.expected_count,
               {.missing = "overflow (>size_t); completeness unavailable"})
        .field("discovered_levels", "Discovered levels", scan.states.size())
        .field("attempts", "Attempts", scan.attempts)
        .field("attempt_budget", "Attempt budget", args.max_attempts.value_or(4000))
        .field("failed_attempts", "Failed attempts", scan.failed_attempts)
        .field("ordering", "Ordering",
               complete ? "numerically complete TL module; count matched, not a rigorous certificate"
                        : "incomplete discoveries; NOT guaranteed lowest levels");
    states = std::move(scan.states);
  }
  auto const ground = context.measure([&] {
    return args.ferromagnetic ? model::ferromagnetic::solve_real<Real>(args.sites, {}, options)
                              : model::ground_state<Real>(args.sites, options);
  });
  complete = complete && ground.reference.converged;
  report.field("tl_through_lines", "TL through-lines", args.q_seed ? states.front().through_lines : ell)
      .field("tl_defects", "TL defects M", (args.sites - (args.q_seed ? states.front().through_lines : ell)) / 2)
      .field("q_system_validation", "Q-system validation",
             args.sites <= 8 ? "within small-chain regression range"
                             : "experimental beyond N=8; no completeness or convergence guarantee")
      .field("residual_convention", "Residual convention",
             "Q-system Wronskian coefficient backward error; not an energy-error bound")
      .field("root_coordinate", "Root coordinate", "x=cosh(2u)=cos(alpha); Bajnok u, alpha=-2iu")
      .field("multiplicity_meaning", "Multiplicity meaning", "physical states per TL eigenvector, not SU(2) multiplets")
      .field("gap_reference", "Gap reference",
             args.ferromagnetic           ? "E-(N-1); exact degenerate ferro ground space"
             : ground.reference.converged ? "E-E0; global singlet ground state"
                                          : "unavailable; ground solve failed")
      .result(complete, complete ? "converged" : "incomplete or unverified");
  std::vector<std::string> names{"states", "reference", "q_coefficients"};
  if (args.roots) names.push_back("roots");
  cli::ResultOutput output(report, args.output, table_names(args, names));
  output.table(
      "states", "Q-system levels",
      [&](auto& table) {
        for (std::size_t i = 0; i < states.size(); ++i)
        {
          auto const& s = states[i];
          auto const& r = s.reference;
          std::optional<Real> gap;
          if (r.converged && ground.reference.converged)
            gap =
                args.ferromagnetic ? s.tl_energy : std::optional<Real>{Real{2} * (*r.energy - ground.reference.energy)};
          table.append(i, s.through_lines, s.multiplicity, s.energy, gap, s.tl_energy, r.energy, r.residual_norm,
                       uni20::isfinite(r.bethe_residual) ? std::optional<Real>{r.bethe_residual} : std::nullopt,
                       uni20::isfinite(r.bethe_residual_bound) ? std::optional<Real>{r.bethe_residual_bound}
                                                               : std::nullopt,
                       r.iterations, r.converged, std::string(status(r.status)));
        }
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("through_lines", "Through-lines"),
      cli::column<std::optional<std::uint64_t>>("multiplicity", "Multiplicity"),
      cli::column<std::optional<Real>>("energy", "Energy"), cli::column<std::optional<Real>>("gap", "E-E0"),
      cli::column<std::optional<Real>>("tl_energy"), cli::column<std::optional<Real>>("reference_energy"),
      cli::column<Real>("residual", "Wronskian residual"),
      cli::column<std::optional<Real>>("bethe_residual", "Bethe residual"),
      cli::column<std::optional<Real>>("bethe_residual_bound", "Bethe uncertainty bound"),
      cli::column<std::size_t>("iterations"), cli::column<bool>("converged"),
      cli::column<std::string>("status", "Status"));
  output.table(
      "reference", "Global ground reference",
      [&](auto& table) {
        table.append(states.size(), ground.energy, ground.reference.residual_norm, ground.reference.converged,
                     std::string(status(ground.reference.status)));
      },
      cli::column<std::size_t>("state_id"), cli::column<Real>("energy", "Energy"), cli::column<Real>("residual"),
      cli::column<bool>("converged"), cli::column<std::string>("status", "Status"));
  output.table(
      "q_coefficients", "Monic Q(x) coefficients (ascending powers; leading 1 omitted)",
      [&](auto& table) {
        for (std::size_t i = 0; i < states.size(); ++i)
          for (std::size_t k = 0; k < states[i].reference.coefficients.size(); ++k)
            table.append(i, k, states[i].reference.coefficients[k]);
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("power"), cli::column<Real>("coefficient"));
  if (args.roots)
    output.table(
        "roots", "Complex Q-system roots (x=cosh(2u))",
        [&](auto& table) {
          for (std::size_t i = 0; i < states.size(); ++i)
            for (std::size_t k = 0; k < states[i].reference.roots.roots.size(); ++k)
            {
              auto const& r = states[i].reference;
              auto x = r.roots.roots[k];
              std::optional<Real> ur, ui;
              if (k < r.rapidities.size())
              {
                ur = r.rapidities[k].real();
                ui = r.rapidities[k].imag();
              }
              auto finite = [](Real value) {
                return uni20::isfinite(value) ? std::optional<Real>{value} : std::nullopt;
              };
              table.append(i, k, finite(x.real()), finite(x.imag()), ur ? finite(*ur) : std::nullopt,
                           ui ? finite(*ui) : std::nullopt);
            }
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("index"),
        cli::column<std::optional<Real>>("x_real"), cli::column<std::optional<Real>>("x_imag"),
        cli::column<std::optional<Real>>("u_real"), cli::column<std::optional<Real>>("u_imag"));
  write_spin_content(output, args, states.size() + 1,
                     [&](std::size_t i) { return i < states.size() ? states[i].through_lines : ground.through_lines; });
  output.finish();
  if (!complete) std::cerr << "Q-system calculation incomplete or unverified; no lowest-level guarantee.\n";
  return complete ? 0 : 2;
}

} // namespace bethe::apps::biquadratic
