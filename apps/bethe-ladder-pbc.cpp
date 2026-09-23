// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "report-common.hpp"
#include <bethe/ladder.hpp>

namespace
{
namespace model = bethe::ladder;
namespace cli = bethe::cli;
struct Arguments
{
    std::size_t rungs = 0, max_iterations = 10000, max_branches = 10000;
    std::optional<std::size_t> singlets;
    std::optional<std::string> rung, tolerance;
    std::string precision = "fp64", format = "auto";
    bool sectors = false, roots = false;
};
auto program_info()
{
  auto info = bethe::cli::program_info("bethe-ladder-pbc",
                                       "Wang's integrable spin-1/2 ladder, L>=2 periodic rungs, zero field.",
                                       bethe::citations::Tool::ladder_pbc);
  info.notes = {"H=sum[S.S_next+T.T_next+4(S.S_next)(T.T_next)]+JR*sum S.T.",
                "This is NOT the ordinary two-leg Heisenberg ladder. JR may have either sign.",
                "Sector minima include compatible SU(4) descendants; triplet populations",
                "are minimized too, not fixed Sz. --sectors --roots prints the best state.",
                "No excitations, fields, arbitrary four-spin couplings or open ends.",
                "An incomplete scan reports only a candidate upper bound, never a minimum.",
                "See docs/ladder.md for normalization, finite-ring labels and scan cost.",
                "Use --references for literature and applicability; see CITATIONS.md."};
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  bethe::cli::option(app, "L", args.rungs, "Number of sites or rungs")->required();
  bethe::cli::option(app, "--rung", args.rung, "Rung coupling, either sign")->required();
  bethe::cli::option(app, "--singlets", args.singlets, "fixed singlet-count sector minimum");
  bethe::cli::option(app, "--sectors", args.sectors,
                     "all singlet-count sector minima (mutually exclusive; default: global minimum)");
  bethe::cli::option(app, "--max-branches", args.max_branches, "total highest-weight sea branches (10000)")
      ->capture_default_str();
  bethe::cli::option(app, "--roots", args.roots, "highest-weight roots/labels for selected state");
  bethe::cli::option(app, "--tolerance", args.tolerance,
                     "max Bethe residual divided by L (default: 32 epsilon; not an energy-error bound)");
  bethe::cli::option(app, "--max-iterations", args.max_iterations, "total attempted Newton corrections (10000)")
      ->capture_default_str();
  bethe::cli::precision_option(app, args.precision);
  app.add_option("--format", args.format, "Stdout layout")
      ->check(CLI::IsMember({"auto", "pretty", "plain"}))
      ->capture_default_str();
}

void validate(Arguments const& args)
{
  if (!args.rung) throw std::invalid_argument("--rung is required");
  if (args.sectors && args.singlets) throw std::invalid_argument("--sectors and --singlets are mutually exclusive");
  if (args.format != "auto" && args.format != "plain" && args.format != "pretty")
    throw std::invalid_argument("unknown output format: " + std::string(args.format));
}
char const* status(model::SolveStatus value)
{
  switch (value)
  {
    case model::SolveStatus::converged:
      return "converged";
    case model::SolveStatus::iteration_limit:
      return "Newton budget exhausted; incomplete scan";
    case model::SolveStatus::branch_limit:
      return "branch budget exhausted; incomplete scan";
    case model::SolveStatus::stalled:
      return "line search stalled; incomplete scan";
    case model::SolveStatus::ill_conditioned:
      return "ill-conditioned Newton system; incomplete scan";
  }
  return "unknown";
}
std::string shape_text(model::detail::Shape const& shape)
{
  return fmt::format("{}, {}, {}, {}", shape[0], shape[1], shape[2], shape[3]);
}

template <uni20::Real Real> int run(Arguments const& args)
{
  Real const rung = uni20::parse_real<Real>(*args.rung);
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  options.max_branches = args.max_branches;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  cli::CpuTimer const timer;
  model::State<Real> state;
  std::optional<model::SectorScan<Real>> scan;
  if (args.sectors)
  {
    scan = model::sector_ground_states(args.rungs, rung, options);
    state = scan->sectors.front();
    for (auto const& s : scan->sectors)
      if (s.energy && (!state.energy || *s.energy < *state.energy)) state = s;
  }
  else if (args.singlets)
    state = model::sector_ground_state(args.rungs, *args.singlets, rung, options);
  else
    state = model::ground_state(args.rungs, rung, options);
  auto const cpu_time = timer.elapsed_text();
  cli::report_builder report("Integrable spin ladder (periodic)");
  report.status(state.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning, status(state.status))
      .field("Hamiltonian", "H=sum[S.S_next+T.T_next+4(S.S_next)(T.T_next)]+J_r*sum S.T")
      .field("Calculation", args.singlets  ? "fixed singlet-count sector"
                            : args.sectors ? "all singlet-count sectors"
                                           : "global ground state")
      .field("Rungs", args.rungs)
      .field("Rung coupling J_r", uni20::format_real(rung))
      .field("Precision", args.precision)
      .field("Residual tolerance", uni20::format_real(options.residual_tolerance))
      .field("Status", status(state.status));
  if (state.energy)
  {
    report.field(state.converged ? "Total energy" : "Candidate energy (upper bound)", uni20::format_real(*state.energy))
        .field("Singlets N_s", state.singlets)
        .field("Populations (s,t+,t0,t-)", shape_text(state.populations))
        .field("Highest-weight rows", shape_text(state.highest_weight.shape))
        .field("SU(4) descendant", state.descendant ? "yes; roots belong to the highest-weight representative"
                                                    : "no (up to color permutation)")
        .field("Permutation energy", uni20::format_real(state.highest_weight.energy))
        .field("Momentum index", state.highest_weight.momentum_index)
        .field("Reflected momentum index", (state.rungs - state.highest_weight.momentum_index) % state.rungs)
        .field("Selected branch residual", uni20::format_real(state.highest_weight.residual));
    if (state.converged) report.field("Energy per rung", uni20::format_real(*state.energy / Real(state.rungs)));
    if (state.analytic) report.field("Exact limit", "rung-singlet product (no root solve)");
  }
  else
    report.field("Energy", "unavailable: no converged branch");
  report.field("Highest weights visited", state.tableaux)
      .field("Sea branches attempted", state.branches)
      .field("Newton corrections", state.iterations)
      .field("CPU time", cpu_time);
  if (scan)
  {
    auto& table = report.table(scan->complete ? "Singlet-sector minima"
                                              : "Singlet-sector candidate upper bounds (scan incomplete)");
    table.header_separator()
        .column("N_s")
        .column("Energy", cli::table_alignment::decimal)
        .column("Highest-weight rows")
        .column("Descendant")
        .column("Momentum index");
    for (auto const& s : scan->sectors)
      table.row(s.singlets, s.energy ? uni20::format_real(*s.energy) : "unavailable",
                s.energy ? shape_text(s.highest_weight.shape) : "-", s.energy ? (s.descendant ? "yes" : "no") : "-",
                s.energy ? std::to_string(s.highest_weight.momentum_index) : "-");
  }
  if (args.roots && state.energy)
    for (std::size_t a = 0; a < state.highest_weight.rapidities.size(); ++a)
    {
      auto& table = report.table(fmt::format("Highest-weight rapidities: level {}", a + 1));
      table.header_separator().column("Index").column("Bethe label").column("Rapidity", cli::table_alignment::decimal);
      for (std::size_t j = 0; j < state.highest_weight.rapidities[a].size(); ++j)
        table.row(j, uni20::to_string_fraction(state.highest_weight.labels[a][j]),
                  uni20::format_real(state.highest_weight.rapidities[a][j]));
    }
  cli::print_report(report, args.format);
  if (!state.converged) std::cerr << "Ladder scan incomplete: candidate energies are not certified sector minima.\n";
  return state.converged ? 0 : 2;
}
} // namespace
int main(int argc, char** argv)
{
  Arguments args;
  return bethe::cli::program_main(
      argc, argv, program_info(), [&](auto& app) { add_options(app, args); },
      [&](auto&) {
        validate(args);
        return bethe::cli::dispatch_precision(args.precision, [&]<uni20::Real Real> { return run<Real>(args); });
      });
}
