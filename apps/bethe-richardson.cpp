// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "report-common.hpp"
#include <bethe/richardson.hpp>

namespace
{
namespace model = bethe::richardson;
namespace cli = bethe::cli;
struct Arguments
{
    std::optional<std::string> levels, g, tolerance;
    std::optional<std::size_t> pairs;
    std::string blocked, precision = "fp64", format = "auto";
    std::size_t max_iterations = 10000, max_stages = 10000;
    bool variables = false;
};
auto program_info()
{
  auto info =
      bethe::cli::program_info("bethe-richardson", "Reduced BCS ground state in a specified blocked-level sector.",
                               bethe::citations::Tool::richardson);
  info.notes = {"H=sum_i epsilon_i*(n_up+n_down)-g*sum_ij b_i^dagger*b_j, including i=j.",
                "Distinct ascending single-particle energies; each level is a time-reversed doublet.",
                "Default Newton budget: 10000; zero budget returns the zero-coupling seed.",
                "Incomplete solves report energy at the REACHED g, not the requested g.",
                "Repeated levels, higher degeneracies, repulsive g<0 and excitations are not implemented.",
                "Pair rapidities are not reconstructed; no lattice momentum or PBC/OBC applies.",
                "See docs/richardson.md for energy shifts, blocking and continuation controls.",
                "Use --references for literature and applicability; see CITATIONS.md."};
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  bethe::cli::option(app, "--levels", args.levels, "required comma-separated single-particle energies")->required();
  bethe::cli::option(app, "--pairs", args.pairs, "required pair count in unblocked levels")->required();
  bethe::cli::option(app, "--g", args.g, "required finite g>=0 (attractive pairing)")->required();
  bethe::cli::option(app, "--blocked", args.blocked, "singly occupied levels, zero-based (default: none)");
  bethe::cli::option(app, "--variables", args.variables, "print regularized eigenvalue variables (not occupations)");
  bethe::cli::option(app, "--max-stages", args.max_stages, "attempted continuation stages (default: 10000)")
      ->capture_default_str();
  bethe::cli::option(app, "--tolerance", args.tolerance, "polynomial backward residual (default: 32 epsilon)");
  bethe::cli::option(app, "--max-iterations", args.max_iterations, "attempted Newton corrections, including retries")
      ->capture_default_str();
  bethe::cli::precision_option(app, args.precision);
  app.add_option("--format", args.format, "Stdout layout")
      ->check(CLI::IsMember({"auto", "pretty", "plain"}))
      ->capture_default_str();
}

void validate(Arguments const& args)
{
  if (!args.levels || !args.pairs || !args.g) throw std::invalid_argument("--levels, --pairs and --g are required");
  if (args.format != "auto" && args.format != "pretty" && args.format != "plain")
    throw std::invalid_argument("unknown output format: " + std::string(args.format));
}
template <typename Function> void each_item(std::string_view list, Function function)
{
  if (list.empty()) return;
  for (;;)
  {
    auto const comma = list.find(',');
    if (list.substr(0, comma).empty()) throw std::invalid_argument("empty list item");
    function(list.substr(0, comma));
    if (comma == std::string_view::npos) return;
    list.remove_prefix(comma + 1);
    if (list.empty()) throw std::invalid_argument("empty list item");
  }
}
char const* status(model::SolveStatus value)
{
  switch (value)
  {
    case model::SolveStatus::converged:
      return "converged";
    case model::SolveStatus::iteration_limit:
      return "Newton budget exhausted; target coupling not reached";
    case model::SolveStatus::stage_limit:
      return "stage budget exhausted; target coupling not reached";
    case model::SolveStatus::stalled:
      return "continuation stalled; target coupling not reached";
    case model::SolveStatus::ill_conditioned:
      return "linear system unresolved; target coupling not reached";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& args)
{
  std::vector<Real> levels;
  each_item(*args.levels, [&](auto value) { levels.push_back(uni20::parse_real<Real>(value)); });
  std::vector<std::size_t> blocked;
  each_item(args.blocked, [&](auto value) { blocked.push_back(cli::parse_size(value)); });
  Real const g = uni20::parse_real<Real>(*args.g);
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  options.max_stages = args.max_stages;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  cli::CpuTimer const timer;
  auto const state = model::ground_state<Real>(levels, *args.pairs, g, blocked, options);
  auto const cpu_time = timer.elapsed_text();
  cli::report_builder report("Richardson reduced BCS pairing");
  report.status(state.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning, status(state.status))
      .field("Hamiltonian", "sum epsilon_i*n_i-g*sum_ij b_i^dagger*b_j (including i=j)")
      .field("Calculation", "lowest state in the specified pair/blocked sector")
      .field("Levels", levels.size())
      .field("Unblocked levels", state.active.size())
      .field("Pairs", state.pairs)
      .field("Blocked levels", state.blocked.size())
      .field("Fermions", 2 * state.pairs + state.blocked.size())
      .field("Precision", args.precision)
      .field("Requested coupling g", uni20::format_real(state.coupling))
      .field("Reached coupling g", uni20::format_real(state.reached_coupling))
      .field("Energy evaluated at g", uni20::format_real(state.reached_coupling))
      .field("Status", status(state.status))
      .field("Total energy", uni20::format_real(state.energy))
      .field("Residual tolerance", uni20::format_real(options.residual_tolerance))
      .field("Reached backward residual", uni20::format_real(state.residual_norm))
      .field("Target backward residual", uni20::format_real(state.target_residual_norm))
      .field("Pair-number error", uni20::format_real(state.particle_number_error))
      .field("Newton corrections", state.iterations)
      .field("Continuation stages", state.stages)
      .field("Rejected stages", state.rejected_stages)
      .field("CPU time", cpu_time);
  if (args.variables)
  {
    auto& table = report.table("Level data and eigenvalue variables (not occupations)");
    table.header_separator()
        .column("Index")
        .column("epsilon", cli::table_alignment::decimal)
        .column("Blocked")
        .column("y", cli::table_alignment::decimal);
    std::size_t active = 0;
    for (std::size_t i = 0; i < levels.size(); ++i)
    {
      bool const is_active = active < state.active.size() && state.active[active] == i;
      table.row(i, uni20::format_real(levels[i]), is_active ? "no" : "yes",
                is_active ? uni20::format_real(state.eigenvalue_variables[active++]) : "-");
    }
  }
  cli::print_report(report, args.format);
  if (!state.converged)
    std::cerr
        << "Richardson solve incomplete: the reported energy is at the reached coupling, not the requested one.\n";
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
