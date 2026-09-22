// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "citation-report.hpp"
#include "report-common.hpp"
#include <bethe/richardson.hpp>

namespace
{
namespace model = bethe::richardson;
namespace cli = bethe::cli;
struct Arguments
{
    std::optional<std::string_view> levels, g, tolerance;
    std::optional<std::size_t> pairs;
    std::string_view blocked, precision = "fp64", format = "auto";
    std::size_t max_iterations = 10000, max_stages = 10000;
    bool variables = false;
};
void usage(std::ostream& out)
{
  out << "Usage: bethe-richardson --levels E0,E1,... --pairs M --g G [options]\n"
      << "Reduced BCS ground state in a specified blocked-level sector.\n"
      << "H=sum_i epsilon_i*(n_up+n_down)-g*sum_ij b_i^dagger*b_j, including i=j.\n"
      << "Distinct ascending single-particle energies; each level is a time-reversed doublet.\n"
      << "  --levels LIST                      required comma-separated single-particle energies\n"
      << "  --pairs COUNT                      required pair count in unblocked levels\n"
      << "  --g VALUE                          required finite g>=0 (attractive pairing)\n"
      << "  --blocked INDICES                  singly occupied levels, zero-based (default: none)\n"
      << "  --variables                        print regularized eigenvalue variables (not occupations)\n"
      << "  --precision fp64|long-double|fp128  (default: fp64; fp128 requires MPLAPACK)\n"
      << "  --tolerance VALUE                  polynomial backward residual (default: 32 epsilon)\n"
      << "  --max-iterations COUNT             attempted Newton corrections, including retries\n"
      << "  --max-stages COUNT                 attempted continuation stages (default: 10000)\n"
      << "  --format auto|pretty|plain         (default: auto)\n"
      << "  --help                             show this help and references\n"
      << "Default Newton budget: 10000; zero budget returns the zero-coupling seed.\n"
      << "Incomplete solves report energy at the REACHED g, not the requested g.\n"
      << "Repeated levels, higher degeneracies, repulsive g<0 and excitations are not implemented.\n"
      << "Pair rapidities are not reconstructed; no lattice momentum or PBC/OBC applies.\n"
      << "See docs/richardson.md for energy shifts, blocking and continuation controls.\n";
  cli::print_citations(out, bethe::citations::Tool::richardson);
}
Arguments parse(int argc, char** argv)
{
  Arguments out;
  for (int i = 1; i < argc; ++i)
  {
    std::string_view const option = argv[i];
    if (option == "--variables")
    {
      out.variables = true;
      continue;
    }
    if (option != "--levels" && option != "--pairs" && option != "--g" && option != "--blocked" &&
        option != "--precision" && option != "--format" && option != "--tolerance" && option != "--max-iterations" &&
        option != "--max-stages")
      throw std::invalid_argument("unknown option: " + std::string(option));
    if (++i == argc) throw std::invalid_argument("missing value for " + std::string(option));
    std::string_view const value = argv[i];
    if (option == "--levels")
      out.levels = value;
    else if (option == "--pairs")
      out.pairs = cli::parse_size(value);
    else if (option == "--g")
      out.g = value;
    else if (option == "--blocked")
      out.blocked = value;
    else if (option == "--precision")
      out.precision = value;
    else if (option == "--format")
      out.format = value;
    else if (option == "--tolerance")
      out.tolerance = value;
    else if (option == "--max-iterations")
      out.max_iterations = cli::parse_size(value);
    else
      out.max_stages = cli::parse_size(value);
  }
  if (!out.levels || !out.pairs || !out.g) throw std::invalid_argument("--levels, --pairs and --g are required");
  if (out.format != "auto" && out.format != "pretty" && out.format != "plain")
    throw std::invalid_argument("unknown output format: " + std::string(out.format));
  return out;
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
  if (argc == 2 && std::string_view(argv[1]) == "--help")
  {
    usage(std::cout);
    return 0;
  }
  if (argc < 2)
  {
    usage(std::cerr);
    return 1;
  }
  try
  {
    auto const args = parse(argc, argv);
    return cli::dispatch_precision(args.precision, [&]<uni20::Real Real> { return run<Real>(args); });
  }
  catch (std::exception const& error)
  {
    std::cerr << "bethe-richardson: " << error.what() << '\n';
    return 1;
  }
}
