// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "report-common.hpp"
#include <bethe/su_fermions.hpp>

namespace
{
namespace model = bethe::su_fermions;
namespace cli = bethe::cli;
struct Arguments
{
    std::optional<std::string> populations, length, c, tolerance;
    std::string precision = "fp64", format = "auto";
    std::size_t max_iterations = 10000, max_stages = 10000;
    bool roots = false;
};
auto program_info()
{
  auto info = bethe::cli::program_info("bethe-sun-fermions-pbc",
                                       "Equal-mass SU(n) continuum fermions on a ring; n is the number of components.",
                                       bethe::citations::Tool::sun_fermions_pbc);
  info.notes = {"H=-sum_j d_j^2+2c sum_(i<j) delta(x_i-x_j), hbar^2/(2m)=1; C>=0.",
                "Interacting ground states require EVERY occupied population odd.",
                "At c=0 or with at most one occupied component, arbitrary counts (including vacuum) are supported.",
                "Component order is preserved in the report; nesting uses descending nonzero populations.",
                "Zero budget yields no interacting energy. Incomplete solves retain only the last converged coupling.",
                "Even free seas choose the positive-current degenerate representative. E=sum k_j^2.",
                "Other periodic shell branches, attraction, excitations, open ends and Bose-Fermi mixtures",
                "are not implemented. See docs/su-fermions.md for conventions and diagnostics.",
                "Use --references for literature and applicability; see CITATIONS.md."};
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  bethe::cli::option(app, "--populations", args.populations, "required comma-separated nonnegative component counts")
      ->required();
  bethe::cli::option(app, "--length", args.length, "required finite positive physical circumference")->required();
  bethe::cli::option(app, "--c", args.c, "required finite repulsion, in inverse length")->required();
  bethe::cli::option(app, "--roots", args.roots, "print all nested rapidities/labels, or free modes");
  bethe::cli::option(app, "--max-stages", args.max_stages, "continuation stages (default: 10000)")
      ->capture_default_str();
  bethe::cli::option(app, "--tolerance", args.tolerance, "scaled equation residual (default: 32 epsilon)");
  bethe::cli::option(app, "--max-iterations", args.max_iterations, "attempted Newton corrections (default: 10000)")
      ->capture_default_str();
  bethe::cli::precision_option(app, args.precision);
  app.add_option("--format", args.format, "Stdout layout")
      ->check(CLI::IsMember({"auto", "pretty", "plain"}))
      ->capture_default_str();
}

void validate(Arguments const& args)
{
  if (!args.populations || !args.length || !args.c)
    throw std::invalid_argument("--populations, --length and --c are required");
  if (args.format != "auto" && args.format != "pretty" && args.format != "plain")
    throw std::invalid_argument("unknown output format: " + std::string(args.format));
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
      return "line search stalled; target coupling not reached";
    case model::SolveStatus::ill_conditioned:
      return "linear system unresolved; target coupling not reached";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& args)
{
  std::vector<std::size_t> pop;
  std::string_view list = *args.populations;
  if (list.empty()) throw std::invalid_argument("at least one component population is required");
  for (;;)
  {
    auto const comma = list.find(',');
    if (list.substr(0, comma).empty()) throw std::invalid_argument("empty population list item");
    pop.push_back(cli::parse_size(list.substr(0, comma)));
    if (comma == std::string_view::npos) break;
    list.remove_prefix(comma + 1);
    if (list.empty()) throw std::invalid_argument("empty population list item");
  }
  Real const length = uni20::parse_real<Real>(*args.length), c = uni20::parse_real<Real>(*args.c);
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  options.max_stages = args.max_stages;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  cli::CpuTimer const timer;
  auto const state = model::ground_state<Real>(pop, length, c, options);
  auto const cpu_time = timer.elapsed_text();
  cli::report_builder report("SU(n) fermion gas (periodic)");
  report.status(state.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning, status(state.status))
      .field("Calculation", state.free ? "exact free-fermion ground state" : "odd-population sector ground state")
      .field("Particles", state.particles)
      .field("Components", pop.size())
      .field("Occupied components", state.component_order.size())
      .field("Length", uni20::format_real(length))
      .field("Requested c", uni20::format_real(c))
      .field("Units", "hbar^2/(2m)=1; interaction 2c delta")
      .field("Precision", args.precision)
      .field("Status", status(state.status));
  if (state.energy)
  {
    report.field("Energy evaluated at c", uni20::format_real(*state.reached_interaction))
        .field("Total energy", uni20::format_real(*state.energy))
        .field("Momentum index", state.momentum_index)
        .field("Momentum P", uni20::format_real(state.momentum));
    if (!state.free)
      report.field("Reached residual", uni20::format_real(state.residual_norm))
          .field("Target residual", uni20::format_real(state.target_residual_norm));
  }
  else
    report.field("Energy", "unavailable: no converged coupling stage");
  report.field("Residual tolerance", uni20::format_real(options.residual_tolerance))
      .field("Newton corrections", state.iterations)
      .field("Continuation stages", state.stages)
      .field("CPU time", cpu_time);
  auto& components = report.table("Physical components (original input order)");
  components.header_separator().column("Component").column("Particles").column("Nesting rank");
  for (std::size_t a = 0; a < pop.size(); ++a)
  {
    auto const it = std::find(state.component_order.begin(), state.component_order.end(), a);
    components.row(a, pop[a],
                   it == state.component_order.end() ? "empty" : std::to_string(it - state.component_order.begin()));
  }
  if (args.roots && state.free)
  {
    Real const pi = Real{4} * std::atan(Real{1});
    for (std::size_t a = 0; a < pop.size(); ++a)
    {
      auto& table = report.table("Free modes: component " + std::to_string(a));
      table.header_separator().column("Mode").column("k", cli::table_alignment::decimal);
      for (auto mode : state.free_modes[a])
        table.row(mode, uni20::format_real(Real{2} * pi * Real(mode) / length));
    }
  }
  else if (args.roots)
    for (std::size_t a = 0; a < state.rapidities.size(); ++a)
    {
      auto& table = report.table(a ? "Spin rapidities: level " + std::to_string(a) : "Charge momenta: level 0");
      table.header_separator()
          .column("Index")
          .column(a ? "J" : "I")
          .column(a ? "lambda" : "k", cli::table_alignment::decimal);
      for (std::size_t j = 0; j < state.rapidities[a].size(); ++j)
        table.row(j, uni20::to_string_fraction(state.quantum_numbers[a][j]),
                  uni20::format_real(state.rapidities[a][j]));
    }
  cli::print_report(report, args.format);
  if (!state.converged) std::cerr << "SU fermion solve incomplete: no energy at the requested coupling is claimed.\n";
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
