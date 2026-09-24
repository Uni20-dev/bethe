// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/su_fermions.hpp>

namespace
{
namespace model = bethe::su_fermions;
namespace cli = bethe::cli;
struct Arguments
{
    std::optional<std::string> populations, length, c, tolerance;
    std::string precision = "fp64";
    cli::DataOutputOptions output;
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
  cli::add_data_output_options(app, args.output, true);
}

void validate(Arguments const& args)
{
  if (!args.populations || !args.length || !args.c)
    throw std::invalid_argument("--populations, --length and --c are required");
  args.output.validate();
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
template <uni20::Real Real> int run(Arguments const& args, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
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
  auto computation = context.computation();
  auto const state = model::ground_state<Real>(pop, length, c, options);
  computation.finish();
  cli::RunReport report(context, "SU(n) fermion gas (periodic)");
  report.status(state.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning, status(state.status))
      .field("calculation", "Calculation",
             state.free ? "exact free-fermion ground state" : "odd-population sector ground state")
      .field("particles", "Particles", state.particles)
      .field("components", "Components", pop.size())
      .field("occupied_components", "Occupied components", state.component_order.size())
      .field("length", "Length", length)
      .field("requested_c", "Requested c", c)
      .field("units", "Units", "hbar^2/(2m)=1; interaction 2c delta")
      .field("precision", "Precision", args.precision)
      .result(state.converged, status(state.status));
  if (state.energy)
  {
    report.field("energy_evaluated_at_c", "Energy evaluated at c", *state.reached_interaction)
        .field("total_energy", "Total energy", *state.energy)
        .field("momentum_index", "Momentum index", state.momentum_index)
        .field("momentum_p", "Momentum P", state.momentum);
    if (!state.free)
      report.field("reached_residual", "Reached residual", state.residual_norm)
          .field("target_residual", "Target residual", state.target_residual_norm);
  }
  else
    report.field("energy", "Energy", state.energy, {.missing = "unavailable: no converged coupling stage"});
  report.field("residual_tolerance", "Residual tolerance", options.residual_tolerance)
      .field("newton_corrections", "Newton corrections", state.iterations)
      .field("continuation_stages", "Continuation stages", state.stages);
  std::vector<std::string> tables{"states", "components"};
  if (args.roots) tables.push_back(state.free ? "free_modes" : "roots");
  cli::ResultOutput output(report, args.output, tables);
  output.table(
      "states", "State",
      [&](auto& table) {
        table.append(0, state.energy, c, state.reached_interaction,
                     state.energy ? std::optional{state.momentum_index} : std::nullopt,
                     state.energy ? std::optional{state.momentum} : std::nullopt,
                     state.energy ? std::optional{state.residual_norm} : std::nullopt,
                     state.energy ? std::optional{state.target_residual_norm} : std::nullopt, state.iterations,
                     state.stages, state.converged, std::string(status(state.status)));
      },
      cli::column<std::size_t>("state_id"), cli::column<std::optional<Real>>("energy", "Energy"),
      cli::column<Real>("requested_c"), cli::column<std::optional<Real>>("reached_c"),
      cli::column<std::optional<std::int64_t>>("momentum_index"), cli::column<std::optional<Real>>("p", "P"),
      cli::column<std::optional<Real>>("residual"), cli::column<std::optional<Real>>("target_residual"),
      cli::column<std::size_t>("iterations"), cli::column<std::size_t>("stages"), cli::column<bool>("converged"),
      cli::column<std::string>("status"));
  output.table(
      "components", "Physical components (original input order)",
      [&](auto& table) {
        for (std::size_t a = 0; a < pop.size(); ++a)
        {
          auto it = std::find(state.component_order.begin(), state.component_order.end(), a);
          table.append(0, a, pop[a],
                       it == state.component_order.end()
                           ? std::nullopt
                           : std::optional<std::size_t>(it - state.component_order.begin()));
        }
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("component", "Component"),
      cli::column<std::size_t>("particles", "Particles"),
      cli::column<std::optional<std::size_t>>("nesting_rank", "Nesting rank"));
  if (args.roots && state.free)
  {
    Real const pi = Real{4} * std::atan(Real{1});
    output.table(
        "free_modes", "Free modes (original component order)",
        [&](auto& table) {
          for (std::size_t a = 0; a < pop.size(); ++a)
            for (auto mode : state.free_modes[a])
              table.append(0, a, mode, Real{2} * pi * Real(mode) / length);
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("component", "Component"),
        cli::column<std::int64_t>("mode", "Mode"), cli::column<Real>("k"));
  }
  else if (args.roots)
    output.table(
        "roots", "Nested rapidities (level 0: charge; higher levels: spin)",
        [&](auto& table) {
          for (std::size_t a = 0; a < state.rapidities.size(); ++a)
            for (std::size_t j = 0; j < state.rapidities[a].size(); ++j)
              table.append(0, a, j, state.quantum_numbers[a][j], state.rapidities[a][j]);
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("level", "Level"),
        cli::column<std::size_t>("index", "Index"), cli::column<uni20::half_int>("quantum_number", "I / J"),
        cli::column<Real>("rapidity", "k / lambda"));
  output.finish();
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
        return bethe::cli::dispatch_precision(args.precision,
                                              [&]<uni20::Real Real> { return run<Real>(args, argc, argv); });
      });
}
