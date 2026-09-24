// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/su3.hpp>

namespace
{
namespace model = bethe::su3;
namespace cli = bethe::cli;
struct Arguments
{
    std::size_t sites = 0, max_iterations = 10000;
    std::optional<std::string> tolerance;
    std::string precision = "fp64";
    cli::DataOutputOptions output;
    bool roots = false;
};

auto program_info()
{
  auto info = bethe::cli::program_info("bethe-su3-pbc",
                                       "Fundamental SU(3) permutation chain, periodic: H=sum_j P_(j,j+1), J=1.",
                                       bethe::citations::Tool::su3_pbc);
  info.notes = {"Balanced singlet ground state only: L>=3 divisible by 3.",
                "Each color has L/3 sites; nested real-root counts are M1=2L/3 and M2=L/3.",
                "E=L-sum_j 1/(lambda_j^2+1/4). Momentum is 0 in the supported singlet.",
                "Roots use conventional lambda,mu, not the XXX front end's z=2*lambda.",
                "Other color sectors, excitations, complex strings, twists and open ends",
                "are not implemented. See docs/su3.md, including the spin-1 ULS mapping.",
                "Use --references for literature and applicability; see CITATIONS.md."};
  info.notes.push_back(
      "Tables: states; first_roots and second_roots with --roots. See docs/output.md for file exports.");
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  bethe::cli::option(app, "L", args.sites, "Number of sites or rungs")->required();
  bethe::cli::option(app, "--roots", args.roots, "print both root families and exact I,J labels");
  bethe::cli::option(app, "--tolerance", args.tolerance,
                     "max equation residual divided by L (default: 32 epsilon; not an energy-error bound)");
  bethe::cli::option(app, "--max-iterations", args.max_iterations, "accepted Newton updates (default: 10000)")
      ->capture_default_str();
  bethe::cli::precision_option(app, args.precision);
  cli::add_data_output_options(app, args.output, true);
}

void validate(Arguments const& args) { args.output.validate(); }

char const* status(model::SolveStatus value)
{
  switch (value)
  {
    case model::SolveStatus::converged:
      return "converged";
    case model::SolveStatus::iteration_limit:
      return "iteration limit reached; unconverged estimate";
    case model::SolveStatus::stalled:
      return "line search or representable precision stalled; unconverged estimate";
  }
  return "unknown";
}

template <uni20::Real Real> int run(Arguments const& args, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  auto computation = context.computation();
  auto const state = model::ground_state<Real>(args.sites, options);
  computation.finish();
  cli::RunReport report(context, "SU(3) permutation chain (periodic)");
  report.status(state.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning, status(state.status))
      .field("hamiltonian", "Hamiltonian", "H=sum_j P_(j,j+1); J=1")
      .field("calculation", "Calculation", "balanced singlet ground state")
      .field("sites", "Sites", state.sites)
      .field("color_populations", "Color populations",
             fmt::format("{}, {}, {}", state.populations[0], state.populations[1], state.populations[2]))
      .field("first_level_roots", "First-level roots", state.rapidities[0].size())
      .field("second_level_roots", "Second-level roots", state.rapidities[1].size())
      .field("precision", "Precision", args.precision)
      .field("residual_tolerance", "Residual tolerance", options.residual_tolerance)
      .result(state.converged, status(state.status))
      .field("total_energy", "Total energy", state.energy)
      .field("energy_per_site", "Energy per site", state.energy / Real(state.sites))
      .field("momentum_index", "Momentum index", state.momentum_index)
      .field("momentum_p", "Momentum P", state.momentum)
      .field("first_level_residual", "First-level residual", state.level_residuals[0])
      .field("second_level_residual", "Second-level residual", state.level_residuals[1])
      .field("residual_norm", "Residual norm", state.residual_norm)
      .field("iterations", "Iterations", state.iterations);
  using cli::column;
  std::vector<std::string> names{"states"};
  if (args.roots)
  {
    names.push_back("first_roots");
    names.push_back("second_roots");
  }
  cli::ResultOutput output(report, args.output, names);
  output.table(
      "states", "State",
      [&](auto& t) {
        t.append(0, state.energy, state.momentum_index, state.momentum, state.level_residuals[0],
                 state.level_residuals[1], state.residual_norm, state.iterations, state.converged,
                 status(state.status));
      },
      column<std::size_t>("state_id"), column<Real>("energy"), column<std::size_t>("momentum_index"), column<Real>("p"),
      column<Real>("first_residual"), column<Real>("second_residual"), column<Real>("residual"),
      column<std::size_t>("iterations"), column<bool>("converged"), column<std::string>("status"));
  if (args.roots)
    for (std::size_t a = 0; a < 2; ++a)
      output.table(
          a ? "second_roots" : "first_roots", a ? "Second-level rapidities" : "First-level rapidities",
          [&](auto& t) {
            for (std::size_t j = 0; j < state.rapidities[a].size(); ++j)
              t.append(0, j, state.quantum_numbers[a][j], state.rapidities[a][j]);
          },
          column<std::size_t>("state_id"), column<std::size_t>("index"),
          column<uni20::half_int>("quantum_number", a ? "J" : "I"), column<Real>("rapidity", a ? "mu" : "lambda"));
  output.finish();
  if (!state.converged) std::cerr << "SU(3) solve incomplete; consider a larger budget or higher precision.\n";
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
