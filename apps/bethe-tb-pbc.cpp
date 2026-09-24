// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/takhtajan_babujian.hpp>

namespace
{
namespace model = bethe::takhtajan_babujian;
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
  auto info = bethe::cli::program_info(
      "bethe-tb-pbc", "Spin-1 Takhtajan-Babujian chain, periodic: H=sum_j [S_j.S_(j+1)-(S_j.S_(j+1))^2].",
      bethe::citations::Tool::tb_pbc);
  info.notes = {"Zero-field singlet ground state; even L>=4. Bilinear coefficient is 1.",
                "Complex two-string roots retain finite-size deviations, not ideal strings.",
                "lambda=x +/- i*(1/2+delta); E=-4 sum_j 1/(1+lambda_j^2); momentum 0.",
                "Odd lengths, other spins/sectors, excitations and open ends are not implemented.",
                "This is not the generic spin-1 Heisenberg chain or the SU(3) ULS point.",
                "See docs/takhtajan-babujian.md for normalization and numerical conventions.",
                "Use --references for literature and applicability; see CITATIONS.md."};
  info.notes.push_back(
      "Tables: states; strings and roots with --roots. Complex roots have separate real/imaginary columns.");
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  bethe::cli::option(app, "L", args.sites, "Number of sites or rungs")->required();
  bethe::cli::option(app, "--roots", args.roots, "print string centers, deviations and complex roots");
  bethe::cli::option(app, "--tolerance", args.tolerance,
                     "max phase/modulus residual divided by L (default: 32 epsilon; not an energy-error bound)");
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
    case model::SolveStatus::ill_conditioned:
      return "Newton system unresolved; unconverged estimate";
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
  cli::RunReport report(context, "Spin-1 Takhtajan-Babujian chain (periodic)");
  report.status(state.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning, status(state.status))
      .field("hamiltonian", "Hamiltonian", "H=sum_j [S.S-(S.S)^2]; bilinear coefficient 1")
      .field("calculation", "Calculation", "zero-field singlet ground state with finite string deviations")
      .field("sites", "Sites", state.sites)
      .field("spin", "Spin", 1)
      .field("total_spin", "Total spin", 0)
      .field("two_strings", "Two-strings", state.centers.size())
      .field("complex_roots", "Complex roots", state.rapidities.size())
      .field("precision", "Precision", args.precision)
      .field("residual_tolerance", "Residual tolerance", options.residual_tolerance)
      .result(state.converged, status(state.status))
      .field("total_energy", "Total energy", state.energy)
      .field("energy_per_site", "Energy per site", state.energy / Real(state.sites))
      .field("momentum_index", "Momentum index", state.momentum_index)
      .field("momentum_p", "Momentum P", state.momentum)
      .field("phase_residual", "Phase residual", state.phase_residual)
      .field("modulus_residual", "Modulus residual", state.modulus_residual)
      .field("residual_norm", "Residual norm", state.residual_norm)
      .field("iterations", "Iterations", state.iterations);
  using cli::column;
  std::vector<std::string> names{"states"};
  if (args.roots)
  {
    names.push_back("strings");
    names.push_back("roots");
  }
  cli::ResultOutput output(report, args.output, names);
  output.table(
      "states", "State",
      [&](auto& t) {
        t.append(0, state.energy, state.momentum_index, state.momentum, state.phase_residual, state.modulus_residual,
                 state.residual_norm, state.iterations, state.converged, status(state.status));
      },
      column<std::size_t>("state_id"), column<Real>("energy"), column<std::size_t>("momentum_index"), column<Real>("p"),
      column<Real>("phase_residual"), column<Real>("modulus_residual"), column<Real>("residual"),
      column<std::size_t>("iterations"), column<bool>("converged"), column<std::string>("status"));
  if (args.roots)
  {
    output.table(
        "strings", "Deviated two-strings",
        [&](auto& t) {
          for (std::size_t j = 0; j < state.centers.size(); ++j)
            t.append(0, j, state.string_quantum_numbers[j], state.centers[j], state.deviations[j]);
        },
        column<std::size_t>("state_id"), column<std::size_t>("index"), column<uni20::half_int>("i", "I"),
        column<Real>("center", "Center x"), column<Real>("deviation", "Deviation delta"));
    output.table(
        "roots", "Complex rapidities",
        [&](auto& t) {
          for (std::size_t j = 0; j < state.rapidities.size(); ++j)
            t.append(0, j, state.rapidities[j].real(), state.rapidities[j].imag());
        },
        column<std::size_t>("state_id"), column<std::size_t>("index"), column<Real>("real", "Re lambda"),
        column<Real>("imag", "Im lambda"));
  }
  output.finish();
  if (!state.converged) std::cerr << "TB solve incomplete; consider a larger budget or higher precision.\n";
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
