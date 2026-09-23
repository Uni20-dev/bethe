// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/central_spin.hpp>

namespace
{
namespace model = bethe::central_spin;
namespace cli = bethe::cli;
struct Arguments
{
    std::optional<std::string> couplings, field, tolerance;
    std::optional<uni20::half_int> sz;
    std::string precision = "fp64";
    cli::DataOutputOptions output;
    std::size_t max_iterations = 10000, max_stages = 10000;
    bool variables = false;
};
auto program_info()
{
  auto info = bethe::cli::program_info("bethe-central-spin",
                                       "Rational Gaudin central-spin ground state in a specified total-Sz sector.",
                                       bethe::citations::Tool::central_spin);
  info.notes = {
      "H=B*S0^z+sum_j A_j*S0.Sj; central and bath spins are all 1/2.",
      "Default Newton budget: 10000. A zero budget leaves the infinite-field seed, with no energy.",
      "Incomplete solves report energy only at the REACHED field, not the requested field.",
      "An empty coupling list describes an isolated central spin. Coupling order is immaterial.",
      "Repeated/zero couplings, bath interactions, excitations and rapidity reconstruction are not implemented.",
      "No PBC/OBC or lattice momentum applies; a sector minimum is not necessarily the global minimum.",
      "See docs/central-spin.md for conventions, state selection and continuation controls.",
      "Use --references for literature and applicability; see CITATIONS.md."};
  info.notes.push_back("Tables: states; variables with --variables. Null energy/reached_field means no finite-field stage was reached.");
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  bethe::cli::option(app, "--couplings", args.couplings, "distinct nonzero real bath couplings (either sign)")
      ->required();
  bethe::cli::option(app, "--field", args.field, "central field, either sign or exactly zero")->required();
  bethe::cli::option(app, "--sz", args.sz, "required total magnetization, including the central spin")->required();
  bethe::cli::option(app, "--variables", args.variables, "print compactified eigenvalue variables (not occupations)");
  bethe::cli::option(app, "--max-stages", args.max_stages, "attempted continuation stages (default: 10000)")
      ->capture_default_str();
  bethe::cli::option(app, "--tolerance", args.tolerance, "polynomial backward residual (default: 32 epsilon)");
  bethe::cli::option(app, "--max-iterations", args.max_iterations, "attempted Newton corrections, including retries")
      ->capture_default_str();
  bethe::cli::precision_option(app, args.precision);
  cli::add_data_output_options(app, args.output, true);
}

void validate(Arguments const& args)
{
  if (!args.couplings || !args.field || !args.sz)
    throw std::invalid_argument("--couplings, --field and --sz are required");
  args.output.validate();
}
char const* status(model::SolveStatus value)
{
  switch (value)
  {
    case model::SolveStatus::converged:
      return "converged";
    case model::SolveStatus::iteration_limit:
      return "Newton budget exhausted; target field not reached";
    case model::SolveStatus::stage_limit:
      return "stage budget exhausted; target field not reached";
    case model::SolveStatus::stalled:
      return "continuation stalled; target field not reached";
    case model::SolveStatus::ill_conditioned:
      return "linear system unresolved; target field not reached";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& args, int argc, char** argv)
{
  std::vector<Real> a;
  std::string_view list = *args.couplings;
  while (!list.empty())
  {
    auto const comma = list.find(',');
    if (list.substr(0, comma).empty()) throw std::invalid_argument("empty coupling list item");
    a.push_back(uni20::parse_real<Real>(list.substr(0, comma)));
    if (comma == std::string_view::npos) break;
    list.remove_prefix(comma + 1);
    if (list.empty()) throw std::invalid_argument("empty coupling list item");
  }
  Real const b = uni20::parse_real<Real>(*args.field);
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  options.max_stages = args.max_stages;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  cli::CpuTimer const timer;
  auto const state = model::sector_ground_state<Real>(a, b, *args.sz, options);
  auto const cpu_time = timer.elapsed_text();
  cli::report_builder report("Rational Gaudin central spin");
  report.status(state.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning, status(state.status))
      .field("Hamiltonian", "B*S0^z+sum_j A_j*S0.Sj (all spins 1/2)")
      .field("Calculation", "lowest state in the specified total-Sz sector")
      .field("Bath spins", a.size())
      .field("Total Sz", uni20::to_string_fraction(state.sz))
      .field("Up spins", state.up_spins)
      .field("Precision", args.precision)
      .field("Spin reversed", state.spin_reversed ? "yes" : "no")
      .field("Requested field B", uni20::format_real(state.field))
      .field("Reached field B", state.reached_field ? uni20::format_real(*state.reached_field) : "infinite-field seed")
      .field("Status", status(state.status));
  if (state.energy)
    report.field("Energy evaluated at B", uni20::format_real(*state.reached_field))
        .field("Total energy", uni20::format_real(*state.energy));
  else
    report.field("Energy", "unavailable: no finite-field stage reached");
  report.field("Residual tolerance", uni20::format_real(options.residual_tolerance))
      .field("Reached backward residual", uni20::format_real(state.residual_norm))
      .field("Number-constraint error", uni20::format_real(state.number_error))
      .field("Newton corrections", state.iterations)
      .field("Continuation stages", state.stages)
      .field("Rejected stages", state.rejected_stages)
      .field("CPU time", cpu_time);
  using cli::column;
  std::vector<std::string> names{"states"};
  if (args.variables) names.push_back("variables");
  cli::ResultOutput output(report, args.output, "bethe-central-spin", argc, argv, names);
  output.table(
      "states", "State",
      [&](auto& t) {
        t.append(0, state.sz, state.energy, state.field, state.reached_field, state.residual_norm, state.number_error,
                 state.iterations, state.stages, state.rejected_stages, state.converged, status(state.status));
      },
      column<std::size_t>("state_id"), column<uni20::half_int>("sz"), column<std::optional<Real>>("energy"),
      column<Real>("requested_field"), column<std::optional<Real>>("reached_field"), column<Real>("residual"),
      column<Real>("number_error"), column<std::size_t>("iterations"), column<std::size_t>("stages"),
      column<std::size_t>("rejected_stages"), column<bool>("converged"), column<std::string>("status"));
  if (args.variables)
    output.table(
        "variables", "Eigenvalue variables (not occupations; spin-reversed frame if B<0)",
        [&](auto& t) {
          for (std::size_t j = 0; j <= a.size(); ++j)
            t.append(0, j, j ? std::optional<Real>{a[j - 1]} : std::nullopt,
                     state.eigenvalue_variables.empty() ? std::nullopt
                                                        : std::optional<Real>{state.eigenvalue_variables[j]});
        },
        column<std::size_t>("state_id"), column<std::size_t>("spin"),
        column<std::optional<Real>>("coupling").description("Null identifies the central spin"),
        column<std::optional<Real>>("v").description("Null for an analytic state; not an occupation"));
  output.finish();
  if (!state.converged) std::cerr << "Central-spin solve incomplete: no energy at the requested field is claimed.\n";
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
