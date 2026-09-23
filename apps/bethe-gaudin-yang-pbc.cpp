// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/gaudin_yang.hpp>

namespace
{
namespace model = bethe::gaudin_yang;
namespace cli = bethe::cli;
struct Arguments
{
    std::size_t particles = 0, max_iterations = 10000;
    std::optional<std::string> length, interaction, tolerance;
    std::optional<uni20::half_int> sz;
    std::string precision = "fp64";
    cli::DataOutputOptions output;
    bool roots = false;
};

auto program_info()
{
  auto info = bethe::cli::program_info("bethe-gaudin-yang-pbc",
                                       "Repulsive spin-1/2 continuum fermions on a ring; finite ELL>0 and C>=0.",
                                       bethe::citations::Tool::gaudin_yang_pbc);
  info.notes = {"H=-sum_j d_j^2 + 2c sum_(i<j) delta(x_i-x_j), hbar^2/(2m)=1.",
                "Interacting mixed-spin ground states require odd N_up AND odd N_down.",
                "At c=0 or full polarization any particle count is supported, including vacuum.",
                "E=sum(k_j^2); P is signed with no Brillouin-zone reduction.",
                "Free even-population seas select the positive-current degenerate representative.",
                "Other interacting shell branches, excitations, attraction and open ends",
                "are not implemented. See docs/gaudin-yang.md for conventions and limits.",
                "Use --references for literature and applicability; see CITATIONS.md."};
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  bethe::cli::option(app, "N", args.particles, "Number of particles or sites")->required();
  bethe::cli::option(app, "--length", args.length, "required physical circumference")->required();
  bethe::cli::option(app, "--c", args.interaction, "required coupling (inverse length)")->required();
  bethe::cli::option(app, "--sz", args.sz, "spin projection (default: 0 even N, 1/2 odd N)");
  bethe::cli::option(app, "--roots", args.roots, "print momenta and spin rapidities, or free modes");
  bethe::cli::option(app, "--tolerance", args.tolerance,
                     "component-scaled equation residual (default: 32 epsilon; not an energy-error bound)");
  bethe::cli::option(app, "--max-iterations", args.max_iterations, "accepted Newton updates (default: 10000)")
      ->capture_default_str();
  bethe::cli::precision_option(app, args.precision);
  cli::add_data_output_options(app, args.output, true);
}

void validate(Arguments const& args)
{
  if (!args.length || !args.interaction) throw std::invalid_argument("--length ELL and --c C are required");
  args.output.validate();
}

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
  if (args.particles > std::size_t(std::numeric_limits<std::int64_t>::max() / 4))
    throw std::invalid_argument("Gaudin-Yang particle count exceeds the quantum-number range");
  auto const n = std::int64_t(args.particles);
  auto const sz = args.sz.value_or(uni20::from_twice(n % 2));
  auto const s = sz.twice();
  if (s < -n || s > n || (n - s) % 2)
    throw std::invalid_argument("Sz must satisfy |2*Sz|<=N and have the same parity as N");
  auto const down = std::size_t((n - s) / 2), up = args.particles - down;
  Real const length = uni20::parse_real<Real>(*args.length), c = uni20::parse_real<Real>(*args.interaction);
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  cli::CpuTimer const timer;
  auto const state = model::ground_state(up, down, length, c, options);
  auto const cpu_time = timer.elapsed_text();
  cli::report_builder report("Gaudin-Yang gas (periodic)");
  report.status(state.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning, status(state.status))
      .field("Calculation", state.free ? "exact free-fermion ground state" : "odd-population sector ground state")
      .field("Particles", args.particles)
      .field("N_up", up)
      .field("N_down", down)
      .field("Sz", uni20::to_string_fraction(sz))
      .field("Length", uni20::format_real(length))
      .field("c", uni20::format_real(c))
      .field("Units", "hbar^2/(2m)=1; interaction 2c delta")
      .field("Precision", args.precision)
      .field("Residual tolerance", uni20::format_real(options.residual_tolerance))
      .field("Status", status(state.status))
      .field("Total energy", uni20::format_real(state.energy))
      .field("Momentum index", state.momentum_index)
      .field("Momentum P", uni20::format_real(state.momentum))
      .field("Charge residual", uni20::format_real(state.charge_residual))
      .field("Spin residual", uni20::format_real(state.spin_residual))
      .field("Residual norm", uni20::format_real(state.residual_norm))
      .field("Root coupling reached", uni20::format_real(state.root_interaction))
      .field("Iterations", state.iterations)
      .field("CPU time", cpu_time);
  if (!state.free)
    report.field("Reference spin", state.spin_reversed ? "down (spin reversed)" : "up")
        .field("Spin roots", state.spin_rapidities.size());
  std::vector<std::string> tables{"states"};
  if (args.roots)
  {
    if (state.free)
      tables.insert(tables.end(), {"free_up", "free_down"});
    else
      tables.insert(tables.end(), {"charge_roots", "spin_roots"});
  }
  cli::ResultOutput output(report, args.output, "bethe-gaudin-yang-pbc", argc, argv, tables);
  output.table(
      "states", "State",
      [&](auto& table) {
        table.append(0, sz, state.energy, state.momentum_index, state.momentum, state.charge_residual,
                     state.spin_residual, state.residual_norm, state.root_interaction, state.iterations,
                     state.converged, std::string(status(state.status)));
      },
      cli::column<std::size_t>("state_id"), cli::column<uni20::half_int>("sz", "Sz"),
      cli::column<Real>("energy", "Energy"), cli::column<std::int64_t>("momentum_index"), cli::column<Real>("p", "P"),
      cli::column<Real>("charge_residual"), cli::column<Real>("spin_residual"), cli::column<Real>("residual"),
      cli::column<Real>("root_c"), cli::column<std::size_t>("iterations"), cli::column<bool>("converged"),
      cli::column<std::string>("status"));
  if (args.roots && state.free)
  {
    Real const pi = Real{4} * std::atan(Real{1});
    for (std::size_t a = 0; a < 2; ++a)
      output.table(
          a ? "free_down" : "free_up", a ? "Free down-spin modes" : "Free up-spin modes",
          [&](auto& table) {
            for (auto mode : state.free_modes[a])
              table.append(0, mode, Real{2} * pi * Real(mode) / length);
          },
          cli::column<std::size_t>("state_id"), cli::column<std::int64_t>("mode", "Mode"), cli::column<Real>("k"));
  }
  else if (args.roots)
  {
    output.table(
        "charge_roots", "Charge momenta",
        [&](auto& table) {
          for (std::size_t j = 0; j < state.momenta.size(); ++j)
            table.append(0, j, state.quantum_numbers.charge[j], state.momenta[j]);
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("index", "Index"),
        cli::column<uni20::half_int>("quantum_number", "I"), cli::column<Real>("k"));
    output.table(
        "spin_roots", "Spin rapidities",
        [&](auto& table) {
          for (std::size_t j = 0; j < state.spin_rapidities.size(); ++j)
            table.append(0, j, state.quantum_numbers.spin[j], state.spin_rapidities[j]);
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("index", "Index"),
        cli::column<uni20::half_int>("quantum_number", "J"), cli::column<Real>("rapidity", "lambda"));
  }
  output.finish();
  if (!state.converged)
    std::cerr << "Gaudin-Yang solve incomplete; residuals use requested c, not the intermediate coupling.\n";
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
