// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "report-common.hpp"
#include <bethe/lieb_liniger.hpp>

namespace
{
namespace model = bethe::lieb_liniger;
namespace cli = bethe::cli;
struct Arguments
{
    std::size_t particles = 0;
    std::optional<std::string> length, interaction, numbers, tolerance;
    std::optional<std::size_t> count, padding, max_candidates;
    std::size_t max_iterations = 10000;
    std::string precision = "fp64", format = "auto";
    bool roots = false;
};

auto program_info()
{
  auto info = bethe::cli::program_info("bethe-lieb-liniger-pbc",
                                       "Repulsive continuum bosons on a ring: N>=0, finite ELL>0 and C>0.",
                                       bethe::citations::Tool::lieb_liniger_pbc);
  info.notes = {"H=-sum_j d_j^2 + 2c sum_(i<j) delta(x_i-x_j), hbar^2/(2m)=1.",
                "ELL is a physical length, not a lattice site count. E=sum(k_j^2).",
                "Scans include the ground state; all means the entire specified window,",
                "not the infinite continuum spectrum. COUNT does not reduce solve count.",
                "Momentum P=2*pi*sum(I)/ELL is signed, with no Brillouin-zone reduction.",
                "Attraction, c=0, c=infinity, hard walls, and thermodynamics are not implemented.",
                "See docs/lieb-liniger.md for equations, residual scaling, and numerical limits.",
                "Use --references for literature and applicability; see CITATIONS.md."};
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  bethe::cli::option(app, "N", args.particles, "Number of particles or sites")->required();
  bethe::cli::option(app, "--length", args.length, "required ring circumference")->required();
  bethe::cli::option(app, "--c", args.interaction, "required repulsive coupling (inverse length)")->required();
  bethe::cli::option(app, "--quantum-numbers", args.numbers,
                     "explicit sorted labels (default: ground state) integers for odd N, half-odd integers for even N "
                     "use none for the N=0 vacuum");
  bethe::cli::all_count_option(app, "--excitations", args.count, "retain lowest converged states in a finite window");
  bethe::cli::option(
      app, "--padding", args.padding,
      "REQUIRED for scans: P extra slots at EACH edge choose N labels from N+2P slots; P=0 is ground only");
  bethe::cli::option(app, "--max-candidates", args.max_candidates,
                     "reject larger windows before solving (default: 10000)");
  bethe::cli::option(app, "--roots", args.roots, "print physical momenta k and labels I");
  bethe::cli::option(app, "--tolerance", args.tolerance,
                     "dimensionless component-scaled residual (default: 32 epsilon; not an energy-error bound)");
  bethe::cli::option(app, "--max-iterations", args.max_iterations, "Newton updates per state (default: 10000)")
      ->capture_default_str();
  bethe::cli::precision_option(app, args.precision);
  app.add_option("--format", args.format, "Stdout layout")
      ->check(CLI::IsMember({"auto", "pretty", "plain"}))
      ->capture_default_str();
}

void validate(Arguments const& args)
{
  if (!args.length || !args.interaction) throw std::invalid_argument("--length ELL and --c C are required");
  if (args.count && !args.padding) throw std::invalid_argument("--excitations requires an explicit --padding P");
  if (!args.count && (args.padding || args.max_candidates))
    throw std::invalid_argument("--padding and --max-candidates require --excitations");
  if (args.count && args.numbers)
    throw std::invalid_argument("--quantum-numbers cannot be combined with --excitations");
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
      return "iteration limit reached; unconverged estimate";
    case model::SolveStatus::stalled:
      return "line search or representable precision stalled; unconverged estimate";
  }
  return "unknown";
}

template <uni20::Real Real>
void add_state(cli::report_builder& report, model::State<Real> const& state, std::string prefix = "")
{
  report.field(prefix + "Status", status(state.status))
      .field(prefix + "Total energy", uni20::format_real(state.energy))
      .field(prefix + "Momentum index", state.momentum_index)
      .field(prefix + "Momentum", uni20::format_real(state.momentum))
      .field(prefix + "Residual norm", uni20::format_real(state.residual_norm))
      .field(prefix + "Iterations", state.iterations);
}

template <uni20::Real Real>
void add_momenta(cli::report_builder& report, model::State<Real> const& state, std::string title)
{
  auto& table = report.table(std::move(title) + (state.particles ? "" : " (vacuum; no roots)"));
  table.header_separator().column("Index").column("I").column("k", cli::table_alignment::decimal);
  for (std::size_t j = 0; j < state.particles; ++j)
    table.row(j, uni20::to_string_fraction(state.quantum_numbers[j]), uni20::format_real(state.momenta[j]));
}

template <uni20::Real Real> int run(Arguments const& args)
{
  Real const length = uni20::parse_real<Real>(*args.length), c = uni20::parse_real<Real>(*args.interaction);
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  cli::report_builder report("Lieb-Liniger (periodic)");
  report.field("Particles", args.particles)
      .field("Length", uni20::format_real(length))
      .field("c", uni20::format_real(c))
      .field("Units", "hbar^2/(2m)=1; interaction 2c delta")
      .field("Precision", args.precision)
      .field("Residual tolerance", uni20::format_real(options.residual_tolerance))
      .field("Momentum convention", "P=2*pi*sum(I)/length; no modular reduction");
  bool converged;
  if (args.count)
  {
    bethe::RealExcitationOptions const enumeration{.count = *args.count,
                                                   .max_candidates = args.max_candidates.value_or(10000)};
    cli::CpuTimer const timer;
    auto const scan = model::real_excitations(args.particles, length, c, *args.padding, enumeration, options);
    auto const cpu_time = timer.elapsed_text();
    converged = scan.converged();
    report.field("Calculation", "finite-window excitations (including ground state)")
        .field("Padding per edge", *args.padding)
        .field("Candidate count", scan.candidate_count)
        .field("Converged count", scan.converged_count)
        .field("Retained count", scan.levels.size())
        .field("Status",
               converged ? "converged within the specified window" : "incomplete scan; failed states excluded")
        .field("CPU time", cpu_time);
    add_state(report, scan.ground_state, "Ground ");
    auto& table = report.table("Converged levels (ranked only within this window)");
    table.header_separator()
        .column("Rank")
        .column("Q")
        .column("P")
        .column("Energy")
        .column("Gap")
        .column("Residual")
        .column("Iterations");
    for (std::size_t j = 0; j < scan.levels.size(); ++j)
    {
      auto const& level = scan.levels[j];
      table.row(j, level.state.momentum_index, uni20::format_real(level.state.momentum),
                uni20::format_real(level.state.energy), level.gap ? uni20::format_real(*level.gap) : "unavailable",
                uni20::format_real(level.state.residual_norm), level.state.iterations);
      if (args.roots) add_momenta(report, level.state, "Momenta: rank " + std::to_string(j));
    }
    if (scan.first_unconverged)
    {
      add_state(report, *scan.first_unconverged, "First failed ");
      if (args.roots) add_momenta(report, *scan.first_unconverged, "First failed momenta (estimate)");
    }
  }
  else
  {
    std::optional<model::QuantumNumbers> numbers;
    if (args.numbers)
    {
      numbers = cli::parse_quantum_numbers(*args.numbers == "none" ? "" : *args.numbers);
      if (numbers->size() != args.particles) throw std::invalid_argument("quantum-number count must equal N");
    }
    cli::CpuTimer const timer;
    auto const state = numbers ? model::solve_real(length, c, *numbers, options)
                               : model::ground_state(args.particles, length, c, options);
    auto const cpu_time = timer.elapsed_text();
    converged = state.converged;
    report.field("Calculation", args.numbers ? "specified Bethe state" : "ground state").field("CPU time", cpu_time);
    add_state(report, state);
    if (args.roots) add_momenta(report, state, "Physical momenta");
  }
  report.status(converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning,
                converged ? "converged" : "unconverged estimates are not ranked");
  cli::print_report(report, args.format);
  if (!converged) std::cerr << "Lieb-Liniger solve incomplete; consider a larger budget or higher precision.\n";
  return converged ? 0 : 2;
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
