// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
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
    std::string precision = "fp64";
    cli::DataOutputOptions output;
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
  cli::add_data_output_options(app, args.output, true);
}

void validate(Arguments const& args)
{
  if (!args.length || !args.interaction) throw std::invalid_argument("--length ELL and --c C are required");
  if (args.count && !args.padding) throw std::invalid_argument("--excitations requires an explicit --padding P");
  if (!args.count && (args.padding || args.max_candidates))
    throw std::invalid_argument("--padding and --max-candidates require --excitations");
  if (args.count && args.numbers)
    throw std::invalid_argument("--quantum-numbers cannot be combined with --excitations");
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

template <uni20::Real Real>
void add_state(cli::RunReport& report, model::State<Real> const& state, std::string id_prefix = "",
               std::string prefix = "")
{
  if (id_prefix.empty())
    report.result(state.converged, status(state.status));
  else
    report.field(id_prefix + "status", prefix + "Status", status(state.status));
  report.field(id_prefix + "total_energy", prefix + "Total energy", state.energy)
      .field(id_prefix + "momentum_index", prefix + "Momentum index", state.momentum_index)
      .field(id_prefix + "momentum", prefix + "Momentum", state.momentum)
      .field(id_prefix + "residual_norm", prefix + "Residual norm", state.residual_norm)
      .field(id_prefix + "iterations", prefix + "Iterations", state.iterations);
}

template <uni20::Real Real>
void write_output(cli::RunReport& report, Arguments const& args, std::vector<model::State<Real> const*> const& states,
                  std::vector<std::optional<Real>> const& gaps, model::State<Real> const* reference = nullptr,
                  model::State<Real> const* failed = nullptr)
{
  std::vector<std::string> names{"states"};
  if (reference) names.push_back("reference");
  if (failed) names.push_back("failed");
  if (args.roots) names.push_back("roots");
  cli::ResultOutput output(report, args.output, names);
  auto write_states = [&](std::string name, std::string title, auto const& rows, std::size_t offset, bool ranked) {
    output.table(
        name, title,
        [&](auto& table) {
          for (std::size_t i = 0; i < rows.size(); ++i)
          {
            auto const& s = *rows[i];
            table.append(offset + i, s.energy, ranked ? gaps[i] : std::nullopt, s.momentum_index, s.momentum,
                         s.residual_norm, s.iterations, s.converged, std::string(status(s.status)));
          }
        },
        cli::column<std::size_t>("state_id"), cli::column<Real>("energy", "Energy"),
        cli::column<std::optional<Real>>("gap", "Gap"), cli::column<std::int64_t>("momentum_index", "Q"),
        cli::column<Real>("p", "P"), cli::column<Real>("residual", "Residual"),
        cli::column<std::size_t>("iterations", "Iterations"), cli::column<bool>("converged"),
        cli::column<std::string>("status", "Status"));
  };
  write_states("states", reference ? "Converged levels (ranked only within this window)" : "State", states, 0, true);
  std::vector<model::State<Real> const*> all = states;
  if (reference)
  {
    write_states("reference", "Ground reference", std::vector{reference}, all.size(), false);
    all.push_back(reference);
  }
  if (failed)
  {
    write_states("failed", "First failed state (unranked estimate)", std::vector{failed}, all.size(), false);
    all.push_back(failed);
  }
  if (args.roots)
    output.table(
        "roots", args.particles ? "Physical momenta" : "Physical momenta (vacuum; no roots)",
        [&](auto& table) {
          for (std::size_t i = 0; i < all.size(); ++i)
            for (std::size_t j = 0; j < all[i]->momenta.size(); ++j)
              table.append(i, j, all[i]->quantum_numbers[j], all[i]->momenta[j]);
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("index", "Index"),
        cli::column<uni20::half_int>("quantum_number", "I"), cli::column<Real>("k"));
  output.finish();
}

template <uni20::Real Real> int run(Arguments const& args, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const length = uni20::parse_real<Real>(*args.length), c = uni20::parse_real<Real>(*args.interaction);
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  cli::RunReport report(context, "Lieb-Liniger (periodic)");
  report.field("particles", "Particles", args.particles)
      .field("length", "Length", length)
      .field("c", "c", c)
      .field("units", "Units", "hbar^2/(2m)=1; interaction 2c delta")
      .field("precision", "Precision", args.precision)
      .field("residual_tolerance", "Residual tolerance", options.residual_tolerance)
      .field("momentum_convention", "Momentum convention", "P=2*pi*sum(I)/length; no modular reduction");
  bool converged;
  if (args.count)
  {
    bethe::RealExcitationOptions const enumeration{.count = *args.count,
                                                   .max_candidates = args.max_candidates.value_or(10000)};
    auto computation = context.computation();
    auto const scan = model::real_excitations(args.particles, length, c, *args.padding, enumeration, options);
    computation.finish();
    converged = scan.converged();
    report.field("calculation", "Calculation", "finite-window excitations (including ground state)")
        .field("padding_per_edge", "Padding per edge", *args.padding)
        .field("candidate_count", "Candidate count", scan.candidate_count)
        .field("converged_count", "Converged count", scan.converged_count)
        .field("retained_count", "Retained count", scan.levels.size())
        .result(converged,
                converged ? "converged within the specified window" : "incomplete scan; failed states excluded");
    add_state(report, scan.ground_state, "ground_", "Ground ");
    if (scan.first_unconverged) add_state(report, *scan.first_unconverged, "first_failed_", "First failed ");
    report.status(converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning,
                  converged ? "converged" : "unconverged estimates are not ranked");
    std::vector<model::State<Real> const*> states;
    std::vector<std::optional<Real>> gaps;
    for (auto const& level : scan.levels)
    {
      states.push_back(&level.state);
      gaps.push_back(level.gap);
    }
    write_output(report, args, states, gaps, &scan.ground_state,
                 scan.first_unconverged ? &*scan.first_unconverged : nullptr);
  }
  else
  {
    std::optional<model::QuantumNumbers> numbers;
    if (args.numbers)
    {
      numbers = cli::parse_quantum_numbers(*args.numbers == "none" ? "" : *args.numbers);
      if (numbers->size() != args.particles) throw std::invalid_argument("quantum-number count must equal N");
    }
    auto computation = context.computation();
    auto const state = numbers ? model::solve_real(length, c, *numbers, options)
                               : model::ground_state(args.particles, length, c, options);
    computation.finish();
    converged = state.converged;
    report.field("calculation", "Calculation", args.numbers ? "specified Bethe state" : "ground state");
    add_state(report, state);
    report.status(converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning,
                  converged ? "converged" : "unconverged estimate");
    write_output<Real>(report, args, {&state}, {std::nullopt});
  }

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
        return bethe::cli::dispatch_precision(args.precision,
                                              [&]<uni20::Real Real> { return run<Real>(args, argc, argv); });
      });
}
