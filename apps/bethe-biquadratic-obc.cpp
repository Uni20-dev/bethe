// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "excitation-report.hpp"
#include "program-options.hpp"
#include <bethe/biquadratic.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::biquadratic;
struct Arguments
{
    std::size_t sites = 0, max_iterations = 10000;
    std::optional<std::size_t> through_lines, excitations, max_candidates;
    std::optional<bethe::xxz::QuantumNumbers> numbers;
    std::optional<std::string> tolerance;
    std::string precision = "fp64", format = "auto";
    bool roots = false, sectors = false;
};
auto program_info()
{
  auto info = bethe::cli::program_info("bethe-biquadratic-obc",
                                       "Spin-1 pure biquadratic chain, free ends: H=-sum_i (S_i.S_(i+1))^2.",
                                       bethe::citations::Tool::biquadratic_obc);
  info.notes = {"Unique singlet ground state; even N>=2, coefficient -1.",
                "TL loop weight 3; reference XXZ Delta=3/2 with opposite end fields.",
                "Real-root scans are NOT complete spectra: complex-root levels are excluded.",
                "Multiplicity counts physical states per TL eigenvector, not SU(2) multiplets.",
                "TL through-lines are not physical spin; no odd chains or lattice momentum.",
                "This is not the TB point, ULS point, or zero-boundary-field XXZ chain.",
                "See docs/biquadratic.md for the TL mapping and representation multiplicities.",
                "Use --references for literature and applicability; see CITATIONS.md."};
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  bethe::cli::option(app, "N", args.sites, "Number of particles or sites")->required();
  bethe::cli::option(app, "--through-lines", args.through_lines, "lowest level in an even TL module, 0<=ELL<=N");
  bethe::cli::option(app, "--sectors", args.sectors, "lowest level in every TL module");
  bethe::cli::all_count_option(
      app, "--excitations", args.excitations,
      "lowest COUNT, or all, supported real-root levels default ELL=2; includes module minimum");
  bethe::cli::option(app, "--max-candidates", args.max_candidates, "exhaustive scan limit (default: 10000)");
  bethe::cli::option(app, "--quantum-numbers", args.numbers, "explicit integer labels; none for vacuum");
  bethe::cli::option(app, "--roots", args.roots, "print reference XXZ rapidities and labels");
  bethe::cli::option(app, "--tolerance", args.tolerance,
                     "max logarithmic residual divided by 2*N (default: 32 epsilon; not an energy-error bound)");
  bethe::cli::option(app, "--max-iterations", args.max_iterations, "accepted Newton updates (default: 10000)")
      ->capture_default_str();
  bethe::cli::precision_option(app, args.precision);
  app.add_option("--format", args.format, "Stdout layout")
      ->check(CLI::IsMember({"auto", "pretty", "plain"}))
      ->capture_default_str();
}

void validate(Arguments const& args)
{
  if (args.format != "auto" && args.format != "plain" && args.format != "pretty")
    throw std::invalid_argument("unknown output format: " + std::string(args.format));
  if (args.sectors && (args.through_lines || args.excitations || args.numbers))
    throw std::invalid_argument(
        "--sectors cannot be combined with --through-lines, --excitations or --quantum-numbers");
  if (args.numbers && (args.through_lines || args.excitations))
    throw std::invalid_argument(
        "--quantum-numbers determines the TL module; cannot combine with --through-lines or --excitations");
  if (args.max_candidates && !args.excitations) throw std::invalid_argument("--max-candidates requires --excitations");
}
char const* status(bethe::xxz::quantum_group::SolveStatus value)
{
  using Status = bethe::xxz::quantum_group::SolveStatus;
  switch (value)
  {
    case Status::converged:
      return "converged";
    case Status::iteration_limit:
      return "iteration limit reached; unconverged estimate";
    case Status::singular_jacobian:
      return "singular or ill-conditioned Jacobian; unconverged estimate";
    case Status::stalled:
      return "line search or representable precision stalled; unconverged estimate";
  }
  return "unknown";
}
std::string multiplicity_text(std::optional<std::uint64_t> value)
{
  return value ? std::to_string(*value) : "overflow (>uint64)";
}

template <typename Reference> void add_roots(cli::report_builder& report, Reference const& reference, std::string title)
{
  auto& table = report.table(std::move(title));
  table.header_separator()
      .column("Index")
      .column("I")
      .column("alpha", cli::table_alignment::decimal)
      .column("x=Theta_1/2", cli::table_alignment::decimal);
  for (std::size_t j = 0; j < reference.rapidities.size(); ++j)
    table.row(j, uni20::to_string(reference.quantum_numbers[j]), uni20::format_real(reference.rapidities[j]),
              uni20::format_real(reference.angles[j]));
}

template <uni20::Real Real> auto preamble(Arguments const& args, bethe::SolverOptions<Real> const& options)
{
  cli::report_builder report("Spin-1 pure biquadratic chain (free ends)");
  report.field("Hamiltonian", "H=-sum_i (S_i.S_(i+1))^2")
      .field("Sites", args.sites)
      .field("Spin", 1)
      .field("TL loop weight", 3)
      .field("XXZ Delta", "1.5")
      .field("XXZ reference", "spin-half exchange 1; +sqrt(5)/4*(sz_1-sz_N)")
      .field("Precision", args.precision)
      .field("Residual tolerance", uni20::format_real(options.residual_tolerance));
  return report;
}

int finish(cli::report_builder const& report, Arguments const& args, bool converged)
{
  cli::print_report(report, args.format);
  std::cout.flush();
  if (!std::cout) throw std::ios_base::failure("output flush failed");
  if (!converged) std::cerr << "Biquadratic solve incomplete; consider a larger budget or higher precision.\n";
  return converged ? 0 : 2;
}

template <uni20::Real Real> int run(Arguments const& args)
{
  bethe::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  auto report = preamble(args, options);
  cli::CpuTimer const timer;
  if (args.excitations)
  {
    auto const ell = args.through_lines.value_or(2);
    auto const scan = model::real_excitations<Real>(
        args.sites, ell, {.count = *args.excitations, .max_candidates = args.max_candidates.value_or(10000)}, options);
    auto const cpu_time = timer.elapsed_text();
    auto const& ground = scan.ground_state;
    report.field("Calculation", "restricted real-root excitations")
        .field("Family", "positive finite roots; complex-root levels excluded; NOT the complete TL spectrum")
        .field("TL through-lines", ell)
        .field("Multiplicity per TL eigenvector",
               multiplicity_text(bethe::temperley_lieb::spin_chain_multiplicity(3, ell)))
        .field("Multiplicity meaning", "physical states, not SU(2) multiplets; not an accidental-degeneracy sum")
        .field("Candidates", scan.candidate_count)
        .field("Converged candidates", scan.converged_count)
        .field("Returned levels", scan.levels.size())
        .field("Ordering",
               scan.family_converged() ? "complete within supported family" : "incomplete; failed candidates excluded")
        .field("Ground energy", uni20::format_real(ground.energy))
        .field("Ground status", status(ground.reference.status))
        .field("Ground residual", uni20::format_real(ground.reference.residual_norm))
        .field("Ground iterations", ground.reference.iterations)
        .field("Gap reference",
               ground.reference.converged ? "E-E0; global singlet ground state" : "unavailable; ground solve failed")
        .field("Status", scan.converged() ? "converged" : "incomplete scan or ground reference")
        .field("CPU time", cpu_time);
    if (scan.first_unconverged)
      report.field("First failed I", cli::quantum_number_text(scan.first_unconverged->reference.quantum_numbers))
          .field("First failed status", status(scan.first_unconverged->reference.status))
          .field("First failed residual", uni20::format_real(scan.first_unconverged->reference.residual_norm));
    auto& levels = report.table("Real-root TL levels (module minimum included)");
    levels.header_separator().column("Level").column("Multiplicity").column("Energy").column("E-E0");
    auto& diagnostics = report.table("Convergence and reference quantum numbers");
    diagnostics.header_separator().column("Level").column("Residual").column("Iterations").column("I");
    for (std::size_t j = 0; j < scan.levels.size(); ++j)
    {
      auto const& level = scan.levels[j];
      auto const& r = level.state.reference;
      levels.row(j + 1, multiplicity_text(level.state.multiplicity), uni20::format_real(level.state.energy),
                 level.gap ? uni20::format_real(*level.gap) : "unavailable");
      diagnostics.row(j + 1, uni20::format_real(r.residual_norm), r.iterations,
                      cli::quantum_number_text(r.quantum_numbers));
    }
    if (args.roots)
      for (std::size_t j = 0; j < scan.levels.size(); ++j)
        add_roots(report, scan.levels[j].state.reference, "Reference XXZ roots: level=" + std::to_string(j + 1));
    return finish(report, args, scan.converged());
  }
  if (args.sectors)
  {
    std::vector<model::State<Real>> states;
    states.push_back(model::ground_state<Real>(args.sites, options));
    bool converged = states.front().reference.converged;
    for (std::size_t ell = 2; ell <= args.sites; ell += 2)
    {
      states.push_back(model::sector_ground_state<Real>(args.sites, ell, options));
      converged = converged && states.back().reference.converged;
    }
    auto const cpu_time = timer.elapsed_text();
    report.field("Calculation", "TL module minima (not physical-spin sectors)")
        .field("Multiplicity meaning", "physical states per TL eigenvector, not SU(2) multiplets")
        .field("Status", converged ? "converged" : "incomplete; unconverged estimates")
        .field("CPU time", cpu_time);
    auto& levels = report.table("TL module minima");
    levels.header_separator().column("Through-lines").column("Multiplicity").column("Energy").column("E-E0");
    auto& diagnostics = report.table("Convergence by TL module");
    diagnostics.header_separator().column("Through-lines").column("Residual").column("Iterations").column("Status");
    for (auto const& state : states)
    {
      auto const& r = state.reference;
      auto const gap = r.converged && states.front().reference.converged
                           ? uni20::format_real(Real{2} * (r.energy - states.front().reference.energy))
                           : "unavailable";
      levels.row(state.through_lines, multiplicity_text(state.multiplicity), uni20::format_real(state.energy), gap);
      diagnostics.row(state.through_lines, uni20::format_real(r.residual_norm), r.iterations, status(r.status));
    }
    if (args.roots)
      for (auto const& state : states)
        add_roots(report, state.reference, "Reference XXZ roots: through-lines=" + std::to_string(state.through_lines));
    return finish(report, args, converged);
  }
  auto const state = args.numbers
                         ? model::solve_real<Real>(args.sites, *args.numbers, options)
                         : model::sector_ground_state<Real>(args.sites, args.through_lines.value_or(0), options);
  auto const cpu_time = timer.elapsed_text();
  auto const& reference = state.reference;
  report
      .status(reference.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning,
              status(reference.status))
      .field("Calculation", args.numbers               ? "specified real-root TL level"
                            : state.through_lines == 0 ? "even-chain singlet ground state"
                                                       : "TL module minimum")
      .field("TL through-lines", state.through_lines)
      .field(state.through_lines == 0 ? "Ground-state multiplicity" : "Multiplicity per TL eigenvector",
             multiplicity_text(state.multiplicity))
      .field("Status", status(reference.status))
      .field("Total energy", uni20::format_real(state.energy))
      .field("Energy per site", uni20::format_real(state.energy / Real(args.sites)))
      .field("TL energy (-sum e_i)", uni20::format_real(state.tl_energy))
      .field("XXZ reference energy", uni20::format_real(reference.energy))
      .field("Residual norm", uni20::format_real(reference.residual_norm))
      .field("Iterations", reference.iterations)
      .field("CPU time", cpu_time);
  if (state.through_lines == 0) report.field("Total spin", 0);
  if (state.through_lines != 0) report.field("Multiplicity meaning", "physical states, not a physical-spin label");
  if (args.roots) add_roots(report, reference, "Reference XXZ roots (not physical spin-1 quantum numbers)");
  return finish(report, args, reference.converged);
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
