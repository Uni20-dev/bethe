// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include <bethe/xxz_excitations.hpp>

#include "citation-report.hpp"
#include "excitation-report.hpp"
#include "xxz-cli.hpp"

namespace
{
namespace model = bethe::xxz::open;
using bethe::cli::finish;
using bethe::cli::parse_size;
using bethe::cli::print_roots;

using Arguments = bethe::cli::XxzArguments;

void usage(std::ostream& out)
{
  out << "Usage: bethe-xxz-obc N --delta VALUE [options]\n"
      << "Open spin-1/2 XXZ chain, free ends, J=1, zero field; 0 <= Delta <= 1.\n"
      << "H=sum_(i=0)^(N-2) (Sx_i Sx_(i+1) + Sy_i Sy_(i+1) + Delta Sz_i Sz_(i+1)).\n"
      << "Default: ground state (one Sz=1/2 representative for odd N).\n"
      << "  --delta VALUE                      required anisotropy in [0,1]\n"
      << "  --sz VALUE                         lowest energy in an Sz sector, e.g. -1/2\n"
      << "  --sectors                          lowest energy in every Sz sector\n"
      << "  --excitations COUNT|all             lowest COUNT, or all, states in a restricted real-root family\n"
      << "                                     --sz selects sector (default: 1 even N, 1/2 odd N)\n"
      << "  --max-candidates COUNT             exhaustive scan limit (default: 10000)\n"
      << "  --quantum-numbers I1,I2,...         positive integer labels; use none for vacuum\n"
      << "  --precision fp64|long-double|fp128  (default: fp64)\n"
      << "  --tolerance VALUE                  normalized equation residual, max|F|/(2N)\n"
      << "  --max-iterations COUNT             update budget (default: 10000)\n"
      << "  --roots                            print scaled rapidities z and quantum numbers\n"
      << "  --format auto|pretty|plain         terminal report or script output (default: auto)\n"
      << "  --help                             show this help\n"
      << "z=tanh(lambda)/tan(gamma/2), Delta=cos(gamma); at Delta=1, z=2*lambda_XXX.\n"
      << "Excitation scans include the sector minimum; NOT a complete Sz spectrum.\n"
      << "The finite-real window depends on Delta; strings and infinite rapidities are excluded.\n"
      << "No boundary fields, lattice momentum, or SU(2) multiplet classification.\n"
      << "fp128 requires a Uni20 build with MPLAPACK enabled.\n";
  bethe::cli::print_citations(out, bethe::citations::Tool::xxz_obc);
}

template <uni20::Real Real> int run(Arguments const& args)
{
  Real const delta = uni20::parse_real<Real>(*args.delta);
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  bool const pretty = args.format == "pretty" || (args.format == "auto" && terminal::is_a_terminal(stdout));
  if (!pretty)
    std::cout << "Sites: " << args.sites << '\n'
              << "Model: open XXZ, free ends, J=1, h=0\n"
              << "Delta: " << uni20::format_real(delta) << '\n'
              << "Precision: " << args.precision << '\n'
              << "Residual tolerance: " << uni20::format_real(options.residual_tolerance) << '\n'
              << "Root coordinate: scaled rapidity z\n";
  auto header = [&](std::string_view mode, std::string_view cpu_time) {
    bethe::cli::report_builder report("Heisenberg XXZ (free ends) - " + std::string(mode));
    report.field("Model", "open spin-1/2, free ends, J=1, h=0")
        .field("Sites", args.sites)
        .field("Delta", uni20::format_real(delta))
        .field("Precision", args.precision)
        .field("Residual tolerance", uni20::format_real(options.residual_tolerance))
        .field("CPU time", cpu_time)
        .field("Root coordinate", "scaled rapidity z");
    return report;
  };
  bethe::cli::CpuTimer const timer;
  if (args.excitation_count)
  {
    auto const sz = args.sz.value_or(uni20::from_twice(std::int64_t{args.sites % 2 == 0 ? 2 : 1}));
    auto const scan = model::real_excitations<Real>(
        args.sites, delta, sz, {.count = *args.excitation_count, .max_candidates = args.max_candidates.value_or(10000)},
        options);
    auto const cpu_time = timer.elapsed_text();
    auto report = header("real-root excitations", cpu_time);
    auto const window = scan.window.slots == 0 ? "empty (polarized vacuum)"
                                               : uni20::to_string_fraction(scan.window.first) + ".." +
                                                     uni20::to_string_fraction(scan.window.last);
    report.field("Quantum-number window", window)
        .field("Available slots", scan.window.slots)
        .field("Spin-reversed reference", sz.twice() < 0 ? "yes" : "no");
    if (!pretty)
      std::cout << "# Quantum-number window: " << window << '\n'
                << "# Available slots: " << scan.window.slots << '\n'
                << "# Spin-reversed reference: " << (sz.twice() < 0) << '\n';
    return finish(
        bethe::cli::print_excitation_report(std::move(report), scan,
                                            {.family = "restricted finite-real XXZ states; NOT a complete Sz spectrum",
                                             .sector_label = "Sz",
                                             .sector = sz,
                                             .multiplet_size = std::nullopt},
                                            args.print_roots, cpu_time, pretty));
  }
  if (args.sectors)
  {
    auto const states = model::sector_ground_states<Real>(args.sites, delta, options);
    auto const cpu_time = timer.elapsed_text();
    if (pretty)
      return finish(bethe::cli::print_sector_report(header("sector minima", cpu_time), states, args.print_roots));
    std::cout << "# CPU time: " << cpu_time << '\n' << "# Sz energy residual iterations converged\n";
    bool converged = true;
    for (auto const& state : states)
    {
      std::cout << state.sz << ' ' << uni20::format_real(state.energy) << ' ' << uni20::format_real(state.residual_norm)
                << ' ' << state.iterations << ' ' << state.converged << '\n';
      if (args.print_roots)
      {
        std::cout << "# Roots: Sz=" << state.sz << '\n';
        print_roots(state);
      }
      converged = converged && state.converged;
    }
    return finish(converged);
  }
  auto const state = args.quantum_numbers ? model::solve_real<Real>(args.sites, delta, *args.quantum_numbers, options)
                     : args.sz            ? model::sector_ground_state<Real>(args.sites, delta, *args.sz, options)
                                          : model::ground_state<Real>(args.sites, delta, options);
  auto const cpu_time = timer.elapsed_text();
  if (pretty)
    return finish(bethe::cli::print_state_report<Real>(header(args.quantum_numbers ? "specified real-root state"
                                                              : args.sz            ? "sector minimum"
                                                                                   : "ground state",
                                                              cpu_time),
                                                       args.sites, state, args.print_roots));
  std::cout << "Sz: " << state.sz << '\n'
            << "Spin-reversed reference: " << state.spin_reversed << '\n'
            << "Status: " << (state.converged ? "converged" : "iteration limit reached") << '\n'
            << "Iterations: " << state.iterations << '\n'
            << "CPU time: " << cpu_time << '\n'
            << "Residual norm: " << uni20::format_real(state.residual_norm) << '\n'
            << "Total energy: " << uni20::format_real(state.energy) << '\n'
            << "Energy per site: " << uni20::format_real(state.energy / static_cast<Real>(args.sites)) << '\n';
  if (args.print_roots) print_roots(state);
  return finish(state.converged);
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
    auto const args = bethe::cli::parse_xxz_arguments(argc, argv);
    return bethe::cli::dispatch_precision(args.precision, [&]<uni20::Real Real> { return run<Real>(args); });
  }
  catch (std::exception const& error)
  {
    std::cerr << "bethe-xxz-obc: " << error.what() << '\n';
    return 1;
  }
}
