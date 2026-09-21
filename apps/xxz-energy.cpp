// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include <bethe/xxz_excitations.hpp>

#include "cli-common.hpp"
#include "excitation-report.hpp"

namespace
{
namespace model = bethe::xxz;
using bethe::cli::finish;
using bethe::cli::parse_size;
using bethe::cli::print_roots;

struct Arguments
{
    std::size_t sites;
    std::optional<std::string_view> delta = std::nullopt;
    std::string_view precision = "fp64";
    std::optional<std::string_view> tolerance = std::nullopt;
    std::size_t max_iterations = 10000;
    std::optional<uni20::half_int> sz = std::nullopt;
    bool sectors = false;
    std::optional<std::size_t> excitation_count = std::nullopt;
    std::optional<std::size_t> max_candidates = std::nullopt;
    std::optional<model::QuantumNumbers> quantum_numbers = std::nullopt;
    bool print_roots = false;
    std::string_view format = "auto";
};

void usage(std::ostream& out)
{
  out << "Usage: xxz-energy N --delta VALUE [options]\n"
      << "Periodic spin-1/2 XXZ chain, J=1, zero field; 0 <= Delta <= 1.\n"
      << "H=sum_i (Sx_i Sx_(i+1) + Sy_i Sy_(i+1) + Delta Sz_i Sz_(i+1)).\n"
      << "Default: ground state (one Sz=1/2 representative for odd N).\n"
      << "  --delta VALUE                      required anisotropy in [0,1]\n"
      << "  --sz VALUE                         lowest energy in an Sz sector, e.g. -1/2\n"
      << "  --sectors                          lowest energy in every Sz sector\n"
      << "  --excitations COUNT|all             lowest COUNT, or all, states in a restricted real-root family\n"
      << "                                     --sz selects sector (default: 1 even N, 1/2 odd N)\n"
      << "  --max-candidates COUNT             exhaustive scan limit (default: 10000)\n"
      << "  --quantum-numbers I1,I2,...         specified real-root state; use none for vacuum\n"
      << "  --precision fp64|long-double|fp128  (default: fp64)\n"
      << "  --tolerance VALUE                  normalized equation residual, max|F|/N\n"
      << "  --max-iterations COUNT             update budget (default: 10000)\n"
      << "  --roots                            print scaled rapidities z and quantum numbers\n"
      << "  --format auto|pretty|plain         terminal report or script output (default: auto)\n"
      << "  --help                             show this help\n"
      << "z=tanh(lambda)/tan(gamma/2), Delta=cos(gamma); at Delta=1, z=2*lambda_XXX.\n"
      << "Excitation scans include the sector minimum; NOT a complete Sz spectrum.\n"
      << "The finite-real window depends on Delta; strings and infinite rapidities are excluded.\n"
      << "No OBC or SU(2) multiplet classification; even-N Sz=0 has only one supported state.\n"
      << "fp128 requires a Uni20 build with MPLAPACK enabled.\n";
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
              << "Model: periodic XXZ, J=1, h=0\n"
              << "Delta: " << uni20::format_real(delta) << '\n'
              << "Precision: " << args.precision << '\n'
              << "Residual tolerance: " << uni20::format_real(options.residual_tolerance) << '\n'
              << "Root coordinate: scaled rapidity z\n";
  auto header = [&](std::string_view mode, std::string_view cpu_time) {
    bethe::cli::report_builder report("Heisenberg XXZ - " + std::string(mode));
    report.field("Model", "periodic spin-1/2, J=1, h=0")
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
    std::cout << "# CPU time: " << cpu_time << '\n' << "# Sz momentum_index P energy residual iterations converged\n";
    bool converged = true;
    for (auto const& state : states)
    {
      std::cout << state.sz << ' ' << state.momentum_index << ' ' << uni20::format_real(state.momentum) << ' '
                << uni20::format_real(state.energy) << ' ' << uni20::format_real(state.residual_norm) << ' '
                << state.iterations << ' ' << state.converged << '\n';
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
            << "Momentum index: " << state.momentum_index << '\n'
            << "Momentum: " << uni20::format_real(state.momentum) << '\n'
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
    Arguments args{.sites = parse_size(argv[1])};
    for (int i = 2; i < argc; ++i)
    {
      std::string_view const option = argv[i];
      if (option == "--roots")
        args.print_roots = true;
      else if (option == "--sectors")
        args.sectors = true;
      else if (option == "--delta" || option == "--sz" || option == "--precision" || option == "--tolerance" ||
               option == "--max-iterations" || option == "--format" || option == "--excitations" ||
               option == "--max-candidates" || option == "--quantum-numbers")
      {
        if (++i == argc) throw std::invalid_argument("missing value for " + std::string(option));
        if (option == "--delta")
          args.delta = argv[i];
        else if (option == "--sz")
          args.sz = uni20::half_int::parse(argv[i]);
        else if (option == "--precision")
          args.precision = argv[i];
        else if (option == "--tolerance")
          args.tolerance = argv[i];
        else if (option == "--format")
          args.format = argv[i];
        else if (option == "--excitations")
          args.excitation_count =
              std::string_view(argv[i]) == "all" ? std::numeric_limits<std::size_t>::max() : parse_size(argv[i]);
        else if (option == "--max-candidates")
          args.max_candidates = parse_size(argv[i]);
        else if (option == "--quantum-numbers")
          args.quantum_numbers = std::string_view(argv[i]) == "none" ? model::QuantumNumbers{}
                                                                     : bethe::cli::parse_quantum_numbers(argv[i]);
        else
          args.max_iterations = parse_size(argv[i]);
      }
      else
        throw std::invalid_argument("unknown option: " + std::string(option));
    }
    if (!args.delta) throw std::invalid_argument("--delta VALUE is required");
    if (args.sectors && args.sz) throw std::invalid_argument("--sz and --sectors are mutually exclusive");
    if (args.excitation_count && (args.sectors || args.quantum_numbers))
      throw std::invalid_argument("--excitations is mutually exclusive with --sectors and --quantum-numbers");
    if (args.quantum_numbers && (args.sectors || args.sz))
      throw std::invalid_argument("--quantum-numbers is mutually exclusive with --sectors and --sz");
    if (args.max_candidates && !args.excitation_count)
      throw std::invalid_argument("--max-candidates requires --excitations COUNT|all");
    if (args.format != "auto" && args.format != "pretty" && args.format != "plain")
      throw std::invalid_argument("unknown output format: " + std::string(args.format));
    return bethe::cli::dispatch_precision(args.precision, [&]<uni20::Real Real> { return run<Real>(args); });
  }
  catch (std::exception const& error)
  {
    std::cerr << "xxz-energy: " << error.what() << '\n';
    return 1;
  }
}
