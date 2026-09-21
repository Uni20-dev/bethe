// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include <bethe/heisenberg_open.hpp>

#include "heisenberg-cli.hpp"
#include "heisenberg-report.hpp"

#include <optional>

namespace
{
using bethe::cli::CpuTimer;
using bethe::cli::finish;
using bethe::cli::parse_quantum_numbers;
using bethe::cli::parse_size;
using bethe::cli::print_roots;
namespace model = bethe::heisenberg::open;

struct Arguments
{
    std::size_t sites;
    std::string_view precision = "fp64";
    std::optional<std::string_view> tolerance = std::nullopt;
    std::size_t max_iterations = 10000;
    bool print_roots = false;
    std::optional<uni20::half_int> sz = std::nullopt;
    std::optional<std::string_view> quantum_numbers = std::nullopt;
    bool sectors = false;
    std::string_view format = "auto";
};

void usage(std::ostream& out)
{
  out << "Usage: heisenberg-open-energy N [options]\n"
      << "Open spin-1/2 Heisenberg chain, free ends, J=1, zero field.\n"
      << "Default: ground state (Sz=0 for even N, Sz=1/2 for odd N). Modes:\n"
      << "  --sz VALUE                         lowest energy in an Sz sector, e.g. 1/2\n"
      << "  --sectors                          lowest energy in every Sz sector\n"
      << "  --quantum-numbers I0,I1,...         distinct integers in [1,N-M], M<=N/2\n"
      << "  --precision fp64|long-double|fp128  (default: fp64)\n"
      << "  --tolerance VALUE                  normalized equation residual, max|F|/(2N)\n"
      << "  --max-iterations COUNT             update budget (default: 10000)\n"
      << "  --roots                            print the positive rapidities\n"
      << "  --format auto|pretty|plain         terminal report or script output (default: auto)\n"
      << "  --help                             show this help\n"
      << "No boundary fields, complex strings, or lattice momentum.\n"
      << "fp128 requires a Uni20 build with MPLAPACK enabled.\n";
}

template <uni20::Real Real> int run(Arguments const& args)
{
  bethe::heisenberg::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  bool const pretty = args.format == "pretty" || (args.format == "auto" && terminal::is_a_terminal(stdout));
  if (!pretty)
    std::cout << "Sites: " << args.sites << '\n'
              << "Boundary: open (free ends)\n"
              << "Precision: " << args.precision << '\n'
              << "Residual tolerance: " << uni20::format_real(options.residual_tolerance) << '\n';
  CpuTimer const timer;
  if (args.sectors)
  {
    auto const states = model::sector_ground_states<Real>(args.sites, options);
    auto const cpu_time = timer.elapsed_text();
    if (pretty)
      return finish(bethe::cli::print_sectors(args.sites, args.precision, options, states, args.print_roots, cpu_time));
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
  auto const state = args.quantum_numbers
                         ? model::solve_real<Real>(args.sites, parse_quantum_numbers(*args.quantum_numbers), options)
                     : args.sz ? model::sector_ground_state<Real>(args.sites, *args.sz, options)
                               : model::ground_state<Real>(args.sites, options);
  auto const cpu_time = timer.elapsed_text();
  if (pretty)
    return finish(bethe::cli::print_state(args.sites, args.precision, options, state,
                                          args.quantum_numbers ? "specified real-root state"
                                          : args.sz            ? "sector minimum"
                                                               : "ground state",
                                          args.print_roots, cpu_time));
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
    Arguments args{.sites = parse_size(argv[1])};
    for (int i = 2; i < argc; ++i)
    {
      std::string_view const option = argv[i];
      if (option == "--roots")
        args.print_roots = true;
      else if (option == "--sectors")
        args.sectors = true;
      else if (option == "--precision" || option == "--tolerance" || option == "--max-iterations" || option == "--sz" ||
               option == "--quantum-numbers" || option == "--format")
      {
        if (++i == argc) throw std::invalid_argument("missing value for " + std::string(option));
        if (option == "--precision")
          args.precision = argv[i];
        else if (option == "--format")
          args.format = argv[i];
        else if (option == "--tolerance")
          args.tolerance = argv[i];
        else if (option == "--sz")
          args.sz = uni20::half_int::parse(argv[i]);
        else if (option == "--quantum-numbers")
          args.quantum_numbers = argv[i];
        else
          args.max_iterations = parse_size(argv[i]);
      }
      else
        throw std::invalid_argument("unknown option: " + std::string(option));
    }
    if (static_cast<int>(args.sz.has_value()) + static_cast<int>(args.quantum_numbers.has_value()) +
            static_cast<int>(args.sectors) >
        1)
      throw std::invalid_argument("--sz, --sectors and --quantum-numbers are mutually exclusive");
    if (args.format != "auto" && args.format != "pretty" && args.format != "plain")
      throw std::invalid_argument("unknown output format: " + std::string(args.format));
    if (args.precision == "fp64") return run<double>(args);
    if (args.precision == "long-double") return run<long double>(args);
    if (args.precision == "fp128")
    {
#if UNI20_HAS_FLOAT128
      return run<uni20::float128>(args);
#else
      throw std::invalid_argument("fp128 is unavailable; configure with -DUNI20_ENABLE_MPLAPACK=ON");
#endif
    }
    throw std::invalid_argument("unknown precision: " + std::string(args.precision));
  }
  catch (std::exception const& error)
  {
    std::cerr << "heisenberg-open-energy: " << error.what() << '\n';
    return 1;
  }
}
