// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include "heisenberg-cli.hpp"
#include "heisenberg-report.hpp"

#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace
{
using bethe::cli::CpuTimer;
using bethe::cli::finish;
using bethe::cli::parse_quantum_numbers;
using bethe::cli::parse_size;
using bethe::cli::print_roots;

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
    bool spinons = false;
    bethe::cli::ExcitationArguments excitations = {};
    std::string_view format = "auto";
};

void usage(std::ostream& out)
{
  out << "Usage: bethe-xxx-pbc N [options]\n"
      << "Periodic spin-1/2 Heisenberg chain, J=1, zero field.\n"
      << "Default: ground state (one representative for odd N). Modes:\n"
      << "  --sz VALUE                         lowest energy in an Sz sector, e.g. 1/2\n"
      << "  --sectors                          lowest energy in every Sz sector\n"
      << "  --spinons                          odd-N one-spinon branch (Sz=1/2)\n"
      << "  --quantum-numbers I0,I1,...         specified finite real-root state\n";
  bethe::cli::excitation_usage(out);
  out << "  --precision fp64|long-double|fp128  (default: fp64)\n"
      << "  --tolerance VALUE                  normalized equation residual\n"
      << "  --max-iterations COUNT             update budget (default: 10000)\n"
      << "  --roots                            print the rapidities\n"
      << "  --format auto|pretty|plain         terminal report or script output (default: auto)\n"
      << "  --help                             show this help\n"
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
              << "Precision: " << args.precision << '\n'
              << "Residual tolerance: " << uni20::format_real(options.residual_tolerance) << '\n';
  CpuTimer const timer;
  if (args.excitations.count)
  {
    auto const scan = bethe::heisenberg::real_excitations<Real>(args.sites, args.excitations.selected_spin(args.sites),
                                                                args.excitations.options(), options);
    auto const cpu_time = timer.elapsed_text();
    return finish(
        bethe::cli::print_excitations(args.sites, args.precision, options, scan, args.print_roots, cpu_time, pretty));
  }
  if (args.sectors)
  {
    auto const states = bethe::heisenberg::sector_ground_states<Real>(args.sites, options);
    auto const cpu_time = timer.elapsed_text();
    if (pretty)
      return finish(bethe::cli::print_sectors(args.sites, args.precision, options, states, args.print_roots, cpu_time));
    std::cout << "# CPU time: " << cpu_time << '\n' << "# Sz momentum_index P energy residual iterations converged\n";
    bool converged = true;
    for (auto const& state : states)
    {
      std::cout << state.sz << ' ' << state.momentum_index << ' ' << uni20::format_real(state.momentum) << ' '
                << uni20::format_real(state.energy) << ' ' << uni20::format_real(state.residual_norm) << ' '
                << state.iterations << ' ' << state.converged << '\n';
      if (args.print_roots) print_roots(state);
      converged = converged && state.converged;
    }
    return finish(converged);
  }
  if (args.spinons)
  {
    auto const branch = bethe::heisenberg::one_spinon_branch<Real>(args.sites, options);
    auto const cpu_time = timer.elapsed_text();
    if (pretty)
      return finish(bethe::cli::print_spinons(args.sites, args.precision, options, branch, args.print_roots, cpu_time));
    std::cout << "# CPU time: " << cpu_time << '\n'
              << "# Sz=1/2; k=pi/2-2*pi*I_h/N; bulk reference e_inf=1/4-log(2)\n"
              << "# hole k momentum_index P energy E_minus_N_e_inf epsilon_inf residual iterations converged\n";
    bool converged = true;
    for (auto const& point : branch)
    {
      auto const& state = point.state;
      std::cout << point.hole << ' ' << uni20::format_real(point.spinon_momentum) << ' ' << state.momentum_index << ' '
                << uni20::format_real(state.momentum) << ' ' << uni20::format_real(state.energy) << ' '
                << uni20::format_real(point.bulk_subtracted_energy) << ' '
                << uni20::format_real(bethe::heisenberg::spinon_energy(point.spinon_momentum)) << ' '
                << uni20::format_real(state.residual_norm) << ' ' << state.iterations << ' ' << state.converged << '\n';
      if (args.print_roots) print_roots(state);
      converged = converged && state.converged;
    }
    return finish(converged);
  }
  auto const result =
      args.quantum_numbers
          ? bethe::heisenberg::solve_real<Real>(args.sites, parse_quantum_numbers(*args.quantum_numbers), options)
      : args.sz ? bethe::heisenberg::sector_ground_state<Real>(args.sites, *args.sz, options)
                : bethe::heisenberg::ground_state<Real>(args.sites, options);
  auto const cpu_time = timer.elapsed_text();
  if (pretty)
    return finish(bethe::cli::print_state(args.sites, args.precision, options, result,
                                          args.quantum_numbers ? "specified real-root state"
                                          : args.sz            ? "sector minimum"
                                                               : "ground state",
                                          args.print_roots, cpu_time));
  std::cout << "Sz: " << result.sz << '\n'
            << "Spin-reversed reference: " << result.spin_reversed << '\n'
            << "Momentum index: " << result.momentum_index << '\n'
            << "Momentum: " << uni20::format_real(result.momentum) << '\n'
            << "Status: " << (result.converged ? "converged" : "iteration limit reached") << '\n'
            << "Iterations: " << result.iterations << '\n'
            << "CPU time: " << cpu_time << '\n'
            << "Residual norm: " << uni20::format_real(result.residual_norm) << '\n'
            << "Total energy: " << uni20::format_real(result.energy) << '\n'
            << "Energy per site: " << uni20::format_real(result.energy / static_cast<Real>(args.sites)) << '\n';
  if (args.print_roots) print_roots(result);
  return finish(result.converged);
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
      else if (option == "--spinons")
        args.spinons = true;
      else if (option == "--precision" || option == "--tolerance" || option == "--max-iterations" || option == "--sz" ||
               option == "--quantum-numbers" || option == "--format" || option == "--excitations" ||
               option == "--spin" || option == "--max-candidates")
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
        else if (!args.excitations.parse(option, argv[i]))
          args.max_iterations = parse_size(argv[i]);
      }
      else
        throw std::invalid_argument("unknown option: " + std::string(option));
    }
    args.excitations.validate();
    if (static_cast<int>(args.sz.has_value()) + static_cast<int>(args.quantum_numbers.has_value()) +
            static_cast<int>(args.sectors) + static_cast<int>(args.spinons) +
            static_cast<int>(args.excitations.count.has_value()) >
        1)
      throw std::invalid_argument(
          "--sz, --sectors, --spinons, --quantum-numbers and --excitations are mutually exclusive");
    if (args.format != "auto" && args.format != "pretty" && args.format != "plain")
      throw std::invalid_argument("unknown output format: " + std::string(args.format));
    return bethe::cli::dispatch_precision(args.precision, [&]<uni20::Real Real> { return run<Real>(args); });
  }
  catch (std::exception const& error)
  {
    std::cerr << "bethe-xxx-pbc: " << error.what() << '\n';
    return 1;
  }
}
