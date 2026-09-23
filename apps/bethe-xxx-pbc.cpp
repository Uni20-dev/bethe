// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include "heisenberg-cli.hpp"
#include "heisenberg-report.hpp"
#include "program-options.hpp"

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
    std::size_t sites = 0;
    std::string precision = "fp64";
    std::optional<std::string> tolerance = std::nullopt;
    std::size_t max_iterations = 10000;
    bool print_roots = false;
    std::optional<uni20::half_int> sz = std::nullopt;
    std::optional<std::string> quantum_numbers = std::nullopt;
    bool sectors = false;
    bool spinons = false;
    bethe::cli::ExcitationArguments excitations = {};
    std::string format = "auto";
};

auto program_info()
{
  auto info = bethe::cli::program_info("bethe-xxx-pbc", "Periodic spin-1/2 Heisenberg chain, J=1, zero field.",
                                       bethe::citations::Tool::xxx_pbc);
  info.notes = {"Default: ground state (one representative for odd N).",
                "Excitations enumerate a restricted finite-real highest-weight family including the sector minimum, "
                "NOT a complete spectrum. "
                "The total spin defaults to 1 for even N, 1/2 for odd N. Strings and infinite roots are excluded.",
                "--spinons selects the odd-N one-spinon branch at Sz=1/2. Lattice momentum is modulo 2*pi.",
                "Use --references for literature and applicability; see docs/xxx.md."};
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  namespace cli = bethe::cli;
  cli::count_option(app, "N", args.sites, "Number of sites")->required();
  cli::option(app, "--sz", args.sz, "Lowest energy in an Sz sector");
  cli::option(app, "--quantum-numbers", args.quantum_numbers, "Explicit finite-real labels; empty list for vacuum");
  cli::option(app, "--sectors", args.sectors, "Lowest energy in every Sz sector");
  cli::option(app, "--spinons", args.spinons, "Odd-N one-spinon branch (Sz=1/2)");
  cli::all_count_option(app, "--excitations", args.excitations.count,
                        "Lowest COUNT, or all, multiplets in the supported family");
  cli::option(app, "--spin", args.excitations.spin, "Total spin for the excitation scan")->needs("--excitations");
  cli::option(app, "--max-candidates", args.excitations.max_candidates, "Exhaustive scan limit (default: 10000)")
      ->needs("--excitations");
  cli::option(app, "--roots", args.print_roots, "Print rapidities and exact labels");
  cli::text_option(app, "--tolerance", args.tolerance, "Normalized equation residual in native precision");
  cli::count_option(app, "--max-iterations", args.max_iterations, "Update budget")->capture_default_str();
  cli::precision_option(app, args.precision);
  app.add_option("--format", args.format, "Stdout layout")
      ->check(CLI::IsMember({"auto", "pretty", "plain"}))
      ->capture_default_str();
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
  Arguments args;
  return bethe::cli::program_main(
      argc, argv, program_info(), [&](auto& app) { add_options(app, args); },
      [&](auto&) {
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
      });
}
