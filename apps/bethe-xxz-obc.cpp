// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include <bethe/xxz_excitations.hpp>

#include "excitation-report.hpp"
#include "program-options.hpp"
#include "xxz-cli.hpp"
#include "xxz-open-report.hpp"

namespace
{
namespace model = bethe::xxz::open;
namespace output = bethe::cli::xxz_open_report;
using bethe::cli::parse_size;
using bethe::cli::print_roots;
using output::finish;
using output::print_roots;

using Arguments = bethe::cli::XxzArguments;

auto program_info()
{
  auto info = bethe::cli::program_info("bethe-xxz-obc",
                                       "Open spin-1/2 XXZ chain, free ends, J=1, zero field; Delta > -1 ground states.",
                                       bethe::citations::Tool::xxz_obc);
  info.notes = {"H=sum_(i=0)^(N-2) (Sx_i Sx_(i+1) + Sy_i Sy_(i+1) + Delta Sz_i Sz_(i+1)).",
                "Default: ground state (one Sz=1/2 representative for odd N).",
                "z=tanh(lambda)/tan(gamma/2), Delta=cos(gamma); at Delta=1, z=2*lambda_XXX.",
                "For -1<Delta<0, --roots also prints lambda; residuals use rank-subtracted equations divided by N*s,",
                "where s=sqrt((1+Delta)/(1-Delta)); this avoids false convergence near Delta=-1.",
                "For Delta>1, bulk z=tan(lambda)/tanh(eta/2), Delta=cosh(eta).",
                "Even zero-Sz massive ground states also carry a boundary root as y=1/z_B^2 and log distance w.",
                "--excitations and --quantum-numbers require 0 <= Delta <= 1.",
                "Excitation scans include the sector minimum; NOT a complete Sz spectrum.",
                "The finite-real window depends on Delta; strings and infinite rapidities are excluded.",
                "No boundary fields, lattice momentum, or SU(2) multiplet classification.",
                "Use --references for literature and applicability; see CITATIONS.md."};
  return info;
}

template <uni20::Real Real> int run(Arguments const& args)
{
  Real const delta = uni20::parse_real<Real>(*args.delta);
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  bool const pretty = args.format == "pretty" || (args.format == "auto" && terminal::is_a_terminal(stdout));
  auto const coordinate = delta > Real{1}   ? "bulk scaled rapidity z; boundary y=1/z_B^2, log distance w"
                          : delta < Real{0} ? "scaled rapidity z=s*tanh(lambda); hyperbolic lambda"
                                            : "scaled rapidity z";
  auto const convention = delta < Real{0}   ? model::GroundResidualConvention::negative_rank_scaled
                          : delta > Real{1} ? model::GroundResidualConvention::massive_regularized
                                            : model::GroundResidualConvention::logarithmic_phase;
  if (!pretty)
    std::cout << "Sites: " << args.sites << '\n'
              << "Model: open XXZ, free ends, J=1, h=0\n"
              << "Delta: " << uni20::format_real(delta) << '\n'
              << "Precision: " << args.precision << '\n'
              << "Residual tolerance: " << uni20::format_real(options.residual_tolerance) << '\n'
              << "Root coordinate: " << coordinate << '\n';
  if (!pretty && (delta > Real{1} || delta < Real{0}))
    std::cout << "Residual convention: " << output::residual_description(convention) << '\n';
  auto header = [&](std::string_view mode, std::string_view cpu_time) {
    bethe::cli::report_builder report("Heisenberg XXZ (free ends) - " + std::string(mode));
    report.field("Model", "open spin-1/2, free ends, J=1, h=0")
        .field("Sites", args.sites)
        .field("Delta", uni20::format_real(delta))
        .field("Precision", args.precision)
        .field("Residual tolerance", uni20::format_real(options.residual_tolerance))
        .field("CPU time", cpu_time)
        .field("Root coordinate", coordinate);
    if (delta > Real{1} || delta < Real{0})
      report.field("Residual convention", output::residual_description(convention));
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
    if (pretty) return finish(output::print_sectors(header("sector minima", cpu_time), states, args.print_roots));
    std::cout << "# CPU time: " << cpu_time << '\n' << "# Sz energy residual iterations converged";
    if (delta > Real{1} || delta < Real{0}) std::cout << " root_delta status";
    std::cout << '\n';
    bool converged = true;
    for (auto const& state : states)
    {
      std::cout << state.sz << ' ' << uni20::format_real(state.energy) << ' ' << uni20::format_real(state.residual_norm)
                << ' ' << state.iterations << ' ' << state.converged;
      if (delta > Real{1} || delta < Real{0})
        std::cout << ' ' << uni20::format_real(state.root_delta) << ' ' << output::status_code(state.status);
      std::cout << '\n';
      if (args.print_roots)
      {
        std::cout << "# Roots: Sz=" << state.sz << '\n';
        print_roots(state);
      }
      converged = converged && state.converged;
    }
    return finish(converged);
  }
  auto report_one = [&](auto const& state) {
    auto const cpu_time = timer.elapsed_text();
    if (pretty)
      return finish(output::print_state<Real>(header(args.quantum_numbers ? "specified real-root state"
                                                     : args.sz            ? "sector minimum"
                                                                          : "ground state",
                                                     cpu_time),
                                              args.sites, state, args.print_roots));
    std::cout << "Sz: " << state.sz << '\n'
              << "Spin-reversed reference: " << state.spin_reversed << '\n'
              << "Status: " << output::status_text(state) << '\n'
              << "Iterations: " << state.iterations << '\n'
              << "CPU time: " << cpu_time << '\n'
              << "Residual norm: " << uni20::format_real(state.residual_norm) << '\n'
              << "Total energy: " << uni20::format_real(state.energy) << '\n'
              << "Energy per site: " << uni20::format_real(state.energy / static_cast<Real>(args.sites)) << '\n';
    if constexpr (requires { state.boundary_root; })
      if (delta > Real{1})
        std::cout << "Root Delta: " << uni20::format_real(state.root_delta) << '\n'
                  << "Boundary root: " << (state.boundary_root ? output::boundary_kind(*state.boundary_root) : "none")
                  << '\n';
    if (args.print_roots) print_roots(state);
    return finish(state.converged);
  };
  if (args.quantum_numbers)
    return report_one(model::solve_real<Real>(args.sites, delta, *args.quantum_numbers, options));
  return report_one(args.sz ? model::sector_ground_state<Real>(args.sites, delta, *args.sz, options)
                            : model::ground_state<Real>(args.sites, delta, options));
}
} // namespace

int main(int argc, char** argv)
{
  Arguments args;
  return bethe::cli::program_main(
      argc, argv, program_info(), [&](auto& app) { bethe::cli::add_xxz_options(app, args); },
      [&](auto&) {
        bethe::cli::validate_xxz_arguments(args);
        return bethe::cli::dispatch_precision(args.precision, [&]<uni20::Real Real> { return run<Real>(args); });
      });
}
