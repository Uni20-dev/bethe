// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include <bethe/xxz_excitations.hpp>

#include "excitation-report.hpp"
#include "program-options.hpp"
#include "xxz-cli.hpp"

namespace
{
namespace model = bethe::xxz;
using bethe::cli::finish;

using Arguments = bethe::cli::XxzArguments;

auto program_info()
{
  auto info = bethe::cli::program_info(
      "bethe-xxz-pbc", "Periodic spin-1/2 XXZ chain, J=1, zero field; ground states at finite Delta >= 0,",
      bethe::citations::Tool::xxz_pbc);
  info.notes = {"and -1 < Delta < 0 on even rings. Negative-Delta odd rings are not yet supported.",
                "H=sum_i (Sx_i Sx_(i+1) + Sy_i Sy_(i+1) + Delta Sz_i Sz_(i+1)).",
                "Default: ground state (one Sz=1/2 representative for odd N).",
                "z=tanh(lambda)/tan(gamma/2), Delta=cos(gamma); at Delta=1, z=2*lambda_XXX.",
                "For -1<Delta<0, --roots also prints lambda; residuals use rank-subtracted equations divided by N*s,",
                "where s=sqrt((1+Delta)/(1-Delta)). Ground states and --sz/--sectors only, even N.",
                "For Delta>1: z=tan(lambda)/tanh(eta/2), Delta=cosh(eta).",
                "Delta>1 supports ground states and --sz/--sectors only, not explicit labels or excitations.",
                "Excitation scans include the sector minimum; NOT a complete Sz spectrum.",
                "The finite-real window depends on Delta; strings and infinite rapidities are excluded.",
                "No OBC or SU(2) multiplet classification; the excitation window has one even-N Sz=0 state.",
                "Use --references for literature and applicability; see CITATIONS.md."};
  return info;
}

template <uni20::Real Real> int run(Arguments const& args, int argc, char** argv)
{
  namespace cli = bethe::cli;
  Real const delta = uni20::parse_real<Real>(*args.delta);
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  auto const coordinate = delta < Real{0} ? "scaled rapidity z=s*tanh(lambda); hyperbolic lambda" : "scaled rapidity z";
  constexpr auto negative_residual = "rank-subtracted equations divided by N*s; s=sqrt((1+Delta)/(1-Delta))";
  auto header = [&](std::string_view mode, std::string_view cpu_time) {
    bethe::cli::report_builder report("Heisenberg XXZ - " + std::string(mode));
    report.field("Model", "periodic spin-1/2, J=1, h=0")
        .field("Sites", args.sites)
        .field("Delta", uni20::format_real(delta))
        .field("Precision", args.precision)
        .field("Residual tolerance", uni20::format_real(options.residual_tolerance))
        .field("CPU time", cpu_time)
        .field("Root coordinate", coordinate);
    if (delta < Real{0}) report.field("Residual convention", negative_residual);
    return report;
  };

  cli::CpuTimer const timer;
  if (args.excitation_count)
  {
    auto const sz = args.sz.value_or(uni20::from_twice(std::int64_t{args.sites % 2 == 0 ? 2 : 1}));
    auto const scan = model::real_excitations<Real>(
        args.sites, delta, sz, {.count = *args.excitation_count, .max_candidates = args.max_candidates.value_or(10000)},
        options);
    auto const cpu_time = timer.elapsed_text();
    auto report = header("real-root excitations", cpu_time);
    auto const window = scan.window.slots == 0
                            ? "empty (polarized vacuum)"
                            : uni20::to_string(scan.window.first) + ".." + uni20::to_string(scan.window.last);
    report.field("Quantum-number window", window)
        .field("Available slots", scan.window.slots)
        .field("Spin-reversed reference", sz.twice() < 0 ? 1 : 0);
    return finish(
        cli::print_excitation_report(std::move(report), scan,
                                     {.family = "restricted finite-real XXZ states; NOT a complete Sz spectrum",
                                      .sector_label = "Sz",
                                      .sector = sz,
                                      .multiplet_size = std::nullopt},
                                     args.print_roots, args.output, "bethe-xxz-pbc", argc, argv));
  }
  if (args.sectors)
  {
    auto const states = model::sector_ground_states<Real>(args.sites, delta, options);
    auto const cpu_time = timer.elapsed_text();
    return finish(cli::spin_sectors<Real>(header("sector minima", cpu_time), states, args.print_roots, args.output,
                                          "bethe-xxz-pbc", argc, argv));
  }
  auto report_one = [&](auto const& state) {
    auto const cpu_time = timer.elapsed_text();
    return finish(cli::spin_state<Real>(header(args.quantum_numbers ? "specified real-root state"
                                               : args.sz            ? "sector minimum"
                                                                    : "ground state",
                                               cpu_time),
                                        args.sites, state, args.print_roots, args.output, "bethe-xxz-pbc", argc, argv));
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
        return bethe::cli::dispatch_precision(args.precision,
                                              [&]<uni20::Real Real> { return run<Real>(args, argc, argv); });
      });
}
