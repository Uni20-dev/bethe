// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include <bethe/heisenberg_open.hpp>

#include "heisenberg-cli.hpp"
#include "heisenberg-report.hpp"
#include "program-options.hpp"

#include <optional>

namespace
{
using bethe::cli::finish;
using bethe::cli::parse_quantum_numbers;
namespace model = bethe::heisenberg::open;

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
    bethe::cli::ExcitationArguments excitations = {};
    bethe::cli::DataOutputOptions output;
};

auto program_info()
{
  auto info = bethe::cli::program_info("bethe-xxx-obc", "Free-end spin-1/2 Heisenberg chain, J=1, zero field.",
                                       bethe::citations::Tool::xxx_obc);
  info.notes = {
      "Default: ground state (one representative for odd N).",
      "Excitations enumerate a restricted finite-real highest-weight family including the sector minimum, NOT a "
      "complete spectrum. "
      "The total spin defaults to 1 for even N, 1/2 for odd N. Strings and infinite roots are excluded.",
      "Positive rapidities; distinct integer labels in [1,N-M], M<=N/2. No boundary fields or lattice momentum.",
      "Use --references for literature and applicability; see docs/open-chains.md."};
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  namespace cli = bethe::cli;
  cli::count_option(app, "N", args.sites, "Number of sites")->required();
  cli::option(app, "--sz", args.sz, "Lowest energy in an Sz sector");
  cli::option(app, "--quantum-numbers", args.quantum_numbers, "Explicit finite-real labels; empty list for vacuum");
  cli::option(app, "--sectors", args.sectors, "Lowest energy in every Sz sector");
  cli::all_count_option(app, "--excitations", args.excitations.count,
                        "Lowest COUNT, or all, multiplets in the supported family");
  cli::option(app, "--spin", args.excitations.spin, "Total spin for the excitation scan")->needs("--excitations");
  cli::option(app, "--max-candidates", args.excitations.max_candidates, "Exhaustive scan limit (default: 10000)")
      ->needs("--excitations");
  cli::option(app, "--roots", args.print_roots, "Print rapidities and exact labels");
  cli::text_option(app, "--tolerance", args.tolerance, "Normalized equation residual in native precision");
  cli::count_option(app, "--max-iterations", args.max_iterations, "Update budget")->capture_default_str();
  cli::precision_option(app, args.precision);
  cli::add_data_output_options(app, args.output, true);
}

template <uni20::Real Real> int run(Arguments const& args, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  namespace cli = bethe::cli;
  namespace model = bethe::heisenberg::open;
  bethe::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  auto computation = context.computation();
  auto header = [&](std::string_view mode) {
    return cli::report_header(args.sites, args.precision, options, mode, context, false);
  };
  if (args.excitations.count)
  {
    auto const scan = model::real_excitations<Real>(args.sites, args.excitations.selected_spin(args.sites),
                                                    args.excitations.options(), options);
    computation.finish();
    return finish(cli::print_excitation_report(
        header("real-root excitations"), scan,
        {.family = "restricted real-root highest-weight multiplets; NOT a complete spectrum",
         .sector_label = "S",
         .sector = scan.spin,
         .multiplet_size = scan.spin.twice() + 1},
        args.print_roots, args.output));
  }
  if (args.sectors)
  {
    auto const states = model::sector_ground_states<Real>(args.sites, options);
    computation.finish();
    return finish(cli::spin_sectors<Real>(header("sector minima"), states, args.print_roots, args.output));
  }
  auto const state = args.quantum_numbers
                         ? model::solve_real<Real>(args.sites, parse_quantum_numbers(*args.quantum_numbers), options)
                     : args.sz ? model::sector_ground_state<Real>(args.sites, *args.sz, options)
                               : model::ground_state<Real>(args.sites, options);
  computation.finish();
  return finish(cli::spin_state<Real>(header(args.quantum_numbers ? "specified real-root state"
                                             : args.sz            ? "sector minimum"
                                                                  : "ground state"),
                                      args.sites, state, args.print_roots, args.output));
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
                static_cast<int>(args.sectors) + static_cast<int>(args.excitations.count.has_value()) >
            1)
          throw std::invalid_argument("--sz, --sectors, --quantum-numbers and --excitations are mutually exclusive");
        args.output.validate();
        return bethe::cli::dispatch_precision(args.precision,
                                              [&]<uni20::Real Real> { return run<Real>(args, argc, argv); });
      });
}
