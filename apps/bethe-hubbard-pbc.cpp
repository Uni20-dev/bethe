// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "hubbard-report.hpp"
#include "program-options.hpp"
#include <bethe/hubbard.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::hubbard;
struct Arguments
{
    std::size_t sites = 0;
    std::optional<std::size_t> particles;
    uni20::half_int sz{0};
    std::optional<std::string> interaction, tolerance;
    std::string precision = "fp64";
    cli::DataOutputOptions output;
    std::size_t max_iterations = 10000;
    bool roots = false;
};
auto program_info()
{
  auto info =
      bethe::cli::program_info("bethe-hubbard-pbc", "Periodic Hubbard sector ground state, t=1, either sign of U.",
                               bethe::citations::Tool::hubbard_pbc);
  info.notes = {"H=-sum_(j,sigma)(c^dagger_(j,sigma)c_(j+1,sigma)+h.c.)+U sum_j n_up n_down.",
                "Even L>=2; supported sector families are listed below.",
                "The interaction is U*n_up*n_down, not the particle-hole-shifted convention.",
                "L=2 counts the periodic hopping bond twice. U=0 uses exact free fermions.",
                "Repulsive root sectors: half filling with any Sz, doped odd N_up and N_down,",
                "or a single spin species. Above half filling uses particle-hole symmetry.",
                "Attractive U uses the Shiba mapping; all balanced even-N sectors are supported.",
                "Other sectors work only if their mapped repulsive sector is supported. U=0 is unrestricted.",
                "Other interacting shell parities, excitations and odd rings are not implemented here.",
                "For free ends and unrestricted sectors, use bethe-hubbard-obc.",
                "Mapped roots and residuals explicitly describe the auxiliary sector, not attractive roots.",
                "Residuals use the final root-sector U and are not energy-error bounds.",
                "Use --references for literature and applicability; see CITATIONS.md."};
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  bethe::cli::option(app, "L", args.sites, "Number of sites or rungs")->required();
  bethe::cli::option(app, "--u", args.interaction, "required finite interaction")->required();
  bethe::cli::option(app, "--particles", args.particles, "0 <= N <= 2L (default: L)");
  bethe::cli::option(app, "--sz", args.sz, "integer or half-integer spin projection (default: 0)");
  bethe::cli::option(app, "--roots", args.roots, "print charge momenta k and spin rapidities Lambda");
  bethe::cli::option(app, "--tolerance", args.tolerance, "max normalized charge/spin residual (default: 32 epsilon)");
  bethe::cli::option(app, "--max-iterations", args.max_iterations, "total Newton update budget, including continuation")
      ->capture_default_str();
  bethe::cli::precision_option(app, args.precision);
  cli::add_data_output_options(app, args.output, true);
}

void validate(Arguments const& args)
{
  if (!args.interaction) throw std::invalid_argument("--u VALUE is required");
  args.output.validate();
}
template <uni20::Real Real> int run(Arguments const& args, int argc, char** argv)
{
  Real const interaction = uni20::parse_real<Real>(*args.interaction);
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  bethe::cli::CpuTimer const timer;
  auto const state =
      model::sector_ground_state<Real>(args.sites, args.particles.value_or(args.sites), args.sz, interaction, options);
  auto const cpu_time = timer.elapsed_text();
  return bethe::cli::print_hubbard_state(state, args.sz, args.precision, args.output, options.residual_tolerance,
                                         cpu_time, args.roots, argc, argv);
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
