// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "hubbard-report.hpp"
#include "program-options.hpp"
#include <bethe/hubbard_open.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::hubbard::open;
struct Arguments
{
    std::size_t sites = 0;
    std::optional<std::size_t> particles;
    std::optional<uni20::half_int> sz;
    std::optional<std::string> interaction, tolerance;
    std::string precision = "fp64";
    cli::DataOutputOptions output;
    std::size_t max_iterations = 10000;
    bool roots = false;
};
auto program_info()
{
  auto info = bethe::cli::program_info("bethe-hubbard-obc",
                                       "Free-end Hubbard sector ground state, t=1, either sign of U, L>=1.",
                                       bethe::citations::Tool::hubbard_obc);
  info.notes = {"H=-sum_(j=0..L-2,sigma)(c^dagger_(j,sigma)c_(j+1,sigma)+h.c.)+U sum_j n_up n_down.",
                "No boundary fields; all physically valid particle/spin sectors are supported.",
                "The interaction is U*n_up*n_down, not the particle-hole-shifted convention.",
                "Odd and even lengths are supported. L=2 is the ordinary single-bond dimer.",
                "There is no conserved lattice momentum; k labels standing waves.",
                "Attractive U uses Shiba mapping; N>L uses particle-hole symmetry.",
                "Mapped roots and residuals explicitly describe the auxiliary sector.",
                "U=0 and single-species root sectors use exact free fermions.",
                "Residuals use the final root-sector U and are not energy-error bounds.",
                "Excitations and boundary fields are not implemented.",
                "Use --references for literature and applicability; see CITATIONS.md."};
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  bethe::cli::option(app, "L", args.sites, "Number of sites or rungs")->required();
  bethe::cli::option(app, "--u", args.interaction, "required finite interaction")->required();
  bethe::cli::option(app, "--particles", args.particles, "0 <= N <= 2L (default: L)");
  bethe::cli::option(app, "--sz", args.sz, "spin projection (default: 0 for even N, 1/2 for odd N)");
  bethe::cli::option(app, "--roots", args.roots, "print standing-wave k and spin rapidities Lambda");
  bethe::cli::option(app, "--tolerance", args.tolerance,
                     "max charge/spin residual divided by 2*(L+1) (default: 32 epsilon)");
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
  auto const particles = args.particles.value_or(args.sites);
  auto const sz = args.sz.value_or(uni20::from_twice(static_cast<std::int64_t>(particles % 2)));
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  bethe::cli::CpuTimer const timer;
  auto const state = model::sector_ground_state<Real>(args.sites, particles, sz, interaction, options);
  auto const cpu_time = timer.elapsed_text();
  return bethe::cli::print_hubbard_state(state, sz, args.precision, args.output, options.residual_tolerance, cpu_time,
                                         args.roots, argc, argv);
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
