// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "citation-report.hpp"
#include "hubbard-report.hpp"
#include <bethe/hubbard_open.hpp>

namespace
{
namespace model = bethe::hubbard::open;
struct Arguments
{
    std::size_t sites;
    std::optional<std::size_t> particles;
    std::optional<uni20::half_int> sz;
    std::optional<std::string_view> interaction, tolerance;
    std::string_view precision = "fp64", format = "auto";
    std::size_t max_iterations = 10000;
    bool roots = false;
};
void usage(std::ostream& out)
{
  out << "Usage: bethe-hubbard-obc L --u VALUE [options]\n"
      << "Free-end Hubbard sector ground state, t=1, either sign of U, L>=1.\n"
      << "H=-sum_(j=0..L-2,sigma)(c^dagger_(j,sigma)c_(j+1,sigma)+h.c.)+U sum_j n_up n_down.\n"
      << "No boundary fields; all physically valid particle/spin sectors are supported.\n"
      << "  --u VALUE                          required finite interaction\n"
      << "  --particles COUNT                  0 <= N <= 2L (default: L)\n"
      << "  --sz VALUE                         spin projection (default: 0 for even N, 1/2 for odd N)\n"
      << "  --precision fp64|long-double|fp128  (default: fp64)\n"
      << "  --tolerance VALUE                  max charge/spin residual divided by 2*(L+1)\n"
      << "                                     (default: 32 epsilon)\n"
      << "  --max-iterations COUNT             total Newton update budget, including continuation\n"
      << "  --roots                            print standing-wave k and spin rapidities Lambda\n"
      << "  --format auto|pretty|plain         terminal report or script output (default: auto)\n"
      << "  --help                             show this help\n"
      << "The interaction is U*n_up*n_down, not the particle-hole-shifted convention.\n"
      << "Odd and even lengths are supported. L=2 is the ordinary single-bond dimer.\n"
      << "There is no conserved lattice momentum; k labels standing waves.\n"
      << "Attractive U uses Shiba mapping; N>L uses particle-hole symmetry.\n"
      << "Mapped roots and residuals explicitly describe the auxiliary sector.\n"
      << "U=0 and single-species root sectors use exact free fermions.\n"
      << "Residuals use the final root-sector U and are not energy-error bounds.\n"
      << "Excitations and boundary fields are not implemented.\n"
      << "fp128 requires a Uni20 build with MPLAPACK enabled.\n";
  bethe::cli::print_citations(out, bethe::citations::Tool::hubbard_obc);
}
Arguments parse(int argc, char** argv)
{
  Arguments result{.sites = bethe::cli::parse_size(argv[1]),
                   .particles = std::nullopt,
                   .sz = std::nullopt,
                   .interaction = std::nullopt,
                   .tolerance = std::nullopt};
  for (int i = 2; i < argc; ++i)
  {
    std::string_view const option = argv[i];
    if (option == "--roots")
      result.roots = true;
    else if (option == "--u" || option == "--particles" || option == "--sz" || option == "--precision" ||
             option == "--format" || option == "--tolerance" || option == "--max-iterations")
    {
      if (++i == argc) throw std::invalid_argument("missing value for " + std::string(option));
      if (option == "--u")
        result.interaction = argv[i];
      else if (option == "--particles")
        result.particles = bethe::cli::parse_size(argv[i]);
      else if (option == "--sz")
        result.sz = uni20::half_int::parse(argv[i]);
      else if (option == "--precision")
        result.precision = argv[i];
      else if (option == "--format")
        result.format = argv[i];
      else if (option == "--tolerance")
        result.tolerance = argv[i];
      else
        result.max_iterations = bethe::cli::parse_size(argv[i]);
    }
    else
      throw std::invalid_argument("unknown option: " + std::string(option));
  }
  if (!result.interaction) throw std::invalid_argument("--u VALUE is required");
  if (result.format != "auto" && result.format != "pretty" && result.format != "plain")
    throw std::invalid_argument("unknown output format: " + std::string(result.format));
  return result;
}
template <uni20::Real Real> int run(Arguments const& args)
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
  return bethe::cli::print_hubbard_state(state, sz, args.precision, args.format, options.residual_tolerance, cpu_time,
                                         args.roots);
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
    auto const args = parse(argc, argv);
    return bethe::cli::dispatch_precision(args.precision, [&]<uni20::Real Real> { return run<Real>(args); });
  }
  catch (std::exception const& error)
  {
    std::cerr << "bethe-hubbard-obc: " << error.what() << '\n';
    return 1;
  }
}
