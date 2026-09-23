// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "citation-report.hpp"
#include "report-common.hpp"
#include <bethe/biquadratic.hpp>

namespace
{
namespace cli = bethe::cli;
struct Arguments
{
    std::size_t sites = 0, max_iterations = 10000;
    std::optional<std::string_view> tolerance;
    std::string_view precision = "fp64", format = "auto";
    bool roots = false;
};
void usage(std::ostream& out)
{
  out << "Usage: bethe-biquadratic-obc N [options]\n"
      << "Spin-1 pure biquadratic chain, free ends: H=-sum_i (S_i.S_(i+1))^2.\n"
      << "Unique singlet ground state; even N>=2, coefficient -1.\n"
      << "TL loop weight 3; reference XXZ Delta=3/2 with opposite end fields.\n"
      << "  --precision fp64|long-double|fp128  (default: fp64; fp128 requires MPLAPACK)\n"
      << "  --tolerance VALUE                  max logarithmic residual divided by 2*N\n"
      << "                                     (default: 32 epsilon; not an energy-error bound)\n"
      << "  --max-iterations COUNT             accepted Newton updates (default: 10000)\n"
      << "  --roots                            print reference XXZ rapidities and labels\n"
      << "  --format auto|pretty|plain         (default: auto)\n"
      << "  --help                             show this help and references\n"
      << "No odd chains, excited states, physical-spin sector scans or lattice momentum.\n"
      << "This is not the TB point, ULS point, or zero-boundary-field XXZ chain.\n"
      << "See docs/biquadratic.md for the TL mapping and representation multiplicities.\n";
  cli::print_citations(out, bethe::citations::Tool::biquadratic_obc);
}
Arguments parse(int argc, char** argv)
{
  Arguments args;
  args.sites = cli::parse_size(argv[1]);
  for (int i = 2; i < argc; ++i)
  {
    std::string_view const option = argv[i];
    if (option == "--roots")
    {
      args.roots = true;
      continue;
    }
    if (option != "--precision" && option != "--format" && option != "--tolerance" && option != "--max-iterations")
      throw std::invalid_argument("unknown option: " + std::string(option));
    if (++i == argc) throw std::invalid_argument("missing value for " + std::string(option));
    std::string_view const value = argv[i];
    if (option == "--precision")
      args.precision = value;
    else if (option == "--format")
      args.format = value;
    else if (option == "--tolerance")
      args.tolerance = value;
    else
      args.max_iterations = cli::parse_size(value);
  }
  if (args.format != "auto" && args.format != "plain" && args.format != "pretty")
    throw std::invalid_argument("unknown output format: " + std::string(args.format));
  return args;
}
char const* status(bethe::xxz::quantum_group::SolveStatus value)
{
  using Status = bethe::xxz::quantum_group::SolveStatus;
  switch (value)
  {
    case Status::converged:
      return "converged";
    case Status::iteration_limit:
      return "iteration limit reached; unconverged estimate";
    case Status::singular_jacobian:
      return "singular or ill-conditioned Jacobian; unconverged estimate";
    case Status::stalled:
      return "line search or representable precision stalled; unconverged estimate";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& args)
{
  bethe::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  cli::CpuTimer const timer;
  auto const state = bethe::biquadratic::ground_state<Real>(args.sites, options);
  auto const cpu_time = timer.elapsed_text();
  auto const& reference = state.reference;
  cli::report_builder report("Spin-1 pure biquadratic chain (free ends)");
  report
      .status(reference.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning,
              status(reference.status))
      .field("Hamiltonian", "H=-sum_i (S_i.S_(i+1))^2")
      .field("Calculation", "even-chain singlet ground state")
      .field("Sites", args.sites)
      .field("Spin", 1)
      .field("Total spin", 0)
      .field("TL loop weight", 3)
      .field("TL through-lines", 0)
      .field("Ground-state multiplicity", *bethe::temperley_lieb::spin_chain_multiplicity(3, 0))
      .field("XXZ Delta", "1.5")
      .field("XXZ reference", "spin-half exchange 1; +sqrt(5)/4*(sz_1-sz_N)")
      .field("Precision", args.precision)
      .field("Residual tolerance", uni20::format_real(options.residual_tolerance))
      .field("Status", status(reference.status))
      .field("Total energy", uni20::format_real(state.energy))
      .field("Energy per site", uni20::format_real(state.energy / Real(args.sites)))
      .field("TL energy (-sum e_i)", uni20::format_real(state.tl_energy))
      .field("XXZ reference energy", uni20::format_real(reference.energy))
      .field("Residual norm", uni20::format_real(reference.residual_norm))
      .field("Iterations", reference.iterations)
      .field("CPU time", cpu_time);
  if (args.roots)
  {
    auto& table = report.table("Reference XXZ roots (not physical spin-1 quantum numbers)");
    table.header_separator()
        .column("Index")
        .column("I")
        .column("alpha", cli::table_alignment::decimal)
        .column("x=Theta_1/2", cli::table_alignment::decimal);
    for (std::size_t j = 0; j < reference.rapidities.size(); ++j)
      table.row(j, uni20::to_string(reference.quantum_numbers[j]), uni20::format_real(reference.rapidities[j]),
                uni20::format_real(reference.angles[j]));
  }
  cli::print_report(report, args.format);
  std::cout.flush();
  if (!std::cout) throw std::ios_base::failure("output flush failed");
  if (!reference.converged)
    std::cerr << "Biquadratic solve incomplete; consider a larger budget or higher precision.\n";
  return reference.converged ? 0 : 2;
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
    return cli::dispatch_precision(args.precision, [&]<uni20::Real Real> { return run<Real>(args); });
  }
  catch (std::exception const& error)
  {
    std::cerr << "bethe-biquadratic-obc: " << error.what() << '\n';
    return 1;
  }
}
