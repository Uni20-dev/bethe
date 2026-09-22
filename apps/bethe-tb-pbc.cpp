// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "citation-report.hpp"
#include "report-common.hpp"
#include <bethe/takhtajan_babujian.hpp>

namespace
{
namespace model = bethe::takhtajan_babujian;
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
  out << "Usage: bethe-tb-pbc L [options]\n"
      << "Spin-1 Takhtajan-Babujian chain, periodic: H=sum_j [S_j.S_(j+1)-(S_j.S_(j+1))^2].\n"
      << "Zero-field singlet ground state; even L>=4. Bilinear coefficient is 1.\n"
      << "Complex two-string roots retain finite-size deviations, not ideal strings.\n"
      << "  --precision fp64|long-double|fp128  (default: fp64; fp128 requires MPLAPACK)\n"
      << "  --tolerance VALUE                  max phase/modulus residual divided by L\n"
      << "                                     (default: 32 epsilon; not an energy-error bound)\n"
      << "  --max-iterations COUNT             accepted Newton updates (default: 10000)\n"
      << "  --roots                            print string centers, deviations and complex roots\n"
      << "  --format auto|pretty|plain         (default: auto)\n"
      << "  --help                             show this help and references\n"
      << "lambda=x +/- i*(1/2+delta); E=-4 sum_j 1/(1+lambda_j^2); momentum 0.\n"
      << "Odd lengths, other spins/sectors, excitations and open ends are not implemented.\n"
      << "This is not the generic spin-1 Heisenberg chain or the SU(3) ULS point.\n"
      << "See docs/takhtajan-babujian.md for normalization and numerical conventions.\n";
  cli::print_citations(out, bethe::citations::Tool::tb_pbc);
}
Arguments parse(int argc, char** argv)
{
  Arguments result;
  result.sites = cli::parse_size(argv[1]);
  for (int i = 2; i < argc; ++i)
  {
    std::string_view const option = argv[i];
    if (option == "--roots")
    {
      result.roots = true;
      continue;
    }
    if (option != "--precision" && option != "--format" && option != "--tolerance" && option != "--max-iterations")
      throw std::invalid_argument("unknown option: " + std::string(option));
    if (++i == argc) throw std::invalid_argument("missing value for " + std::string(option));
    std::string_view const value = argv[i];
    if (option == "--precision")
      result.precision = value;
    else if (option == "--format")
      result.format = value;
    else if (option == "--tolerance")
      result.tolerance = value;
    else
      result.max_iterations = cli::parse_size(value);
  }
  if (result.format != "auto" && result.format != "plain" && result.format != "pretty")
    throw std::invalid_argument("unknown output format: " + std::string(result.format));
  return result;
}
char const* status(model::SolveStatus value)
{
  switch (value)
  {
    case model::SolveStatus::converged:
      return "converged";
    case model::SolveStatus::iteration_limit:
      return "iteration limit reached; unconverged estimate";
    case model::SolveStatus::stalled:
      return "line search or representable precision stalled; unconverged estimate";
    case model::SolveStatus::ill_conditioned:
      return "Newton system unresolved; unconverged estimate";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& args)
{
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  cli::CpuTimer const timer;
  auto const state = model::ground_state<Real>(args.sites, options);
  auto const cpu_time = timer.elapsed_text();
  cli::report_builder report("Spin-1 Takhtajan-Babujian chain (periodic)");
  report.status(state.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning, status(state.status))
      .field("Hamiltonian", "H=sum_j [S.S-(S.S)^2]; bilinear coefficient 1")
      .field("Calculation", "zero-field singlet ground state with finite string deviations")
      .field("Sites", state.sites)
      .field("Spin", 1)
      .field("Total spin", 0)
      .field("Two-strings", state.centers.size())
      .field("Complex roots", state.rapidities.size())
      .field("Precision", args.precision)
      .field("Residual tolerance", uni20::format_real(options.residual_tolerance))
      .field("Status", status(state.status))
      .field("Total energy", uni20::format_real(state.energy))
      .field("Energy per site", uni20::format_real(state.energy / Real(state.sites)))
      .field("Momentum index", state.momentum_index)
      .field("Momentum P", uni20::format_real(state.momentum))
      .field("Phase residual", uni20::format_real(state.phase_residual))
      .field("Modulus residual", uni20::format_real(state.modulus_residual))
      .field("Residual norm", uni20::format_real(state.residual_norm))
      .field("Iterations", state.iterations)
      .field("CPU time", cpu_time);
  if (args.roots)
  {
    auto& strings = report.table("Deviated two-strings");
    strings.header_separator()
        .column("Index")
        .column("I")
        .column("Center x", cli::table_alignment::decimal)
        .column("Deviation delta", cli::table_alignment::decimal);
    for (std::size_t j = 0; j < state.centers.size(); ++j)
      strings.row(j, uni20::to_string_fraction(state.string_quantum_numbers[j]), uni20::format_real(state.centers[j]),
                  uni20::format_real(state.deviations[j]));
    auto& roots = report.table("Complex rapidities");
    roots.header_separator()
        .column("Index")
        .column("Re lambda", cli::table_alignment::decimal)
        .column("Im lambda", cli::table_alignment::decimal);
    for (std::size_t j = 0; j < state.rapidities.size(); ++j)
      roots.row(j, uni20::format_real(state.rapidities[j].real()), uni20::format_real(state.rapidities[j].imag()));
  }
  cli::print_report(report, args.format);
  if (!state.converged) std::cerr << "TB solve incomplete; consider a larger budget or higher precision.\n";
  return state.converged ? 0 : 2;
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
    std::cerr << "bethe-tb-pbc: " << error.what() << '\n';
    return 1;
  }
}
