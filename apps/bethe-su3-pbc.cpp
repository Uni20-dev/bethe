// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "citation-report.hpp"
#include "report-common.hpp"
#include <bethe/su3.hpp>

namespace
{
namespace model = bethe::su3;
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
  out << "Usage: bethe-su3-pbc L [options]\n"
      << "Fundamental SU(3) permutation chain, periodic: H=sum_j P_(j,j+1), J=1.\n"
      << "Balanced singlet ground state only: L>=3 divisible by 3.\n"
      << "Each color has L/3 sites; nested real-root counts are M1=2L/3 and M2=L/3.\n"
      << "  --precision fp64|long-double|fp128  (default: fp64; fp128 requires MPLAPACK)\n"
      << "  --tolerance VALUE                  max equation residual divided by L\n"
      << "                                     (default: 32 epsilon; not an energy-error bound)\n"
      << "  --max-iterations COUNT             accepted Newton updates (default: 10000)\n"
      << "  --roots                            print both root families and exact I,J labels\n"
      << "  --format auto|pretty|plain         (default: auto)\n"
      << "  --help                             show this help and references\n"
      << "E=L-sum_j 1/(lambda_j^2+1/4). Momentum is 0 in the supported singlet.\n"
      << "Roots use conventional lambda,mu, not the XXX front end's z=2*lambda.\n"
      << "Other color sectors, excitations, complex strings, twists and open ends\n"
      << "are not implemented. See docs/su3.md, including the spin-1 ULS mapping.\n";
  cli::print_citations(out, bethe::citations::Tool::su3_pbc);
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
  cli::report_builder report("SU(3) permutation chain (periodic)");
  report.status(state.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning, status(state.status))
      .field("Hamiltonian", "H=sum_j P_(j,j+1); J=1")
      .field("Calculation", "balanced singlet ground state")
      .field("Sites", state.sites)
      .field("Color populations",
             fmt::format("{}, {}, {}", state.populations[0], state.populations[1], state.populations[2]))
      .field("First-level roots", state.rapidities[0].size())
      .field("Second-level roots", state.rapidities[1].size())
      .field("Precision", args.precision)
      .field("Residual tolerance", uni20::format_real(options.residual_tolerance))
      .field("Status", status(state.status))
      .field("Total energy", uni20::format_real(state.energy))
      .field("Energy per site", uni20::format_real(state.energy / Real(state.sites)))
      .field("Momentum index", state.momentum_index)
      .field("Momentum P", uni20::format_real(state.momentum))
      .field("First-level residual", uni20::format_real(state.level_residuals[0]))
      .field("Second-level residual", uni20::format_real(state.level_residuals[1]))
      .field("Residual norm", uni20::format_real(state.residual_norm))
      .field("Iterations", state.iterations)
      .field("CPU time", cpu_time);
  if (args.roots)
  {
    for (std::size_t a = 0; a < 2; ++a)
    {
      auto& table = report.table(a == 0 ? "First-level rapidities" : "Second-level rapidities");
      table.header_separator()
          .column("Index")
          .column(a == 0 ? "I" : "J")
          .column(a == 0 ? "lambda" : "mu", cli::table_alignment::decimal);
      for (std::size_t j = 0; j < state.rapidities[a].size(); ++j)
        table.row(j, uni20::to_string_fraction(state.quantum_numbers[a][j]),
                  uni20::format_real(state.rapidities[a][j]));
    }
  }
  cli::print_report(report, args.format);
  if (!state.converged) std::cerr << "SU(3) solve incomplete; consider a larger budget or higher precision.\n";
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
    std::cerr << "bethe-su3-pbc: " << error.what() << '\n';
    return 1;
  }
}
