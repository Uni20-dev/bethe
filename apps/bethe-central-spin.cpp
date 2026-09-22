// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "citation-report.hpp"
#include "report-common.hpp"
#include <bethe/central_spin.hpp>

namespace
{
namespace model = bethe::central_spin;
namespace cli = bethe::cli;
struct Arguments
{
    std::optional<std::string_view> couplings, field, tolerance;
    std::optional<uni20::half_int> sz;
    std::string_view precision = "fp64", format = "auto";
    std::size_t max_iterations = 10000, max_stages = 10000;
    bool variables = false;
};
void usage(std::ostream& out)
{
  out << "Usage: bethe-central-spin --couplings A1,A2,... --field B --sz SZ [options]\n"
      << "Rational Gaudin central-spin ground state in a specified total-Sz sector.\n"
      << "H=B*S0^z+sum_j A_j*S0.Sj; central and bath spins are all 1/2.\n"
      << "  --couplings LIST                   distinct nonzero real bath couplings (either sign)\n"
      << "  --field VALUE                      central field, either sign or exactly zero\n"
      << "  --sz VALUE                         required total magnetization, including the central spin\n"
      << "  --variables                        print compactified eigenvalue variables (not occupations)\n"
      << "  --precision fp64|long-double|fp128  (default: fp64; fp128 requires MPLAPACK)\n"
      << "  --tolerance VALUE                  polynomial backward residual (default: 32 epsilon)\n"
      << "  --max-iterations COUNT             attempted Newton corrections, including retries\n"
      << "  --max-stages COUNT                 attempted continuation stages (default: 10000)\n"
      << "  --format auto|pretty|plain         (default: auto)\n"
      << "  --help                             show this help and references\n"
      << "Default Newton budget: 10000. A zero budget leaves the infinite-field seed, with no energy.\n"
      << "Incomplete solves report energy only at the REACHED field, not the requested field.\n"
      << "An empty coupling list describes an isolated central spin. Coupling order is immaterial.\n"
      << "Repeated/zero couplings, bath interactions, excitations and rapidity reconstruction are not implemented.\n"
      << "No PBC/OBC or lattice momentum applies; a sector minimum is not necessarily the global minimum.\n"
      << "See docs/central-spin.md for conventions, state selection and continuation controls.\n";
  cli::print_citations(out, bethe::citations::Tool::central_spin);
}
Arguments parse(int argc, char** argv)
{
  Arguments out;
  for (int i = 1; i < argc; ++i)
  {
    std::string_view const option = argv[i];
    if (option == "--variables")
    {
      out.variables = true;
      continue;
    }
    if (option != "--couplings" && option != "--field" && option != "--sz" && option != "--precision" &&
        option != "--format" && option != "--tolerance" && option != "--max-iterations" && option != "--max-stages")
      throw std::invalid_argument("unknown option: " + std::string(option));
    if (++i == argc) throw std::invalid_argument("missing value for " + std::string(option));
    std::string_view const value = argv[i];
    if (option == "--couplings")
      out.couplings = value;
    else if (option == "--field")
      out.field = value;
    else if (option == "--sz")
      out.sz = uni20::half_int::parse(value);
    else if (option == "--precision")
      out.precision = value;
    else if (option == "--format")
      out.format = value;
    else if (option == "--tolerance")
      out.tolerance = value;
    else if (option == "--max-iterations")
      out.max_iterations = cli::parse_size(value);
    else
      out.max_stages = cli::parse_size(value);
  }
  if (!out.couplings || !out.field || !out.sz)
    throw std::invalid_argument("--couplings, --field and --sz are required");
  if (out.format != "auto" && out.format != "pretty" && out.format != "plain")
    throw std::invalid_argument("unknown output format: " + std::string(out.format));
  return out;
}
char const* status(model::SolveStatus value)
{
  switch (value)
  {
    case model::SolveStatus::converged:
      return "converged";
    case model::SolveStatus::iteration_limit:
      return "Newton budget exhausted; target field not reached";
    case model::SolveStatus::stage_limit:
      return "stage budget exhausted; target field not reached";
    case model::SolveStatus::stalled:
      return "continuation stalled; target field not reached";
    case model::SolveStatus::ill_conditioned:
      return "linear system unresolved; target field not reached";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& args)
{
  std::vector<Real> a;
  auto list = *args.couplings;
  while (!list.empty())
  {
    auto const comma = list.find(',');
    if (list.substr(0, comma).empty()) throw std::invalid_argument("empty coupling list item");
    a.push_back(uni20::parse_real<Real>(list.substr(0, comma)));
    if (comma == std::string_view::npos) break;
    list.remove_prefix(comma + 1);
    if (list.empty()) throw std::invalid_argument("empty coupling list item");
  }
  Real const b = uni20::parse_real<Real>(*args.field);
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  options.max_stages = args.max_stages;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  cli::CpuTimer const timer;
  auto const state = model::sector_ground_state<Real>(a, b, *args.sz, options);
  auto const cpu_time = timer.elapsed_text();
  cli::report_builder report("Rational Gaudin central spin");
  report.status(state.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning, status(state.status))
      .field("Hamiltonian", "B*S0^z+sum_j A_j*S0.Sj (all spins 1/2)")
      .field("Calculation", "lowest state in the specified total-Sz sector")
      .field("Bath spins", a.size())
      .field("Total Sz", uni20::to_string_fraction(state.sz))
      .field("Up spins", state.up_spins)
      .field("Precision", args.precision)
      .field("Spin reversed", state.spin_reversed ? "yes" : "no")
      .field("Requested field B", uni20::format_real(state.field))
      .field("Reached field B", state.reached_field ? uni20::format_real(*state.reached_field) : "infinite-field seed")
      .field("Status", status(state.status));
  if (state.energy)
    report.field("Energy evaluated at B", uni20::format_real(*state.reached_field))
        .field("Total energy", uni20::format_real(*state.energy));
  else
    report.field("Energy", "unavailable: no finite-field stage reached");
  report.field("Residual tolerance", uni20::format_real(options.residual_tolerance))
      .field("Reached backward residual", uni20::format_real(state.residual_norm))
      .field("Number-constraint error", uni20::format_real(state.number_error))
      .field("Newton corrections", state.iterations)
      .field("Continuation stages", state.stages)
      .field("Rejected stages", state.rejected_stages)
      .field("CPU time", cpu_time);
  if (args.variables)
  {
    auto& table = report.table("Eigenvalue variables (not occupations; spin-reversed frame if B<0)");
    table.header_separator()
        .column("Spin")
        .column("Coupling A", cli::table_alignment::decimal)
        .column("v", cli::table_alignment::decimal);
    for (std::size_t j = 0; j <= a.size(); ++j)
      table.row(j, j ? uni20::format_real(a[j - 1]) : "central",
                state.eigenvalue_variables.empty() ? "analytic" : uni20::format_real(state.eigenvalue_variables[j]));
  }
  cli::print_report(report, args.format);
  if (!state.converged) std::cerr << "Central-spin solve incomplete: no energy at the requested field is claimed.\n";
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
    std::cerr << "bethe-central-spin: " << error.what() << '\n';
    return 1;
  }
}
