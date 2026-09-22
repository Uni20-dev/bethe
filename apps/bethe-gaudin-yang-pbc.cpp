// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "citation-report.hpp"
#include "report-common.hpp"
#include <bethe/gaudin_yang.hpp>

namespace
{
namespace model = bethe::gaudin_yang;
namespace cli = bethe::cli;
struct Arguments
{
    std::size_t particles = 0, max_iterations = 10000;
    std::optional<std::string_view> length, interaction, tolerance;
    std::optional<uni20::half_int> sz;
    std::string_view precision = "fp64", format = "auto";
    bool roots = false;
};

void usage(std::ostream& out)
{
  out << "Usage: bethe-gaudin-yang-pbc N --length ELL --c C [options]\n"
      << "Repulsive spin-1/2 continuum fermions on a ring; finite ELL>0 and C>=0.\n"
      << "H=-sum_j d_j^2 + 2c sum_(i<j) delta(x_i-x_j), hbar^2/(2m)=1.\n"
      << "Interacting mixed-spin ground states require odd N_up AND odd N_down.\n"
      << "At c=0 or full polarization any particle count is supported, including vacuum.\n"
      << "  --length ELL                       required physical circumference\n"
      << "  --c C                              required coupling (inverse length)\n"
      << "  --sz VALUE                         spin projection (default: 0 even N, 1/2 odd N)\n"
      << "  --precision fp64|long-double|fp128  (default: fp64; fp128 requires MPLAPACK)\n"
      << "  --tolerance VALUE                  component-scaled equation residual\n"
      << "                                     (default: 32 epsilon; not an energy-error bound)\n"
      << "  --max-iterations COUNT             accepted Newton updates (default: 10000)\n"
      << "  --roots                            print momenta and spin rapidities, or free modes\n"
      << "  --format auto|pretty|plain         (default: auto)\n"
      << "  --help                             show this help and references\n"
      << "E=sum(k_j^2); P is signed with no Brillouin-zone reduction.\n"
      << "Free even-population seas select the positive-current degenerate representative.\n"
      << "Other interacting shell branches, excitations, attraction and open ends\n"
      << "are not implemented. See docs/gaudin-yang.md for conventions and limits.\n";
  cli::print_citations(out, bethe::citations::Tool::gaudin_yang_pbc);
}

Arguments parse(int argc, char** argv)
{
  Arguments result;
  result.particles = cli::parse_size(argv[1]);
  for (int i = 2; i < argc; ++i)
  {
    std::string_view const option = argv[i];
    if (option == "--roots")
    {
      result.roots = true;
      continue;
    }
    if (option != "--length" && option != "--c" && option != "--sz" && option != "--precision" &&
        option != "--format" && option != "--tolerance" && option != "--max-iterations")
      throw std::invalid_argument("unknown option: " + std::string(option));
    if (++i == argc) throw std::invalid_argument("missing value for " + std::string(option));
    std::string_view const value = argv[i];
    if (option == "--length")
      result.length = value;
    else if (option == "--c")
      result.interaction = value;
    else if (option == "--sz")
      result.sz = uni20::half_int::parse(value);
    else if (option == "--precision")
      result.precision = value;
    else if (option == "--format")
      result.format = value;
    else if (option == "--tolerance")
      result.tolerance = value;
    else
      result.max_iterations = cli::parse_size(value);
  }
  if (!result.length || !result.interaction) throw std::invalid_argument("--length ELL and --c C are required");
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
  if (args.particles > std::size_t(std::numeric_limits<std::int64_t>::max() / 4))
    throw std::invalid_argument("Gaudin-Yang particle count exceeds the quantum-number range");
  auto const n = std::int64_t(args.particles);
  auto const sz = args.sz.value_or(uni20::from_twice(n % 2));
  auto const s = sz.twice();
  if (s < -n || s > n || (n - s) % 2)
    throw std::invalid_argument("Sz must satisfy |2*Sz|<=N and have the same parity as N");
  auto const down = std::size_t((n - s) / 2), up = args.particles - down;
  Real const length = uni20::parse_real<Real>(*args.length), c = uni20::parse_real<Real>(*args.interaction);
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  cli::CpuTimer const timer;
  auto const state = model::ground_state(up, down, length, c, options);
  auto const cpu_time = timer.elapsed_text();
  cli::report_builder report("Gaudin-Yang gas (periodic)");
  report.status(state.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning, status(state.status))
      .field("Calculation", state.free ? "exact free-fermion ground state" : "odd-population sector ground state")
      .field("Particles", args.particles)
      .field("N_up", up)
      .field("N_down", down)
      .field("Sz", uni20::to_string_fraction(sz))
      .field("Length", uni20::format_real(length))
      .field("c", uni20::format_real(c))
      .field("Units", "hbar^2/(2m)=1; interaction 2c delta")
      .field("Precision", args.precision)
      .field("Residual tolerance", uni20::format_real(options.residual_tolerance))
      .field("Status", status(state.status))
      .field("Total energy", uni20::format_real(state.energy))
      .field("Momentum index", state.momentum_index)
      .field("Momentum P", uni20::format_real(state.momentum))
      .field("Charge residual", uni20::format_real(state.charge_residual))
      .field("Spin residual", uni20::format_real(state.spin_residual))
      .field("Residual norm", uni20::format_real(state.residual_norm))
      .field("Root coupling reached", uni20::format_real(state.root_interaction))
      .field("Iterations", state.iterations)
      .field("CPU time", cpu_time);
  if (!state.free)
    report.field("Reference spin", state.spin_reversed ? "down (spin reversed)" : "up")
        .field("Spin roots", state.spin_rapidities.size());
  if (args.roots && state.free)
  {
    Real const pi = Real{4} * std::atan(Real{1});
    for (std::size_t species = 0; species < 2; ++species)
    {
      auto& table = report.table(species ? "Free down-spin modes" : "Free up-spin modes");
      table.header_separator().column("Mode").column("k", cli::table_alignment::decimal);
      for (auto mode : state.free_modes[species])
        table.row(mode, uni20::format_real(Real{2} * pi * Real(mode) / length));
    }
  }
  else if (args.roots)
  {
    auto& charges = report.table("Charge momenta");
    charges.header_separator().column("Index").column("I").column("k", cli::table_alignment::decimal);
    for (std::size_t j = 0; j < state.momenta.size(); ++j)
      charges.row(j, uni20::to_string_fraction(state.quantum_numbers.charge[j]), uni20::format_real(state.momenta[j]));
    auto& spins = report.table("Spin rapidities");
    spins.header_separator().column("Index").column("J").column("lambda", cli::table_alignment::decimal);
    for (std::size_t a = 0; a < state.spin_rapidities.size(); ++a)
      spins.row(a, uni20::to_string_fraction(state.quantum_numbers.spin[a]),
                uni20::format_real(state.spin_rapidities[a]));
  }
  cli::print_report(report, args.format);
  if (!state.converged)
    std::cerr << "Gaudin-Yang solve incomplete; residuals use requested c, not the intermediate coupling.\n";
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
    std::cerr << "bethe-gaudin-yang-pbc: " << error.what() << '\n';
    return 1;
  }
}
