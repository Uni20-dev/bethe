// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "citation-report.hpp"
#include "report-common.hpp"
#include <bethe/hubbard.hpp>

namespace
{
namespace model = bethe::hubbard;
struct Arguments
{
    std::size_t sites;
    std::optional<std::string_view> interaction, tolerance;
    std::string_view precision = "fp64", format = "auto";
    std::size_t max_iterations = 10000;
    bool roots = false;
};
void usage(std::ostream& out)
{
  out << "Usage: bethe-hubbard-pbc L --u VALUE [options]\n"
      << "Periodic Hubbard ground state, t=1, repulsive U>=0.\n"
      << "H=-sum_(j,sigma)(c^dagger_(j,sigma)c_(j+1,sigma)+h.c.)+U sum_j n_up n_down.\n"
      << "First implementation: even L>=2, half filling (N=L), Sz=0 only.\n"
      << "  --u VALUE                          required interaction, finite U>=0\n"
      << "  --particles COUNT                  currently must equal L (default: L)\n"
      << "  --sz VALUE                         currently must be 0 (default: 0)\n"
      << "  --precision fp64|long-double|fp128  (default: fp64)\n"
      << "  --tolerance VALUE                  max normalized charge/spin residual (default: 32 epsilon)\n"
      << "  --max-iterations COUNT             total Newton update budget, including continuation\n"
      << "  --roots                            print charge momenta k and spin rapidities Lambda\n"
      << "  --format auto|pretty|plain         terminal report or script output (default: auto)\n"
      << "  --help                             show this help\n"
      << "The interaction is U*n_up*n_down, not the particle-hole-shifted convention.\n"
      << "L=2 counts the periodic hopping bond twice. U=0 uses exact free fermions.\n"
      << "At U>0 roots have separate charge I and spin J quantum numbers.\n"
      << "Residuals are evaluated at the requested U and are not energy-error bounds.\n"
      << "Doping, polarized sectors, excitations, attractive U, and OBC are not yet supported.\n"
      << "fp128 requires a Uni20 build with MPLAPACK enabled.\n";
  bethe::cli::print_citations(out, bethe::citations::Tool::hubbard_pbc);
}
Arguments parse(int argc, char** argv)
{
  Arguments result{.sites = bethe::cli::parse_size(argv[1]), .interaction = std::nullopt, .tolerance = std::nullopt};
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
      else if (option == "--precision")
        result.precision = argv[i];
      else if (option == "--format")
        result.format = argv[i];
      else if (option == "--tolerance")
        result.tolerance = argv[i];
      else if (option == "--max-iterations")
        result.max_iterations = bethe::cli::parse_size(argv[i]);
      else if (option == "--particles")
      {
        if (bethe::cli::parse_size(argv[i]) != result.sites)
          throw std::invalid_argument("currently only half filling is supported: --particles must equal L");
      }
      else if (uni20::half_int::parse(argv[i]) != uni20::half_int{0})
        throw std::invalid_argument("currently only Sz=0 is supported");
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
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  bethe::cli::CpuTimer const timer;
  auto const state = model::ground_state<Real>(args.sites, interaction, options);
  auto const cpu_time = timer.elapsed_text();
  auto const status = state.converged                               ? "converged"
                      : state.status == model::SolveStatus::stalled ? "line search stalled; unconverged estimate"
                                                                    : "iteration limit reached; unconverged estimate";
  auto const method = state.free_fermion ? "free fermions (U=0)" : "Lieb-Wu (damped Newton + continuation)";
  bool const pretty = args.format == "pretty" || (args.format == "auto" && terminal::is_a_terminal(stdout));
  if (pretty)
  {
    using bethe::cli::semantic_glyph;
    using bethe::cli::table_alignment;
    bethe::cli::report_builder report("Hubbard (periodic) - half-filled ground state");
    report.status(state.converged ? semantic_glyph::success : semantic_glyph::warning, status)
        .field("Model", "periodic Hubbard, t=1, U*n_up*n_down")
        .field("Sites", state.sites)
        .field("Particles", state.particles)
        .field("Down spins", state.down_spins)
        .field("Sz", "0")
        .field("U", uni20::format_real(interaction))
        .field("Precision", args.precision)
        .field("Method", method)
        .field("Residual tolerance", uni20::format_real(options.residual_tolerance))
        .field("CPU time", cpu_time)
        .field("Iterations", state.iterations)
        .field("Completed continuation stages", state.continuation_steps)
        .field("Charge residual", uni20::format_real(state.charge_residual))
        .field("Spin residual", uni20::format_real(state.spin_residual))
        .field("Residual norm", uni20::format_real(state.residual_norm))
        .field("Momentum index", state.momentum_index)
        .field("Momentum P", uni20::format_real(state.momentum))
        .field("Total energy", uni20::format_real(state.energy))
        .field("Energy per site", uni20::format_real(state.energy / Real(state.sites)));
    if (args.roots)
    {
      auto& charge =
          report.table(state.free_fermion ? "Free-fermion occupied momenta (no Bethe labels)" : "Charge momenta");
      charge.header_separator().column("Index").column("k", table_alignment::decimal);
      if (!state.free_fermion) charge.column("I");
      for (std::size_t j = 0; j < state.charge_momenta.size(); ++j)
        if (state.free_fermion)
          charge.row(j, uni20::format_real(state.charge_momenta[j]));
        else
          charge.row(j, uni20::format_real(state.charge_momenta[j]),
                     uni20::to_string_fraction(state.quantum_numbers.charge[j]));
      if (!state.free_fermion)
      {
        auto& spin = report.table("Spin rapidities (conventional Lambda)");
        spin.header_separator().column("Index").column("Lambda", table_alignment::decimal).column("J");
        for (std::size_t a = 0; a < state.spin_rapidities.size(); ++a)
          spin.row(a, uni20::format_real(state.spin_rapidities[a]),
                   uni20::to_string_fraction(state.quantum_numbers.spin[a]));
      }
    }
    bethe::cli::print_report(report);
  }
  else
  {
    std::cout << "Sites: " << state.sites << "\nModel: periodic Hubbard, t=1, U*n_up*n_down\n"
              << "Particles: " << state.particles << "\nDown spins: " << state.down_spins << "\nSz: 0\n"
              << "U: " << uni20::format_real(interaction) << "\nPrecision: " << args.precision << '\n'
              << "Method: " << method << "\nResidual tolerance: " << uni20::format_real(options.residual_tolerance)
              << '\n'
              << "Status: " << status << "\nIterations: " << state.iterations << "\nCPU time: " << cpu_time << '\n'
              << "Completed continuation stages: " << state.continuation_steps << '\n'
              << "Charge residual: " << uni20::format_real(state.charge_residual) << '\n'
              << "Spin residual: " << uni20::format_real(state.spin_residual) << '\n'
              << "Residual norm: " << uni20::format_real(state.residual_norm) << '\n'
              << "Momentum index: " << state.momentum_index << "\nMomentum: " << uni20::format_real(state.momentum)
              << '\n'
              << "Total energy: " << uni20::format_real(state.energy) << '\n'
              << "Energy per site: " << uni20::format_real(state.energy / Real(state.sites)) << '\n';
    if (args.roots)
    {
      std::cout << (state.free_fermion ? "# index k (free-fermion occupations; no Bethe labels)\n" : "# index k I\n");
      for (std::size_t j = 0; j < state.charge_momenta.size(); ++j)
      {
        std::cout << j << ' ' << uni20::format_real(state.charge_momenta[j]);
        if (!state.free_fermion) std::cout << ' ' << state.quantum_numbers.charge[j];
        std::cout << '\n';
      }
      if (!state.free_fermion)
      {
        std::cout << "# index Lambda J\n";
        for (std::size_t a = 0; a < state.spin_rapidities.size(); ++a)
          std::cout << a << ' ' << uni20::format_real(state.spin_rapidities[a]) << ' ' << state.quantum_numbers.spin[a]
                    << '\n';
      }
    }
  }
  if (!state.converged) std::cerr << "Hubbard solve: " << status << "; consider a larger budget or higher precision.\n";
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
    return bethe::cli::dispatch_precision(args.precision, [&]<uni20::Real Real> { return run<Real>(args); });
  }
  catch (std::exception const& error)
  {
    std::cerr << "bethe-hubbard-pbc: " << error.what() << '\n';
    return 1;
  }
}
