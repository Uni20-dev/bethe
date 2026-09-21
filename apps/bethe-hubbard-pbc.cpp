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
    std::optional<std::size_t> particles;
    uni20::half_int sz{0};
    std::optional<std::string_view> interaction, tolerance;
    std::string_view precision = "fp64", format = "auto";
    std::size_t max_iterations = 10000;
    bool roots = false;
};
void usage(std::ostream& out)
{
  out << "Usage: bethe-hubbard-pbc L --u VALUE [options]\n"
      << "Periodic Hubbard sector ground state, t=1, either sign of U.\n"
      << "H=-sum_(j,sigma)(c^dagger_(j,sigma)c_(j+1,sigma)+h.c.)+U sum_j n_up n_down.\n"
      << "Even L>=2; supported sector families are listed below.\n"
      << "  --u VALUE                          required finite interaction\n"
      << "  --particles COUNT                  0 <= N <= 2L (default: L)\n"
      << "  --sz VALUE                         integer or half-integer spin projection (default: 0)\n"
      << "  --precision fp64|long-double|fp128  (default: fp64)\n"
      << "  --tolerance VALUE                  max normalized charge/spin residual (default: 32 epsilon)\n"
      << "  --max-iterations COUNT             total Newton update budget, including continuation\n"
      << "  --roots                            print charge momenta k and spin rapidities Lambda\n"
      << "  --format auto|pretty|plain         terminal report or script output (default: auto)\n"
      << "  --help                             show this help\n"
      << "The interaction is U*n_up*n_down, not the particle-hole-shifted convention.\n"
      << "L=2 counts the periodic hopping bond twice. U=0 uses exact free fermions.\n"
      << "Repulsive root sectors: half filling with any Sz, doped odd N_up and N_down,\n"
      << "or a single spin species. Above half filling uses particle-hole symmetry.\n"
      << "Attractive U uses the Shiba mapping; all balanced even-N sectors are supported.\n"
      << "Other sectors work only if their mapped repulsive sector is supported. U=0 is unrestricted.\n"
      << "Other interacting shell parities, excitations, odd rings and OBC are not implemented.\n"
      << "Mapped roots and residuals explicitly describe the auxiliary sector, not attractive roots.\n"
      << "Residuals use the final root-sector U and are not energy-error bounds.\n"
      << "fp128 requires a Uni20 build with MPLAPACK enabled.\n";
  bethe::cli::print_citations(out, bethe::citations::Tool::hubbard_pbc);
}
Arguments parse(int argc, char** argv)
{
  Arguments result{.sites = bethe::cli::parse_size(argv[1]),
                   .particles = std::nullopt,
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
      else if (option == "--precision")
        result.precision = argv[i];
      else if (option == "--format")
        result.format = argv[i];
      else if (option == "--tolerance")
        result.tolerance = argv[i];
      else if (option == "--max-iterations")
        result.max_iterations = bethe::cli::parse_size(argv[i]);
      else if (option == "--particles")
        result.particles = bethe::cli::parse_size(argv[i]);
      else
        result.sz = uni20::half_int::parse(argv[i]);
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
  auto const state =
      model::sector_ground_state<Real>(args.sites, args.particles.value_or(args.sites), args.sz, interaction, options);
  auto const cpu_time = timer.elapsed_text();
  auto const status = state.converged                               ? "converged"
                      : state.status == model::SolveStatus::stalled ? "line search stalled; unconverged estimate"
                                                                    : "iteration limit reached; unconverged estimate";
  auto const method = state.free_fermion ? "exact free fermions" : "Lieb-Wu (damped Newton + continuation)";
  std::string mapping;
  if (state.shiba_transformed) mapping += "Shiba (down-spin particle-hole); ";
  if (state.particle_hole_transformed) mapping += "full particle-hole; ";
  if (state.spin_reversed) mapping += "spin reversal; ";
  if (mapping.empty())
    mapping = "none";
  else
    mapping.resize(mapping.size() - 2);
  bool const pretty = args.format == "pretty" || (args.format == "auto" && terminal::is_a_terminal(stdout));
  if (pretty)
  {
    using bethe::cli::semantic_glyph;
    using bethe::cli::table_alignment;
    bethe::cli::report_builder report("Hubbard (periodic) - sector ground state");
    report.status(state.converged ? semantic_glyph::success : semantic_glyph::warning, status)
        .field("Model", "periodic Hubbard, t=1, U*n_up*n_down")
        .field("Sites", state.sites)
        .field("Particles", state.particles)
        .field("Down spins", state.down_spins)
        .field("Sz", uni20::to_string_fraction(args.sz))
        .field("U", uni20::format_real(interaction))
        .field("Precision", args.precision)
        .field("Method", method)
        .field("Symmetry mapping", mapping)
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
    if (state.auxiliary_roots())
      report.field("Roots and residuals", "auxiliary sector (not physical-sector Bethe roots)")
          .field("Root particles", state.root_particles)
          .field("Root down spins", state.root_down_spins)
          .field("Root U", uni20::format_real(state.root_interaction))
          .field("Energy offset", uni20::format_real(state.energy_offset))
          .field("Momentum index offset", state.momentum_offset);
    if (args.roots)
    {
      auto& charge =
          report.table(std::string(state.auxiliary_roots() ? "Auxiliary: " : "") +
                       (state.free_fermion ? "Free-fermion occupied momenta (no Bethe labels)" : "Charge momenta"));
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
        auto& spin = report.table(std::string(state.auxiliary_roots() ? "Auxiliary: " : "") +
                                  "Spin rapidities (conventional Lambda)");
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
              << "Particles: " << state.particles << "\nDown spins: " << state.down_spins
              << "\nSz: " << uni20::to_string_fraction(args.sz) << '\n'
              << "U: " << uni20::format_real(interaction) << "\nPrecision: " << args.precision << '\n'
              << "Method: " << method << "\nSymmetry mapping: " << mapping
              << "\nResidual tolerance: " << uni20::format_real(options.residual_tolerance) << '\n'
              << "Status: " << status << "\nIterations: " << state.iterations << "\nCPU time: " << cpu_time << '\n'
              << "Completed continuation stages: " << state.continuation_steps << '\n'
              << "Charge residual: " << uni20::format_real(state.charge_residual) << '\n'
              << "Spin residual: " << uni20::format_real(state.spin_residual) << '\n'
              << "Residual norm: " << uni20::format_real(state.residual_norm) << '\n'
              << "Momentum index: " << state.momentum_index << "\nMomentum: " << uni20::format_real(state.momentum)
              << '\n'
              << "Total energy: " << uni20::format_real(state.energy) << '\n'
              << "Energy per site: " << uni20::format_real(state.energy / Real(state.sites)) << '\n';
    if (state.auxiliary_roots())
      std::cout << "Roots and residuals: auxiliary sector (not physical-sector Bethe roots)\n"
                << "Root particles: " << state.root_particles << "\nRoot down spins: " << state.root_down_spins
                << "\nRoot U: " << uni20::format_real(state.root_interaction)
                << "\nEnergy offset: " << uni20::format_real(state.energy_offset)
                << "\nMomentum index offset: " << state.momentum_offset << '\n';
    if (args.roots)
    {
      if (state.auxiliary_roots()) std::cout << "# Auxiliary-sector roots/occupations follow\n";
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
