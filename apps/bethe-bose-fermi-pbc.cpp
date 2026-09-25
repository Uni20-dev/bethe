// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/bose_fermi.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::bose_fermi;
struct Arguments
{
    std::size_t bosons = 0, fermions = 0, iterations = 10000;
    std::optional<std::string> length, coupling, tolerance;
    std::string precision = "fp64";
    bool roots = false;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-bose-fermi-pbc", "Equal-coupling Bose-Fermi ground states on a ring.",
                                bethe::citations::Tool::bose_fermi_pbc);
  info.examples = {
      {"bethe-bose-fermi-pbc --bosons 2 --fermions 3 --length 5 --c 1 --roots", "Mixed periodic ground state"},
      {"bethe-bose-fermi-pbc --bosons 4 --fermions 0 --length 4 --c 2", "Pure-boson Lieb-Liniger limit"}};
  info.notes = {
      "Equal masses and equal BB/BF repulsion: H=-sum_j d_j^2+2c sum_(i<j) delta(x_i-x_j), hbar^2/(2m)=1.",
      "Interacting mixed sectors require odd N_f. Pure species, vacuum and c=0 allow arbitrary counts.",
      "Free even fermion shells select the positive-momentum degenerate representative. Momentum is not folded.",
      "No auxiliary-auxiliary scattering; this is not an extra color of the SU(n) Fermi gas.",
      "Default tolerance: 32 epsilon. Failed energies/roots are missing, with exit status 2.",
      "Tables: states; --roots adds charge_roots and auxiliary_roots, or free_modes in free limits.",
      "See docs/bose-fermi.md for branches and conventions; --references for literature. fp128 requires MPLAPACK."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::count_option(app, "--bosons", a.bosons, "Boson population N_b>=0")->required();
  cli::count_option(app, "--fermions", a.fermions, "Spinless fermion population N_f>=0")->required();
  cli::text_option(app, "--length", a.length, "Physical circumference ell>0")->required()->type_name("REAL");
  cli::text_option(app, "--c", a.coupling, "Equal contact parameter c>=0 (inverse length)")
      ->required()
      ->type_name("REAL");
  cli::text_option(app, "--tolerance", a.tolerance, "Scaled equation residual; default 32 epsilon")->type_name("REAL");
  cli::count_option(app, "--max-iterations", a.iterations, "Accepted Newton updates; zero checks the seed")
      ->capture_default_str();
  app.add_flag("--roots", a.roots, "Export charge/auxiliary roots or exact free modes");
  cli::precision_option(app, a.precision);
  cli::add_data_output_options(app, a.output, true);
}
char const* name(model::Status status)
{
  switch (status)
  {
    case model::Status::converged:
      return "converged";
    case model::Status::iteration_limit:
      return "iteration_limit";
    case model::Status::stalled:
      return "stalled";
    case model::Status::precision_limit:
      return "precision_limit";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const length = uni20::parse_real<Real>(*a.length), c = uni20::parse_real<Real>(*a.coupling);
  model::SolverOptions<Real> controls;
  controls.max_iterations = a.iterations;
  if (a.tolerance) controls.residual_tolerance = uni20::parse_real<Real>(*a.tolerance);
  // Validate and solve before opening output files, including --force targets.
  auto const state = context.measure([&] { return model::ground_state(a.bosons, a.fermions, length, c, controls); });
  bool const free = c == Real{0} || a.bosons == 0;
  using Optional = std::optional<Real>;
  cli::RunReport report(context, "Bose-Fermi mixture (periodic)");
  report.field("bosons", "Bosons N_b", a.bosons)
      .field("fermions", "Fermions N_f", a.fermions)
      .field("length", "Circumference", length)
      .field("c", "Contact parameter c", c)
      .field("hamiltonian", "Hamiltonian", "-sum d_j^2+2c sum delta; equal masses, equal BB/BF couplings")
      .field("units", "Units", "hbar^2/(2m)=1; physical circumference, not lattice sites")
      .field("calculation", "Calculation",
             free         ? "exact free ground state"
             : a.fermions ? "odd-fermion mixed ground state"
                          : "Lieb-Liniger reduction")
      .field("precision", "Precision", a.precision)
      .field("tolerance", "Residual tolerance", controls.residual_tolerance)
      .field("max_iterations", "Max Newton updates", controls.max_iterations)
      .field("energy", "Energy", state.energy)
      .field("momentum", "Momentum", state.converged ? Optional(state.momentum) : std::nullopt)
      .field("momentum_convention", "Momentum convention",
             "P=2*pi*index/ell; no folding; positive representative for free even N_f")
      .field("residual", "Scaled equation residual", state.residual_norm)
      .field("iterations", "Newton updates", state.iterations)
      .result(state.converged, name(state.status));
  std::vector<std::string> tables{"states"};
  if (a.roots)
  {
    if (free)
      tables.push_back("free_modes");
    else
      tables.insert(tables.end(), {"charge_roots", "auxiliary_roots"});
  }
  cli::ResultOutput output(report, a.output, tables);
  using cli::column;
  output.table(
      "states", "Ground state",
      [&](auto& table) {
        table.append(0, state.energy, state.momentum_index, state.converged ? Optional(state.momentum) : std::nullopt,
                     state.residual_norm, state.iterations, state.converged, std::string(name(state.status)));
      },
      column<std::size_t>("state_id"), column<Optional>("energy"), column<std::int64_t>("momentum_index"),
      column<Optional>("momentum"), column<Real>("residual"), column<std::size_t>("iterations"),
      column<bool>("converged"), column<std::string>("status"));
  if (a.roots && free)
    output.table(
        "free_modes", "Free occupation modes",
        [&](auto& table) {
          for (std::size_t j = 0; j < a.bosons; ++j)
            table.append(0, std::string("boson"), j, std::int64_t{0},
                         state.converged ? Optional(Real{0}) : std::nullopt);
          for (std::size_t j = 0; j < a.fermions; ++j)
            table.append(0, std::string("fermion"), j, state.free_fermion_modes[j],
                         state.converged ? Optional(state.momenta[a.bosons + j]) : std::nullopt);
        },
        column<std::size_t>("state_id"), column<std::string>("species"), column<std::size_t>("index"),
        column<std::int64_t>("mode"), column<Optional>("k"));
  else if (a.roots)
    for (bool auxiliary : {false, true})
    {
      auto const& numbers = auxiliary ? state.auxiliary_numbers : state.charge_numbers;
      auto const& roots = auxiliary ? state.auxiliary : state.momenta;
      output.table(
          auxiliary ? "auxiliary_roots" : "charge_roots", auxiliary ? "Auxiliary rapidities" : "Charge momenta",
          [&](auto& table) {
            for (std::size_t j = 0; j < numbers.size(); ++j)
              table.append(0, j, numbers[j], state.converged ? Optional(roots[j]) : std::nullopt);
          },
          column<std::size_t>("state_id"), column<std::size_t>("index"), column<uni20::half_int>(auxiliary ? "J" : "I"),
          column<Optional>(auxiliary ? "lambda" : "k"));
    }
  output.finish();
  if (!state.converged) std::cerr << "Bose-Fermi solve incomplete; no energy or converged roots published.\n";
  return state.converged ? 0 : 2;
}
} // namespace
int main(int argc, char** argv)
{
  Arguments a;
  return cli::program_main(
      argc, argv, program_info(), [&](auto& app) { add_options(app, a); },
      [&](auto&) {
        a.output.validate();
        return cli::dispatch_precision(a.precision, [&]<typename Real>() { return run<Real>(a, argc, argv); });
      });
}
