// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/xyz.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::xyz;
struct Arguments
{
    std::size_t sites = 0, iterations = 10000;
    std::optional<std::string> eta, t, tolerance;
    std::string precision = "fp64";
    bool roots = false;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-xyz-pbc", "Even periodic XYZ ground states.", bethe::citations::Tool::xyz_pbc);
  info.examples = {{"bethe-xyz-pbc 16 --eta 0.4 --t 0.7 --roots", "Regular ground branch"},
                   {"bethe-xyz-pbc 32 --eta 0.7 --t 2 --precision long-double --json xyz.json", "Negative Jz"}};
  info.notes = {
      "H=sum_j(Jx Sx_j Sx_(j+1)+Jy Sy_j Sy_(j+1)+Jz Sz_j Sz_(j+1)), S=sigma/2; L=2 includes both bonds.",
      "Jx=theta4(eta)/theta4(0), Jy=theta3(eta)/theta3(0), Jz=theta2(eta)/theta2(0); theta(u)=vartheta(pi*u|i*t).",
      "Even L, 0<eta<1, t>0 only. Symmetric imaginary regular roots; no excited states or conserved Sz sectors.",
      "t->infinity gives XXZ with Delta=cos(pi*eta); eta=1/2 gives the anisotropic XY chain.",
      "Residual=max|F|/[L*(1-eta)], default tolerance 32 epsilon. Extreme parameters may stall or exceed range.",
      "Failed energies/roots are missing, with exit status 2. Tables: states; --roots adds roots.",
      "See docs/xyz.md for conventions and limits; --references for literature. fp128 requires MPLAPACK."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::count_option(app, "L", a.sites, "Even periodic sites, L>=2")->required();
  cli::text_option(app, "--eta", a.eta, "Real elliptic crossing parameter, 0<eta<1")->required()->type_name("REAL");
  cli::text_option(app, "--t", a.t, "Imaginary period tau=i*t, t>0")->required()->type_name("REAL");
  cli::text_option(app, "--tolerance", a.tolerance, "Scaled residual tolerance; default 32 epsilon")->type_name("REAL");
  cli::count_option(app, "--max-iterations", a.iterations, "Accepted Newton updates; zero checks the seed")
      ->capture_default_str();
  app.add_flag("--roots", a.roots, "Export labels and converged complex roots as real/imaginary columns");
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
  Real const eta = uni20::parse_real<Real>(*a.eta), t = uni20::parse_real<Real>(*a.t);
  model::SolverOptions<Real> controls;
  controls.max_iterations = a.iterations;
  if (a.tolerance) controls.residual_tolerance = uni20::parse_real<Real>(*a.tolerance);
  // Validate and solve before opening files, including --force targets.
  auto const state = context.measure([&] { return model::ground_state(a.sites, eta, t, controls); });
  using Optional = std::optional<Real>;
  // Coupling construction assigns all three together; x remains zero on failure.
  bool const have_exchange = state.exchange.x > Real{0};
  Optional const momentum = state.converged ? Optional(state.momentum) : std::nullopt;
  cli::RunReport report(context, "XYZ ground state (periodic)");
  report.field("sites", "Sites", a.sites)
      .field("eta", "Crossing parameter eta", eta)
      .field("t", "Imaginary period t", t)
      .field("jx", "Jx", have_exchange ? Optional(state.exchange.x) : std::nullopt)
      .field("jy", "Jy", have_exchange ? Optional(state.exchange.y) : std::nullopt)
      .field("jz", "Jz", have_exchange ? Optional(state.exchange.z) : std::nullopt)
      .field("hamiltonian", "Hamiltonian", "sum(Jx SxSx+Jy SySy+Jz SzSz); S=sigma/2; periodic")
      .field("calculation", "Calculation", "even-chain symmetric regular ground branch; xi=0; sum(lambda)=0")
      .field("precision", "Precision", a.precision)
      .field("tolerance", "Residual tolerance", controls.residual_tolerance)
      .field("max_iterations", "Max Newton updates", controls.max_iterations)
      .field("energy", "Energy", state.energy)
      .field("momentum", "Momentum", momentum)
      .field("momentum_convention", "Momentum convention", "lattice spacing one; P in (-pi,pi]")
      .field("residual", "Scaled equation residual", state.residual_norm)
      .field("iterations", "Newton updates", state.iterations)
      .result(state.converged, name(state.status));
  std::vector<std::string> tables{"states"};
  if (a.roots) tables.push_back("roots");
  cli::ResultOutput output(report, a.output, tables);
  using cli::column;
  output.table(
      "states", "Ground state",
      [&](auto& table) {
        table.append(0, state.energy, state.energy ? Optional(*state.energy / Real(a.sites)) : std::nullopt, momentum,
                     state.residual_norm, state.iterations, state.converged, std::string(name(state.status)));
      },
      column<std::size_t>("state_id"), column<Optional>("energy"), column<Optional>("energy_per_site"),
      column<Optional>("momentum"), column<Real>("residual"), column<std::size_t>("iterations"),
      column<bool>("converged"), column<std::string>("status"));
  if (a.roots)
    output.table(
        "roots", "Regular Bethe roots",
        [&](auto& table) {
          for (std::size_t j = 0; j < state.quantum_numbers.size(); ++j)
            table.append(0, j, state.quantum_numbers[j],
                         state.converged ? Optional(state.roots[j].real()) : std::nullopt,
                         state.converged ? Optional(state.roots[j].imag()) : std::nullopt);
        },
        column<std::size_t>("state_id"), column<std::size_t>("index"), column<uni20::half_int>("I"),
        column<Optional>("lambda_real"), column<Optional>("lambda_imag"));
  output.finish();
  if (!state.converged) std::cerr << "XYZ solve incomplete; no energy or converged roots published.\n";
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
