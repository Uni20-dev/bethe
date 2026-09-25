// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "exclusion-output.hpp"
#include "program-options.hpp"
#include <bethe/tasep.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::tasep;
struct Arguments
{
    std::size_t sites = 0, particles = 0, iterations = 100, seed_iterations = 1000, max_sites = 256;
    std::optional<std::string> rate, tolerance;
    std::string precision = "fp64";
    bool roots = false;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-tasep-pbc", "Periodic totally asymmetric exclusion relaxation gaps.",
                                bethe::citations::Tool::tasep_pbc);
  info.examples = {{"bethe-tasep-pbc 32 --particles 8 --roots", "Leading complex decay mode"},
                   {"bethe-tasep-pbc 64 --particles 32 --rate 2 --precision long-double --json gap.json",
                    "Half-filled relaxation gap"}};
  info.notes = {
      "Each particle jumps right at rate r if the next site is empty. dP/dt=M P; these are rates, not energies.",
      "The stationary eigenvalue is zero. Gap=-Re(lambda)>0; frequency=Im(lambda)>=0 selects one conjugate partner.",
      "Empty/full sectors have no relaxation mode: eigenvalue, gap and frequency are missing, not zero.",
      "Particle-hole reduction uses min(N,L-N) roots; no physical momentum is assigned to this representative.",
      "Default tolerance: 128 epsilon. Separate seed/Newton budgets; one particle or hole is analytic.",
      "Tables: relaxation; --roots adds roots. Failed observables and roots are missing, with exit status 2.",
      "Only totally asymmetric PBC gaps: no backward hopping, reservoirs or full-spectrum scan.",
      "See docs/tasep.md; --references for literature. fp128 requires MPLAPACK."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::count_option(app, "L", a.sites, "Periodic sites, L>=2")->required();
  cli::count_option(app, "--particles", a.particles, "Particle count, 0<=N<=L")->required();
  cli::text_option(app, "--rate", a.rate, "Right-hop rate r>0 (default 1)")->type_name("REAL");
  cli::text_option(app, "--tolerance", a.tolerance, "Scaled logarithmic residual; default 128 epsilon")
      ->type_name("REAL");
  cli::count_option(app, "--max-iterations", a.iterations, "Accepted coupled Newton updates")->capture_default_str();
  cli::count_option(app, "--max-seed-iterations", a.seed_iterations, "All-root seed sweeps")->capture_default_str();
  cli::count_option(app, "--max-sites", a.max_sites, "Explicit site/work budget")->capture_default_str();
  app.add_flag("--roots", a.roots, "Export converged reduced-filling fugacity roots Z");
  cli::precision_option(app, a.precision);
  cli::add_data_output_options(app, a.output, true);
}
char const* name(model::Status s)
{
  switch (s)
  {
    case model::Status::converged:
      return "converged";
    case model::Status::stationary_only:
      return "stationary_only";
    case model::Status::seed_limit:
      return "seed_limit";
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
  Real const rate = a.rate ? uni20::parse_real<Real>(*a.rate) : Real{1};
  model::Options<Real> controls;
  controls.max_iterations = a.iterations;
  controls.max_seed_iterations = a.seed_iterations;
  controls.max_sites = a.max_sites;
  if (a.tolerance) controls.tolerance = uni20::parse_real<Real>(*a.tolerance);
  // Invalid inputs must not touch output targets, including --force files.
  auto const state = context.measure([&] { return model::relaxation_gap(a.sites, a.particles, rate, controls); });
  using Optional = std::optional<Real>;
  Optional const real = state.eigenvalue ? Optional(state.eigenvalue->real()) : std::nullopt;
  Optional const imag = state.eigenvalue ? Optional(state.eigenvalue->imag()) : std::nullopt;
  bool const has_mode = state.effective_particles != 0;
  cli::RunReport report(context, "TASEP relaxation (periodic)");
  report.field("sites", "Sites", a.sites)
      .field("particles", "Particles", a.particles)
      .field("effective_particles", "Reduced root count", state.effective_particles)
      .field("rate", "Right-hop rate", rate)
      .field("generator", "Generator", "dP/dt=M P; right hops only; columns sum to zero")
      .field("units", "Units", "inverse time; not quantum energies")
      .field("stationary_eigenvalue", "Stationary eigenvalue", Real{0})
      .field("calculation", "Calculation", "leading relaxation pair; representative Im(lambda)>=0")
      .field("root_convention", "Root convention", "Z=2/z-1; particle-hole reduced filling min(N,L-N)")
      .field("has_relaxation_mode", "Has relaxation mode", has_mode)
      .field("precision", "Precision", a.precision)
      .field("tolerance", "Residual tolerance", controls.tolerance)
      .field("max_iterations", "Max Newton updates", controls.max_iterations)
      .field("max_seed_iterations", "Max seed sweeps", controls.max_seed_iterations)
      .field("max_sites", "Site budget", controls.max_sites)
      .field("lambda_real", "Eigenvalue real part", real)
      .field("lambda_imag", "Eigenvalue imaginary part", imag)
      .field("gap", "Relaxation gap", state.gap)
      .field("frequency", "Angular frequency", imag)
      .field("residual", "Scaled equation residual", state.residual_norm)
      .field("iterations", "Newton updates", state.iterations)
      .field("seed_iterations", "Seed sweeps", state.seed_iterations)
      .result(state.converged, name(state.status));
  std::vector<std::string> tables{"relaxation"};
  if (a.roots) tables.push_back("roots");
  cli::ResultOutput output(report, a.output, tables);
  using cli::column;
  cli::relaxation_table(output, state, name(state.status));
  if (a.roots)
    output.table(
        "roots", "Reduced-filling fugacity roots",
        [&](auto& table) {
          for (std::size_t j = 0; j < state.effective_particles; ++j)
            table.append(j, state.converged ? Optional(state.roots[j].real()) : std::nullopt,
                         state.converged ? Optional(state.roots[j].imag()) : std::nullopt);
        },
        column<std::size_t>("index"), column<Optional>("Z_real"), column<Optional>("Z_imag"));
  output.finish();
  if (!state.converged) std::cerr << "TASEP solve incomplete; no decay eigenvalue, gap or converged roots published.\n";
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
