// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "exclusion-output.hpp"
#include "program-options.hpp"
#include <bethe/asep.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::asep;
struct Arguments
{
    std::size_t sites = 0, particles = 0, iterations = 100, steps = 256;
    std::size_t seed_iterations = 1000, seed_newton = 100, max_sites = 256;
    std::optional<std::string> right, left, tolerance, seed_tolerance;
    std::string precision = "fp64";
    bool roots = false;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-asep-pbc", "Periodic asymmetric exclusion relaxation gaps.",
                                bethe::citations::Tool::asep_pbc);
  info.examples = {
      {"bethe-asep-pbc 32 --particles 8 --right-rate 1 --left-rate 0.5 --roots", "Bidirectional relaxation pair"},
      {"bethe-asep-pbc 64 --particles 32 --left-rate 1 --precision long-double --json gap.json",
       "Symmetric exclusion gap"}};
  info.notes = {
      "Particles hop to empty neighbours at nonnegative right/left rates; at least one rate must be positive.",
      "dP/dt=M P. Gap=-Re(lambda)>0; frequency=Im(lambda)>=0 selects one conjugate partner. These are rates, not "
      "energies.",
      "Empty/full sectors have no relaxation mode: observables are missing, not zero.",
      "Reflection/particle-hole reduction uses min(N,L-N) roots; no physical momentum assignment.",
      "TASEP seed followed by bounded continuation. Equal rates and one particle/hole are analytic, with no root rows.",
      "Tables: relaxation; --roots adds scaled roots v, not ordinary z. See metadata for reconstruction.",
      "Failed observables/root coordinates are missing; exit status 2. No open reservoirs or spectrum scan.",
      "See docs/asep.md; --references for literature. fp128 requires MPLAPACK."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::count_option(app, "L", a.sites, "Periodic sites, L>=2")->required();
  cli::count_option(app, "--particles", a.particles, "Particle count, 0<=N<=L")->required();
  cli::text_option(app, "--right-rate", a.right, "Right-hop rate (default 1)")->type_name("REAL");
  cli::text_option(app, "--left-rate", a.left, "Left-hop rate (default 0)")->type_name("REAL");
  cli::text_option(app, "--tolerance", a.tolerance, "Rescaled residual tolerance; default 128 epsilon")
      ->type_name("REAL");
  cli::text_option(app, "--seed-tolerance", a.seed_tolerance, "TASEP residual tolerance; default 128 epsilon")
      ->type_name("REAL");
  cli::count_option(app, "--max-iterations", a.iterations, "Newton updates per continuation attempt")
      ->capture_default_str();
  cli::count_option(app, "--max-continuation-steps", a.steps, "Continuation attempts, including rejected steps")
      ->capture_default_str();
  cli::count_option(app, "--max-seed-iterations", a.seed_iterations, "TASEP all-root seed sweeps")
      ->capture_default_str();
  cli::count_option(app, "--max-seed-newton-iterations", a.seed_newton, "TASEP coupled Newton updates")
      ->capture_default_str();
  cli::count_option(app, "--max-sites", a.max_sites, "Explicit site/work budget")->capture_default_str();
  app.add_flag("--roots", a.roots, "Export scaled reduced-filling roots v; analytic cases have no root rows");
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
    case model::Status::continuation_limit:
      return "continuation_limit";
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
  Real const right = a.right ? uni20::parse_real<Real>(*a.right) : Real{1};
  Real const left = a.left ? uni20::parse_real<Real>(*a.left) : Real{0};
  model::Options<Real> controls;
  controls.max_iterations = a.iterations;
  controls.max_continuation_steps = a.steps;
  controls.seed_options.max_iterations = a.seed_newton;
  controls.seed_options.max_seed_iterations = a.seed_iterations;
  controls.seed_options.max_sites = a.max_sites;
  if (a.tolerance) controls.tolerance = uni20::parse_real<Real>(*a.tolerance);
  if (a.seed_tolerance) controls.seed_options.tolerance = uni20::parse_real<Real>(*a.seed_tolerance);
  // Validate even unused analytic-branch controls before opening output files.
  if (!uni20::isfinite(controls.seed_options.tolerance) || !(controls.seed_options.tolerance > Real{0}))
    throw std::invalid_argument("ASEP seed tolerance must be finite and positive");
  auto const state =
      context.measure([&] { return model::relaxation_gap(a.sites, a.particles, right, left, controls); });
  using Optional = std::optional<Real>;
  Optional const real = state.eigenvalue ? Optional(state.eigenvalue->real()) : std::nullopt;
  Optional const imag = state.eigenvalue ? Optional(state.eigenvalue->imag()) : std::nullopt;
  bool const regular = state.converged && !state.analytic;
  cli::RunReport report(context, "ASEP relaxation (periodic)");
  report.field("sites", "Sites", a.sites)
      .field("particles", "Particles", a.particles)
      .field("effective_particles", "Reduced root count", state.effective_particles)
      .field("right_rate", "Right-hop rate", right)
      .field("left_rate", "Left-hop rate", left)
      .field("generator", "Generator", "dP/dt=M P; right/left hops; columns sum to zero")
      .field("units", "Units", "inverse time; not quantum energies")
      .field("stationary_eigenvalue", "Stationary eigenvalue", Real{0})
      .field("calculation", "Calculation", "leading relaxation pair; representative Im(lambda)>=0")
      .field("root_convention", "Root convention",
             "n=min(N,L-N); x=min(r,s)/max(r,s); delta=(max(r,s)-min(r,s))/max(r,s); "
             "z_j=1+delta*v_j+[j=wave_index]*wave_base")
      .field("wave_index", "Wave root index", regular ? std::optional<std::size_t>(state.wave_index) : std::nullopt)
      .field("wave_base_real", "Wave base real", regular ? Optional(state.wave_base.real()) : std::nullopt)
      .field("wave_base_imag", "Wave base imaginary", regular ? Optional(state.wave_base.imag()) : std::nullopt)
      .field("has_relaxation_mode", "Has relaxation mode", state.effective_particles != 0)
      .field("analytic", "Analytic result", state.analytic)
      .field("precision", "Precision", a.precision)
      .field("tolerance", "Residual tolerance", controls.tolerance)
      .field("seed_tolerance", "Seed residual tolerance", controls.seed_options.tolerance)
      .field("max_iterations", "Max Newton updates per stage", controls.max_iterations)
      .field("max_continuation_steps", "Max continuation attempts", controls.max_continuation_steps)
      .field("max_seed_iterations", "Max seed sweeps", controls.seed_options.max_seed_iterations)
      .field("max_seed_newton_iterations", "Max seed Newton updates", controls.seed_options.max_iterations)
      .field("max_sites", "Site budget", controls.seed_options.max_sites)
      .field("lambda_real", "Eigenvalue real part", real)
      .field("lambda_imag", "Eigenvalue imaginary part", imag)
      .field("gap", "Relaxation gap", state.gap)
      .field("frequency", "Angular frequency", imag)
      .field("residual", "Last reached residual", state.residual_norm)
      .field("reached_ratio", "Last reached rate ratio", state.reached_ratio)
      .field("iterations", "Continuation Newton updates", state.iterations)
      .field("continuation_steps", "Continuation attempts", state.continuation_steps)
      .field("seed_iterations", "Seed sweeps", state.seed_iterations)
      .field("seed_newton_iterations", "Seed Newton updates", state.seed_newton_iterations)
      .result(state.converged, name(state.status));
  std::vector<std::string> tables{"relaxation"};
  if (a.roots) tables.push_back("roots");
  cli::ResultOutput output(report, a.output, tables);
  cli::relaxation_table(output, state, name(state.status));
  if (a.roots)
    output.table(
        "roots", "Scaled reduced-filling roots",
        [&](auto& table) {
          if (!state.analytic)
            for (std::size_t j = 0; j < state.effective_particles; ++j)
              table.append(j, regular ? Optional(state.scaled_roots[j].real()) : std::nullopt,
                           regular ? Optional(state.scaled_roots[j].imag()) : std::nullopt);
        },
        cli::column<std::size_t>("index"), cli::column<Optional>("v_real"), cli::column<Optional>("v_imag"));
  output.finish();
  if (!state.converged) std::cerr << "ASEP solve incomplete; no decay eigenvalue, gap or converged roots published.\n";
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
