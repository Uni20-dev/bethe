// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/q_boson.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::q_boson;
struct Arguments
{
    std::size_t sites = 0, particles = 0, iterations = 10000;
    std::optional<std::string> eta, tolerance;
    std::optional<std::size_t> count, max_candidates;
    std::string precision = "fp64";
    bool phase = false, roots = false;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-q-boson-pbc", "Periodic q-boson fixed-particle-number states.",
                                bethe::citations::Tool::q_boson_pbc);
  info.examples = {{"bethe-q-boson-pbc 16 --particles 8 --eta 0.5 --roots", "Finite deformation q=exp(eta)"},
                   {"bethe-q-boson-pbc 16 --particles 8 --phase --json phase.json", "Exact phase-model limit"},
                   {"bethe-q-boson-pbc 5 --particles 3 --eta 0.5 --excitations all --roots",
                    "Scan every canonical state in a small sector"}};
  info.notes = {
      "H=-sum_j(B_j^dagger B_{j+1}+h.c.-2N_j), hopping and lattice spacing one. The +2N shift is included.",
      "B|n>=sqrt([n]_q)|n-1>, [n]_q=(1-exp(-2*eta*n))/(1-exp(-2*eta)); q=exp(eta). "
      "This is deformed hopping, not the Bose-Hubbard Hamiltonian.",
      "eta=0 gives free bosons; --phase gives eta=+infinity, with multiple occupancy still allowed. "
      "For L=2 the periodic sum includes both bonds, doubling the hopping between the sites.",
      "Default: fixed-N ground state. --excitations COUNT|all scans all binomial(L+N-1,N) canonical mode patterns. "
      "COUNT limits retained levels, not solves. --max-candidates defaults to 10000.",
      "The default residual tolerance is 32 epsilon in the selected precision. Failed energies/roots are omitted "
      "(JSON null, CSV/TSV empty), with exit status 2. Tables: states, plus roots with --roots.",
      "See docs/q-boson.md for equations and limits. Use --references for literature. fp128 requires MPLAPACK."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::count_option(app, "L", a.sites, "Periodic sites, L>=2")->required();
  cli::count_option(app, "--particles", a.particles, "Particle number N>=0 (not restricted by L)")->required();
  auto* deformation = app.add_option_group("Deformation", "Choose finite eta or the phase limit")->require_option(1, 1);
  cli::text_option(*deformation, "--eta", a.eta, "Finite log(q)>=0")->type_name("REAL");
  deformation->add_flag("--phase", a.phase, "Exact eta=+infinity limit");
  app.add_flag("--roots", a.roots, "Export Bethe labels and converged momenta");
  cli::all_count_option(app, "--excitations", a.count, "Retain lowest converged levels, including ground state");
  cli::count_option(app, "--max-candidates", a.max_candidates, "Scan budget (default 10000); requires --excitations");
  cli::text_option(app, "--tolerance", a.tolerance, "Relative equation residual (default 32 epsilon)")
      ->type_name("REAL");
  cli::count_option(app, "--max-iterations", a.iterations, "Accepted Newton updates; zero only checks the seed")
      ->capture_default_str();
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
  Real const eta = a.phase ? uni20::numeric_limits<Real>::infinity() : uni20::parse_real<Real>(*a.eta);
  if (!a.phase && !uni20::isfinite(eta)) throw std::invalid_argument("--eta must be finite; use --phase for the limit");
  model::SolverOptions<Real> controls;
  controls.max_iterations = a.iterations;
  if (a.tolerance) controls.residual_tolerance = uni20::parse_real<Real>(*a.tolerance);
  // Validation and calculation precede opening any output files.
  std::optional<model::ExcitationScan<Real>> scan;
  model::State<Real> single;
  context.measure([&] {
    if (a.count)
      scan = model::real_excitations(a.sites, a.particles, eta,
                                     {.count = *a.count, .max_candidates = a.max_candidates.value_or(10000)}, controls);
    else
      single = model::ground_state(a.sites, a.particles, eta, controls);
  });
  auto const& state = scan ? scan->ground_state : single;
  bool const converged = scan ? scan->converged() : state.converged;
  cli::RunReport report(context, "Periodic q-boson states");
  report.field("sites", "Sites L", state.sites)
      .field("particles", "Particles N", state.particles)
      .field("deformation", "Deformation",
             a.phase          ? "phase limit"
             : eta == Real{0} ? "free bosons"
                              : "finite eta")
      .field("hamiltonian", "Hamiltonian", "-sum_j(B_j^dagger B_{j+1}+h.c.-2N_j); hopping=1, lattice spacing=1")
      .field("energy_reference", "Energy reference", "+2N diagonal shift included; subtract 2N for pure hopping")
      .field("calculation", "Calculation", scan ? "canonical fixed-N excitation scan" : "fixed-N ground state")
      .field("momentum_convention", "Momentum convention", "(-pi,pi]; index=sum(m) modulo L; lifted roots")
      .field("precision", "Precision", a.precision)
      .field("tolerance", "Residual tolerance", controls.residual_tolerance)
      .field("max_iterations", "Max Newton updates", controls.max_iterations)
      .field("energy", "Ground energy", state.energy)
      .field("residual", "Scaled equation residual", state.residual_norm)
      .field("iterations", "Newton updates", state.iterations)
      .result(converged, scan ? (converged ? "converged canonical family" : "incomplete scan; failed states excluded")
                              : name(state.status));
  if (scan)
    report.field("candidate_count", "Candidate count", scan->candidate_count)
        .field("converged_count", "Converged count", scan->converged_count)
        .field("retained_count", "Retained count", scan->levels.size())
        .field("max_candidates", "Candidate budget", a.max_candidates.value_or(10000));
  if (!a.phase) report.field("eta", "eta=log(q)", eta);
  using Optional = std::optional<Real>;
  using cli::column;
  std::vector<std::string> names{"states"};
  if (scan) names.push_back("reference");
  if (scan && scan->first_unconverged) names.push_back("failed");
  if (a.roots) names.push_back("roots");
  cli::ResultOutput output(report, a.output, names);
  std::vector<model::State<Real> const*> all;
  std::vector<Optional> gaps;
  if (scan)
    for (auto const& level : scan->levels)
    {
      all.push_back(&level.state);
      gaps.push_back(level.gap);
    }
  else
  {
    all.push_back(&state);
    gaps.push_back(state.converged ? Optional(Real{0}) : std::nullopt);
  }
  auto write_states = [&](std::string id, std::string title, std::size_t begin, std::size_t end, bool ranked) {
    output.table(
        id, title,
        [&](auto& table) {
          for (auto i = begin; i < end; ++i)
          {
            auto const& s = *all[i];
            table.append(i, s.energy, s.energy ? Optional(*s.energy / Real(s.sites)) : std::nullopt,
                         ranked ? gaps[i] : std::nullopt, s.momentum_index,
                         s.converged ? Optional(s.momentum) : std::nullopt, s.residual_norm, s.iterations, s.converged,
                         std::string(name(s.status)));
          }
        },
        column<std::size_t>("state_id"), column<Optional>("energy"), column<Optional>("energy_per_site"),
        column<Optional>("gap"), column<std::size_t>("momentum_index"), column<Optional>("momentum"),
        column<Real>("residual"), column<std::size_t>("iterations"), column<bool>("converged"),
        column<std::string>("status"));
  };
  write_states("states", scan ? "Lowest converged levels" : "Ground state", 0, all.size(), true);
  if (scan)
  {
    all.push_back(&state);
    write_states("reference", "Ground reference", all.size() - 1, all.size(), false);
    if (scan->first_unconverged)
    {
      all.push_back(&*scan->first_unconverged);
      write_states("failed", "First failed candidate (unranked)", all.size() - 1, all.size(), false);
    }
  }
  if (a.roots)
    output.table(
        "roots", "Lifted Bethe momenta: k=2*pi*m/L+deviation",
        [&](auto& table) {
          for (std::size_t i = 0; i < all.size(); ++i)
            for (std::size_t j = 0; j < all[i]->particles; ++j)
            {
              auto const& s = *all[i];
              table.append(i, j, s.quantum_numbers[j], s.modes[j], s.converged ? Optional(s.momenta[j]) : std::nullopt,
                           s.converged ? Optional(s.deviations[j]) : std::nullopt);
            }
        },
        column<std::size_t>("state_id"), column<std::size_t>("root_index"), column<uni20::half_int>("I"),
        column<std::size_t>("m"), column<Optional>("k").description("Lifted root; null on failure"),
        column<Optional>("deviation").description("Retains weak cluster splitting when rounded k values coincide"));
  output.finish();
  if (!converged) std::cerr << "q-boson calculation incomplete; failed energies and roots are unavailable.\n";
  return converged ? 0 : 2;
}
} // namespace
int main(int argc, char** argv)
{
  Arguments a;
  return cli::program_main(
      argc, argv, program_info(), [&](auto& app) { add_options(app, a); },
      [&](auto&) {
        a.output.validate();
        if (a.phase == a.eta.has_value()) throw std::invalid_argument("choose --eta or --phase");
        if (a.max_candidates && !a.count) throw std::invalid_argument("--max-candidates requires --excitations");
        return cli::dispatch_precision(a.precision, [&]<typename Real>() { return run<Real>(a, argc, argv); });
      });
}
