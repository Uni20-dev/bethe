// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/xxx_structure_factor.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::heisenberg;
struct Arguments
{
    std::size_t sites = 0, max_candidates = 10000, max_iterations = 10000;
    std::string precision = "fp64", channel = "zz";
    std::optional<std::string> tolerance;
    bool roots = false, diagnostics = false;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-xxx-structure-factor", "Finite-ring XXX two-spinon spectral weights.",
                                bethe::citations::Tool::xxx_structure_factor);
  info.examples = {{"bethe-xxx-structure-factor 64 --csv lines.csv", "Zero-field two-spinon Szz lines"},
                   {"bethe-xxx-structure-factor 32 --channel raising --precision long-double --json spectrum.json",
                    "S-+ channel with native-precision export"}};
  info.notes = {
      "Even periodic N; H=sum S_j.S_(j+1), J=1, zero field; singlet ground reference.",
      "Exhaustive conventional real-root S=1 two-spinon family, NOT the full DSF.",
      "S(q,w)=2*pi*sum weight*delta(w-gap); S_q^+=sum exp(-iqj)S_j^+/sqrt(N).",
      "raising denotes S-+; zz weights are half the raising weights by SU(2).",
      "q is momentum transfer (P_ground-P_excited mod 2*pi), not the absolute state momentum.",
      "No broadening or renormalization to full sum rules. Failed weights are omitted, with exit status 2.",
      "Root tolerance does not certify form-factor accuracy; diagnostics include a pivot rejection indicator.",
      "See docs/xxx-structure-factor.md; --references for literature."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::count_option(app, "N", a.sites, "Even number of sites")->required();
  cli::option(app, "--channel", a.channel, "zz or raising (S-+)")
      ->check(CLI::IsMember({"zz", "raising"}))
      ->capture_default_str();
  cli::count_option(app, "--max-candidates", a.max_candidates, "Reject larger two-spinon families before solving")
      ->capture_default_str();
  cli::count_option(app, "--max-iterations", a.max_iterations, "Root iteration budget per state")
      ->capture_default_str();
  cli::text_option(app, "--tolerance", a.tolerance, "Normalized root residual tolerance")->type_name("REAL");
  cli::precision_option(app, a.precision);
  app.add_flag("--roots", a.roots, "Show ground and converged-state z=2*lambda roots on screen");
  app.add_flag("--diagnostics", a.diagnostics, "Show state convergence and form-factor diagnostics on screen");
  cli::add_data_output_options(
      app, a.output,
      {{.name = "spectrum", .description = "Accepted discrete (q,gap,weight) lines", .primary = true},
       {.name = "moments", .description = "Partial zeroth/first moments and full first-moment sum rule"},
       {.name = "diagnostics",
        .description = "Ground, converged root states and first failed root example",
        .screen_option = "--diagnostics"},
       {.name = "roots", .description = "Ground and converged-state roots (z=2*lambda)", .screen_option = "--roots"}});
}
std::string status(model::FormFactorStatus s)
{
  switch (s)
  {
    case model::FormFactorStatus::converged:
      return "converged";
    case model::FormFactorStatus::roots_unconverged:
      return "roots_unconverged";
    case model::FormFactorStatus::precision_limit:
      return "precision_limit";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  model::SolverOptions<Real> controls;
  controls.max_iterations = a.max_iterations;
  if (a.tolerance) controls.residual_tolerance = uni20::parse_real<Real>(*a.tolerance);
  auto const result =
      context.measure([&] { return model::two_spinon_structure_factor<Real>(a.sites, controls, a.max_candidates); });
  auto const& scan = result.scan;
  Real const scale = a.channel == "zz" ? Real{0.5} : Real{1};
  using Optional = std::optional<Real>;
  cli::RunReport report(context, "XXX two-spinon structure factor");
  report.field("sites", "Sites", a.sites)
      .field("hamiltonian", "Hamiltonian", "H=sum S_j.S_(j+1); periodic; J=1; zero field")
      .field("channel", "Channel", a.channel)
      .field("fourier_convention", "Fourier convention", "S_q^+=sum exp(-iqj)S_j^+/sqrt(N); q=P0-Pn")
      .field("weight_convention", "Weight convention", "S(q,w)=2*pi*sum weight*delta(w-gap); zz=raising/2")
      .field("family", "Family", "Two-spinon real-root S=1; partial DSF; no rescaling")
      .field("precision", "Precision", a.precision)
      .field("tolerance", "Root residual tolerance", controls.residual_tolerance)
      .field("max_iterations", "Max root iterations", controls.max_iterations)
      .field("max_candidates", "Max candidates", a.max_candidates)
      .field("ground_energy", "Ground energy",
             scan.ground_state.converged ? Optional(scan.ground_state.energy) : Optional{})
      .field("ground_momentum_index", "Ground momentum index", scan.ground_state.momentum_index)
      .field("candidates", "Candidates", scan.candidate_count)
      .field("root_converged", "Converged root states", scan.converged_count)
      .field("accepted_lines", "Accepted lines", result.lines.size())
      .field("integrated_weight", "Integrated weight (sum_q/N)", scale * result.integrated_weight)
      .field("full_weight_sum_rule", "Full integrated sum rule", scale / Real{2})
      .field("weight_fraction", "Integrated sum-rule fraction", result.weight_fraction)
      .field("first_moment_fraction", "Integrated first-moment fraction", result.first_moment_fraction)
      .result(result.converged(), result.converged() ? "converged" : "incomplete_two_spinon_family");
  cli::ResultOutput output(report, a.output, {"spectrum", "moments", "diagnostics", "roots"});
  using cli::column;
  output.table(
      "spectrum", "Accepted two-spinon lines",
      [&](auto& table) {
        for (auto const& line : result.lines)
          table.append(line.state_id, line.momentum_index, line.momentum, line.gap, scale * line.weight);
      },
      column<std::size_t>("state_id"), column<std::size_t>("momentum_index"), column<Real>("q"), column<Real>("gap"),
      column<Real>("weight"));
  using std::atan;
  Real const two_pi = Real{8} * atan(Real{1});
  output.table(
      "moments", "Partial moments and full first-moment sum rule",
      [&](auto& table) {
        for (std::size_t q = 0; q < a.sites; ++q)
        {
          Optional full, fraction;
          if (scan.ground_state.converged)
          {
            full = scale * model::raising_first_moment(a.sites, q, scan.ground_state.energy);
            if (*full > Real{0}) fraction = scale * result.moments[q].first_moment / *full;
          }
          table.append(q, two_pi * Real(q) / Real(a.sites), scale * result.moments[q].weight,
                       scale * result.moments[q].first_moment, full, fraction);
        }
      },
      column<std::size_t>("momentum_index"), column<Real>("q"), column<Real>("weight"), column<Real>("first_moment"),
      column<Optional>("full_first_moment"), column<Optional>("first_moment_fraction"));
  output.table(
      "diagnostics", "State diagnostics (failed root: first example only)",
      [&](auto& table) {
        auto append = [&](std::optional<std::size_t> id, std::string role, auto const& state, Optional pivot,
                          std::string ff_status) {
          table.append(id, role, state.momentum_index, state.energy, state.residual_norm, state.iterations,
                       state.converged, pivot, ff_status);
        };
        append(0, "ground", scan.ground_state, {}, "not_applicable");
        for (std::size_t i = 0; i < scan.levels.size(); ++i)
          append(i + 1, "excited", scan.levels[i].state, result.form_factors[i].pivot_margin,
                 status(result.form_factors[i].status));
        if (scan.first_unconverged) append({}, "failed_example", *scan.first_unconverged, {}, "roots_unconverged");
      },
      column<std::optional<std::size_t>>("state_id"), column<std::string>("role"),
      column<std::size_t>("state_momentum_index"), column<Real>("iterate_energy"), column<Real>("residual_norm"),
      column<std::size_t>("iterations"), column<bool>("roots_converged"), column<Optional>("pivot_margin"),
      column<std::string>("form_factor_status"));
  output.table(
      "roots", "Bethe roots (ground=0; z=2*lambda)",
      [&](auto& table) {
        auto append = [&](std::size_t id, auto const& state) {
          for (std::size_t j = 0; j < state.rapidities.size(); ++j)
            table.append(id, j, state.quantum_numbers[j], state.rapidities[j], state.converged);
        };
        append(0, scan.ground_state);
        for (std::size_t i = 0; i < scan.levels.size(); ++i)
          append(i + 1, scan.levels[i].state);
      },
      column<std::size_t>("state_id"), column<std::size_t>("root_index"), column<uni20::half_int>("I"),
      column<Real>("z"), column<bool>("converged"));
  output.finish();
  if (!result.converged())
    std::cerr << "Incomplete two-spinon calculation: inspect counts and --diagnostics; no sum-rule rescaling.\n";
  return result.converged() ? 0 : 2;
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
