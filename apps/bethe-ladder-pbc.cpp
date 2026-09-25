// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/ladder.hpp>

namespace
{
namespace model = bethe::ladder;
namespace cli = bethe::cli;
struct Arguments
{
    std::size_t rungs = 0, max_iterations = 10000, max_branches = 10000;
    std::optional<std::size_t> singlets;
    std::optional<std::string> rung, field, tolerance;
    std::string precision = "fp64";
    cli::DataOutputOptions output;
    bool sectors = false, roots = false;
};
auto program_info()
{
  auto info = bethe::cli::program_info("bethe-ladder-pbc",
                                       "Wang's integrable spin-1/2 ladder, L>=2 periodic rungs, longitudinal field.",
                                       bethe::citations::Tool::ladder_pbc);
  info.notes = {"H=sum[S.S_next+T.T_next+4(S.S_next)(T.T_next)]+JR*sum S.T-h*sum(Sz+Tz).",
                "This is NOT the ordinary two-leg Heisenberg ladder. JR may have either sign.",
                "Sector minima include compatible SU(4) descendants; triplet populations",
                "are minimized too, not fixed Sz. --sectors --roots prints the best state.",
                "--field h uses energy units (g*mu_B absorbed); default zero. M=N_t+ - N_t- is total Sz.",
                "At zero field, populations are a balanced representative, not a unique magnetization.",
                "No excitations, fixed Sz, arbitrary four-spin couplings or open ends.",
                "An incomplete scan reports only a candidate upper bound, never a minimum.",
                "See docs/ladder.md for normalization, finite-ring labels and scan cost.",
                "Use --references for literature and applicability; see CITATIONS.md."};
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  bethe::cli::option(app, "L", args.rungs, "Number of sites or rungs")->required();
  bethe::cli::option(app, "--rung", args.rung, "Rung coupling, either sign")->required();
  bethe::cli::option(app, "--field", args.field, "Longitudinal Zeeman coefficient h, either sign; default 0");
  bethe::cli::option(app, "--singlets", args.singlets, "fixed singlet-count sector minimum");
  bethe::cli::option(app, "--sectors", args.sectors,
                     "all singlet-count sector minima (mutually exclusive; default: global minimum)");
  bethe::cli::option(app, "--max-branches", args.max_branches, "total highest-weight sea branches (10000)")
      ->capture_default_str();
  bethe::cli::option(app, "--roots", args.roots, "highest-weight roots/labels for selected state");
  bethe::cli::option(app, "--tolerance", args.tolerance,
                     "max Bethe residual divided by L (default: 32 epsilon; not an energy-error bound)");
  bethe::cli::option(app, "--max-iterations", args.max_iterations, "total attempted Newton corrections (10000)")
      ->capture_default_str();
  bethe::cli::precision_option(app, args.precision);
  cli::add_data_output_options(app, args.output, true);
}

void validate(Arguments const& args)
{
  if (!args.rung) throw std::invalid_argument("--rung is required");
  if (args.sectors && args.singlets) throw std::invalid_argument("--sectors and --singlets are mutually exclusive");
  args.output.validate();
}
char const* status(model::SolveStatus value)
{
  switch (value)
  {
    case model::SolveStatus::converged:
      return "converged";
    case model::SolveStatus::iteration_limit:
      return "Newton budget exhausted; incomplete scan";
    case model::SolveStatus::branch_limit:
      return "branch budget exhausted; incomplete scan";
    case model::SolveStatus::stalled:
      return "line search stalled; incomplete scan";
    case model::SolveStatus::ill_conditioned:
      return "ill-conditioned Newton system; incomplete scan";
  }
  return "unknown";
}
std::string shape_text(model::detail::Shape const& shape)
{
  return fmt::format("{}, {}, {}, {}", shape[0], shape[1], shape[2], shape[3]);
}

template <uni20::Real Real> int run(Arguments const& args, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const rung = uni20::parse_real<Real>(*args.rung);
  Real const field = args.field ? uni20::parse_real<Real>(*args.field) : Real{0};
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  options.max_branches = args.max_branches;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  auto computation = context.computation();
  model::State<Real> state;
  std::optional<model::SectorScan<Real>> scan;
  if (args.sectors)
  {
    scan = model::sector_ground_states(args.rungs, rung, field, options);
    state = model::minimum_candidate(*scan);
  }
  else if (args.singlets)
    state = model::sector_ground_state(args.rungs, *args.singlets, rung, field, options);
  else
    state = model::ground_state(args.rungs, rung, field, options);
  computation.finish();
  cli::RunReport report(context, "Integrable spin ladder (periodic)");
  report.status(state.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning, status(state.status))
      .field("hamiltonian", "Hamiltonian", "H=sum[S.S_next+T.T_next+4(S.S_next)(T.T_next)]+J_r*sum S.T-h*sum(Sz+Tz)")
      .field("calculation", "Calculation",
             args.singlets  ? "fixed singlet-count sector"
             : args.sectors ? "all singlet-count sectors"
                            : "global ground state")
      .field("rungs", "Rungs", args.rungs)
      .field("rung_coupling_j_r", "Rung coupling J_r", rung)
      .field("magnetic_field", "Magnetic field h", field)
      .field("magnetization_convention", "Magnetization convention",
             "M=N_t+ - N_t-; total Sz, one minimizing representative, not a degeneracy average")
      .field("precision", "Precision", args.precision)
      .field("residual_tolerance", "Residual tolerance", options.residual_tolerance)
      .result(state.converged, status(state.status));
  if (state.energy)
  {
    report.field("total_energy", state.converged ? "Total energy" : "Candidate energy (upper bound)", *state.energy)
        .field("singlets_n_s", "Singlets N_s", state.singlets)
        .field("magnetization", "Total magnetization M", state.magnetization)
        .field("populations_s_t_t0_t", "Populations (s,t+,t0,t-)", shape_text(state.populations))
        .field("highest_weight_rows", "Highest-weight rows", shape_text(state.highest_weight.shape))
        .field("su_4_descendant", "SU(4) descendant",
               state.descendant ? "yes; roots belong to the highest-weight representative"
                                : "no (up to color permutation)")
        .field("permutation_energy", "Permutation energy", state.highest_weight.energy)
        .field("momentum_index", "Momentum index", state.highest_weight.momentum_index)
        .field("reflected_momentum_index", "Reflected momentum index",
               (state.rungs - state.highest_weight.momentum_index) % state.rungs)
        .field("selected_branch_residual", "Selected branch residual", state.highest_weight.residual);
    if (state.converged) report.field("energy_per_rung", "Energy per rung", *state.energy / Real(state.rungs));
    if (state.analytic)
      report.field("exact_limit", "Exact limit",
                   state.singlets == state.rungs ? "rung-singlet product (no root solve)"
                                                 : "fully polarized triplet product (no root solve)");
  }
  else
    report.field("energy", "Energy", state.energy, {.missing = "unavailable: no converged branch"});
  report.field("highest_weights_visited", "Highest weights visited", state.tableaux)
      .field("sea_branches_attempted", "Sea branches attempted", state.branches)
      .field("newton_corrections", "Newton corrections", state.iterations);
  std::vector<std::string> tables{"states", "representations"};
  if (args.roots) tables.push_back("roots");
  cli::ResultOutput output(report, args.output, tables);
  // Scan rows retain their sector identity; roots below belong only to the selected row.
  auto const selected_id = scan ? state.singlets : std::size_t{0};
  auto for_states = [&](auto&& append) {
    if (scan)
      for (std::size_t j = 0; j < scan->sectors.size(); ++j)
        append(j, scan->sectors[j]);
    else
      append(std::size_t{0}, state);
  };
  output.table(
      "states",
      scan ? (scan->complete ? "Singlet-sector minima" : "Singlet-sector candidate upper bounds (scan incomplete)")
           : "State",
      [&](auto& table) {
        for_states([&](std::size_t id, auto const& s) {
          table.append(id, s.singlets, s.energy, s.energy ? std::optional{s.magnetization} : std::nullopt,
                       id == selected_id, s.energy ? std::optional{s.descendant} : std::nullopt,
                       s.energy ? std::optional{s.highest_weight.momentum_index} : std::nullopt,
                       s.energy ? std::optional{s.highest_weight.residual} : std::nullopt, s.iterations, s.branches,
                       s.tableaux, s.converged, std::string(status(s.status)));
        });
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("singlets", "N_s"),
      cli::column<std::optional<Real>>("energy", "Energy"),
      cli::column<std::optional<std::int64_t>>("magnetization", "M"), cli::column<bool>("selected"),
      cli::column<std::optional<bool>>("descendant", "Descendant"),
      cli::column<std::optional<std::size_t>>("momentum_index", "Momentum index"),
      cli::column<std::optional<Real>>("residual"), cli::column<std::size_t>("iterations"),
      cli::column<std::size_t>("branches"), cli::column<std::size_t>("tableaux"), cli::column<bool>("converged"),
      cli::column<std::string>("status"));
  output.table(
      "representations", "Physical populations and highest-weight rows",
      [&](auto& table) {
        for_states([&](std::size_t id, auto const& s) {
          if (s.energy)
            for (std::size_t a = 0; a < 4; ++a)
              table.append(id, a, s.populations[a], s.highest_weight.shape[a]);
        });
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("component", "Component (s,t+,t0,t-)"),
      cli::column<std::size_t>("population", "Physical population"),
      cli::column<std::size_t>("highest_weight", "Highest-weight row"));
  if (args.roots)
    output.table(
        "roots", "Highest-weight rapidities (selected state)",
        [&](auto& table) {
          if (state.energy)
            for (std::size_t a = 0; a < state.highest_weight.rapidities.size(); ++a)
              for (std::size_t j = 0; j < state.highest_weight.rapidities[a].size(); ++j)
                table.append(selected_id, a + 1, j, state.highest_weight.labels[a][j],
                             state.highest_weight.rapidities[a][j]);
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("level", "Level"),
        cli::column<std::size_t>("index", "Index"), cli::column<uni20::half_int>("quantum_number", "Bethe label"),
        cli::column<Real>("rapidity", "Rapidity"));
  output.finish();
  if (!state.converged) std::cerr << "Ladder scan incomplete: candidate energies are not certified sector minima.\n";
  return state.converged ? 0 : 2;
}
} // namespace
int main(int argc, char** argv)
{
  Arguments args;
  return bethe::cli::program_main(
      argc, argv, program_info(), [&](auto& app) { add_options(app, args); },
      [&](auto&) {
        validate(args);
        return bethe::cli::dispatch_precision(args.precision,
                                              [&]<uni20::Real Real> { return run<Real>(args, argc, argv); });
      });
}
