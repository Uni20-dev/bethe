// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "execution-options.hpp"
#include "result-output.hpp"
#include "solver-precision.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/potts.hpp>

namespace
{
namespace cli = bethe::cli;
struct Arguments
{
    std::size_t sites = 0, iterations = 10000, max_sites = 512;
    std::optional<std::size_t> momentum_index;
    std::optional<std::string> tolerance;
    std::string branch = "all", charge = "both", precision = "fp64";
    int threads = 1;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-potts-pbc", "Critical three-state Potts vacuum and charged one-hole levels.",
                                bethe::citations::Tool::potts_pbc);
  info.examples = {{"bethe-potts-pbc 16 --json potts.json", "Vacuum and one charged level per momentum and charge"},
                   {"bethe-potts-pbc 64 --momentum-index 1 --precision long-double",
                    "First charged descendants, with finite-size gaps"}};
  info.notes = {
      "H=-sum(X+X^dagger+Z_j Z_{j+1}^dagger+Z_j^dagger Z_{j+1}), X^3=Z^3=1, periodic, J=1. "
      "At L=2 the bond is counted twice. Critical ferromagnet only; no off-critical/chiral coupling.",
      "Charge q=0,+1,-1 labels the global clock rotation exp(2*pi*i*q/3); p=2*pi*k/L is physical "
      "one-site momentum. The auxiliary XXZ chain has 2*L sites, not L.",
      "all means vacuum plus the selected charged one-hole family, NOT the full spectrum or all excitations. "
      "Neutral excited levels and other charged families are not implemented.",
      "--momentum-index and --charge select charged rows; the q=0,k=0 vacuum remains the gap reference. "
      "--branch ground outputs only that vacuum; --branch charged omits its row, not its calculation.",
      "x_scaled=L*(E-E0)/(2*pi*v), v=3*sqrt(3)/2. This is a finite-size estimator, not an exact CFT dimension. "
      "The k=0 charged level tends to x=2/15; k=1,L-1 tend to x=17/15.",
      "Residual=max|F|/(2*L) for the auxiliary logarithmic equations; default tolerance=32 epsilon. "
      "Unconverged energies and gaps are missing, with exit status 2. No spectral weights are calculated.",
      "bits-N/digits-N use explicit MPFR precision; no MPLAPACK provider is needed.",
      "--threads overlaps independent charged momenta after one vacuum solve. Charge partners share a solve.",
      "Rows stream in momentum/charge order; --no-retain discards delivered rows. One worker is the default.",
      "See docs/potts.md for physical-root selection, state counting and validation; --references for literature."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::count_option(app, "L", a.sites, "Periodic clock sites, L>=2")->required();
  app.add_option("--branch", a.branch, "Selected family, not completeness")
      ->check(CLI::IsMember({"ground", "charged", "all"}))
      ->capture_default_str();
  cli::count_option(app, "--momentum-index", a.momentum_index, "One charged momentum k, 0<=k<L (default: all)");
  app.add_option("--charge", a.charge, "Charged rows, +1 and -1 are degenerate")
      ->check(CLI::IsMember({"1", "+1", "-1", "both"}))
      ->capture_default_str();
  cli::count_option(app, "--max-iterations", a.iterations, "Sweeps per Bethe solve")->capture_default_str();
  cli::count_option(app, "--max-sites", a.max_sites, "Explicit work budget, default 512")->capture_default_str();
  cli::text_option(app, "--tolerance", a.tolerance, "Normalized residual tolerance (default 32 epsilon)")
      ->type_name("REAL");
  cli::precision_option(app, a.precision, cli::ArbitraryPrecision::real);
  cli::threads_option(app, a.threads)->description("Scheduler concurrency limit for independent charged momenta");
  cli::add_data_output_options(app, a.output,
                               {{.name = "levels", .description = "Selected energy levels", .primary = true}});
}

template <uni20::Real Real> int run(Arguments const& a, bethe::SolverOptions<Real> controls, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  bethe::potts::detail::validate_sites(a.sites);
  if (a.sites > a.max_sites) throw std::invalid_argument("sites exceed --max-sites work budget");
  if (a.momentum_index && *a.momentum_index >= a.sites) throw std::invalid_argument("momentum index requires 0<=k<L");
  if (a.branch == "ground" && (a.momentum_index || a.charge != "both"))
    throw std::invalid_argument("--momentum-index/--charge select the charged branch, not the vacuum");
  if (!uni20::isfinite(controls.residual_tolerance) || controls.residual_tolerance <= Real{0})
    throw std::invalid_argument("tolerance must be finite and positive");
  Real const pi = bethe::detail::pi_like(bethe::detail::solver_zero(controls));
  Real const velocity = bethe::potts::velocity<Real>(controls.precision());
  using Optional = std::optional<Real>;
  using State = bethe::potts::State<Real>;
  struct Solve
  {
      std::optional<State> state;
      std::string status;
  };
  auto solve = [sites = a.sites, controls](int charge, std::size_t k) {
    Solve result;
    try
    {
      result.state = charge == 0 ? bethe::potts::ground_state<Real>(sites, controls)
                                 : bethe::potts::charged_one_hole<Real>(sites, charge, k, controls);
      result.status = result.state->converged ? "converged" : "iteration_limit";
    }
    catch (std::runtime_error const&)
    {
      result.status = "precision_limit";
    }
    return result;
  };
  auto const vacuum = context.measure([&] { return solve(0, 0); });
  Optional const reference = vacuum.state && vacuum.state->converged ? Optional(vacuum.state->energy) : Optional{};
  std::size_t const first = a.momentum_index.value_or(0);
  std::size_t const count = a.branch == "ground" ? 0 : a.momentum_index ? 1 : a.sites;
  bool const parallel = a.threads > 1 && count > 1;
  std::size_t const window = std::max(std::size_t{1}, std::min(count, std::size_t{128}));
  int const concurrency = parallel ? static_cast<int>(std::min(window, std::size_t(a.threads))) : 1;
  std::optional<uni20::async::TbbScheduler> scheduler;
  if (parallel) scheduler.emplace(concurrency);
  bethe::ExecutionOptions const execution{.scheduler = scheduler ? &*scheduler : nullptr, .max_in_flight = window};
  cli::RunReport report(context, "Critical three-state Potts (periodic)");
  report.field("sites", "Sites", a.sites)
      .field("hamiltonian", "Hamiltonian", "H=-sum(X+X^dagger+ZZ^dagger+Z^dagger Z); critical, J=1, PBC")
      .field("selection", "State selection", "vacuum and charged one-hole family; not a complete spectrum")
      .field("momentum_convention", "Momentum convention", "p=2*pi*k/L; physical one-clock-site translation")
      .field("charge_convention", "Charge convention", "global clock rotation exp(2*pi*i*q/3); q=0,+1,-1")
      .field("root_convention", "Root convention", "auxiliary 2L-site XXZ at Delta=sqrt(3)/2; twists pi/3,pi")
      .field("residual_convention", "Residual convention", "max|F|/(2L), auxiliary logarithmic equations")
      .precision(a.precision)
      .field("branches", "Branches", a.branch)
      .field("charge", "Selected charged rows", a.charge)
      .field("tolerance", "Residual tolerance", controls.residual_tolerance)
      .field("max_iterations", "Max sweeps per solve", a.iterations)
      .field("max_sites", "Site budget", a.max_sites)
      .field("threads", "Requested concurrency limit", a.threads)
      .field("concurrency", "Effective concurrency limit", concurrency)
      .field("execution_mode", "Execution mode", parallel ? "bounded Async momentum grid" : "direct serial")
      .field("max_in_flight", "Maximum charged momenta in flight", parallel ? window : std::size_t{1})
      .field("cpu_timing_scope", "CPU timing scope", "vacuum solve and ordered level grid including row output")
      .field("velocity", "Velocity", velocity)
      .field("bulk_energy", "Bulk energy density", bethe::potts::bulk_energy_density<Real>(controls.precision()))
      .field("reference_status", "Vacuum status", vacuum.status);
  if (reference) report.field("ground_energy", "Vacuum energy", *reference);
  if (a.momentum_index) report.field("momentum_index", "Requested charged momentum index", *a.momentum_index);
  namespace data = cli::data;
  auto table = data::make_data_table(
      "Selected Potts levels",
      {.retain = a.output.retain ? data::retention::all : data::retention::none, .metadata = report.metadata()},
      cli::column<std::string>("branch"), cli::column<int>("charge"), cli::column<std::size_t>("k"),
      cli::column<Real>("momentum").unit("radians/site"), cli::column<Optional>("energy"), cli::column<Optional>("gap"),
      cli::column<Optional>("x_scaled"), cli::column<Optional>("residual"),
      cli::column<std::optional<std::size_t>>("iterations"), cli::column<std::string>("status"));
  bool const complete = cli::stream_result_table(
      report, context, a.output, "levels", table,
      [&](auto& table) {
        bool success = reference.has_value();
        context.measure([&] {
          auto emit = [&](Solve const& solved, int charge, std::size_t k) {
            Optional energy, gap, scaled, residual;
            std::optional<std::size_t> iterations;
            std::string status = solved.status;
            if (solved.state)
            {
              residual = solved.state->residual_norm;
              iterations = solved.state->iterations;
              if (solved.state->converged)
              {
                energy = solved.state->energy;
                if (reference)
                {
                  gap = *energy - *reference;
                  scaled = Real(a.sites) * *gap / (Real{2} * pi * velocity);
                }
                else
                  status = "reference_unavailable";
              }
            }
            success = success && energy.has_value() && gap.has_value();
            table.append(charge == 0 ? "ground" : "charged-one-hole", charge, k,
                         Real{2} * pi * (Real(k) / Real(a.sites)), energy, gap, scaled, residual, iterations, status);
          };
          if (a.branch != "charged") emit(vacuum, 0, 0);
          auto evaluate = [solve, first](std::size_t i) { return solve(1, first + i); };
          auto collect = [&](std::size_t i, Solve const& state) {
            auto const k = first + i;
            // Charge conjugation, not a second root solve.
            if (a.charge != "-1") emit(state, 1, k);
            if (a.charge == "both" || a.charge == "-1") emit(state, -1, k);
          };
          if (parallel)
            bethe::evaluate_observable_grid(count, evaluate, collect, execution);
          else
            for (std::size_t i = 0; i < count; ++i)
              collect(i, evaluate(i));
        });
        return success;
      },
      "incomplete; unavailable quantities omitted");
  if (!complete) std::cerr << "Some Potts energies or reference gaps are unavailable; inspect row statuses.\n";
  return complete ? 0 : 2;
}
} // namespace
int main(int argc, char** argv)
{
  Arguments a;
  return cli::program_main(
      argc, argv, program_info(), [&](auto& app) { add_options(app, a); },
      [&](auto&) {
        a.output.validate();
        return cli::dispatch_solver_precision<cli::ArbitraryPrecision::real>(
            a.precision, a.tolerance, a.iterations,
            [&](auto controls) { return run(a, std::move(controls), argc, argv); });
      });
}
