// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/lieb_liniger_thermal.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::lieb_liniger::thermal;
struct Arguments
{
    std::string c, temperature, precision = "fp64";
    std::optional<std::string> mu, density, end, tolerance, density_tolerance, cutoff;
    std::size_t points = 17, initial_nodes = 16, max_nodes = 256, iterations = 512, cutoffs = 4, evaluations = 128;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-lieb-liniger-thermal", "Repulsive Lieb-Liniger Yang-Yang thermodynamics.",
                                bethe::citations::Tool::lieb_liniger_thermal);
  info.examples = {{"bethe-lieb-liniger-thermal --c 4 --temperature 1 --mu -1", "Grand-canonical equilibrium"},
                   {"bethe-lieb-liniger-thermal --c 4 --temperature 0.5 --temperature-end 2 --points 9 --density 1 "
                    "--max-nodes 512 --csv thermal.csv",
                    "Temperature scan at fixed density"}};
  info.notes = {
      "H=-sum d_i^2+2c sum delta(x_i-x_j), hbar=2m=k_B=1; c>0 and T>0. Bulk equilibrium, not a finite-size spectrum.",
      "Choose exactly one of --mu (grand canonical) and --density (solve for mu). A temperature scan keeps that "
      "quantity and c fixed. Pressure, energy and entropy are per length; grand potential/length is -pressure.",
      "Tolerance defaults to 16384 epsilon; density tolerance to 65536 epsilon. Nonlinear, mesh, cutoff and "
      "density inversion budgets are independent. Weak coupling or low temperature can exceed the mesh budget.",
      "Failed observables are omitted (JSON null, CSV/TSV empty), with exit status 2. fp128 requires MPLAPACK. "
      "Exports stream during computation; --no-retain discards delivered rows. See docs/lieb-liniger-thermal.md. "
      "Use --references for literature."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  app.add_option("--c", a.c, "Repulsive interaction c>0")->required()->type_name("REAL");
  app.add_option("--temperature", a.temperature, "Temperature T>0, or scan start")->required()->type_name("REAL");
  auto* ensemble = app.add_option_group("Ensemble", "Choose chemical potential or density")->require_option(1, 1);
  cli::text_option(*ensemble, "--mu", a.mu, "Chemical potential, any finite real")->type_name("REAL");
  cli::text_option(*ensemble, "--density", a.density, "Fixed particle density n>0")->type_name("REAL");
  auto* end =
      cli::text_option(app, "--temperature-end", a.end, "Positive scan endpoint (inclusive)")->type_name("REAL");
  cli::count_option(app, "--points", a.points, "Temperature samples, 2..1000000")->capture_default_str()->needs(end);
  cli::precision_option(app, a.precision);
  cli::text_option(app, "--tolerance", a.tolerance, "Inner relative tolerance (default 16384 epsilon)");
  cli::text_option(app, "--density-tolerance", a.density_tolerance,
                   "Relative density tolerance (default 65536 epsilon)");
  cli::text_option(app, "--initial-cutoff", a.cutoff, "Initial physical rapidity cutoff (otherwise automatic)");
  cli::count_option(app, "--initial-nodes", a.initial_nodes, "Initial positive-half quadrature order")
      ->capture_default_str();
  cli::count_option(app, "--max-nodes", a.max_nodes, "Mesh budget, <=512")->capture_default_str();
  cli::count_option(app, "--max-iterations", a.iterations, "Accepted Newton updates per equilibrium solve")
      ->capture_default_str();
  cli::count_option(app, "--max-cutoffs", a.cutoffs, "Cutoff refinements per equilibrium solve, <=32")
      ->capture_default_str();
  cli::count_option(app, "--max-evaluations", a.evaluations, "Equilibrium solves per fixed-density point")
      ->capture_default_str();
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
    case model::Status::mesh_limit:
      return "mesh_limit";
    case model::Status::cutoff_limit:
      return "cutoff_limit";
    case model::Status::density_limit:
      return "density_limit";
    case model::Status::precision_limit:
      return "precision_limit";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const c = uni20::parse_real<Real>(a.c), start = uni20::parse_real<Real>(a.temperature);
  Real const end = a.end ? uni20::parse_real<Real>(*a.end) : start;
  Real const fixed = uni20::parse_real<Real>(a.mu ? *a.mu : *a.density);
  model::DensityOptions<Real> controls;
  auto& inner = controls.equilibrium;
  inner.initial_nodes = a.initial_nodes;
  inner.max_nodes = a.max_nodes;
  inner.max_iterations = a.iterations;
  inner.max_cutoffs = a.cutoffs;
  controls.max_evaluations = a.evaluations;
  if (a.tolerance) inner.tolerance = uni20::parse_real<Real>(*a.tolerance);
  if (a.density_tolerance) controls.tolerance = uni20::parse_real<Real>(*a.density_tolerance);
  if (a.cutoff) inner.initial_cutoff = uni20::parse_real<Real>(*a.cutoff);
  // Validate before opening any output file, including when --force is used.
  model::detail::validate(c, start, a.mu ? fixed : Real{0}, inner);
  if (!(end > Real{0}) || !uni20::isfinite(end) || !uni20::isfinite(fixed) || (a.density && !(fixed > Real{0})) ||
      !(controls.tolerance > Real{0}) || !(controls.tolerance < Real{1}) || !uni20::isfinite(controls.tolerance))
    throw std::invalid_argument("require finite T>0, density>0 and 0<density tolerance<1");
  std::size_t const count = a.end ? a.points : 1;
  cli::RunReport report(context, "Lieb-Liniger thermal equilibrium");
  report.field("interaction", "c", c)
      .field("ensemble", "Ensemble", a.mu ? "fixed chemical potential" : "fixed density")
      .field("hamiltonian", "Hamiltonian", "H=-sum d_i^2+2c sum delta(x_i-x_j); hbar=2m=k_B=1")
      .field("temperature", "Initial temperature", start)
      .field("temperature_end", "Final temperature", end)
      .field("points", "Points", count)
      .field("precision", "Precision", a.precision)
      .field("tolerance", "Requested inner tolerance", inner.tolerance)
      .field("initial_nodes", "Initial positive-half nodes", inner.initial_nodes)
      .field("max_nodes", "Max nodes", inner.max_nodes)
      .field("max_iterations", "Max Newton updates per solve", inner.max_iterations)
      .field("max_cutoffs", "Max cutoffs per solve", inner.max_cutoffs);
  if (inner.initial_cutoff) report.field("initial_cutoff", "Initial physical cutoff", *inner.initial_cutoff);
  if (a.mu)
    report.field("chemical_potential", "Fixed chemical potential", fixed);
  else
    report.field("density", "Requested density", fixed)
        .field("density_tolerance", "Relative density tolerance", controls.tolerance)
        .field("effective_tolerance", "Effective inner tolerance",
               std::min(inner.tolerance, controls.tolerance / Real{8}))
        .field("max_evaluations", "Max equilibrium evaluations per point", controls.max_evaluations);
  namespace data = cli::data;
  using Optional = std::optional<Real>;
  using Count = std::optional<std::size_t>;
  auto table = data::make_data_table(
      "Lieb-Liniger thermodynamics",
      {.retain = a.output.retain ? data::retention::all : data::retention::none, .metadata = report.metadata()},
      cli::column<Real>("temperature").unit("1/length^2"),
      cli::column<Optional>("chemical_potential").unit("1/length^2"), cli::column<Optional>("density").unit("1/length"),
      cli::column<Optional>("pressure").unit("1/length^3"),
      cli::column<Optional>("energy_per_length").unit("1/length^3"),
      cli::column<Optional>("entropy_per_length").unit("1/length"), cli::column<Count>("nodes"),
      cli::column<Optional>("cutoff").unit("1/length"), cli::column<Optional>("mesh_error"),
      cli::column<Optional>("cutoff_error"), cli::column<Optional>("density_error"),
      cli::column<std::size_t>("evaluations"), cli::column<std::size_t>("iterations"),
      cli::column<std::string>("status"));
  cli::DataOutput output(a.output, {"thermodynamics"});
  bool complete = true;
  try
  {
    output.attach(table, "thermodynamics");
    for (std::size_t i = 0; i < count; ++i)
    {
      Real const fraction = count == 1 ? Real{0} : Real(i) / Real(count - 1);
      Real const temperature = i == count - 1 ? end : (Real{1} - fraction) * start + fraction * end;
      std::optional<model::State<Real>> state;
      Optional density_error;
      std::size_t evaluations = 1, iterations = 0;
      model::Status status;
      if (a.mu)
      {
        state = context.measure([&] { return model::equilibrium(c, temperature, fixed, inner); });
        status = state->status;
        iterations = state->iterations;
      }
      else
      {
        auto result = context.measure([&] { return model::at_density(c, temperature, fixed, controls); });
        status = result.status;
        evaluations = result.evaluations;
        iterations = result.iterations;
        if (result.converged) density_error = result.density_error;
        state = std::move(result.state);
      }
      bool const converged = status == model::Status::converged;
      complete = complete && converged;
      table.append(temperature,
                   a.mu    ? Optional(fixed)
                   : state ? Optional(state->chemical_potential)
                           : std::nullopt,
                   state ? state->density : std::nullopt, state ? state->pressure : std::nullopt,
                   state ? state->energy_per_length : std::nullopt, state ? state->entropy_per_length : std::nullopt,
                   state ? Count(state->nodes) : std::nullopt, state ? Optional(state->cutoff) : std::nullopt,
                   state ? Optional(state->mesh_error) : std::nullopt,
                   state ? Optional(state->cutoff_error) : std::nullopt, density_error, evaluations, iterations,
                   std::string(name(status)));
    }
    report.result(complete, complete ? "converged" : "incomplete; failed thermodynamics omitted");
    auto const summary = report.finish();
    output.overview(report.overview(summary));
    output.finish(table, summary);
    output.finish_document();
  }
  catch (...)
  {
    data::table_metadata aborted{{"Status", "aborted"}};
    if (!context.finished()) try
      {
        aborted = cli::run_summary(context.finish(uni20::run_outcome::failed));
      }
      catch (...)
      {}
    output.abort(table, std::move(aborted));
    throw;
  }
  if (!complete) std::cerr << "Some thermal points did not converge; failed observables are unavailable.\n";
  return complete ? 0 : 2;
}
} // namespace
int main(int argc, char** argv)
{
  Arguments a;
  return cli::program_main(
      argc, argv, program_info(), [&](auto& app) { add_options(app, a); },
      [&](auto& app) {
        a.output.validate();
        if (a.end && (a.points < 2 || a.points > 1000000)) throw std::invalid_argument("require 2<=points<=1000000");
        if (a.mu && (a.density_tolerance || app.count("--max-evaluations")))
          throw std::invalid_argument("density inversion controls require --density");
        return cli::dispatch_precision(a.precision, [&]<typename Real>() { return run<Real>(a, argc, argv); });
      });
}
