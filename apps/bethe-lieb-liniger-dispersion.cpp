// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/lieb_liniger_thermo.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::lieb_liniger::thermo;
struct Arguments
{
    std::string c, density = "1", branch = "all", precision = "fp64";
    std::optional<std::string> momentum, pmax, tolerance;
    std::size_t points = 33, initial_nodes = 16, max_nodes = 256, background_iterations = 128, iterations = 160;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-lieb-liniger-dispersion",
                                "Repulsive Lieb-Liniger bulk ground state and type-I/type-II excitation curves.",
                                bethe::citations::Tool::lieb_liniger_dispersion);
  info.examples = {
      {"bethe-lieb-liniger-dispersion --c 4 --points 33 --csv curves.csv", "Both elementary branches"},
      {"bethe-lieb-liniger-dispersion --c 1 --density 0.5 --branch type-i --p-max 8 --precision long-double",
       "Particle branch over a larger physical momentum interval"}};
  info.notes = {
      "H=-sum d_i^2+2c sum delta(x_i-x_j), hbar=2m=1; c>0, n=N/L>0, zero temperature.",
      "Type I: p>=0; type II: 0<=p<=2*pi*n. p is physical momentum (inverse length), not rapidity or p modulo 2*pi. "
      "Energies are fixed-N excitation gaps. Negative-momentum branches follow by parity.",
      "The default grid ends at 2*pi*n on both branches. --p-max changes only the type-I grid. "
      "These are elementary lines, not matrix elements, spectral weights, or a complete multiparticle spectrum.",
      "Tolerance defaults to 4096 epsilon in the selected precision. Independent mesh refinement and momentum "
      "inversion are required; very weak coupling can exceed the node budget. fp128 requires MPLAPACK.",
      "Failed energies are omitted (JSON null, CSV/TSV empty) and exit status is 2. Exports stream during computation; "
      "--no-retain discards delivered rows. See docs/lieb-liniger-thermo.md and docs/output.md. "
      "Use --references for literature and applicability."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  app.add_option("--c", a.c, "Repulsive coupling c>0")->required()->type_name("REAL");
  app.add_option("--density", a.density, "Particle density n>0")->capture_default_str()->type_name("REAL");
  app.add_option("--branch", a.branch, "Elementary branch selection")
      ->check(CLI::IsMember({"type-i", "type-ii", "all"}))
      ->capture_default_str();
  auto* points = cli::count_option(app, "--points", a.points, "Uniform momentum samples per branch, 2..1000000")
                     ->capture_default_str();
  auto* momentum = cli::text_option(app, "--momentum", a.momentum, "One physical momentum instead of a grid")
                       ->type_name("REAL")
                       ->excludes(points);
  cli::text_option(app, "--p-max", a.pmax, "Type-I grid upper momentum (default 2*pi*n)")
      ->type_name("REAL")
      ->excludes(momentum);
  cli::precision_option(app, a.precision);
  cli::text_option(app, "--tolerance", a.tolerance, "Relative numerical tolerance (default 4096 epsilon)")
      ->type_name("REAL");
  cli::count_option(app, "--initial-nodes", a.initial_nodes, "Initial full-interval Gauss-Legendre order")
      ->capture_default_str();
  cli::count_option(app, "--max-nodes", a.max_nodes, "Mesh budget, <=512")->capture_default_str();
  cli::count_option(app, "--max-background-iterations", a.background_iterations, "Density iterations across meshes")
      ->capture_default_str();
  cli::count_option(app, "--max-iterations", a.iterations, "Momentum inversion updates per point, both meshes")
      ->capture_default_str();
  cli::add_data_output_options(app, a.output, true);
}
char const* name(model::Branch branch) { return branch == model::Branch::type_i ? "type-i" : "type-ii"; }
char const* name(model::Status status)
{
  switch (status)
  {
    case model::Status::converged:
      return "converged";
    case model::Status::mesh_limit:
      return "mesh_limit";
    case model::Status::density_limit:
      return "density_limit";
    case model::Status::momentum_limit:
      return "momentum_limit";
    case model::Status::precision_limit:
      return "precision_limit";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const c = uni20::parse_real<Real>(a.c), density = uni20::parse_real<Real>(a.density);
  Real const pi = Real{4} * std::atan(Real{1}), end = Real{2} * pi * density;
  Real const pmax = a.pmax ? uni20::parse_real<Real>(*a.pmax) : end;
  std::optional<Real> momentum;
  if (a.momentum) momentum = uni20::parse_real<Real>(*a.momentum);
  if (!uni20::isfinite(pmax) || !(pmax > Real{0}) || !uni20::isfinite(end) || !(end > Real{0}))
    throw std::invalid_argument("grid endpoints and density must be finite and positive");
  if (momentum && (!uni20::isfinite(*momentum) || *momentum < Real{0} || (a.branch != "type-i" && *momentum > end)))
    throw std::invalid_argument("momentum outside selected branch range");
  model::Options<Real> controls;
  controls.initial_nodes = a.initial_nodes;
  controls.max_nodes = a.max_nodes;
  controls.max_iterations = a.background_iterations;
  controls.max_momentum_iterations = a.iterations;
  if (a.tolerance) controls.tolerance = uni20::parse_real<Real>(*a.tolerance);
  auto solver = context.measure([&] { return model::Solver<Real>(c, density, controls); });
  auto const& bg = solver.background();
  cli::RunReport report(context, "Lieb-Liniger thermodynamic dispersions");
  report.field("interaction", "c", c)
      .field("density", "Density N/L", density)
      .field("hamiltonian", "Hamiltonian", "H=-sum d_i^2+2c sum delta(x_i-x_j); hbar=2m=1")
      .field("energy_reference", "Energy reference", "fixed-N excitation gap above the bulk ground state")
      .field("momentum_convention", "Momentum convention", "physical inverse length; not reduced modulo 2*pi")
      .field("branches", "Branches", a.branch)
      .field("precision", "Precision", a.precision)
      .field("tolerance", "Relative tolerance", controls.tolerance)
      .field("initial_nodes", "Initial nodes", controls.initial_nodes)
      .field("max_nodes", "Max nodes", controls.max_nodes)
      .field("max_background_iterations", "Max background iterations", controls.max_iterations)
      .field("max_iterations", "Max momentum iterations", controls.max_momentum_iterations)
      .field("background_status", "Background status", name(bg.status))
      .field("nodes", "Nodes", bg.nodes)
      .field("background_iterations", "Background iterations", bg.iterations)
      .field("mesh_error", "Background relative mesh error", bg.mesh_error)
      .field("density_error", "Relative density residual", bg.density_error)
      .field("fermi_rapidity", "Fermi rapidity Q", bg.fermi_rapidity)
      .field("energy_per_length", "Ground energy/length", bg.energy_per_length)
      .field("chemical_potential", "Chemical potential", bg.chemical_potential);
  if (momentum)
    report.field("momentum", "Requested momentum", *momentum);
  else
  {
    report.field("points", "Points per branch", a.points);
    if (a.branch != "type-ii") report.field("type_i_max", "Type-I grid maximum", pmax);
  }
  namespace data = cli::data;
  using Optional = std::optional<Real>;
  auto table = data::make_data_table(
      "Lieb-Liniger elementary branches",
      {.retain = a.output.retain ? data::retention::all : data::retention::none, .metadata = report.metadata()},
      cli::column<std::string>("branch"), cli::column<Real>("p").unit("1/length"), cli::column<Real>("p_over_pi_n"),
      cli::column<Optional>("energy").unit("1/length^2"), cli::column<Optional>("rapidity").unit("1/length"),
      cli::column<Optional>("edge_distance").unit("1/length").description("Distance from the nearest Fermi edge"),
      cli::column<Optional>("energy_error").unit("1/length^2"),
      cli::column<Optional>("momentum_error").unit("1/length"), cli::column<std::size_t>("iterations"),
      cli::column<std::string>("status"));
  cli::DataOutput output(a.output, {"dispersion"});
  bool complete = bg.converged;
  try
  {
    output.attach(table, "dispersion");
    for (auto branch : {model::Branch::type_i, model::Branch::type_ii})
    {
      if (a.branch != "all" && a.branch != name(branch)) continue;
      std::size_t const count = momentum ? 1 : a.points;
      Real const high = branch == model::Branch::type_i ? pmax : end;
      for (std::size_t i = 0; i < count; ++i)
      {
        Real const p = momentum ? *momentum : i == count - 1 ? high : high * (Real(i) / Real(count - 1));
        auto const point = context.measure([&] { return solver.at_momentum(branch, p); });
        complete = complete && point.converged;
        table.append(std::string(name(branch)), p, p / (pi * density), point.energy, point.rapidity,
                     point.edge_distance, point.converged ? Optional(point.energy_error) : std::nullopt,
                     point.converged ? Optional(point.momentum_error) : std::nullopt, point.iterations,
                     std::string(name(point.status)));
      }
    }
    report.result(complete, complete ? "converged" : "incomplete; failed energies omitted");
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
  if (!complete) std::cerr << "Some Lieb-Liniger results did not converge; failed energies are unavailable.\n";
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
        if (!a.momentum && (a.points < 2 || a.points > 1000000))
          throw std::invalid_argument("require 2<=points<=1000000");
        if (a.pmax && a.branch == "type-ii") throw std::invalid_argument("--p-max changes only the type-I grid");
        return cli::dispatch_precision(a.precision, [&]<typename Real>() { return run<Real>(a, argc, argv); });
      });
}
