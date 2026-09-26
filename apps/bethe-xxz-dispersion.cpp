// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/xxz_bulk.hpp>
#include <bethe/xxz_dispersion.hpp>

namespace
{
namespace cli = bethe::cli;
struct Arguments
{
    std::string delta, exchange = "1", branch = "all", precision = "fp64";
    std::optional<std::string> momentum, tolerance;
    std::size_t points = 65, evaluations = 200000, levels = 18;
    bool folded = false;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-xxz-dispersion",
                                "Zero-field XXZ bulk energy, spinon dispersion and two-spinon continuum edges.",
                                bethe::citations::Tool::xxz_dispersion);
  info.examples = {
      {"bethe-xxz-dispersion --delta 2 --branch spinon --csv spinons.csv", "Massive topological spinons"},
      {"bethe-xxz-dispersion --delta 2 --branch two-spinon --folded --tsv continuum.tsv", "Two-site MPS continuum"},
      {"bethe-xxz-dispersion --delta 1.001 --precision long-double --json curves.json", "Near-isotropic tiny gap"}};
  info.notes = {
      "H=J sum (SxSx+SySy+Delta SzSz), spin 1/2, J>0, Delta>-1, infinite chain, zero field and temperature.",
      "Spinon: p in [0,pi], Sz=+/-1/2. At Delta>1 it connects different Neel vacua. "
      "A fixed-vacuum local excitation needs an even number of spinons; the two-spinon sectors have Sz=-1,0,+1.",
      "Two-spinon: Q in [0,2*pi], unfolded by default (p1+p2=Q). --folded takes the union with its pi-shifted "
      "image, as appropriate to a two-site unit cell. cell_momentum=2*p modulo 2*pi is the two-site translation phase.",
      "Continuum edges are kinematic energy bounds, not spectral intensity, a full excited spectrum or finite-size "
      "levels.",
      "Bulk tolerance defaults to 256 epsilon in units J*max(1,abs(Delta)). Numerical failures leave values empty "
      "(JSON null) and exit 2. A massive gap below the normal scalar range is rejected, never reported as zero.",
      "See docs/xxz-dispersion.md for conventions and limits. Use --references for literature and applicability."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  app.add_option("--delta", a.delta, "Anisotropy Delta>-1")->required()->type_name("REAL");
  app.add_option("--exchange", a.exchange, "Antiferromagnetic exchange J>0")->capture_default_str()->type_name("REAL");
  app.add_option("--branch", a.branch, "Output branch selection")
      ->check(CLI::IsMember({"spinon", "two-spinon", "all"}))
      ->capture_default_str();
  app.add_flag("--folded", a.folded, "Fold the two-spinon continuum modulo pi (two-site unit cell)");
  auto* points =
      cli::count_option(app, "--points", a.points, "Uniform samples per branch, 2..1000000")->capture_default_str();
  cli::text_option(app, "--momentum", a.momentum, "One momentum instead of a grid (radians/site)")
      ->type_name("REAL")
      ->excludes(points);
  cli::precision_option(app, a.precision);
  cli::text_option(app, "--tolerance", a.tolerance, "Bulk-energy tolerance (default 256 epsilon)")->type_name("REAL");
  cli::count_option(app, "--max-evaluations", a.evaluations, "Bulk integral/series evaluation budget")
      ->capture_default_str();
  cli::count_option(app, "--max-levels", a.levels, "Bulk quadrature refinement budget, <=30")->capture_default_str();
  cli::add_data_output_options(app, a.output, true);
}
char const* name(bethe::xxz::BulkStatus status)
{
  switch (status)
  {
    case bethe::xxz::BulkStatus::converged:
      return "converged";
    case bethe::xxz::BulkStatus::evaluation_limit:
      return "evaluation_limit";
    case bethe::xxz::BulkStatus::mesh_limit:
      return "mesh_limit";
    case bethe::xxz::BulkStatus::precision_limit:
      return "precision_limit";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const delta = uni20::parse_real<Real>(a.delta), exchange = uni20::parse_real<Real>(a.exchange);
  Real const pi = bethe::detail::pi<Real>();
  std::optional<Real> momentum;
  if (a.momentum) momentum = uni20::parse_real<Real>(*a.momentum);
  if (momentum && (!uni20::isfinite(*momentum) || *momentum < Real{0} ||
                   *momentum > (a.branch == "two-spinon" ? Real{2} * pi : pi)))
    throw std::invalid_argument("momentum outside the selected branch range");
  bethe::xxz::BulkOptions<Real> controls;
  controls.max_evaluations = a.evaluations;
  controls.max_levels = a.levels;
  if (a.tolerance) controls.tolerance = uni20::parse_real<Real>(*a.tolerance);
  auto const bulk = context.measure([&] { return bethe::xxz::bulk_energy_density(delta, exchange, controls); });
  std::optional<bethe::xxz::SpinonDispersion<Real>> band;
  std::string band_status = "converged", failure;
  try
  {
    band = context.measure([&] { return bethe::xxz::SpinonDispersion<Real>(delta, exchange); });
  }
  catch (std::underflow_error const& e)
  {
    band_status = "precision_limit";
    failure = e.what();
  }
  catch (std::overflow_error const& e)
  {
    band_status = "precision_limit";
    failure = e.what();
  }
  cli::RunReport report(context, "XXZ thermodynamic spinons");
  report.field("delta", "Delta", delta)
      .field("exchange", "Exchange J", exchange)
      .field("hamiltonian", "Hamiltonian", "H=J sum (SxSx+SySy+Delta SzSz); spin 1/2")
      .field("energy_reference", "Energy reference",
             "excitation energy above the zero-field thermodynamic ground state")
      .field("spinon_sector", "Spinon sector", "Sz=+/-1/2; for Delta>1 a domain wall between different Neel vacua")
      .field("two_spinon_sectors", "Two-spinon sectors", "Sz=-1,0,+1; same asymptotic vacuum")
      .field("momentum_convention", "Continuum momentum",
             a.folded ? "Q modulo pi; union of both translation branches" : "p1+p2=Q modulo 2*pi; unfolded")
      .field("cell_momentum", "Two-site translation phase", "2*p modulo 2*pi; p is in radians/site")
      .field("branches", "Branches", a.branch)
      .field("precision", "Precision", a.precision)
      .field("bulk_tolerance", "Bulk tolerance", controls.tolerance)
      .field("max_evaluations", "Max bulk evaluations", controls.max_evaluations)
      .field("max_levels", "Max bulk levels", controls.max_levels)
      .field("evaluations", "Bulk evaluations", bulk.evaluations)
      .field("bulk_status", "Bulk status", name(bulk.status))
      .field("band_status", "Spinon status", band_status);
  if (bulk.energy)
    report.field("bulk_energy", "Bulk energy/site", *bulk.energy)
        .field("bulk_error", "Bulk energy error estimate", *bulk.error);
  if (band)
    report.field("spinon_gap", "Single-spinon gap", band->gap())
        .field("band_maximum", "Spinon maximum energy", band->maximum_energy());
  if (!failure.empty()) report.field("failure", "Numerical failure", failure);
  if (momentum)
    report.field("momentum", "Requested momentum", *momentum);
  else
    report.field("points", "Points per branch", a.points);

  namespace data = cli::data;
  using Optional = std::optional<Real>;
  auto table = data::make_data_table(
      "XXZ spinon lines and continuum edges",
      {.retain = a.output.retain ? data::retention::all : data::retention::none, .metadata = report.metadata()},
      cli::column<std::string>("branch"), cli::column<Real>("p").unit("radians/site"), cli::column<Real>("p_over_pi"),
      cli::column<Real>("cell_momentum").unit("radians/cell"), cli::column<Optional>("energy"),
      cli::column<Optional>("lower"), cli::column<Optional>("upper"), cli::column<std::string>("status"));
  cli::DataOutput output(a.output, {"dispersion"});
  bool complete = bulk.converged && band.has_value();
  try
  {
    output.attach(table, "dispersion");
    for (std::string const branch : {"spinon", "two-spinon"})
    {
      if (a.branch != "all" && a.branch != branch) continue;
      Real const end = branch == "spinon" ? pi : Real{2} * pi;
      std::size_t const count = momentum ? 1 : a.points;
      for (std::size_t i = 0; i < count; ++i)
      {
        Real const p = momentum ? *momentum : i == count - 1 ? end : end * (Real(i) / Real(count - 1));
        Optional energy, lower, upper;
        std::string status = band_status;
        if (band)
        {
          try
          {
            context.measure([&] {
              if (branch == "spinon")
                energy = band->energy(p);
              else
              {
                auto const edges =
                    band->continuum(p, a.folded ? bethe::SpinonMomentum::folded : bethe::SpinonMomentum::unfolded);
                lower = edges.lower;
                upper = edges.upper;
              }
            });
          }
          catch (std::overflow_error const&)
          {
            status = "precision_limit";
            complete = false;
          }
        }
        Real const cell = Real{2} * (p < pi ? p : p < Real{2} * pi ? p - pi : Real{0});
        table.append(branch, p, p / pi, cell, energy, lower, upper, status);
      }
    }
    report.result(complete, complete ? "converged" : "incomplete; unavailable quantities omitted");
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
  if (!complete)
    std::cerr << "Some XXZ thermodynamic quantities are unavailable; inspect statuses and numerical controls.\n";
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
        if (a.folded && a.branch == "spinon")
          throw std::invalid_argument("--folded changes only the two-spinon continuum");
        return cli::dispatch_precision(a.precision, [&]<typename Real>() { return run<Real>(a, argc, argv); });
      });
}
