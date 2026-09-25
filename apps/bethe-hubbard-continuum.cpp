// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/hubbard_continuum.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::hubbard::thermo;
struct Arguments
{
    std::string u, precision = "fp64";
    std::optional<std::string> momentum, tolerance;
    std::size_t points = 33, evaluations = 1000000, levels = 12, iterations = 160;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-hubbard-continuum", "Half-filled Hubbard two-spinon continuum edges.",
                                bethe::citations::Tool::hubbard_continuum);
  info.examples = {{"bethe-hubbard-continuum --u 4 --points 33 --csv edges.csv", "Uniform total-momentum grid"},
                   {"bethe-hubbard-continuum --u 4 --momentum -1 --precision long-double", "One momentum"}};
  info.notes = {"Infinite half-filled chain, zero magnetic field, U>0, hopping t=1. Two-spinon family only.",
                "Total momentum P in [-pi,pi], modulo 2*pi. Constituent spinon momenta are in [0,pi].",
                "DeltaN=0: symmetric/unshifted Hamiltonian and Fermi energy references coincide.",
                "The upper edge bounds two spinons, not all excitations. No spectral weights or finite-ring levels.",
                "Default tolerance 256 epsilon. Energy quadrature errors exclude momentum-inversion errors.",
                "Each momentum shares one evaluation budget across both edges. Failed energies are missing (exit 2).",
                "Table: two_spinon. Exports stream during computation; --no-retain discards delivered rows.",
                "See docs/hubbard-continuum.md; --references for literature."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  app.add_option("--u", a.u, "Repulsive interaction U>0 (t=1)")->required()->type_name("REAL");
  auto* points =
      cli::count_option(app, "--points", a.points, "Uniform grid on [-pi,pi], 2..1000000")->capture_default_str();
  cli::text_option(app, "--momentum", a.momentum, "One total momentum in radians")->type_name("REAL")->excludes(points);
  cli::text_option(app, "--tolerance", a.tolerance, "Relative constituent tolerance; default 256 epsilon");
  cli::count_option(app, "--max-evaluations", a.evaluations, "Quadrature budget per pair of edges")
      ->capture_default_str();
  cli::count_option(app, "--max-levels", a.levels, "Quadrature refinement levels, <=24")->capture_default_str();
  cli::count_option(app, "--max-iterations", a.iterations, "Momentum inversion updates per constituent")
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
    case model::Status::quadrature_limit:
      return "quadrature_limit";
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
  Real const u = uni20::parse_real<Real>(a.u), pi = Real{4} * std::atan(Real{1});
  std::optional<Real> momentum;
  if (a.momentum) momentum = uni20::parse_real<Real>(*a.momentum);
  model::Options<Real> controls;
  if (a.tolerance) controls.relative_tolerance = uni20::parse_real<Real>(*a.tolerance);
  controls.max_evaluations = a.evaluations;
  controls.max_levels = a.levels;
  controls.max_iterations = a.iterations;
  // Validate before opening any output files, including --force targets.
  model::detail::validate(u, controls, model::Convention::symmetric);
  if (momentum && (!uni20::isfinite(*momentum) || std::abs(*momentum) > pi))
    throw std::invalid_argument("total momentum must be in [-pi,pi]");
  std::size_t const count = momentum ? 1 : a.points;
  cli::RunReport report(context, "Hubbard two-spinon continuum");
  report.field("interaction", "U (t=1)", u)
      .field("background", "Background", "half filling, zero field, infinite chain")
      .field("family", "Scattering family", "two spinons only; no spectral weights")
      .field("energy_convention", "Energy convention",
             "DeltaN=0; symmetric/unshifted and Hamiltonian/Fermi references coincide")
      .field("momentum_convention", "Momentum convention", "one-site radians, total P modulo 2*pi")
      .field("precision", "Precision", a.precision)
      .field("points", "Points", count)
      .field("tolerance", "Relative tolerance", controls.relative_tolerance)
      .field("max_evaluations", "Max evaluations per edge pair", controls.max_evaluations)
      .field("max_levels", "Max quadrature levels", controls.max_levels)
      .field("max_iterations", "Max iterations per constituent", controls.max_iterations);
  if (momentum) report.field("momentum", "Requested momentum", *momentum);
  namespace data = cli::data;
  using Optional = std::optional<Real>;
  using cli::column;
  auto table = data::make_data_table(
      "Two-spinon edges",
      {.retain = a.output.retain ? data::retention::all : data::retention::none, .metadata = report.metadata()},
      column<Real>("momentum"), column<Real>("momentum_over_pi"), column<Optional>("lower_energy"),
      column<Optional>("upper_energy"), column<Real>("lower_p1"), column<Real>("lower_p2"), column<Real>("upper_p1"),
      column<Real>("upper_p2"), column<Real>("lower_energy_quad_error"), column<Real>("upper_energy_quad_error"),
      column<Real>("lower_momentum_error"), column<Real>("upper_momentum_error"), column<std::size_t>("evaluations"),
      column<std::size_t>("iterations"), column<bool>("converged"), column<std::string>("status"));
  cli::DataOutput output(a.output, {"two_spinon"});
  bool complete = true;
  try
  {
    output.attach(table, "two_spinon");
    for (std::size_t i = 0; i < count; ++i)
    {
      Real const p = momentum ? *momentum : pi * (Real{2} * Real(i) / Real(count - 1) - Real{1});
      auto const s = context.measure([&] { return model::two_spinon_continuum(u, p, controls); });
      complete = complete && s.converged;
      table.append(p, p / pi, s.lower_energy, s.upper_energy, s.lower_momenta[0], s.lower_momenta[1],
                   s.upper_momenta[0], s.upper_momenta[1], s.lower_energy_quad_error, s.upper_energy_quad_error,
                   s.lower_momentum_error, s.upper_momentum_error, s.evaluations, s.iterations, s.converged,
                   std::string(name(s.status)));
    }
    report.result(complete, complete ? "converged" : "incomplete; failed edges omitted");
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
  if (!complete) std::cerr << "Some continuum points did not converge; failed edges are unavailable.\n";
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
        return cli::dispatch_precision(a.precision, [&]<typename Real>() { return run<Real>(a, argc, argv); });
      });
}
