// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/hubbard_charge_continuum.hpp>
#include <bethe/hubbard_continuum.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::hubbard::thermo;
struct Arguments
{
    std::string u, precision = "fp64";
    std::string channel = "two-spinon", convention = "symmetric", reference = "hamiltonian";
    std::optional<std::string> momentum, tolerance;
    std::optional<std::string> search_tolerance, position_tolerance;
    std::size_t points = 33, evaluations = 1000000, levels = 12, iterations = 160;
    std::size_t total_evaluations = 200000000, search_evaluations = 20000, search_iterations = 256;
    std::size_t initial_intervals = 16, max_intervals = 128;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-hubbard-continuum", "Half-filled Hubbard two-particle continuum edges.",
                                bethe::citations::Tool::hubbard_continuum);
  info.examples = {{"bethe-hubbard-continuum --u 4 --points 33 --csv edges.csv", "Uniform total-momentum grid"},
                   {"bethe-hubbard-continuum --u 4 --momentum -1 --precision long-double", "One momentum"}};
  info.notes = {
      "Infinite half-filled chain, zero magnetic field, U>0, hopping t=1.",
      "Total momentum P in [-pi,pi], modulo 2*pi. Constituent spinon momenta are in [0,pi].",
      "Unshifted Hamiltonian energies add U*DeltaN/2; the half-filled Fermi reference equals symmetric energies.",
      "Edges bound only the selected family, not all excitations. No spectral weights or finite-ring levels.",
      "Default tolerance 256 epsilon. Energy quadrature errors exclude momentum-inversion errors.",
      "Charge channels use numerical extrema searches, not certified global bounds; native fp128 can be slow.",
      "Failed energies are missing (exit 2). Charge search tolerance applies to E/max(1,U).",
      "Tables: two_spinon or charge_continuum. Exports stream; --no-retain discards delivered rows.",
      "See docs/hubbard-continuum.md; --references for literature."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  app.add_option("--u", a.u, "Repulsive interaction U>0 (t=1)")->required()->type_name("REAL");
  app.add_option("--channel", a.channel, "Scattering family")
      ->check(CLI::IsMember({"two-spinon", "spinon-holon", "spinon-antiholon", "holon-antiholon"}))
      ->capture_default_str();
  app.add_option("--convention", a.convention, "Interaction convention")
      ->check(CLI::IsMember({"symmetric", "unshifted"}))
      ->capture_default_str();
  app.add_option("--reference", a.reference, "Hamiltonian DeltaE or Fermi DeltaE-mu*DeltaN")
      ->check(CLI::IsMember({"hamiltonian", "fermi"}))
      ->capture_default_str();
  auto* points =
      cli::count_option(app, "--points", a.points, "Uniform grid on [-pi,pi], 2..1000000")->capture_default_str();
  cli::text_option(app, "--momentum", a.momentum, "One total momentum in radians")->type_name("REAL")->excludes(points);
  cli::text_option(app, "--tolerance", a.tolerance, "Relative constituent tolerance; default 256 epsilon");
  cli::count_option(app, "--max-evaluations", a.evaluations,
                    "Quadrature budget per constituent (two spinons: per edge pair)")
      ->capture_default_str();
  auto* search = app.add_option_group("Charge-channel search");
  cli::text_option(*search, "--search-tolerance", a.search_tolerance, "Tolerance on E/max(1,U); default 65536 epsilon");
  cli::text_option(*search, "--position-tolerance", a.position_tolerance,
                   "Relative bracket tolerance; default sqrt(epsilon)");
  cli::count_option(*search, "--max-total-evaluations", a.total_evaluations, "Total quadrature budget per momentum")
      ->capture_default_str();
  cli::count_option(*search, "--max-search-evaluations", a.search_evaluations, "Objective budget per momentum")
      ->capture_default_str();
  cli::count_option(*search, "--max-search-iterations", a.search_iterations, "Updates per local extremum bracket")
      ->capture_default_str();
  cli::count_option(*search, "--initial-intervals", a.initial_intervals, "Initial search mesh intervals, >=4")
      ->capture_default_str();
  cli::count_option(*search, "--max-intervals", a.max_intervals, "Maximum mesh intervals, <=4096")
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
// Both families use the same incremental export and run-lifecycle handling.
template <typename Table, typename Append>
int stream_rows(Arguments const& a, uni20::run_context& context, cli::RunReport& report, Table& table,
                std::string const& id, std::size_t count, Append append)
{
  cli::DataOutput output(a.output, {id});
  bool complete = true;
  try
  {
    output.attach(table, id);
    for (std::size_t i = 0; i < count; ++i)
      complete = append(i) && complete;
    report.result(complete, complete ? "converged" : "incomplete; failed edges omitted");
    auto const summary = report.finish();
    output.overview(report.overview(summary));
    output.finish(table, summary);
    output.finish_document();
  }
  catch (...)
  {
    cli::data::table_metadata aborted{{"Status", "aborted"}};
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
char const* name(bethe::detail::ExtremaStatus status)
{
  using S = bethe::detail::ExtremaStatus;
  switch (status)
  {
    case S::converged:
      return "converged";
    case S::objective_failure:
      return "objective_failure";
    case S::evaluation_limit:
      return "evaluation_limit";
    case S::iteration_limit:
      return "iteration_limit";
    case S::mesh_limit:
      return "mesh_limit";
    case S::precision_limit:
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
  auto momentum_at = [&](std::size_t i) {
    return momentum ? *momentum : pi * (Real{2} * Real(i) / Real(count - 1) - Real{1});
  };
  if (a.channel != "two-spinon")
  {
    auto const channel = a.channel == "spinon-holon"       ? model::ChargeChannel::spinon_holon
                         : a.channel == "spinon-antiholon" ? model::ChargeChannel::spinon_antiholon
                                                           : model::ChargeChannel::holon_antiholon;
    auto const convention = a.convention == "symmetric" || a.reference == "fermi" ? model::Convention::symmetric
                                                                                  : model::Convention::unshifted;
    model::ChargeContinuumOptions<Real> options;
    options.constituent = controls;
    options.max_quadrature_evaluations = a.total_evaluations;
    if (a.search_tolerance) options.search.tolerance = uni20::parse_real<Real>(*a.search_tolerance);
    if (a.position_tolerance) options.search.position_tolerance = uni20::parse_real<Real>(*a.position_tolerance);
    options.search.initial_intervals = a.initial_intervals;
    options.search.max_intervals = a.max_intervals;
    options.search.max_evaluations = a.search_evaluations;
    options.search.max_iterations = a.search_iterations;
    // Reuse library validation without evaluating an objective or opening files.
    auto validation = options.search;
    validation.max_evaluations = 0;
    (void)bethe::detail::bounded_extrema([](Real) -> std::optional<bethe::detail::ObjectiveSample<Real>> { return {}; },
                                         channel == model::ChargeChannel::holon_antiholon ? -pi : Real{0}, pi,
                                         validation);
    int const delta = a.channel == "spinon-holon" ? -1 : a.channel == "spinon-antiholon" ? 1 : 0;
    auto const spin = delta == 0 ? uni20::half_int{0} : uni20::from_twice(1);
    cli::RunReport report(context, "Hubbard charge-containing continuum");
    report.field("interaction", "U (t=1)", u)
        .field("background", "Background", "half filling, zero field, infinite chain")
        .field("channel", "Scattering family", a.channel)
        .field("delta_particles", "DeltaN", delta)
        .field("spin", "Spin", spin)
        .field("convention", "Energy convention", a.convention)
        .field("reference", "Energy reference", a.reference)
        .field("energy_offset", "Offset from symmetric energy",
               convention == model::Convention::unshifted ? u / Real{2} * Real(delta) : Real{0})
        .field("momentum_convention", "Momentum convention", "one-site radians, total P modulo 2*pi")
        .field("precision", "Precision", a.precision)
        .field("points", "Points", count)
        .field("tolerance", "Relative constituent tolerance", controls.relative_tolerance)
        .field("search_tolerance", "Search tolerance on E/max(1,U)", options.search.tolerance)
        .field("position_tolerance", "Relative bracket tolerance", options.search.position_tolerance)
        .field("max_evaluations", "Max evaluations per constituent", controls.max_evaluations)
        .field("max_levels", "Max quadrature levels", controls.max_levels)
        .field("max_iterations", "Max iterations per constituent", controls.max_iterations)
        .field("max_total_evaluations", "Max total quadrature evaluations", a.total_evaluations)
        .field("max_search_evaluations", "Max objective evaluations", a.search_evaluations)
        .field("max_search_iterations", "Max updates per local bracket", a.search_iterations)
        .field("initial_intervals", "Initial mesh intervals", a.initial_intervals)
        .field("max_intervals", "Max mesh intervals", a.max_intervals)
        .field("error_contract", "Error estimates",
               "Heuristic extrema; vertical errors exclude momentum-inversion propagation")
        .field("scope", "Scope", "Selected two-particle family only; no spectral weights or global certification");
    if (momentum) report.field("momentum", "Requested momentum", *momentum);
    using Optional = std::optional<Real>;
    using cli::column;
    auto table = cli::data::make_data_table(
        "Charge-containing edges",
        {.retain = a.output.retain ? cli::data::retention::all : cli::data::retention::none,
         .metadata = report.metadata()},
        column<Real>("momentum"), column<Real>("momentum_over_pi"), column<Optional>("lower_energy"),
        column<Optional>("upper_energy"), column<Optional>("lower_symmetric_energy"),
        column<Optional>("upper_symmetric_energy"), column<Optional>("lower_p1"), column<Optional>("lower_p2"),
        column<Optional>("upper_p1"), column<Optional>("upper_p2"), column<Optional>("lower_error"),
        column<Optional>("upper_error"), column<Optional>("lower_momentum_error"),
        column<Optional>("upper_momentum_error"), column<std::size_t>("quadrature_evaluations"),
        column<std::size_t>("constituent_iterations"), column<std::size_t>("objective_evaluations"),
        column<std::size_t>("search_iterations"), column<std::size_t>("intervals"), column<std::size_t>("meshes"),
        column<bool>("converged"), column<std::string>("search_status"), column<std::string>("constituent_status"));
    return stream_rows(a, context, report, table, "charge_continuum", count, [&](std::size_t i) {
      Real const p = momentum_at(i);
      auto const s = context.measure([&] { return model::charge_continuum(channel, u, p, convention, options); });
      auto available = [&](Real value) -> Optional { return s.converged ? Optional(value) : std::nullopt; };
      table.append(p, p / pi, s.lower_energy, s.upper_energy, s.lower_symmetric_energy, s.upper_symmetric_energy,
                   available(s.lower_momenta[0]), available(s.lower_momenta[1]), available(s.upper_momenta[0]),
                   available(s.upper_momenta[1]), available(s.lower_error), available(s.upper_error),
                   available(s.lower_momentum_error), available(s.upper_momentum_error), s.quadrature_evaluations,
                   s.constituent_iterations, s.objective_evaluations, s.search_iterations, s.intervals, s.meshes,
                   s.converged, std::string(name(s.search_status)), std::string(name(s.constituent_status)));
      return s.converged;
    });
  }
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
  return stream_rows(a, context, report, table, "two_spinon", count, [&](std::size_t i) {
    Real const p = momentum_at(i);
    auto const s = context.measure([&] { return model::two_spinon_continuum(u, p, controls); });
    table.append(p, p / pi, s.lower_energy, s.upper_energy, s.lower_momenta[0], s.lower_momenta[1], s.upper_momenta[0],
                 s.upper_momenta[1], s.lower_energy_quad_error, s.upper_energy_quad_error, s.lower_momentum_error,
                 s.upper_momentum_error, s.evaluations, s.iterations, s.converged, std::string(name(s.status)));
    return s.converged;
  });
}
} // namespace
int main(int argc, char** argv)
{
  Arguments a;
  return cli::program_main(
      argc, argv, program_info(), [&](auto& app) { add_options(app, a); },
      [&](auto& app) {
        a.output.validate();
        if (a.channel == "two-spinon")
          for (auto const* option :
               {"--search-tolerance", "--position-tolerance", "--max-total-evaluations", "--max-search-evaluations",
                "--max-search-iterations", "--initial-intervals", "--max-intervals"})
            if (app.count(option))
              throw std::invalid_argument(std::string(option) + " requires a charge-containing channel");
        if (!a.momentum && (a.points < 2 || a.points > 1000000))
          throw std::invalid_argument("require 2<=points<=1000000");
        return cli::dispatch_precision(a.precision, [&]<typename Real>() { return run<Real>(a, argc, argv); });
      });
}
