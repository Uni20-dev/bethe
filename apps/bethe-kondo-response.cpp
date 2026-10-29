// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "execution-options.hpp"
#include "result-output.hpp"
#include "solver-precision.hpp"
#include <bethe/kondo.hpp>
#include <bethe/observable_grid.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::kondo;
struct Arguments
{
    std::optional<std::string> field, end, scale, tolerance;
    std::optional<std::size_t> terms, lobes, evaluations;
    std::size_t levels = 12;
    std::size_t points = 17;
    int threads = 1;
    std::string precision = "fp64";
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-kondo-response", "Universal zero-temperature Kondo impurity response.",
                                bethe::citations::Tool::kondo_response);
  info.examples = {
      {"bethe-kondo-response --field 2 --scale 1", "Uniform-field impurity response"},
      {"bethe-kondo-response --field -4 --field-end 4 --points 129 --scale 1 --threads 4 --csv response.csv",
       "Ordered field-response curve"},
      {"bethe-kondo-response --field -0.1 --scale 1 --precision long-double --json response.json",
       "Negative field with native-precision export"}};
  info.notes = {"Spin-1/2, single-channel antiferromagnetic isotropic Kondo scaling limit, T=0.",
                "b is the full Zeeman splitting: H_field=-b*(S_imp^z+S_host^z), with equal g factors.",
                "Positive scale T_B=2*T1 of Barcza et al.; chi0=1/(sqrt(2*pi*e)*T_B).",
                "Delta E_imp=E_imp(b)-E_imp(0); M_imp=-d Delta E_imp/db includes the host response change.",
                "No finite-band absolute energy, bare J/bandwidth matching, or finite-size excited states.",
                "Tolerance is absolute in M_imp and Delta E_imp/|b|; default 1048576 epsilon.",
                "Field scans include both endpoints, in either direction; --threads overlaps independent fields.",
                "Exports stream in field order; --no-retain discards delivered rows. Single points stay serial.",
                "Failed rows retain diagnostics with missing observables, and cause exit status 2.",
                "See docs/kondo.md; --references for literature."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::text_option(app, "--field", a.field, "Full Zeeman splitting b, or scan start, either sign")
      ->required()
      ->type_name("REAL");
  auto* end = cli::text_option(app, "--field-end", a.end, "Inclusive field scan endpoint")->type_name("REAL");
  cli::count_option(app, "--points", a.points, "Field samples, 2..1000000")->capture_default_str()->needs(end);
  cli::threads_option(app, a.threads)->description("Scheduler concurrency limit for independent fields (1 = serial)");
  cli::text_option(app, "--scale", a.scale, "Positive universal energy scale T_B")->required()->type_name("REAL");
  cli::text_option(app, "--tolerance", a.tolerance, "Absolute tolerance in M and Delta E/|b|")->type_name("REAL");
  cli::count_option(app, "--max-series-terms", a.terms,
                    "Low-field series budget, <=4096; default 256 native, 1024 MPFR")
      ->capture_default_str();
  cli::count_option(app, "--max-lobes", a.lobes,
                    "High-field Fourier lobe budget, <=4096; default 256 native, 1024 MPFR")
      ->capture_default_str();
  cli::count_option(app, "--max-evaluations", a.evaluations,
                    "Total high-field evaluations; default 200000 native, 2000000 MPFR")
      ->capture_default_str();
  cli::count_option(app, "--max-quadrature-levels", a.levels, "Quadrature refinement levels per lobe")
      ->capture_default_str();
  cli::precision_option(app, a.precision, cli::ArbitraryPrecision::real);
  cli::add_data_output_options(
      app, a.output, {{.name = "response", .description = "Impurity response, one row per field", .primary = true}});
}
// Inclusive linear grid, avoiding overflow of end-start for opposite signs.
// Exact endpoints and a monotone clamp also protect very short/equal intervals.
template <class Real> Real field_at(Real start, Real end, std::size_t i, std::size_t count)
{
  if (i == 0 || start == end) return start;
  if (i == count - 1) return end;
  Real const fraction = uni20::scalar_like(start, i) / Real(count - 1);
  Real const value = (start < Real{0}) != (end < Real{0}) ? (Real{1} - fraction) * start + fraction * end
                                                          : start + fraction * (end - start);
  return std::clamp(value, std::min(start, end), std::max(start, end));
}
char const* name(model::Status status)
{
  switch (status)
  {
    case model::Status::converged:
      return "converged";
    case model::Status::series_limit:
      return "series_limit";
    case model::Status::quadrature_limit:
      return "quadrature_limit";
    case model::Status::tail_limit:
      return "tail_limit";
    case model::Status::precision_limit:
      return "precision_limit";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& a, model::Options<Real> controls, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const field = cli::parse_solver_real<Real>(*a.field, controls),
             scale = cli::parse_solver_real<Real>(*a.scale, controls);
  Real const end = a.end ? cli::parse_solver_real<Real>(*a.end, controls) : field;
  if (a.tolerance) controls.tolerance = cli::parse_solver_real<Real>(*a.tolerance, controls);
  if (a.terms) controls.max_series_terms = *a.terms;
  if (a.lobes) controls.max_lobes = *a.lobes;
  if (a.evaluations) controls.max_evaluations = *a.evaluations;
  controls.max_quadrature_levels = a.levels;
  // Validate before opening files, including --force targets, without a solve.
  model::detail::validate(field, scale, controls);
  model::detail::validate(end, scale, controls);
  std::size_t const count = a.end ? a.points : 1;
  bool const parallel = a.threads > 1 && count > 1 && (field != Real{0} || end != Real{0});
  std::size_t const window = std::min(count, std::size_t{128});
  int const concurrency = parallel ? static_cast<int>(std::min(window, std::size_t(a.threads))) : 1;
  std::optional<uni20::async::TbbScheduler> scheduler;
  if (parallel) scheduler.emplace(concurrency);
  bethe::ExecutionOptions const execution{.scheduler = scheduler ? &*scheduler : nullptr, .max_in_flight = window};
  cli::RunReport report(context, "Kondo universal impurity response");
  report.field("field", "Full Zeeman splitting", field)
      .field("field_end", "Final field", end)
      .field("points", "Points", count)
      .field("threads", "Requested concurrency limit", a.threads)
      .field("concurrency", "Effective concurrency limit", concurrency)
      .field("execution_mode", "Execution mode", parallel ? "bounded Async field grid" : "direct serial")
      .field("max_in_flight", "Maximum fields in flight", parallel ? window : std::size_t{1})
      .field("cpu_timing_scope", "CPU timing scope", "field grid including ordered row output")
      .field("scale", "Universal scale T_B", scale)
      .field("scale_convention", "Scale convention", "T_B=2*T1 (Barcza et al.); chi0=1/(sqrt(2*pi*e)*T_B)")
      .field("field_convention", "Field convention", "H_field=-b*(S_imp^z+S_host^z); equal g factors")
      .field("energy_convention", "Energy convention", "Delta E_imp=E_imp(b)-E_imp(0); E_imp=E_with-E_clean_host")
      .field("response_convention", "Response convention",
             "Impurity-induced total magnetization, including host change")
      .field("calculation", "Calculation", "T=0; isotropic antiferromagnetic spin-1/2 single-channel scaling limit")
      .precision(a.precision)
      .field("tolerance", "Absolute M and Delta E/|b| tolerance", controls.tolerance)
      .field("max_series_terms", "Max series terms", controls.max_series_terms)
      .field("max_lobes", "Max Fourier lobes", controls.max_lobes)
      .field("max_evaluations", "Max quadrature evaluations", controls.max_evaluations)
      .field("max_quadrature_levels", "Max quadrature levels", controls.max_quadrature_levels);
  using cli::column;
  using Optional = std::optional<Real>;
  auto table = cli::data::make_data_table(
      "Universal impurity response",
      {.retain = a.output.retain ? cli::data::retention::all : cli::data::retention::none,
       .metadata = report.metadata()},
      column<Real>("field"), column<Optional>("energy_change"), column<Optional>("magnetization"),
      column<Optional>("zero_field_susceptibility"), column<Real>("magnetization_error"),
      column<Real>("scaled_energy_error"), column<std::size_t>("series_terms"), column<std::size_t>("lobes"),
      column<std::size_t>("evaluations"), column<bool>("converged"), column<std::string>("status"));
  std::size_t failures = 0;
  bool const complete = cli::stream_result_table(
      report, context, a.output, "response", table,
      [&](auto& table) {
        context.measure([&] {
          auto evaluate = [field, end, scale, controls, count](std::size_t i) {
            return model::ground_response(field_at(field, end, i, count), scale, controls);
          };
          auto collect = [&](std::size_t, model::State<Real> const& state) {
            failures += !state.converged;
            table.append(state.field, state.energy_change, state.magnetization, state.zero_field_susceptibility,
                         state.magnetization_error, state.scaled_energy_error, state.series_terms, state.lobes,
                         state.evaluations, state.converged, std::string(name(state.status)));
          };
          if (parallel)
            bethe::evaluate_observable_grid(count, evaluate, collect, execution);
          else
            for (std::size_t i = 0; i < count; ++i)
              collect(i, evaluate(i));
        });
        return failures == 0;
      },
      "incomplete; failed observables omitted");
  if (!complete) std::cerr << failures << " Kondo response points incomplete; failed observables are unavailable.\n";
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
        if (a.end && (a.points < 2 || a.points > 1000000)) throw std::invalid_argument("require 2<=points<=1000000");
        return cli::dispatch_options_precision<cli::ArbitraryPrecision::real, model::Options>(
            a.precision, [&](auto controls) { return run(a, std::move(controls), argc, argv); });
      });
}
