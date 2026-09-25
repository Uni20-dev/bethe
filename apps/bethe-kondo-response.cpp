// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/kondo.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::kondo;
struct Arguments
{
    std::optional<std::string> field, scale, tolerance;
    std::size_t terms = 256, lobes = 256, evaluations = 200000, levels = 12;
    std::string precision = "fp64";
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-kondo-response", "Universal zero-temperature Kondo impurity response.",
                                bethe::citations::Tool::kondo_response);
  info.examples = {{"bethe-kondo-response --field 2 --scale 1", "Uniform-field impurity response"},
                   {"bethe-kondo-response --field -0.1 --scale 1 --precision long-double --json response.json",
                    "Negative field with native-precision export"}};
  info.notes = {"Spin-1/2, single-channel antiferromagnetic isotropic Kondo scaling limit, T=0.",
                "b is the full Zeeman splitting: H_field=-b*(S_imp^z+S_host^z), with equal g factors.",
                "Positive scale T_B=2*T1 of Barcza et al.; chi0=1/(sqrt(2*pi*e)*T_B).",
                "Delta E_imp=E_imp(b)-E_imp(0); M_imp=-d Delta E_imp/db includes the host response change.",
                "No finite-band absolute energy, bare J/bandwidth matching, or finite-size excited states.",
                "Tolerance is absolute in M_imp and Delta E_imp/|b|; default 1048576 epsilon.",
                "Failed observables are missing, with exit status 2. Table: response.",
                "See docs/kondo.md; --references for literature."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::text_option(app, "--field", a.field, "Full Zeeman splitting b, either sign")->required()->type_name("REAL");
  cli::text_option(app, "--scale", a.scale, "Positive universal energy scale T_B")->required()->type_name("REAL");
  cli::text_option(app, "--tolerance", a.tolerance, "Absolute tolerance in M and Delta E/|b|")->type_name("REAL");
  cli::count_option(app, "--max-series-terms", a.terms, "Low-field series budget, <=4096")->capture_default_str();
  cli::count_option(app, "--max-lobes", a.lobes, "High-field Fourier lobe budget, <=4096")->capture_default_str();
  cli::count_option(app, "--max-evaluations", a.evaluations, "Total high-field quadrature evaluations")
      ->capture_default_str();
  cli::count_option(app, "--max-quadrature-levels", a.levels, "Quadrature refinement levels per lobe")
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
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const field = uni20::parse_real<Real>(*a.field), scale = uni20::parse_real<Real>(*a.scale);
  model::Options<Real> controls;
  if (a.tolerance) controls.tolerance = uni20::parse_real<Real>(*a.tolerance);
  controls.max_series_terms = a.terms;
  controls.max_lobes = a.lobes;
  controls.max_evaluations = a.evaluations;
  controls.max_quadrature_levels = a.levels;
  // Validate and solve before opening files, including --force targets.
  auto const state = context.measure([&] { return model::ground_response(field, scale, controls); });
  cli::RunReport report(context, "Kondo universal impurity response");
  report.field("field", "Full Zeeman splitting", field)
      .field("scale", "Universal scale T_B", scale)
      .field("scale_convention", "Scale convention", "T_B=2*T1 (Barcza et al.); chi0=1/(sqrt(2*pi*e)*T_B)")
      .field("field_convention", "Field convention", "H_field=-b*(S_imp^z+S_host^z); equal g factors")
      .field("energy_convention", "Energy convention", "Delta E_imp=E_imp(b)-E_imp(0); E_imp=E_with-E_clean_host")
      .field("response_convention", "Response convention",
             "Impurity-induced total magnetization, including host change")
      .field("calculation", "Calculation", "T=0; isotropic antiferromagnetic spin-1/2 single-channel scaling limit")
      .field("precision", "Precision", a.precision)
      .field("tolerance", "Absolute M and Delta E/|b| tolerance", controls.tolerance)
      .field("max_series_terms", "Max series terms", controls.max_series_terms)
      .field("max_lobes", "Max Fourier lobes", controls.max_lobes)
      .field("max_evaluations", "Max quadrature evaluations", controls.max_evaluations)
      .field("max_quadrature_levels", "Max quadrature levels", controls.max_quadrature_levels)
      .result(state.converged, name(state.status));
  cli::ResultOutput output(report, a.output, {"response"});
  using cli::column;
  using Optional = std::optional<Real>;
  output.table(
      "response", "Universal impurity response",
      [&](auto& table) {
        table.append(state.energy_change, state.magnetization, state.zero_field_susceptibility,
                     state.magnetization_error, state.scaled_energy_error, state.series_terms, state.lobes,
                     state.evaluations, state.converged, std::string(name(state.status)));
      },
      column<Optional>("energy_change"), column<Optional>("magnetization"),
      column<Optional>("zero_field_susceptibility"), column<Real>("magnetization_error"),
      column<Real>("scaled_energy_error"), column<std::size_t>("series_terms"), column<std::size_t>("lobes"),
      column<std::size_t>("evaluations"), column<bool>("converged"), column<std::string>("status"));
  output.finish();
  if (!state.converged) std::cerr << "Kondo response incomplete; no observables published.\n";
  return state.converged ? 0 : 2;
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
