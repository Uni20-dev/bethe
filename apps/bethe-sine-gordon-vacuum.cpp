// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/sine_gordon_vacuum.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::sine_gordon;
struct Arguments
{
    std::optional<std::string> mass, length, coupling, tolerance, contour, cutoff;
    std::size_t initial_intervals = 64, max_intervals = 2048, iterations = 10000, cutoffs = 3;
    std::size_t kernel_evaluations = 100000, kernel_levels = 16, fourier_cutoffs = 64;
    std::string precision = "fp64";
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-sine-gordon-vacuum", "Finite-volume sine-Gordon vacuum energy.",
                                bethe::citations::Tool::sine_gordon_vacuum);
  info.examples = {{"bethe-sine-gordon-vacuum --length 1 --p 2", "Repulsive vacuum, soliton mass one"},
                   {"bethe-sine-gordon-vacuum --mass 2 --length 0.5 --p 0.5 --json vacuum.json",
                    "Attractive vacuum with file export"}};
  info.notes = {"M is the soliton mass, L the circumference, u=M*L; velocity=hbar=1.",
                "p=beta^2/(8*pi-beta^2) for kinetic term (1/2)(partial phi)^2; p=1 is free Dirac.",
                "Bulk-subtracted energy E_C=E0-L*e_bulk; Y=L*E_C; c_eff=-6*Y/pi. No absolute bulk energy.",
                "Untwisted zero-topological-charge vacuum only; no excited states or Bethe-Yang approximation.",
                "Tolerance is absolute in Y, default 262144 epsilon. Two contours are independently resolved.",
                "Failed observables are missing, with exit status 2. Table: vacuum. fp128 solves can take minutes.",
                "See docs/sine-gordon.md for numerical controls; --references for literature."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::text_option(app, "--mass", a.mass, "Positive soliton mass M; default 1")->type_name("REAL");
  cli::text_option(app, "--length", a.length, "Positive circumference L")->required()->type_name("REAL");
  cli::text_option(app, "--p", a.coupling, "Positive sine-Gordon coupling p")->required()->type_name("REAL");
  cli::text_option(app, "--tolerance", a.tolerance, "Absolute tolerance in Y; default 262144 epsilon")
      ->type_name("REAL");
  cli::text_option(app, "--contour-shift", a.contour, "eta in (0,pi*min(1,p)/2); default pi*min(1,p)/4")
      ->type_name("REAL");
  cli::text_option(app, "--initial-cutoff", a.cutoff, "Initial rapidity cutoff; default estimated from tolerance and u")
      ->type_name("REAL");
  cli::count_option(app, "--initial-intervals", a.initial_intervals, "Initial even rapidity intervals, >=8")
      ->capture_default_str();
  cli::count_option(app, "--max-intervals", a.max_intervals, "Maximum rapidity intervals, <=16384")
      ->capture_default_str();
  cli::count_option(app, "--max-iterations", a.iterations, "Total nonlinear updates across both contours")
      ->capture_default_str();
  cli::count_option(app, "--max-cutoffs", a.cutoffs, "Rapidity cutoff attempts per contour")->capture_default_str();
  cli::count_option(app, "--max-kernel-evaluations", a.kernel_evaluations,
                    "Fourier quadrature evaluations per kernel table")
      ->capture_default_str();
  cli::count_option(app, "--max-kernel-levels", a.kernel_levels, "Fourier quadrature refinement levels per table")
      ->capture_default_str();
  cli::count_option(app, "--max-fourier-cutoffs", a.fourier_cutoffs, "Fourier cutoff attempts per table")
      ->capture_default_str();
  cli::precision_option(app, a.precision);
  cli::add_data_output_options(app, a.output, true);
}
char const* name(model::VacuumStatus status)
{
  switch (status)
  {
    case model::VacuumStatus::converged:
      return "converged";
    case model::VacuumStatus::kernel_limit:
      return "kernel_limit";
    case model::VacuumStatus::iteration_limit:
      return "iteration_limit";
    case model::VacuumStatus::mesh_limit:
      return "mesh_limit";
    case model::VacuumStatus::cutoff_limit:
      return "cutoff_limit";
    case model::VacuumStatus::contour_limit:
      return "contour_limit";
    case model::VacuumStatus::precision_limit:
      return "precision_limit";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const mass = a.mass ? uni20::parse_real<Real>(*a.mass) : Real{1};
  Real const length = uni20::parse_real<Real>(*a.length), p = uni20::parse_real<Real>(*a.coupling);
  model::VacuumOptions<Real> controls;
  if (a.tolerance) controls.tolerance = uni20::parse_real<Real>(*a.tolerance);
  if (a.contour) controls.contour_shift = uni20::parse_real<Real>(*a.contour);
  if (a.cutoff) controls.initial_cutoff = uni20::parse_real<Real>(*a.cutoff);
  controls.initial_intervals = a.initial_intervals;
  controls.max_intervals = a.max_intervals;
  controls.max_iterations = a.iterations;
  controls.max_cutoffs = a.cutoffs;
  controls.kernel.max_evaluations = a.kernel_evaluations;
  controls.kernel.max_levels = a.kernel_levels;
  controls.kernel.max_cutoffs = a.fourier_cutoffs;
  // Validate and solve before opening files, including --force targets.
  auto const state = context.measure([&] { return model::vacuum_energy(mass, length, p, controls); });
  cli::RunReport report(context, "Sine-Gordon vacuum (bulk-subtracted)");
  report.field("mass", "Soliton mass", mass)
      .field("length", "Circumference", length)
      .field("p", "Coupling p", p)
      .field("scaled_length", "Scaled length u", state.scaled_length)
      .field("units", "Units", "velocity=hbar=1; M is the soliton mass")
      .field("coupling_convention", "Coupling convention", "p=beta^2/(8*pi-beta^2); kinetic term (1/2)(partial phi)^2")
      .field("subtraction", "Energy convention", "bulk-subtracted E_C=E0-L*e_bulk; Y=L*E_C; c_eff=-6*Y/pi")
      .field("calculation", "Calculation", "untwisted zero-topological-charge vacuum NLIE")
      .field("precision", "Precision", a.precision)
      .field("tolerance", "Absolute Y tolerance", controls.tolerance)
      .field("contour_shift", "Contour shift", state.contour_shift)
      .field("verification_contour_shift", "Verification contour shift", state.verification_contour_shift)
      .field("initial_cutoff", "Requested initial rapidity cutoff", controls.initial_cutoff)
      .field("initial_intervals", "Initial intervals", controls.initial_intervals)
      .field("max_intervals", "Max intervals", controls.max_intervals)
      .field("max_iterations", "Max total nonlinear updates", controls.max_iterations)
      .field("max_cutoffs", "Max rapidity cutoffs per contour", controls.max_cutoffs)
      .field("max_kernel_evaluations", "Max evaluations per kernel table", controls.kernel.max_evaluations)
      .field("max_kernel_levels", "Max kernel levels", controls.kernel.max_levels)
      .field("max_fourier_cutoffs", "Max Fourier cutoffs", controls.kernel.max_cutoffs)
      .result(state.converged, name(state.status));
  cli::ResultOutput output(report, a.output, {"vacuum"});
  using cli::column;
  using Optional = std::optional<Real>;
  output.table(
      "vacuum", "Bulk-subtracted vacuum",
      [&](auto& table) {
        table.append(state.casimir_energy, state.scaling_function, state.effective_central_charge,
                     state.nonlinear_residual, state.kernel_error, state.mesh_error, state.cutoff_error,
                     state.contour_error, state.cutoff, state.intervals, state.iterations, state.kernel_evaluations,
                     state.cutoffs, state.converged, std::string(name(state.status)));
      },
      column<Optional>("casimir_energy"), column<Optional>("scaling_function"),
      column<Optional>("effective_central_charge"), column<Real>("nonlinear_residual"), column<Real>("kernel_error"),
      column<Real>("mesh_error"), column<Real>("cutoff_error"), column<Real>("contour_error"), column<Real>("cutoff"),
      column<std::size_t>("intervals"), column<std::size_t>("iterations"), column<std::size_t>("kernel_evaluations"),
      column<std::size_t>("cutoffs"), column<bool>("converged"), column<std::string>("status"));
  output.finish();
  if (!state.converged) std::cerr << "Sine-Gordon vacuum solve incomplete; no observables published.\n";
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
