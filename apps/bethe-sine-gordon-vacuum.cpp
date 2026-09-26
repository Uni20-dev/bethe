// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "sine-gordon-options.hpp"

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::sine_gordon;
using Arguments = cli::sine_gordon::Arguments;
using cli::sine_gordon::name;
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
  cli::sine_gordon::add_options(app, a, "Positive sine-Gordon coupling p",
                                "Absolute tolerance in Y; default 262144 epsilon");
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const mass = a.mass ? uni20::parse_real<Real>(*a.mass) : Real{1};
  Real const length = uni20::parse_real<Real>(*a.length), p = uni20::parse_real<Real>(*a.coupling);
  model::VacuumOptions<Real> controls;
  cli::sine_gordon::configure<Real>(a, controls);
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
      .result(state.converged, name(state.status));
  cli::sine_gordon::report_controls<Real>(report, controls, state);
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
