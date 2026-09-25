// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "lee-yang-common.hpp"
#include "result-output.hpp"
#include <bethe/lee_yang.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::lee_yang;
using cli::lee_yang::Arguments;
using cli::lee_yang::name;
auto program_info()
{
  auto info = cli::program_info("bethe-lee-yang-vacuum", "Periodic scaling Lee-Yang ground-state TBA.",
                                bethe::citations::Tool::lee_yang_vacuum);
  info.notes = {
      "Positive particle mass m and circumference L; r=mL; velocity=hbar=1.",
      "Bulk-subtracted E_C=E0-L*e_bulk; Y=L*E_C. No absolute bulk energy, excited states, boundaries or defects.",
      "Computed c_eff(r)=-6Y/pi approaches 2/5 in the UV; it is not the theory's c=-22/5 (h_min=-1/5).",
      "Tolerance is absolute in Y, default 8192 epsilon. Nonlinear, mesh and cutoff checks are separate.",
      "Failed observables and unavailable diagnostics are missing; incomplete solves exit 2. Table: vacuum.",
      "See docs/lee-yang.md; --references for literature and conventions."};
  info.examples = {
      {"bethe-lee-yang-vacuum --length 1", "Unit-mass periodic ground state"},
      {"bethe-lee-yang-vacuum --mass 2 --length 0.5 --precision fp128 --json vacuum.json", "Native fp128 export"}};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::lee_yang::add_options(app, a, "Absolute Y tolerance; default 8192 epsilon");
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const mass = a.mass ? uni20::parse_real<Real>(*a.mass) : Real{1};
  Real const length = uni20::parse_real<Real>(*a.length);
  auto const options = cli::lee_yang::options_from<Real>(a, model::Options<Real>{});
  // Validation and the solve precede opening any --force targets.
  auto const s = context.measure([&] { return model::ground_state(mass, length, options); });
  cli::RunReport report(context, "Scaling Lee-Yang ground state (bulk-subtracted)");
  report.field("mass", "Particle mass", mass)
      .field("length", "Circumference", length)
      .field("scaled_length", "Scaled length r", s.scaled_length)
      .field("units", "Units", "velocity=hbar=1; r=mL")
      .field("scope", "Scope", "periodic source-free ground-state TBA; no excited states, boundaries or defects")
      .field("energy_convention", "Energy convention",
             "bulk-subtracted E_C=E0-L*e_bulk; Y=L*E_C; no absolute bulk constant")
      .field("effective_charge", "Effective charge convention", "computed c_eff(r)=-6Y/pi; not the CFT central charge")
      .field("theory_c", "CFT central charge (theory)", "-22/5")
      .field("theory_h_min", "Lowest conformal weight (theory)", "-1/5")
      .field("theory_c_eff_uv", "UV effective charge (theory)", "2/5")
      .field("precision", "Precision", a.precision)
      .field("tolerance", "Absolute Y tolerance", options.tolerance)
      .field("initial_cutoff", "Requested initial cutoff", options.initial_cutoff)
      .field("initial_intervals", "Initial intervals", options.initial_intervals)
      .field("max_intervals", "Max intervals", options.max_intervals)
      .field("max_iterations", "Max updates per grid", options.max_iterations)
      .field("max_cutoffs", "Max cutoff trials", options.max_cutoffs)
      .field("max_kernel_products", "Max kernel products", options.max_kernel_products)
      .field("error_convention", "Error convention",
             "errors target absolute Y; residual is in pseudoenergy units; refinement estimates, not certificates")
      .result(s.converged, name(s.status));
  cli::ResultOutput output(report, a.output, {"vacuum"});
  using Optional = std::optional<Real>;
  using cli::column;
  auto finite = [](Real x) -> Optional { return uni20::isfinite(x) ? Optional{x} : std::nullopt; };
  output.table(
      "vacuum", "Bulk-subtracted periodic ground state",
      [&](auto& table) {
        table.append(s.casimir_energy, s.scaling_function, s.effective_central_charge, finite(s.nonlinear_residual),
                     finite(s.nonlinear_error), finite(s.mesh_error), finite(s.cutoff_error),
                     finite(s.direct_tail_bound), s.cutoff, s.intervals, s.iterations, s.cutoffs, s.kernel_products,
                     s.converged, std::string(name(s.status)));
      },
      column<Optional>("casimir_energy"), column<Optional>("scaling_function"),
      column<Optional>("effective_central_charge"), column<Optional>("nonlinear_residual"),
      column<Optional>("nonlinear_error"), column<Optional>("mesh_error"), column<Optional>("cutoff_error"),
      column<Optional>("direct_tail_bound"), column<Real>("cutoff"), column<std::size_t>("intervals"),
      column<std::size_t>("iterations"), column<std::size_t>("cutoffs"), column<std::size_t>("kernel_products"),
      column<bool>("converged"), column<std::string>("status"));
  output.finish();
  if (!s.converged) std::cerr << "Lee-Yang ground-state solve incomplete; no observables published.\n";
  return s.converged ? 0 : 2;
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
