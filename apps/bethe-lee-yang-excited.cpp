// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "lee-yang-common.hpp"
#include "result-output.hpp"
#include <bethe/lee_yang_excited.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::lee_yang;
using cli::lee_yang::name;
struct Arguments : cli::lee_yang::Arguments
{
    std::size_t root_iterations = 256;
    Arguments() { products = 1000000000; }
};
auto program_info()
{
  auto info = cli::program_info("bethe-lee-yang-excited", "Periodic Lee-Yang spin-zero one-particle TBA.",
                                bethe::citations::Tool::lee_yang_excited);
  info.notes = {
      "Positive mass m and circumference L, with 5<=r=mL<=30; velocity=hbar=1.",
      "Regular one-particle branch only: no UV continuation, moving particles, higher levels or boundaries.",
      "E1_C and E0_C are bulk-subtracted levels; the excitation gap is E1_C-E0_C, not E1_C.",
      "Tolerance is absolute in each Y=L*E_C, default 65536 epsilon; work limits apply separately to each state.",
      "Tables: levels (both states), source (one-particle quantization), gap (requires both states converged).",
      "Failures exit 2 with missing affected observables. See docs/lee-yang.md; --references for literature."};
  info.examples = {{"bethe-lee-yang-excited --length 5", "Unit-mass zero-momentum one-particle gap"},
                   {"bethe-lee-yang-excited --length 10 --precision fp128 --json levels.json", "Native fp128 export"}};
  return info;
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const mass = a.mass ? uni20::parse_real<Real>(*a.mass) : Real{1};
  Real const length = uni20::parse_real<Real>(*a.length);
  auto options = cli::lee_yang::options_from<Real>(a, model::OneParticleOptions<Real>{});
  options.max_root_iterations = a.root_iterations;
  // Both validations/solves precede opening output files, including --force.
  auto const excited = context.measure([&] { return model::one_particle(mass, length, options); });
  auto const vacuum = context.measure(
      [&] { return model::ground_state(mass, length, static_cast<model::Options<Real> const&>(options)); });
  using Optional = std::optional<Real>;
  auto finite = [](Real x) -> Optional { return uni20::isfinite(x) ? Optional{x} : std::nullopt; };
  Optional gap, scaled_gap, gap_error;
  bool converged = excited.converged && vacuum.converged;
  if (converged)
  {
    // Subtract dimensionless levels before dividing by L; each state has
    // separately passed mesh/cutoff verification. Include subtraction rounding.
    Real const y = *excited.scaling_function - *vacuum.scaling_function;
    Real const rounding = Real{8} * uni20::numeric_limits<Real>::epsilon() *
                          (std::abs(*excited.scaling_function) + std::abs(*vacuum.scaling_function));
    scaled_gap = finite(y);
    gap = finite(y / length);
    gap_error = finite((excited.cutoff_error + vacuum.cutoff_error + rounding) / length);
    converged = scaled_gap && gap && gap_error;
    if (!converged)
    {
      scaled_gap.reset();
      gap.reset();
      gap_error.reset();
    }
  }
  std::string const status = converged ? "converged"
                                       : (!excited.converged  ? std::string("excited_") + name(excited.status)
                                          : !vacuum.converged ? std::string("vacuum_") + name(vacuum.status)
                                                              : "gap_precision_limit");
  cli::RunReport report(context, "Scaling Lee-Yang one-particle level and gap");
  report.field("mass", "Particle mass", mass)
      .field("length", "Circumference", length)
      .field("scaled_length", "Scaled length r", excited.scaled_length)
      .field("units", "Units", "velocity=hbar=1; r=mL")
      .field("scope", "Scope", "periodic spin-zero one-particle regular branch; 5<=mL<=30; no UV continuation")
      .field("energy_convention", "Energy convention", "bulk-subtracted E_j_C=E_j-L*e_bulk; Y_j=L*E_j_C")
      .field("gap_convention", "Gap convention", "gap=E1_C-E0_C; scaled_gap=L*gap; requires both states converged")
      .field("precision", "Precision", a.precision)
      .field("tolerance", "Absolute Y tolerance per state", options.tolerance)
      .field("initial_cutoff", "Requested initial cutoff", options.initial_cutoff)
      .field("initial_intervals", "Initial intervals", options.initial_intervals)
      .field("max_intervals", "Max intervals", options.max_intervals)
      .field("max_iterations", "Max updates per source evaluation", options.max_iterations)
      .field("max_root_iterations", "Max source iterations per grid", options.max_root_iterations)
      .field("max_cutoffs", "Max cutoff trials per state", options.max_cutoffs)
      .field("max_kernel_products", "Max kernel products per state", options.max_kernel_products)
      .field("error_convention", "Error convention",
             "level errors in Y; source_error in Y1; gap_error in energy units sums both cutoff errors and rounding; "
             "estimates, not certificates")
      .result(converged, status);
  cli::ResultOutput output(report, a.output, {"levels", "source", "gap"});
  using cli::column;
  output.table(
      "levels", "Bulk-subtracted finite-volume levels",
      [&](auto& table) {
        auto append = [&](std::string label, auto const& s) {
          table.append(label, s.casimir_energy, s.scaling_function, finite(s.nonlinear_residual),
                       finite(s.nonlinear_error), finite(s.mesh_error), finite(s.cutoff_error),
                       finite(s.direct_tail_bound), s.cutoff, s.intervals, s.iterations, s.cutoffs, s.kernel_products,
                       s.converged, std::string(name(s.status)));
        };
        append("vacuum", vacuum);
        append("one_particle", excited);
      },
      column<std::string>("state"), column<Optional>("casimir_energy"), column<Optional>("scaling_function"),
      column<Optional>("nonlinear_residual"), column<Optional>("nonlinear_error"), column<Optional>("mesh_error"),
      column<Optional>("cutoff_error"), column<Optional>("direct_tail_bound"), column<Real>("cutoff"),
      column<std::size_t>("intervals"), column<std::size_t>("iterations"), column<std::size_t>("cutoffs"),
      column<std::size_t>("kernel_products"), column<bool>("converged"), column<std::string>("status"));
  output.table(
      "source", "One-particle source quantization",
      [&](auto& table) {
        table.append(excited.beta, excited.pole_displacement, finite(excited.quantization_residual),
                     finite(excited.source_error), excited.root_iterations, excited.converged,
                     std::string(name(excited.status)));
      },
      column<Optional>("beta"), column<Optional>("pole_displacement"), column<Optional>("quantization_residual"),
      column<Optional>("source_error"), column<std::size_t>("root_iterations"), column<bool>("converged"),
      column<std::string>("status"));
  output.table(
      "gap", "Vacuum-relative zero-momentum excitation",
      [&](auto& table) { table.append(gap, scaled_gap, gap_error, converged, status); }, column<Optional>("gap"),
      column<Optional>("scaled_gap"), column<Optional>("gap_error"), column<bool>("converged"),
      column<std::string>("status"));
  output.finish();
  if (!converged) std::cerr << "Lee-Yang excitation incomplete; gap unavailable (" << status << ").\n";
  return converged ? 0 : 2;
}
} // namespace
int main(int argc, char** argv)
{
  Arguments a;
  return cli::program_main(
      argc, argv, program_info(),
      [&](auto& app) {
        cli::lee_yang::add_options(app, a, "Absolute Y tolerance per state; default 65536 epsilon");
        cli::count_option(app, "--max-root-iterations", a.root_iterations, "Source-root updates per grid")
            ->capture_default_str();
      },
      [&](auto&) {
        a.output.validate();
        return cli::dispatch_precision(a.precision, [&]<typename Real>() { return run<Real>(a, argc, argv); });
      });
}
