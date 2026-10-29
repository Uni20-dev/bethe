// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "execution-options.hpp"
#include "independent-solves.hpp"
#include "sine-gordon-options.hpp"
#include <bethe/sine_gordon_excited.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::sine_gordon;
using cli::sine_gordon::name;
struct Arguments : cli::sine_gordon::Arguments
{
    Arguments() { iterations = model::TwoSolitonOptions<double>{}.max_iterations; }
    std::string number = "0.5";
    int charge = 2;
    int threads = 1;
    std::size_t root_iterations = 1024;
};
auto program_info()
{
  auto info = cli::program_info("bethe-sine-gordon-excited", "Exact finite-volume sine-Gordon two-soliton levels.",
                                bethe::citations::Tool::sine_gordon_excited);
  info.examples = {{"bethe-sine-gordon-excited --length 1 --p 2", "Lowest same-charge pair and vacuum-relative gap"},
                   {"bethe-sine-gordon-excited --length 0.1 --p 2 --number 1.5 --json levels.json",
                    "Selected descendant towards the ultraviolet limit"}};
  info.notes = {
      "--threads overlaps the independent vacuum and excited solves (at most two); each solve stays sequential.",
      "M>0 is the soliton mass, L>0 the circumference, u=M*L; velocity=hbar=1, p>=1.",
      "Selected symmetric two-hole states only: delta=0, I=0.5 or 1.5, roots +/-H, total momentum zero.",
      "Winding charge +/-2 describes two solitons or two antisolitons, NOT a neutral soliton-antisoliton pair.",
      "Exact continuum NLIE, including sea/wrapping terms; not the asymptotic Bethe-Yang equation.",
      "E_C=E-L*e_bulk is bulk-subtracted. Gap subtracts a separately converged vacuum; no absolute bulk energy.",
      "Y=L*E_C; scaled_gap=L*gap/(2*pi), tending to Delta+ + Delta- in the ultraviolet.",
      "UV weights: Delta+=Delta-=(p+1)/(4*p)+I-0.5. No complete spectrum or spectral weights.",
      "Tolerance is absolute in Y=L*E_C per state, default 16777216 epsilon. fp128 may take many minutes.",
      "bits-N/digits-N selects MPFR; nonlinear, mesh, hole and Fourier budgets remain independent.",
      "Failed observables are missing; exit 2 if either solve is incomplete.",
      "See docs/sine-gordon-excited.md for scope and controls; --references prints literature."};
  return info;
}
template <uni20::Real Real> int run(Arguments const& a, model::TwoSolitonOptions<Real> controls, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const mass = a.mass ? cli::parse_solver_real<Real>(*a.mass, controls)
                           : uni20::scalar_at_precision<Real>(1, controls.precision());
  Real const length = cli::parse_solver_real<Real>(*a.length, controls),
             p = cli::parse_solver_real<Real>(*a.coupling, controls);
  auto const number = uni20::half_int::parse(a.number);
  cli::sine_gordon::configure<Real>(a, controls);
  controls.max_root_iterations = a.root_iterations;
  // Reject invalid requests before admitting either independent solve.
  (void)model::detail::validate_two_soliton(mass, length, p, number, controls);
  // Both solves finish before opening exports, including --force targets.
  auto const [excited, vacuum] = context.measure([&] {
    return cli::independent_solves(
        [=] { return model::two_soliton_level(mass, length, p, number, controls); },
        [=] { return model::vacuum_energy(mass, length, p, static_cast<model::VacuumOptions<Real> const&>(controls)); },
        a.threads);
  });
  bool const converged = excited.converged && vacuum.converged;
  using Optional = std::optional<Real>;
  auto finite = [](Real x) -> Optional { return uni20::isfinite(x) ? Optional(x) : std::nullopt; };
  Optional gap, scaled_gap;
  if (converged)
  {
    gap = finite(*excited.casimir_energy - *vacuum.casimir_energy);
    scaled_gap = finite((*excited.scaling_function - *vacuum.scaling_function) /
                        (Real{8} * uni20::math::atan(uni20::scalar_at_precision<Real>(1, controls.precision()))));
  }
  bool const complete = converged && gap && scaled_gap;
  std::string const status = complete             ? "converged"
                             : !excited.converged ? "excited:" + std::string(name(excited.status))
                             : !vacuum.converged  ? "vacuum:" + std::string(name(vacuum.status))
                                                  : "precision_limit";
  cli::RunReport report(context, "Sine-Gordon exact finite-volume two-soliton level");
  report.field("mass", "Soliton mass", mass)
      .field("length", "Circumference", length)
      .field("p", "Coupling p", p)
      .field("scaled_length", "Scaled length u", excited.scaled_length)
      .field("units", "Units", "velocity=hbar=1; M is the soliton mass")
      .field("coupling_convention", "Coupling convention", "p=beta^2/(8*pi-beta^2); kinetic term (1/2)(partial phi)^2")
      .field("subtraction", "Energy convention", "bulk-subtracted E_C=E-L*e_bulk; gap=E_C(excited)-E_C(vacuum)")
      .field("scaling", "Dimensionless observables", "Y=L*E_C; scaled_gap=L*gap/(2*pi)")
      .field("calculation", "Calculation", "exact continuum two-hole NLIE; not Bethe-Yang")
      .field("sector", "Sector", "same-charge pair, delta=0, momentum zero; no complex/special roots")
      .field("number", "Positive Bethe number", number)
      .field("charge", "Topological winding charge", a.charge)
      .field("uv_weight", "UV weight Delta+=Delta-",
             (Real{1} + Real{1} / p) / Real{4} + Real(number.twice() - 1) / Real{2})
      .field("threads", "Requested threads", a.threads)
      .field("solver_concurrency", "Solver concurrency", std::min(a.threads, 2))
      .field("execution", "Execution", a.threads == 1 ? "serial" : "independent vacuum/excited solves")
      .precision(a.precision)
      .field("max_root_iterations", "Max hole trials across meshes/contours", controls.max_root_iterations)
      .result(complete, status);
  cli::sine_gordon::report_controls<Real>(report, controls, excited);
  cli::ResultOutput output(report, a.output, {"levels", "source", "gap"});
  using cli::column;
  output.table(
      "levels", "Bulk-subtracted finite-volume levels",
      [&](auto& table) {
        auto append = [&](std::string state, model::VacuumState<Real> const& s) {
          table.append(state, s.casimir_energy, s.scaling_function, finite(s.nonlinear_residual),
                       finite(s.kernel_error), finite(s.mesh_error), finite(s.cutoff_error), finite(s.contour_error),
                       s.cutoff, s.intervals, s.iterations, s.kernel_evaluations, s.cutoffs, s.converged,
                       std::string(name(s.status)));
        };
        append("vacuum", vacuum);
        append("two_soliton", excited);
      },
      column<std::string>("state"), column<Optional>("casimir_energy"), column<Optional>("scaling_function"),
      column<Optional>("nonlinear_residual"), column<Optional>("kernel_error"), column<Optional>("mesh_error"),
      column<Optional>("cutoff_error"), column<Optional>("contour_error"), column<Real>("cutoff"),
      column<std::size_t>("intervals"), column<std::size_t>("iterations"), column<std::size_t>("kernel_evaluations"),
      column<std::size_t>("cutoffs"), column<bool>("converged"), column<std::string>("status"));
  output.table(
      "source", "Two-hole quantization",
      [&](auto& table) {
        table.append(number, a.charge, excited.rapidity, finite(excited.quantization_residual),
                     finite(excited.hole_error), finite(excited.source_quadrature_error),
                     finite(excited.source_tail_bound), excited.root_iterations, excited.source_evaluations,
                     excited.converged, std::string(name(excited.status)));
      },
      column<uni20::half_int>("number"), column<int>("charge"), column<Optional>("rapidity"),
      column<Optional>("quantization_residual"), column<Optional>("hole_error"),
      column<Optional>("source_quadrature_error"), column<Optional>("source_tail_bound"),
      column<std::size_t>("root_iterations"), column<std::size_t>("source_evaluations"), column<bool>("converged"),
      column<std::string>("status"));
  output.table(
      "gap", "Vacuum-relative zero-momentum excitation",
      [&](auto& table) { table.append(gap, scaled_gap, complete, status); }, column<Optional>("gap"),
      column<Optional>("scaled_gap"), column<bool>("converged"), column<std::string>("status"));
  output.finish();
  if (!complete) std::cerr << "Sine-Gordon excitation incomplete; gap unavailable (" << status << ").\n";
  return complete ? 0 : 2;
}
} // namespace
int main(int argc, char** argv)
{
  Arguments a;
  return cli::program_main(
      argc, argv, program_info(),
      [&](auto& app) {
        cli::sine_gordon::add_options(app, a, "Sine-Gordon coupling p>=1",
                                      "Absolute Y tolerance per state; default 16777216 epsilon",
                                      cli::ArbitraryPrecision::real);
        bethe::cli::add_data_output_options(
            app, a.output,
            {{.name = "levels", .description = "Vacuum and excited levels", .primary = true},
             {.name = "source", .description = "Excited-source quantization and diagnostics"},
             {.name = "gap", .description = "Vacuum-relative gap; values missing if either solve fails"}});
        app.add_option("--number", a.number, "Positive half-odd Bethe number: 0.5 or 1.5")->capture_default_str();
        app.add_option("--charge", a.charge, "Pair winding charge")
            ->check(CLI::IsMember({-2, 2}))
            ->capture_default_str();
        cli::threads_option(app, a.threads)
            ->description("Concurrency for vacuum and excited solves (at most 2; default 1)");
        cli::count_option(app, "--max-root-iterations", a.root_iterations,
                          "Total hole trials across meshes and contours")
            ->capture_default_str();
      },
      [&](auto&) {
        a.output.validate();
        return cli::dispatch_options_precision<cli::ArbitraryPrecision::real, model::TwoSolitonOptions>(
            a.precision, [&](auto controls) { return run(a, std::move(controls), argc, argv); });
      });
}
