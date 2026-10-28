// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "execution-options.hpp"
#include "result-output.hpp"
#include "solver-precision.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/sine_gordon_bethe_yang.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::sine_gordon;
struct Arguments
{
    std::string mass = "1", number = "0.5", precision = "fp64";
    std::optional<std::string> coupling, length, tolerance;
    std::size_t levels = 1, iterations = 256, evaluations = 100000, quadrature_levels = 16, cutoffs = 64;
    int charge = 2;
    int threads = 1;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info =
      cli::program_info("bethe-sine-gordon-bethe-yang", "Asymptotic sine-Gordon same-charge two-soliton levels.",
                        bethe::citations::Tool::sine_gordon_bethe_yang);
  info.examples = {{"bethe-sine-gordon-bethe-yang --p 2 --length 10 --levels 4 --csv pairs.csv",
                    "Lowest opposite-rapidity soliton pairs"}};
  info.notes = {
      "Repulsive/free p>=1; soliton mass M>0, circumference L>0, velocity=hbar=1.",
      "Same-charge pair, total momentum zero, topological winding charge +/-2. Not a neutral pair.",
      "I=0.5,1.5,...; roots +/-theta solve M*L*sinh(theta)+chi(2*theta)=2*pi*I, S_ss=-exp(i*chi).",
      "Energy=2*M*cosh(theta) above the vacuum at Bethe-Yang order. Wrapping corrections are omitted.",
      "These are NOT exact finite-volume energies. Small residual certifies the equation, not the large-volume "
      "approximation.",
      "Quadrature/cutoff errors and counting residual are separate; no error bound on omitted wrapping terms.",
      "--threads schedules independent levels; each root search and phase quadrature stays sequential. Rows keep "
      "Bethe-number order, including failed levels. --no-retain permits bounded-memory streaming.",
      "One worker, one level, or the analytic free point p=1 uses direct serial evaluation.",
      "See docs/sine-gordon-excitations.md; --references prints literature."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::text_option(app, "--p", a.coupling, "Coupling p>=1")->required()->type_name("REAL");
  cli::text_option(app, "--length", a.length, "Circumference L>0")->required()->type_name("REAL");
  app.add_option("--mass", a.mass, "Soliton mass M>0")->capture_default_str()->type_name("REAL");
  app.add_option("--number", a.number, "First positive half-odd Bethe number I")->capture_default_str();
  app.add_option("--charge", a.charge, "Pair winding charge")->check(CLI::IsMember({-2, 2}))->capture_default_str();
  cli::count_option(app, "--levels", a.levels, "Consecutive levels, 1..10000")->capture_default_str();
  cli::text_option(app, "--tolerance", a.tolerance, "Absolute counting residual; default 16384 epsilon")
      ->type_name("REAL");
  cli::count_option(app, "--max-iterations", a.iterations, "Bisection updates per level")->capture_default_str();
  cli::count_option(app, "--max-phase-evaluations", a.evaluations, "Fourier evaluations per phase integral")
      ->capture_default_str();
  cli::count_option(app, "--max-phase-levels", a.quadrature_levels, "Quadrature refinement levels")
      ->capture_default_str();
  cli::count_option(app, "--max-fourier-cutoffs", a.cutoffs, "Fourier cutoff attempts")->capture_default_str();
  cli::precision_option(app, a.precision, cli::ArbitraryPrecision::real);
  cli::threads_option(app, a.threads)->description("Scheduler concurrency limit for independent Bethe-Yang levels");
  cli::add_data_output_options(app, a.output,
                               {{.name = "levels", .description = "Selected energy levels", .primary = true}});
}
char const* name(model::BetheYangStatus status)
{
  switch (status)
  {
    case model::BetheYangStatus::converged:
      return "converged";
    case model::BetheYangStatus::phase_limit:
      return "phase_limit";
    case model::BetheYangStatus::iteration_limit:
      return "iteration_limit";
    case model::BetheYangStatus::precision_limit:
      return "precision_limit";
  }
  return "unknown";
}
char const* name(model::KernelStatus status)
{
  switch (status)
  {
    case model::KernelStatus::converged:
      return "converged";
    case model::KernelStatus::quadrature_limit:
      return "quadrature_limit";
    case model::KernelStatus::cutoff_limit:
      return "cutoff_limit";
    case model::KernelStatus::precision_limit:
      return "precision_limit";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv, model::BetheYangOptions<Real> options)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const mass = cli::parse_solver_real<Real>(a.mass, options),
             p = cli::parse_solver_real<Real>(*a.coupling, options),
             length = cli::parse_solver_real<Real>(*a.length, options);
  auto const first = uni20::half_int::parse(a.number);
  if (first.twice() > std::numeric_limits<std::int64_t>::max() - 2 * std::int64_t(a.levels))
    throw std::invalid_argument("Bethe label overflow");
  if (a.tolerance) options.tolerance = cli::parse_solver_real<Real>(*a.tolerance, options);
  options.max_iterations = a.iterations;
  options.phase.max_evaluations = a.evaluations;
  options.phase.max_levels = a.quadrature_levels;
  options.phase.max_cutoffs = a.cutoffs;
  // Validate before opening/truncating exports, without spending a solve.
  auto validation = options;
  validation.max_iterations = 0;
  (void)model::same_charge_pair(mass, length, p, first, validation);
  bool const parallel = a.threads > 1 && a.levels > 1 && p != Real{1};
  std::optional<uni20::async::TbbScheduler> scheduler;
  if (parallel) scheduler.emplace(a.threads);
  bethe::ExecutionOptions const execution{.scheduler = scheduler ? &*scheduler : nullptr};
  cli::RunReport report(context, "Sine-Gordon asymptotic Bethe-Yang pairs");
  report.field("mass", "Soliton mass M", mass)
      .field("length", "Circumference L", length)
      .field("p", "Coupling p", p)
      .field("approximation", "Approximation",
             "Bethe-Yang; wrapping corrections omitted; not exact finite-volume levels")
      .field("energy_reference", "Energy reference", "excitation energy above vacuum at Bethe-Yang order")
      .field("sector", "Sector", "same-charge pair; total momentum zero; winding charge +/-2")
      .field("charge", "Topological charge", a.charge)
      .field("first", "First Bethe number", first)
      .field("levels", "Levels", a.levels)
      .field("threads", "Scheduler concurrency limit", a.threads)
      .field("execution_mode", "Execution mode", parallel ? "bounded Async level grid" : "direct serial")
      .field("max_in_flight", "Maximum levels in flight", parallel ? execution.max_in_flight : std::size_t{1})
      .field("cpu_timing_scope", "CPU timing scope", "level grid including ordered row output")
      .field("tolerance", "Counting tolerance", options.tolerance)
      .precision(a.precision);
  using Optional = std::optional<Real>;
  auto table = cli::data::make_data_table(
      "Sine-Gordon Bethe-Yang levels",
      {.retain = a.output.retain ? cli::data::retention::all : cli::data::retention::none,
       .metadata = report.metadata()},
      cli::column<uni20::half_int>("number"), cli::column<int>("charge"), cli::column<Optional>("rapidity"),
      cli::column<Optional>("energy"), cli::column<Real>("residual"), cli::column<Real>("quadrature_error"),
      cli::column<Real>("tail_bound"), cli::column<std::size_t>("iterations"), cli::column<std::size_t>("evaluations"),
      cli::column<std::string>("phase_status"), cli::column<std::string>("status"));
  bool const complete = cli::stream_result_table(
      report, context, a.output, "levels", table,
      [&](auto& table) {
        bool success = true;
        context.measure([&] {
          auto evaluate = [mass, length, p, first, options](std::size_t j) {
            auto const number = uni20::from_twice(first.twice() + 2 * std::int64_t(j));
            return model::same_charge_pair(mass, length, p, number, options);
          };
          auto collect = [&](std::size_t j, model::BetheYangResult<Real> const& state) {
            auto const number = uni20::from_twice(first.twice() + 2 * std::int64_t(j));
            success = success && state.converged;
            table.append(number, a.charge, state.rapidity, state.energy, state.residual, state.quadrature_error,
                         state.tail_bound, state.iterations, state.evaluations, name(state.phase_status),
                         name(state.status));
          };
          if (parallel)
            bethe::evaluate_observable_grid(a.levels, evaluate, collect, execution);
          else
            for (std::size_t j = 0; j < a.levels; ++j)
              collect(j, evaluate(j));
        });
        return success;
      },
      "incomplete; unavailable quantities omitted");
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
        if (a.levels < 1 || a.levels > 10000) throw std::invalid_argument("require 1<=levels<=10000");
        return cli::dispatch_options_precision<cli::ArbitraryPrecision::real, model::BetheYangOptions>(
            a.precision, [&]<typename Real>(model::BetheYangOptions<Real> options) {
              return run<Real>(a, argc, argv, std::move(options));
            });
      });
}
