// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "execution-options.hpp"
#include "program-options.hpp"
#include "result-output.hpp"
#include "solver-precision.hpp"
#include <bethe/experimental/xxx_boundary_dynamics.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::experimental::xxx::open;
struct Arguments
{
    std::size_t sites = 0, max_candidates = 10000, time_points = 101;
    std::size_t max_iterations = bethe::SolverOptions<double>{}.max_iterations;
    int threads = 1;
    std::string precision = "fp64", time_max = "100";
    std::optional<std::string> tolerance;
    bool spectrum = false, diagnostics = false, roots = false, density = false;
    cli::DataOutputOptions output;
};

auto program_info()
{
  auto info = cli::program_info("bethe-experimental-xxx-boundary-dynamics",
                                "Free-end XXX local and end-to-end time correlations.",
                                bethe::citations::Tool::xxx_boundary_dynamics);
  info.examples = {{"bethe-experimental-xxx-boundary-dynamics 16 --threads 2 --csv time.csv",
                    "Boundary correlations at times 0 to 100 from the finite-real triplets"},
                   {"bethe-experimental-xxx-boundary-dynamics 8 --precision bits-160 --spectrum --diagnostics "
                    "--roots --density --json result.json",
                    "Explicit-precision calculation and all declared tables"}};
  info.notes = {
      "Experimental, even N>=2, H=sum S_j.S_(j+1), free ends, J=1, zero field.",
      "C11(t)=<S_1^z(t) S_1^z(0)> and C1N(t)=<S_1^z(t) S_N^z(0)>, phase exp(-i gap t).",
      "Only finite positive-real-root S=1 triplets are included; this is a partial spectral sum.",
      "Local weights are raising/2 with no Fourier N or 2*pi factor. No missing-weight rescaling.",
      "The full local norm is 1/4. The full first moment comes independently from the ground-state boundary bond.",
      "Reflection parities fix end-to-end signs; squared local weights alone do not fix general amplitudes.",
      "Exit 2 means incomplete roots/weights, a rejected bond solve or inconsistent sum rules.",
      "See docs/experimental-xxx-boundary-dynamics.md; --references for literature."};
  return info;
}

void add_options(CLI::App& app, Arguments& a)
{
  cli::count_option(app, "N", a.sites, "Even chain length, at least two")->required();
  cli::count_option(app, "--max-candidates", a.max_candidates, "Reject larger families before solving")
      ->capture_default_str();
  cli::option(app, "--time-max", a.time_max, "Last time in inverse exchange units, finite and nonnegative")
      ->capture_default_str();
  cli::count_option(app, "--time-points", a.time_points, "Uniform grid including both endpoints, 2 to 1000000")
      ->capture_default_str();
  cli::count_option(app, "--max-iterations", a.max_iterations, "Iteration limit per root solve")->capture_default_str();
  cli::text_option(app, "--tolerance", a.tolerance, "Positive root residual tolerance at the selected precision");
  cli::precision_option(app, a.precision, cli::ArbitraryPrecision::real_lu);
  cli::threads_option(app, a.threads);
  app.add_flag("--spectrum", a.spectrum, "Show accepted local weights and signed endpoint coefficients");
  app.add_flag("--diagnostics", a.diagnostics, "Show every requested label, including failures");
  app.add_flag("--roots", a.roots, "Show ground and attempted excited roots in z=2*lambda coordinates");
  app.add_flag("--density", a.density, "Show the complete ground-state two-site boundary density matrix");
  cli::add_data_output_options(
      app, a.output,
      {{.name = "time", .description = "Local and end-to-end complex correlations", .primary = true},
       {.name = "moments", .description = "Retained moments, full local sum rules and unscaled omissions"},
       {.name = "spectrum", .description = "Accepted weights, gaps and endpoint signs", .screen_option = "--spectrum"},
       {.name = "diagnostics", .description = "All requested labels and failures", .screen_option = "--diagnostics"},
       {.name = "roots", .description = "Ground and attempted excited z roots", .screen_option = "--roots"},
       {.name = "density", .description = "Full ground-root two-site density matrix", .screen_option = "--density"}});
}

std::string status(model::BoundaryFormFactorStatus value)
{
  switch (value)
  {
    case model::BoundaryFormFactorStatus::converged:
      return "converged";
    case model::BoundaryFormFactorStatus::roots_unconverged:
      return "roots_unconverged";
    case model::BoundaryFormFactorStatus::precision_limit:
      return "precision_limit";
  }
  throw std::logic_error("unknown boundary form-factor status");
}

std::string status(model::BoundaryCorrelationStatus value)
{
  switch (value)
  {
    case model::BoundaryCorrelationStatus::converged:
      return "converged";
    case model::BoundaryCorrelationStatus::roots_unconverged:
      return "roots_unconverged";
    case model::BoundaryCorrelationStatus::linear_solve_failed:
      return "linear_solve_failed";
    case model::BoundaryCorrelationStatus::precision_limit:
      return "precision_limit";
  }
  throw std::logic_error("unknown boundary correlation status");
}

template <uni20::Real R> int run(Arguments const& a, int argc, char** argv, bethe::SolverOptions<R> const& controls)
{
  using Optional = std::optional<R>;
  R const zero = bethe::detail::solver_zero(controls), one = uni20::scalar_like(zero, 1);
  R const last = cli::parse_solver_real<R>(a.time_max, controls);
  if (!uni20::isfinite(last) || last < zero) throw std::invalid_argument("time-max must be finite and nonnegative");
  if (a.time_points < 2 || a.time_points > 1000000) throw std::invalid_argument("time-points must be in [2,1000000]");
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  uni20::async::TbbScheduler scheduler(a.threads);
  // Await the final owned result. No scheduler-wide run_all or intermediate batches.
  struct TimeRow
  {
      R time;
      uni20::complex<R> local, cross;
  };
  auto calculation = context.measure([&] {
    auto spectrum = model::boundary_real_triplet_spectrum<R>(a.sites, controls, a.max_candidates, {&scheduler});
    std::vector<TimeRow> times;
    if (spectrum.ground_status == model::BoundaryFormFactorStatus::converged)
    {
      times.reserve(a.time_points);
      for (std::size_t i = 0; i < a.time_points; ++i)
      {
        R const time = last * (uni20::scalar_like(zero, i) / uni20::scalar_like(zero, a.time_points - 1));
        times.push_back({time, bethe::local_time_correlation<R>(spectrum.lines, time),
                         bethe::cross_time_correlation<R>(spectrum.end_to_end_lines, time)});
      }
    }
    return std::pair(std::move(spectrum), std::move(times));
  });
  auto const& result = calculation.first;
  bool const complete = result.converged();
  bool const ground_ready = result.ground_state.converged;
  auto const retained = result.retained_moments();
  auto const cross = result.retained_end_to_end_moments();
  cli::RunReport report(context, "Experimental free-end XXX boundary dynamics");
  report.field("sites", "Sites", a.sites)
      .field("hamiltonian", "Hamiltonian", "H=sum S_j.S_(j+1); free ends; J=1; zero field")
      .field("family", "Family", "Positive finite-real-root S=1 triplets; partial spectrum")
      .field("channel", "Channel", "zz")
      .field("phase", "Time phase", "exp(-i gap t)")
      .field("weight_convention", "Weight convention", "local raising/2; no Fourier N or 2*pi factor")
      .field("root_coordinate", "Root coordinate", "z=2*lambda_ABA")
      .precision(a.precision)
      .field("threads", "Scheduler concurrency limit", a.threads)
      .field("max_iterations", "Root iteration limit", controls.max_iterations)
      .field("root_tolerance", "Root residual tolerance", controls.residual_tolerance)
      .field("max_candidates", "Max candidates", a.max_candidates)
      .field("ground_iterations", "Ground iterations", result.ground_state.iterations)
      .field("ground_solver_residual", "Ground solver residual", result.ground_state.residual_norm)
      .field("ground_energy", "Ground energy", ground_ready ? Optional(result.ground_state.energy) : Optional{})
      .field("ground_solver_converged", "Ground root solver converged", ground_ready)
      .field("ground_status", "Spectral ground preparation", status(result.ground_status))
      .field("ground_residual", "Ground equation residual", result.ground_residual)
      .field("bond_status", "Boundary bond status", status(result.density.status))
      .field("bond", "Full boundary bond", result.density.bond_correlation)
      .field("bond_backward_error", "Bond solve backward error", result.density.solve_backward_error)
      .field("sum_rules_consistent", "Sum rules consistent at working precision", result.sum_rules_consistent())
      .field("requested", "Requested labels", result.records.size())
      .field("accepted", "Accepted weights", result.lines.size())
      .field("time_max", "Last time", last)
      .field("time_points", "Time points", a.time_points)
      .result(complete, complete ? "requested_real_family_obtained" : "incomplete_or_inconsistent_boundary_result");
  cli::ResultOutput output(report, a.output, {"time", "moments", "spectrum", "diagnostics", "roots", "density"});
  using cli::column;
  output.table(
      "time", "Retained boundary correlations",
      [&](auto& table) {
        for (auto const& row : calculation.second)
          table.append(row.time, row.local.real(), row.local.imag(), uni20::math::abs(row.local), row.cross.real(),
                       row.cross.imag(), uni20::math::abs(row.cross));
      },
      column<R>("time"), column<R>("local_real"), column<R>("local_imag"), column<R>("local_abs"),
      column<R>("end_to_end_real"), column<R>("end_to_end_imag"), column<R>("end_to_end_abs"));
  output.table(
      "moments", "Retained moments and complete local sum rules",
      [&](auto& table) {
        if (!retained && !result.density.converged()) return;
        auto const zz = result.density.zz();
        table.append(retained ? Optional(retained->weight) : Optional{}, one / 4, result.omitted_weight(),
                     retained ? Optional(retained->first_moment) : Optional{}, zz ? Optional(-*zz) : Optional{},
                     result.omitted_first_moment(), cross ? Optional(cross->value.real()) : Optional{},
                     cross ? Optional(cross->first_moment.real()) : Optional{});
      },
      column<Optional>("retained_weight"), column<R>("full_local_weight"), column<Optional>("omitted_weight"),
      column<Optional>("retained_first_moment"), column<Optional>("full_local_first_moment"),
      column<Optional>("omitted_first_moment"), column<Optional>("end_to_end_static"),
      column<Optional>("end_to_end_first_moment"));
  output.table(
      "spectrum", "Accepted local weights and signed endpoint coefficients",
      [&](auto& table) {
        for (auto const& row : result.records)
          if (row.line)
            table.append(row.candidate_id, row.line->gap, row.line->weight, 2 * row.line->weight, row.reflection_parity,
                         result.ground_reflection_parity * row.reflection_parity * row.line->weight);
      },
      column<std::string>("candidate_id"), column<R>("gap"), column<R>("zz_weight"), column<R>("raising_weight"),
      column<int>("reflection_parity"), column<R>("end_to_end_coefficient"));
  output.table(
      "diagnostics", "Every requested real-root triplet",
      [&](auto& table) {
        for (auto const& row : result.records)
        {
          std::string labels;
          for (auto q : row.quantum_numbers)
          {
            if (!labels.empty()) labels += ',';
            labels += std::to_string(q.twice() / 2);
          }
          table.append(row.candidate_id, labels,
                       row.state ? (row.state->converged ? "converged" : "roots_unconverged") : "not_attempted",
                       row.factor ? status(row.factor->status) : "not_attempted", row.line.has_value(),
                       row.state ? Optional(row.state->residual_norm) : Optional{},
                       row.factor ? row.factor->pivot_margin : Optional{}, row.diagnostic);
        }
      },
      column<std::string>("candidate_id"), column<std::string>("quantum_numbers"), column<std::string>("root_status"),
      column<std::string>("factor_status"), column<bool>("accepted"), column<Optional>("residual"),
      column<Optional>("pivot_margin"), column<std::string>("diagnostic"));
  output.table(
      "roots", "Ground and attempted excited roots in z=2*lambda coordinates",
      [&](auto& table) {
        auto append = [&](std::string const& id, auto const& state) {
          for (std::size_t i = 0; i < state.rapidities.size(); ++i)
            table.append(id, i, state.quantum_numbers[i], state.rapidities[i], state.converged);
        };
        append("ground", result.ground_state);
        for (auto const& row : result.records)
          if (row.state) append(row.candidate_id, *row.state);
      },
      column<std::string>("candidate_id"), column<std::size_t>("root_index"), column<uni20::half_int>("quantum_number"),
      column<R>("z"), column<bool>("roots_converged"));
  output.table(
      "density", "Complete ground-state boundary two-site density matrix",
      [&](auto& table) {
        if (auto rho = result.density.matrix())
          for (std::size_t i = 0; i < 16; ++i)
            table.append(i / 4, i % 4, (*rho)[i]);
      },
      column<std::size_t>("row"), column<std::size_t>("column"), column<R>("value"));
  output.finish();
  if (!complete) std::cerr << "Boundary result is incomplete or inconsistent; inspect diagnostics and summary.\n";
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
        return cli::dispatch_solver_precision<cli::ArbitraryPrecision::real_lu>(
            a.precision, a.tolerance, a.max_iterations,
            [&](auto const& controls) { return run(a, argc, argv, controls); });
      });
}
