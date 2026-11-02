// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "execution-options.hpp"
#include "program-options.hpp"
#include "result-output.hpp"
#include "solver-precision.hpp"
#include <bethe/experimental/xxx_boundary_magnetization.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::experimental::xxx::open;
struct Arguments
{
    std::size_t sites = 0, max_pairs = 250000, time_points = 101;
    std::size_t max_iterations = bethe::SolverOptions<double>{}.max_iterations;
    int threads = 1;
    std::string precision = "fp64", time_max = "100";
    std::optional<std::string> tolerance;
    bool states = false, pairs = false, roots = false;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-experimental-xxx-boundary-magnetization",
                                "Endpoint magnetization after a free-end XXX boundary spin flip.",
                                bethe::citations::Tool::xxx_boundary_magnetization);
  info.examples = {{"bethe-experimental-xxx-boundary-magnetization 16 --threads 2 --csv time.csv",
                    "Both endpoints after sqrt(2) S_1^+ acting on the singlet"},
                   {"bethe-experimental-xxx-boundary-magnetization 8 --precision bits-160 --json result.json",
                    "Explicit precision, signed pair coefficients and state diagnostics"}};
  info.notes = {
      "Experimental, even N>=2, H=sum S_j.S_(j+1), free ends, J=1, zero field.",
      "Initial full state sqrt(2) S_1^+|ground> has norm one and total Sz=1.",
      "Retain only the positive finite-real-root highest triplets, without renormalizing the projection.",
      "Pair coefficients are signed interference amplitudes, not DSF weights or probabilities.",
      "Only the two endpoints are calculated; no bulk magnetization profile.",
      "The omitted-norm bound excludes numerical error; all pairs must resolve before a time curve is emitted.",
      "Exit 2 indicates an unresolved or inconsistent result. See docs/experimental-xxx-boundary-magnetization.md."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::count_option(app, "N", a.sites, "Even chain length, at least two")->required();
  cli::count_option(app, "--max-pairs", a.max_pairs, "Upper-triangle pair budget, checked before roots")
      ->capture_default_str();
  cli::option(app, "--time-max", a.time_max, "Last time in inverse exchange units, finite and nonnegative")
      ->capture_default_str();
  cli::count_option(app, "--time-points", a.time_points, "Uniform grid including both endpoints, 2 to 1000000")
      ->capture_default_str();
  cli::count_option(app, "--max-iterations", a.max_iterations, "Iteration limit per root solve")->capture_default_str();
  cli::text_option(app, "--tolerance", a.tolerance, "Positive root residual tolerance at the selected precision");
  cli::precision_option(app, a.precision, cli::ArbitraryPrecision::real_lu);
  cli::threads_option(app, a.threads);
  app.add_flag("--states", a.states, "Show every requested state, its amplitude and preparation status");
  app.add_flag("--pairs", a.pairs, "Show every signed matrix element and cosine coefficient, including failures");
  app.add_flag("--roots", a.roots, "Show ground and attempted excited roots in z=2*lambda coordinates");
  cli::add_data_output_options(
      app, a.output,
      {{.name = "time", .description = "Retained magnetization at both endpoints", .primary = true},
       {.name = "summary", .description = "Projected norm and uniform omitted-state bound"},
       {.name = "states",
        .description = "Every requested state and preparation diagnostic",
        .screen_option = "--states"},
       {.name = "pairs",
        .description = "Every requested signed matrix pair and cosine term",
        .screen_option = "--pairs"},
       {.name = "roots", .description = "Ground and attempted excited z roots", .screen_option = "--roots"}});
}
std::string status(model::BoundarySpinStatus s)
{
  switch (s)
  {
    case model::BoundarySpinStatus::converged:
      return "converged";
    case model::BoundarySpinStatus::roots_unconverged:
      return "roots_unconverged";
    case model::BoundarySpinStatus::linear_solve_failed:
      return "linear_solve_failed";
    case model::BoundarySpinStatus::precision_limit:
      return "precision_limit";
  }
  throw std::logic_error("unknown boundary spin status");
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
  struct TimeRow
  {
      R time, left, right;
  };
  auto calculation = context.measure([&] {
    auto result = model::boundary_real_triplet_magnetization<R>(a.sites, controls, a.max_pairs, {&scheduler});
    std::vector<TimeRow> times;
    if (result.converged())
      for (std::size_t i = 0; i < a.time_points; ++i)
      {
        R const time = last * (uni20::scalar_like(zero, i) / uni20::scalar_like(zero, a.time_points - 1));
        times.push_back({time, result.at(time), result.at(time, true)});
      }
    return std::pair(std::move(result), std::move(times));
  });
  auto const& result = calculation.first;
  auto const& spectrum = result.spectrum;
  bool const complete = result.converged();
  auto const resolved =
      std::count_if(result.pairs.begin(), result.pairs.end(), [](auto const& p) { return p.term.has_value(); });
  cli::RunReport report(context, "Experimental XXX magnetization after a boundary spin flip");
  report.field("sites", "Sites", a.sites)
      .field("hamiltonian", "Hamiltonian", "H=sum S_j.S_(j+1); free ends; J=1; zero field")
      .field("family", "Family", "Positive finite-real-root S=1 triplets; partial spectrum")
      .field("initial_state", "Initial state", "sqrt(2) S_1^+|singlet>, projected without renormalization")
      .field("observable", "Observable", "<psi_A(t)|S_1^z|psi_A(t)> and <psi_A(t)|S_N^z|psi_A(t)>")
      .field("coefficient_convention", "Pair convention",
             "Signed cosine terms; positive raising gauge; no Fourier factor")
      .field("root_coordinate", "Root coordinate", "z=2*lambda_ABA")
      .precision(a.precision)
      .field("threads", "Scheduler concurrency limit", a.threads)
      .field("max_iterations", "Root iteration limit", controls.max_iterations)
      .field("root_tolerance", "Root residual tolerance", controls.residual_tolerance)
      .field("max_pairs", "Max pairs", a.max_pairs)
      .field("requested_states", "Requested states", spectrum.records.size())
      .field("requested_pairs", "Requested pairs", result.pairs.size())
      .field("resolved_pairs", "Resolved pairs", resolved)
      .field("time_max", "Last time", last)
      .field("time_points", "Time points", a.time_points)
      .result(complete, complete ? "requested_real_family_and_pairs_obtained"
                                 : "unresolved_or_inconsistent_boundary_magnetization");
  cli::ResultOutput output(report, a.output, {"time", "summary", "states", "pairs", "roots"});
  using cli::column;
  output.table(
      "time", "Retained endpoint magnetization",
      [&](auto& table) {
        for (auto const& row : calculation.second)
          table.append(row.time, row.left, row.right);
      },
      column<R>("time"), column<R>("left_magnetization"), column<R>("right_magnetization"));
  output.table(
      "summary", "Norm and omitted-state bound",
      [&](auto& table) {
        table.append(result.projected_norm, one, result.omitted_norm(), result.truncation_bound(), result.matrix_ready);
      },
      column<Optional>("projected_norm"), column<R>("full_norm"), column<Optional>("omitted_norm"),
      column<Optional>("uniform_truncation_bound"), column<bool>("matrix_ready"));
  output.table(
      "states", "Every requested triplet",
      [&](auto& table) {
        for (std::size_t i = 0; i < spectrum.records.size(); ++i)
        {
          auto const& row = spectrum.records[i];
          auto const& s = result.states[i];
          table.append(i + 1, row.candidate_id, row.line ? Optional(row.line->gap) : Optional{}, row.reflection_parity,
                       row.state ? (row.state->converged ? "converged" : "roots_unconverged") : "not_attempted",
                       s.attempted() ? status(s.status()) : "not_attempted", s.raising_amplitude(), s.diagonal(),
                       s.solve_backward_error(), s.pivot_margin(), row.diagnostic);
        }
      },
      column<std::size_t>("state_id"), column<std::string>("candidate_id"), column<Optional>("gap"),
      column<int>("reflection_parity"), column<std::string>("root_status"), column<std::string>("spin_status"),
      column<Optional>("raising_amplitude"), column<Optional>("diagonal_sz"), column<Optional>("solve_backward_error"),
      column<Optional>("pivot_margin"), column<std::string>("diagnostic"));
  output.table(
      "pairs", "Signed upper-triangle matrix and cosine terms",
      [&](auto& table) {
        for (auto const& p : result.pairs)
          table.append(p.left_state, p.right_state, spectrum.records[p.left_state - 1].candidate_id,
                       spectrum.records[p.right_state - 1].candidate_id, p.element ? p.element->value : Optional{},
                       p.term ? Optional(p.term->frequency) : Optional{},
                       p.term ? Optional(p.term->coefficient) : Optional{},
                       p.term ? std::optional<int>(p.term->opposite_end_parity) : std::optional<int>{},
                       p.element ? status(p.element->status) : "not_attempted",
                       p.element ? p.element->pivot_margin : Optional{}, p.diagnostic);
      },
      column<std::size_t>("left_state"), column<std::size_t>("right_state"), column<std::string>("left_candidate"),
      column<std::string>("right_candidate"), column<Optional>("matrix_element"), column<Optional>("frequency"),
      column<Optional>("coefficient"), column<std::optional<int>>("opposite_end_parity"), column<std::string>("status"),
      column<Optional>("pivot_margin"), column<std::string>("diagnostic"));
  output.table(
      "roots", "Ground and attempted excited roots",
      [&](auto& table) {
        auto append = [&](std::string const& id, auto const& state) {
          for (std::size_t i = 0; i < state.rapidities.size(); ++i)
            table.append(id, i, state.quantum_numbers[i], state.rapidities[i], state.converged);
        };
        append("ground", spectrum.ground_state);
        for (auto const& row : spectrum.records)
          if (row.state) append(row.candidate_id, *row.state);
      },
      column<std::string>("candidate_id"), column<std::size_t>("root_index"), column<uni20::half_int>("quantum_number"),
      column<R>("z"), column<bool>("roots_converged"));
  output.finish();
  if (!complete) std::cerr << "Magnetization is unresolved; no time curve was produced. Inspect states and pairs.\n";
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
