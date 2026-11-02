// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "execution-options.hpp"
#include "program-options.hpp"
#include "result-output.hpp"
#include "solver-precision.hpp"
#include <algorithm>
#include <bethe/experimental/xxx_four_spinon_scan.hpp>
#include <bethe/xxx_structure_factor.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::experimental::xxx;
struct Arguments
{
    int sites = 0, threads = 1;
    std::size_t max_candidates = 10000;
    std::optional<int> momentum, branch;
    std::vector<int> holes;
    std::string precision = "fp64", channel = "zz";
    bool diagnostics = false, roots = false;
    cli::DataOutputOptions output;
};

auto program_info()
{
  auto info = cli::program_info("bethe-experimental-xxx-four-spinon", "Finite-ring XXX four-spinon spectral lines.",
                                bethe::citations::Tool::xxx_structure_factor);
  info.examples = {{"bethe-experimental-xxx-four-spinon 12 --holes=-3,-1,0,2 --branch 1 --csv line.csv",
                    "One resolved triplet at momentum transfer pi/2"},
                   {"bethe-experimental-xxx-four-spinon 16 --momentum-index 4 --threads 2 --json slice.json",
                    "Attempt every four-hole request at q=pi/2, retaining failures"}};
  info.notes = {
      "Experimental finite-deviation one-two-string family; N>=8 divisible by four, J=1, periodic, zero field.",
      "Default: attempt all four-hole/branch labels. Singular/collapsed and unconverged states can be rejected.",
      "This is a partial DSF, not a complete finite four-spinon spectrum. No missing-weight estimate or rescaling.",
      "S(q,w)=2*pi*sum weight*delta(w-gap); S_q^+=sum exp(-iqj)S_j^+/sqrt(N).",
      "q=P0-Pn; --momentum-index j selects q=2*pi*j/N. Each explicit request has multiplicity one.",
      "raising is S-+; zz weights are raising/2. Root output uses standard ABA lambda, not z=2*lambda.",
      "The signed pair deviation is retained separately; rounded complex roots can lose it.",
      "Exit status 2 means some requested weights were not obtained; inspect diagnostics.",
      "See docs/experimental-xxx-four-spinon-scan.md; --references for literature."};
  return info;
}

void add_options(CLI::App& app, Arguments& a)
{
  app.add_option("N", a.sites, "Sites, divisible by four and at least eight")->required();
  app.add_option("--momentum-index", a.momentum, "Fourier transfer index in [0,N)");
  app.add_option("--holes", a.holes, "One selected request: four distinct integer holes, comma separated")
      ->delimiter(',')
      ->expected(4);
  app.add_option("--branch", a.branch, "Selected two-string branch: -1, 0 or 1");
  cli::count_option(app, "--max-candidates", a.max_candidates, "Reject larger catalogues before solving")
      ->capture_default_str();
  cli::option(app, "--channel", a.channel, "zz or raising (S-+)")
      ->check(CLI::IsMember({"zz", "raising"}))
      ->capture_default_str();
  cli::precision_option(app, a.precision, cli::ArbitraryPrecision::complex_linalg);
  cli::threads_option(app, a.threads);
  app.add_flag("--diagnostics", a.diagnostics, "Show every requested state and rejection reason on screen");
  app.add_flag("--roots", a.roots, "Show standard ABA roots on screen; see separate pair deviations");
  cli::add_data_output_options(
      app, a.output,
      {{.name = "spectrum", .description = "Accepted finite-four spectral lines; unit multiplicity", .primary = true},
       {.name = "moments", .description = "Retained zeroth/first moments and full first-moment rule"},
       {.name = "diagnostics",
        .description = "Every requested label, solve/factor status and signed deviation",
        .screen_option = "--diagnostics"},
       {.name = "roots",
        .description = "Prepared ground and accepted-solve standard ABA roots",
        .screen_option = "--roots"}});
}

std::vector<model::FourSpinonRequest> requests(Arguments const& a)
{
  if (a.holes.empty())
  {
    if (a.branch) throw std::invalid_argument("--branch requires --holes");
    return model::four_spinon_catalogue(a.sites, a.momentum, a.max_candidates);
  }
  if (!a.branch || a.holes.size() != 4) throw std::invalid_argument("--holes requires four values and --branch");
  if (!a.max_candidates) throw std::invalid_argument("max_candidates must be positive");
  model::FourSpinonRequest request{a.sites, {a.holes[0], a.holes[1], a.holes[2], a.holes[3]}, *a.branch};
  std::sort(request.holes.begin(), request.holes.end());
  (void)model::four_spinon_id(request); // Validate before ground solve or opening exports.
  if (a.momentum && *a.momentum != model::four_spinon_transfer_index(request))
    throw std::invalid_argument("selected holes have a different momentum transfer");
  return {request};
}

std::string status(model::FourSpinonStatus s)
{
  switch (s)
  {
    case model::FourSpinonStatus::converged:
      return "converged";
    case model::FourSpinonStatus::unsupported_geometry:
      return "unsupported_geometry";
    case model::FourSpinonStatus::unconverged:
      return "unconverged";
    case model::FourSpinonStatus::precision_limit:
      return "precision_limit";
  }
  throw std::logic_error("unknown four-spinon solve status");
}
std::string status(bethe::heisenberg::FormFactorStatus s)
{
  switch (s)
  {
    case bethe::heisenberg::FormFactorStatus::converged:
      return "converged";
    case bethe::heisenberg::FormFactorStatus::roots_unconverged:
      return "roots_unconverged";
    case bethe::heisenberg::FormFactorStatus::precision_limit:
      return "precision_limit";
  }
  throw std::logic_error("unknown form factor status");
}

template <uni20::Real R> int run(Arguments const& a, int argc, char** argv, bethe::SolverOptions<R> const& controls)
{
  auto const candidates = requests(a);
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  uni20::async::TbbScheduler scheduler(a.threads);
  auto result = context.measure([&] {
    auto state = bethe::heisenberg::sector_ground_state<R>(a.sites, uni20::half_int{0}, controls);
    model::PreparedSinglet<R> ground(a.sites, state);
    std::optional<model::FourSpinonScan<R>> scan;
    if (ground.ready()) scan = model::scan_four_spinons(ground, candidates, {&scheduler}, a.max_candidates);
    return std::pair(std::move(ground), std::move(scan));
  });
  auto const& ground = result.first;
  auto const& scan = result.second;
  bool const complete = scan && scan->converged();
  R const zero = bethe::detail::solver_zero(controls), one = uni20::scalar_like(zero, 1);
  R const scale = one / (a.channel == "zz" ? 2 : 1), pi = bethe::detail::pi_like(zero);
  using Optional = std::optional<R>;
  auto const selected_momentum =
      a.holes.empty() ? a.momentum : std::optional<int>(model::four_spinon_transfer_index(candidates.front()));
  cli::RunReport report(context, "Experimental finite XXX four-spinon spectrum");
  report.field("sites", "Sites", a.sites)
      .field("hamiltonian", "Hamiltonian", "H=sum S_j.S_(j+1); periodic; J=1; zero field")
      .field("family", "Family", "Finite four-hole/one-two-string S=1 requests; partial DSF")
      .field("channel", "Channel", a.channel)
      .field("fourier_convention", "Fourier convention", "S_q^+=sum exp(-iqj)S_j^+/sqrt(N); q=P0-Pn")
      .field("weight_convention", "Weight convention", "S(q,w)=2*pi*sum weight*delta(w-gap); zz=raising/2")
      .field("multiplicity", "Multiplicity per label", 1)
      .field("selection", "Selection",
             !a.holes.empty() ? "single label"
             : a.momentum     ? "fixed-momentum catalogue"
                              : "all-momentum catalogue")
      .field("momentum_index", "Requested transfer index (absent=all)", selected_momentum)
      .precision(a.precision)
      .field("threads", "Scheduler concurrency limit", a.threads)
      .field("max_candidates", "Max candidates", a.max_candidates)
      .field("ground_energy", "Ground energy", ground.energy())
      .field("requested", "Requested labels", candidates.size())
      .field("accepted", "Accepted lines", scan ? scan->lines.size() : 0)
      .field("rejected", "Missing requested weights", candidates.size() - (scan ? scan->lines.size() : 0))
      .result(complete, complete ? "requested_weights_obtained" : "incomplete_requested_family");
  cli::ResultOutput output(report, a.output, {"spectrum", "moments", "diagnostics", "roots"});
  using cli::column;
  output.table(
      "spectrum", "Accepted four-spinon lines",
      [&](auto& table) {
        if (scan)
          for (auto const& row : scan->records)
            if (row.line)
              table.append(row.candidate_id, row.line->momentum_index, row.line->momentum, row.line->gap,
                           scale * row.line->weight, 1);
      },
      column<std::string>("candidate_id"), column<std::size_t>("momentum_index"), column<R>("q"), column<R>("gap"),
      column<R>("weight"), column<int>("multiplicity"));
  output.table(
      "diagnostics", "Every requested label, including rejected states",
      [&](auto& table) {
        auto append = [&](model::FourSpinonRequest const& r, std::string solve_status, std::string factor_status,
                          bool accepted, Optional delta, Optional logarithm, Optional residual, Optional pivot,
                          std::string diagnostic) {
          table.append(model::four_spinon_id(r), r.holes[0], r.holes[1], r.holes[2], r.holes[3], r.branch,
                       model::four_spinon_transfer_index(r), solve_status, factor_status, accepted, delta, logarithm,
                       residual, pivot, diagnostic);
        };
        if (scan)
          for (auto const& row : scan->records)
          {
            auto const& state = row.solved.state;
            append(row.solved.request, status(row.solved.status),
                   row.factor ? status(row.factor->status) : "not_attempted", row.line.has_value(),
                   state ? Optional(state->deviation()) : Optional{},
                   state ? Optional(state->log_abs_deviation()) : Optional{},
                   state ? Optional(state->rational_residual()) : Optional{},
                   row.factor ? row.factor->pivot_margin : Optional{}, row.diagnostic);
          }
        else
          for (auto const& r : candidates)
            append(r, "not_attempted", "not_attempted", false, {}, {}, {}, {}, "singlet preparation failed");
      },
      column<std::string>("candidate_id"), column<int>("hole_0"), column<int>("hole_1"), column<int>("hole_2"),
      column<int>("hole_3"), column<int>("branch"), column<int>("momentum_index"), column<std::string>("root_status"),
      column<std::string>("factor_status"), column<bool>("accepted"), column<Optional>("deviation"),
      column<Optional>("log_abs_deviation"), column<Optional>("rational_residual"), column<Optional>("pivot_margin"),
      column<std::string>("diagnostic"));
  output.table(
      "moments", "Retained moments and full first-moment rule",
      [&](auto& table) {
        if (scan)
          for (int q = 0; q < a.sites; ++q)
          {
            R const full = scale * bethe::heisenberg::raising_first_moment(a.sites, q, scan->ground_energy);
            Optional fraction;
            if (full > zero) fraction = scale * scan->moments[q].first_moment / full;
            table.append(q, 2 * pi * q / a.sites, scale * scan->moments[q].weight,
                         scale * scan->moments[q].first_moment, full, fraction);
          }
      },
      column<int>("momentum_index"), column<R>("q"), column<R>("retained_weight"), column<R>("retained_first_moment"),
      column<R>("full_first_moment"), column<Optional>("first_moment_fraction"));
  output.table(
      "roots", "Standard ABA roots; tiny deviations are in diagnostics",
      [&](auto& table) {
        std::size_t j = 0;
        for (auto const& root : ground.aba_rapidities())
          table.append(std::string("ground"), j++, root, zero);
        if (scan)
          for (auto const& row : scan->records)
            if (row.solved.state)
            {
              auto const roots = row.solved.state->roots();
              for (std::size_t i = 0; i < roots.size(); ++i)
                table.append(row.candidate_id, i, roots[i].real(), roots[i].imag());
            }
      },
      column<std::string>("candidate_id"), column<std::size_t>("root_index"), column<R>("lambda_real"),
      column<R>("lambda_imag"));
  output.finish();
  if (!complete)
    std::cerr << "Some requested weights are unavailable; inspect diagnostics. No missing-weight rescaling.\n";
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
        return cli::dispatch_options_precision<cli::ArbitraryPrecision::complex_linalg>(
            a.precision, [&](auto const& controls) { return run(a, argc, argv, controls); });
      });
}
