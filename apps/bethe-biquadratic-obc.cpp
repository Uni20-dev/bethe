// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "excitation-report.hpp"
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/biquadratic.hpp>
#include <bethe/biquadratic_qsystem.hpp>
#include <bethe/biquadratic_two_string.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::biquadratic;
struct Arguments
{
    std::size_t sites = 0, max_iterations = 10000;
    std::optional<std::size_t> through_lines, excitations, max_candidates;
    std::optional<bethe::xxz::QuantumNumbers> numbers;
    std::optional<std::string> tolerance;
    std::optional<std::string> q_seed;
    std::optional<std::size_t> max_attempts;
    std::string precision = "fp64";
    cli::DataOutputOptions output;
    bool roots = false, sectors = false;
    bool q_spectrum = false, singlet_excitation = false;
};
auto program_info()
{
  auto info = bethe::cli::program_info("bethe-biquadratic-obc",
                                       "Spin-1 pure biquadratic chain, free ends: H=-sum_i (S_i.S_(i+1))^2.",
                                       bethe::citations::Tool::biquadratic_obc);
  info.notes = {"Unique singlet ground state; even N>=2, coefficient -1.",
                "TL loop weight 3; reference XXZ Delta=3/2 with opposite end fields.",
                "Real-root scans are NOT complete spectra: complex-root levels are excluded.",
                "--q-spectrum searches real and complex levels of one TL module; validated through N=8.",
                "Q-system modes have no site cutoff; larger sizes are experimental and may be costly or unresolved.",
                "--q-seed selects a polynomial branch, not necessarily a low-lying state.",
                "--singlet-excitation targets one two-string singlet above a real sea, including on long chains.",
                "Multiplicity counts physical states per TL eigenvector, not SU(2) multiplets.",
                "TL through-lines are not physical spin; no odd chains or lattice momentum.",
                "This is not the TB point, ULS point, or zero-boundary-field XXZ chain.",
                "See docs/biquadratic.md for the TL mapping and representation multiplicities.",
                "Use --references for literature and applicability; see CITATIONS.md."};
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  bethe::cli::option(app, "N", args.sites, "Number of particles or sites")->required();
  bethe::cli::option(app, "--through-lines", args.through_lines, "lowest level in an even TL module, 0<=ELL<=N");
  bethe::cli::option(app, "--sectors", args.sectors, "lowest level in every TL module");
  bethe::cli::all_count_option(
      app, "--excitations", args.excitations,
      "lowest COUNT, or all, supported real-root levels default ELL=2; includes module minimum");
  bethe::cli::option(app, "--max-candidates", args.max_candidates, "exhaustive scan limit (default: 10000)");
  bethe::cli::option(app, "--quantum-numbers", args.numbers, "explicit integer labels; none for vacuum");
  bethe::cli::option(app, "--q-spectrum", args.q_spectrum,
                     "budgeted Q-system search including complex roots; default ELL=0; larger N is experimental");
  bethe::cli::option(app, "--q-seed", args.q_seed,
                     "selected Q-system branch: c0,c1,... for monic Q(x), x=cosh(2u); none for vacuum");
  bethe::cli::option(app, "--max-attempts", args.max_attempts, "Q-system seed budget (default: 4000)");
  bethe::cli::option(app, "--singlet-excitation", args.singlet_excitation,
                     "target the low-lying complex-root singlet (one two-string); even N>=4; not a spectrum scan");
  bethe::cli::option(app, "--roots", args.roots, "print reference roots (real labels or complex coordinates)");
  bethe::cli::option(
      app, "--tolerance", args.tolerance,
      "equation tolerance (default: 32 epsilon); Q-system uses coefficient backward error, not phase error");
  bethe::cli::option(app, "--max-iterations", args.max_iterations, "accepted Newton updates (default: 10000)")
      ->capture_default_str();
  bethe::cli::precision_option(app, args.precision);
  cli::add_data_output_options(app, args.output, true);
}

void validate(Arguments const& args)
{
  args.output.validate();
  if (args.singlet_excitation && (args.q_seed || args.q_spectrum || args.sectors || args.excitations || args.numbers ||
                                  args.through_lines || args.max_candidates || args.max_attempts))
    throw std::invalid_argument("--singlet-excitation cannot combine with other state selections or scan budgets");
  if ((args.q_seed || args.q_spectrum) && (args.sectors || args.excitations || args.numbers || args.max_candidates))
    throw std::invalid_argument("Q-system modes cannot be combined with real-root selection or --sectors");
  if (args.q_seed && (args.q_spectrum || args.through_lines))
    throw std::invalid_argument("--q-seed determines the module; cannot combine with --q-spectrum or --through-lines");
  if (args.max_attempts && !args.q_spectrum) throw std::invalid_argument("--max-attempts requires --q-spectrum");
  if (args.sectors && (args.through_lines || args.excitations || args.numbers))
    throw std::invalid_argument(
        "--sectors cannot be combined with --through-lines, --excitations or --quantum-numbers");
  if (args.numbers && (args.through_lines || args.excitations))
    throw std::invalid_argument(
        "--quantum-numbers determines the TL module; cannot combine with --through-lines or --excitations");
  if (args.max_candidates && !args.excitations) throw std::invalid_argument("--max-candidates requires --excitations");
}
char const* status(bethe::xxz::quantum_group::SolveStatus value)
{
  using Status = bethe::xxz::quantum_group::SolveStatus;
  switch (value)
  {
    case Status::converged:
      return "converged";
    case Status::iteration_limit:
      return "iteration limit reached; unconverged estimate";
    case Status::singular_jacobian:
      return "singular or ill-conditioned Jacobian; unconverged estimate";
    case Status::stalled:
      return "line search or representable precision stalled; unconverged estimate";
  }
  return "unknown";
}
std::string multiplicity_text(std::optional<std::uint64_t> value)
{
  return value ? std::to_string(*value) : "overflow (>uint64)";
}

char const* status(bethe::xxz::quantum_group::qsystem::Status value)
{
  using S = bethe::xxz::quantum_group::qsystem::Status;
  switch (value)
  {
    case S::converged:
      return "converged; numerical admissibility checks passed";
    case S::iteration_limit:
      return "iteration limit reached; unconverged estimate";
    case S::singular_jacobian:
      return "singular Jacobian; unconverged estimate";
    case S::stalled:
      return "line search stalled; unconverged estimate";
    case S::unresolved_roots:
      return "unresolved polynomial roots; unverified estimate";
    case S::inadmissible:
      return "root admissibility or original Bethe-equation check failed";
  }
  return "unknown";
}

template <uni20::Real Real> std::vector<Real> parse_q_seed(std::string_view text)
{
  if (text == "none") return {};
  std::vector<Real> result;
  for (;;)
  {
    auto const comma = text.find(',');
    result.push_back(uni20::parse_real<Real>(std::string(text.substr(0, comma))));
    if (comma == std::string_view::npos) return result;
    text.remove_prefix(comma + 1);
  }
}

template <uni20::Real Real>
int run_qsystem(Arguments const& args, int argc, char** argv, bethe::SolverOptions<Real> const& options,
                cli::report_builder report)
{
  cli::CpuTimer const timer;
  if (args.sites > 8)
    std::cerr << "Warning: Q-system spectrum validation covers N<=8. Larger chains are experimental; "
                 "memory, search cost and string conditioning may prevent convergence or completeness.\n";
  std::vector<model::qsystem::State<Real>> states;
  bool complete = false;
  if (args.q_seed)
  {
    states.push_back(model::qsystem::solve<Real>(args.sites, parse_q_seed<Real>(*args.q_seed), options));
    complete = states.front().reference.converged;
    report.field("Calculation", "selected Q-system branch; not necessarily a lowest level");
  }
  else
  {
    auto scan = model::qsystem::spectrum<Real>(args.sites, args.through_lines.value_or(0),
                                               {.max_attempts = args.max_attempts.value_or(4000)}, options);
    complete = scan.complete();
    report.field("Calculation", "Q-system spectrum search (real and complex roots)")
        .field("Expected module dimension", scan.expected_count ? std::to_string(*scan.expected_count)
                                                                : "overflow (>size_t); completeness unavailable")
        .field("Discovered levels", scan.states.size())
        .field("Attempts", scan.attempts)
        .field("Attempt budget", args.max_attempts.value_or(4000))
        .field("Failed attempts", scan.failed_attempts)
        .field("Ordering", complete ? "numerically complete TL module; count matched, not a rigorous certificate"
                                    : "incomplete discoveries; NOT guaranteed lowest levels");
    states = std::move(scan.states);
  }
  auto const ground = model::ground_state<Real>(args.sites, options);
  complete = complete && ground.reference.converged;
  report.field("TL through-lines", args.q_seed ? states.front().through_lines : args.through_lines.value_or(0))
      .field("Q-system validation", args.sites <= 8
                                        ? "within small-chain regression range"
                                        : "experimental beyond N=8; no completeness or convergence guarantee")
      .field("Residual convention", "Q-system Wronskian coefficient backward error; not an energy-error bound")
      .field("Root coordinate", "x=cosh(2u)=cos(alpha); Bajnok u, alpha=-2iu")
      .field("Multiplicity meaning", "physical states per TL eigenvector, not SU(2) multiplets")
      .field("Gap reference",
             ground.reference.converged ? "E-E0; global singlet ground state" : "unavailable; ground solve failed")
      .field("Status", complete ? "converged" : "incomplete or unverified")
      .field("CPU time", timer.elapsed_text());
  std::vector<std::string> names{"states", "reference", "q_coefficients"};
  if (args.roots) names.push_back("roots");
  cli::ResultOutput output(report, args.output, "bethe-biquadratic-obc", argc, argv, names);
  output.table(
      "states", "Q-system levels",
      [&](auto& table) {
        for (std::size_t i = 0; i < states.size(); ++i)
        {
          auto const& s = states[i];
          auto const& r = s.reference;
          std::optional<Real> gap;
          if (r.converged && ground.reference.converged) gap = Real{2} * (*r.energy - ground.reference.energy);
          table.append(i, s.through_lines, s.multiplicity, s.energy, gap, s.tl_energy, r.energy, r.residual_norm,
                       uni20::isfinite(r.bethe_residual) ? std::optional<Real>{r.bethe_residual} : std::nullopt,
                       uni20::isfinite(r.bethe_residual_bound) ? std::optional<Real>{r.bethe_residual_bound}
                                                               : std::nullopt,
                       r.iterations, r.converged, std::string(status(r.status)));
        }
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("through_lines", "Through-lines"),
      cli::column<std::optional<std::uint64_t>>("multiplicity", "Multiplicity"),
      cli::column<std::optional<Real>>("energy", "Energy"), cli::column<std::optional<Real>>("gap", "E-E0"),
      cli::column<std::optional<Real>>("tl_energy"), cli::column<std::optional<Real>>("reference_energy"),
      cli::column<Real>("residual", "Wronskian residual"),
      cli::column<std::optional<Real>>("bethe_residual", "Bethe residual"),
      cli::column<std::optional<Real>>("bethe_residual_bound", "Bethe uncertainty bound"),
      cli::column<std::size_t>("iterations"), cli::column<bool>("converged"),
      cli::column<std::string>("status", "Status"));
  output.table(
      "reference", "Global ground reference",
      [&](auto& table) {
        table.append(states.size(), ground.energy, ground.reference.residual_norm, ground.reference.converged,
                     std::string(status(ground.reference.status)));
      },
      cli::column<std::size_t>("state_id"), cli::column<Real>("energy", "Energy"), cli::column<Real>("residual"),
      cli::column<bool>("converged"), cli::column<std::string>("status", "Status"));
  output.table(
      "q_coefficients", "Monic Q(x) coefficients (ascending powers; leading 1 omitted)",
      [&](auto& table) {
        for (std::size_t i = 0; i < states.size(); ++i)
          for (std::size_t k = 0; k < states[i].reference.coefficients.size(); ++k)
            table.append(i, k, states[i].reference.coefficients[k]);
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("power"), cli::column<Real>("coefficient"));
  if (args.roots)
    output.table(
        "roots", "Complex Q-system roots (x=cosh(2u))",
        [&](auto& table) {
          for (std::size_t i = 0; i < states.size(); ++i)
            for (std::size_t k = 0; k < states[i].reference.roots.roots.size(); ++k)
            {
              auto const& r = states[i].reference;
              auto x = r.roots.roots[k];
              std::optional<Real> ur, ui;
              if (k < r.rapidities.size())
              {
                ur = r.rapidities[k].real();
                ui = r.rapidities[k].imag();
              }
              auto finite = [](Real value) {
                return uni20::isfinite(value) ? std::optional<Real>{value} : std::nullopt;
              };
              table.append(i, k, finite(x.real()), finite(x.imag()), ur ? finite(*ur) : std::nullopt,
                           ui ? finite(*ui) : std::nullopt);
            }
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("index"),
        cli::column<std::optional<Real>>("x_real"), cli::column<std::optional<Real>>("x_imag"),
        cli::column<std::optional<Real>>("u_real"), cli::column<std::optional<Real>>("u_imag"));
  output.finish();
  if (!complete) std::cerr << "Q-system calculation incomplete or unverified; no lowest-level guarantee.\n";
  return complete ? 0 : 2;
}

template <uni20::Real Real>
int run_singlet(Arguments const& args, int argc, char** argv, bethe::SolverOptions<Real> const& options,
                cli::report_builder report)
{
  cli::CpuTimer const timer;
  auto const state = model::two_string::singlet<Real>(args.sites, options);
  auto const ground = model::ground_state<Real>(args.sites, options);
  auto const& r = state.reference;
  bool const complete = r.converged && ground.reference.converged;
  Real const eta = std::acosh(r.delta), d = std::exp(-r.log_deviation);
  std::optional<Real> gap;
  if (complete) gap = Real{2} * (r.energy - ground.reference.energy);
  report.field("Calculation", "selected complex-root singlet excitation")
      .field("Family", "one positive-deviation two-string; real I=1,...,N/2-2; string label 1")
      .field("Ordering", "targeted branch, not an exhaustive search or a global first-excitation guarantee")
      .field("TL through-lines", 0)
      .field("Multiplicity meaning", "one physical singlet per TL eigenvector")
      .field("Residual convention", "max phase/log-modulus equation residual divided by 2N; not an energy-error bound")
      .field("String coordinate", "u=(eta+d)/2 +/- i*a/2; d=exp(-L)>0; L is authoritative")
      .field("String deviation", eta + d == eta ? "unresolved in rounded u; retained by L=-log(d)"
                                                : "resolved in rounded u; L=-log(d) also retained")
      .field("Gap reference",
             ground.reference.converged ? "E-E0; global singlet ground state" : "unavailable; ground solve failed")
      .field("Status", complete ? "converged" : "incomplete; unconverged estimate")
      .field("CPU time", timer.elapsed_text());
  std::vector<std::string> names{"states", "reference", "string"};
  if (args.roots) names.push_back("roots");
  cli::ResultOutput output(report, args.output, "bethe-biquadratic-obc", argc, argv, names);
  output.table(
      "states", "Selected two-string singlet",
      [&](auto& table) {
        table.append(std::size_t{0}, state.through_lines, state.multiplicity, state.energy, gap, state.tl_energy,
                     r.energy, r.residual_norm, r.phase_residual, r.modulus_residual, r.iterations, r.converged,
                     std::string(status(r.status)));
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("through_lines"),
      cli::column<std::uint64_t>("multiplicity"), cli::column<Real>("energy", "Energy"),
      cli::column<std::optional<Real>>("gap", "E-E0"), cli::column<Real>("tl_energy"),
      cli::column<Real>("reference_energy"), cli::column<Real>("residual"), cli::column<Real>("phase_residual"),
      cli::column<Real>("modulus_residual"), cli::column<std::size_t>("iterations"), cli::column<bool>("converged"),
      cli::column<std::string>("status", "Status"));
  output.table(
      "reference", "Global ground reference",
      [&](auto& table) {
        table.append(std::size_t{1}, ground.energy, ground.reference.residual_norm, ground.reference.converged,
                     std::string(status(ground.reference.status)));
      },
      cli::column<std::size_t>("state_id"), cli::column<Real>("energy", "Energy"), cli::column<Real>("residual"),
      cli::column<bool>("converged"), cli::column<std::string>("status", "Status"));
  output.table(
      "string", "Two-string parameters (L remains valid below root-coordinate resolution)",
      [&](auto& table) {
        table.append(std::size_t{0}, r.center, r.log_deviation, d > Real{0} ? std::optional<Real>{d} : std::nullopt);
      },
      cli::column<std::size_t>("state_id"), cli::column<Real>("center", "a"),
      cli::column<Real>("log_deviation", "L=-log(d)"),
      cli::column<std::optional<Real>>("deviation", "d (null if underflow)"));
  if (args.roots)
    output.table(
        "roots", "Reference roots (rounded u; use L for the string deviation)",
        [&](auto& table) {
          for (std::size_t j = 0; j < r.rapidities.size(); ++j)
            table.append(std::size_t{0}, j, std::string("real sea"), Real{0}, r.rapidities[j] / Real{2});
          for (std::size_t j = 0; j < 2; ++j)
            table.append(std::size_t{0}, r.rapidities.size() + j, std::string("two-string"), (eta + d) / Real{2},
                         (j ? -r.center : r.center) / Real{2});
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("index"), cli::column<std::string>("kind"),
        cli::column<Real>("u_real"), cli::column<Real>("u_imag"));
  output.finish();
  if (!complete) std::cerr << "Two-string singlet or ground reference unconverged; no verified gap.\n";
  return complete ? 0 : 2;
}

template <uni20::Real Real>
void write_output(cli::report_builder const& report, Arguments const& args, int argc, char** argv,
                  std::vector<model::State<Real> const*> const& states, std::vector<std::optional<Real>> const& gaps,
                  model::State<Real> const* reference = nullptr, model::State<Real> const* failed = nullptr)
{
  std::vector<std::string> names{"states", "quantum_numbers"};
  if (reference) names.push_back("reference");
  if (failed) names.push_back("failed");
  if (args.roots) names.push_back("roots");
  cli::ResultOutput output(report, args.output, "bethe-biquadratic-obc", argc, argv, names);
  auto write_states = [&](std::string name, std::string title, auto const& rows, std::size_t offset, bool ranked) {
    output.table(
        name, title,
        [&](auto& table) {
          for (std::size_t i = 0; i < rows.size(); ++i)
          {
            auto const& s = *rows[i];
            auto const& r = s.reference;
            table.append(offset + i, s.through_lines, s.multiplicity, s.energy, ranked ? gaps[i] : std::nullopt,
                         s.tl_energy, r.energy, r.residual_norm, r.iterations, r.converged,
                         std::string(status(r.status)));
          }
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("through_lines", "Through-lines"),
        cli::column<std::optional<std::uint64_t>>("multiplicity", "Multiplicity"),
        cli::column<Real>("energy", "Energy"), cli::column<std::optional<Real>>("gap", "E-E0"),
        cli::column<Real>("tl_energy"), cli::column<Real>("reference_energy"),
        cli::column<Real>("residual", "Residual"), cli::column<std::size_t>("iterations", "Iterations"),
        cli::column<bool>("converged"), cli::column<std::string>("status", "Status"));
  };
  write_states("states",
               args.excitations ? "Real-root TL levels (module minimum included)"
               : args.sectors   ? "TL module minima"
                                : "State",
               states, 0, true);
  auto all = states;
  if (reference)
  {
    write_states("reference", "Ground reference", std::vector{reference}, all.size(), false);
    all.push_back(reference);
  }
  if (failed)
  {
    write_states("failed", "First failed state (unranked estimate)", std::vector{failed}, all.size(), false);
    all.push_back(failed);
  }
  output.table(
      "quantum_numbers", "Reference quantum numbers (not physical spin-1 labels)",
      [&](auto& table) {
        for (std::size_t i = 0; i < all.size(); ++i)
          for (std::size_t j = 0; j < all[i]->reference.quantum_numbers.size(); ++j)
            table.append(i, j, all[i]->reference.quantum_numbers[j]);
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("index", "Index"),
      cli::column<uni20::half_int>("quantum_number", "I"));
  if (args.roots)
    output.table(
        "roots", "Reference XXZ roots (not physical spin-1 quantum numbers)",
        [&](auto& table) {
          for (std::size_t i = 0; i < all.size(); ++i)
            for (std::size_t j = 0; j < all[i]->reference.rapidities.size(); ++j)
              table.append(i, j, all[i]->reference.quantum_numbers[j], all[i]->reference.rapidities[j],
                           all[i]->reference.angles[j]);
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("index", "Index"),
        cli::column<uni20::half_int>("quantum_number", "I"), cli::column<Real>("rapidity", "alpha"),
        cli::column<Real>("angle", "x=Theta_1/2"));
  output.finish();
}

template <uni20::Real Real> auto preamble(Arguments const& args, bethe::SolverOptions<Real> const& options)
{
  cli::report_builder report("Spin-1 pure biquadratic chain (free ends)");
  report.field("Hamiltonian", "H=-sum_i (S_i.S_(i+1))^2")
      .field("Sites", args.sites)
      .field("Spin", 1)
      .field("TL loop weight", 3)
      .field("XXZ Delta", "1.5")
      .field("XXZ reference", "spin-half exchange 1; +sqrt(5)/4*(sz_1-sz_N)")
      .field("Precision", args.precision)
      .field("Residual tolerance", uni20::format_real(options.residual_tolerance));
  return report;
}

int finish(bool converged)
{
  if (!converged) std::cerr << "Biquadratic solve incomplete; consider a larger budget or higher precision.\n";
  return converged ? 0 : 2;
}

template <uni20::Real Real> int run(Arguments const& args, int argc, char** argv)
{
  bethe::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  auto report = preamble(args, options);
  if (args.q_seed || args.q_spectrum) return run_qsystem(args, argc, argv, options, std::move(report));
  if (args.singlet_excitation) return run_singlet(args, argc, argv, options, std::move(report));
  cli::CpuTimer const timer;
  if (args.excitations)
  {
    auto const ell = args.through_lines.value_or(2);
    auto const scan = model::real_excitations<Real>(
        args.sites, ell, {.count = *args.excitations, .max_candidates = args.max_candidates.value_or(10000)}, options);
    auto const cpu_time = timer.elapsed_text();
    auto const& ground = scan.ground_state;
    report.field("Calculation", "restricted real-root excitations")
        .field("Family", "positive finite roots; complex-root levels excluded; NOT the complete TL spectrum")
        .field("TL through-lines", ell)
        .field("Multiplicity per TL eigenvector",
               multiplicity_text(bethe::temperley_lieb::spin_chain_multiplicity(3, ell)))
        .field("Multiplicity meaning", "physical states, not SU(2) multiplets; not an accidental-degeneracy sum")
        .field("Candidates", scan.candidate_count)
        .field("Converged candidates", scan.converged_count)
        .field("Returned levels", scan.levels.size())
        .field("Ordering",
               scan.family_converged() ? "complete within supported family" : "incomplete; failed candidates excluded")
        .field("Ground energy", uni20::format_real(ground.energy))
        .field("Ground status", status(ground.reference.status))
        .field("Ground residual", uni20::format_real(ground.reference.residual_norm))
        .field("Ground iterations", ground.reference.iterations)
        .field("Gap reference",
               ground.reference.converged ? "E-E0; global singlet ground state" : "unavailable; ground solve failed")
        .field("Status", scan.converged() ? "converged" : "incomplete scan or ground reference")
        .field("CPU time", cpu_time);
    if (scan.first_unconverged)
      report.field("First failed I", cli::quantum_number_text(scan.first_unconverged->reference.quantum_numbers))
          .field("First failed status", status(scan.first_unconverged->reference.status))
          .field("First failed residual", uni20::format_real(scan.first_unconverged->reference.residual_norm));
    std::vector<model::State<Real> const*> states;
    std::vector<std::optional<Real>> gaps;
    for (auto const& level : scan.levels)
    {
      states.push_back(&level.state);
      gaps.push_back(level.gap);
    }
    write_output(report, args, argc, argv, states, gaps, &ground,
                 scan.first_unconverged ? &*scan.first_unconverged : nullptr);
    return finish(scan.converged());
  }
  if (args.sectors)
  {
    std::vector<model::State<Real>> states;
    states.push_back(model::ground_state<Real>(args.sites, options));
    bool converged = states.front().reference.converged;
    for (std::size_t ell = 2; ell <= args.sites; ell += 2)
    {
      states.push_back(model::sector_ground_state<Real>(args.sites, ell, options));
      converged = converged && states.back().reference.converged;
    }
    auto const cpu_time = timer.elapsed_text();
    report.field("Calculation", "TL module minima (not physical-spin sectors)")
        .field("Multiplicity meaning", "physical states per TL eigenvector, not SU(2) multiplets")
        .field("Status", converged ? "converged" : "incomplete; unconverged estimates")
        .field("CPU time", cpu_time);
    std::vector<model::State<Real> const*> rows;
    std::vector<std::optional<Real>> gaps;
    for (auto const& state : states)
    {
      rows.push_back(&state);
      gaps.push_back(state.reference.converged && states.front().reference.converged
                         ? std::optional<Real>{Real{2} * (state.reference.energy - states.front().reference.energy)}
                         : std::nullopt);
    }
    write_output(report, args, argc, argv, rows, gaps);
    return finish(converged);
  }
  auto const state = args.numbers
                         ? model::solve_real<Real>(args.sites, *args.numbers, options)
                         : model::sector_ground_state<Real>(args.sites, args.through_lines.value_or(0), options);
  auto const cpu_time = timer.elapsed_text();
  auto const& reference = state.reference;
  report
      .status(reference.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning,
              status(reference.status))
      .field("Calculation", args.numbers               ? "specified real-root TL level"
                            : state.through_lines == 0 ? "even-chain singlet ground state"
                                                       : "TL module minimum")
      .field("TL through-lines", state.through_lines)
      .field(state.through_lines == 0 ? "Ground-state multiplicity" : "Multiplicity per TL eigenvector",
             multiplicity_text(state.multiplicity))
      .field("Status", status(reference.status))
      .field("Total energy", uni20::format_real(state.energy))
      .field("Energy per site", uni20::format_real(state.energy / Real(args.sites)))
      .field("TL energy (-sum e_i)", uni20::format_real(state.tl_energy))
      .field("XXZ reference energy", uni20::format_real(reference.energy))
      .field("Residual norm", uni20::format_real(reference.residual_norm))
      .field("Iterations", reference.iterations)
      .field("CPU time", cpu_time);
  if (state.through_lines == 0) report.field("Total spin", 0);
  if (state.through_lines != 0) report.field("Multiplicity meaning", "physical states, not a physical-spin label");
  write_output<Real>(report, args, argc, argv, {&state}, {std::nullopt});
  return finish(reference.converged);
}
} // namespace
int main(int argc, char** argv)
{
  Arguments args;
  return bethe::cli::program_main(
      argc, argv, program_info(), [&](auto& app) { add_options(app, args); },
      [&](auto&) {
        validate(args);
        return bethe::cli::dispatch_precision(args.precision,
                                              [&]<uni20::Real Real> { return run<Real>(args, argc, argv); });
      });
}
