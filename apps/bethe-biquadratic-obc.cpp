// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "excitation-report.hpp"
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/biquadratic.hpp>
#include <bethe/biquadratic_ferromagnetic.hpp>
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
    bool ferromagnetic = false, one_defect = false;
};
auto program_info()
{
  auto info = bethe::cli::program_info("bethe-biquadratic-obc",
                                       "Spin-1 pure biquadratic chain, free ends: H=-sum_i (S_i.S_(i+1))^2 by default.",
                                       bethe::citations::Tool::biquadratic_obc);
  info.notes = {"Default: unique singlet ground state; even N>=2, coefficient -1.",
                "--ferromagnetic reverses the sign: H=+sum (S.S)^2; gaps use the exact degenerate ground space.",
                "--ferromagnetic --one-defect gives the complete ell=N-2 band analytically, including odd N.",
                "Ferromagnetic real-root scans exclude complex-root levels; NOT general module minima.",
                "TL loop weight 3; reference XXZ Delta=3/2 with opposite end fields.",
                "Real-root scans are NOT complete spectra: complex-root levels are excluded.",
                "--q-spectrum searches real and complex levels of one TL module; validated through N=8.",
                "Q-system modes have no site cutoff; larger sizes are experimental and may be costly or unresolved.",
                "--q-seed selects a polynomial branch, not necessarily a low-lying state.",
                "--singlet-excitation targets one two-string singlet above a real sea, including on long chains.",
                "Multiplicity counts physical states per TL eigenvector, not SU(2) multiplets.",
                "TL through-lines are not physical spin; numerical root modes require even N; no lattice momentum.",
                "This is not the TB point, ULS point, or zero-boundary-field XXZ chain.",
                "See docs/biquadratic.md for the TL mapping and representation multiplicities.",
                "Use --references for literature and applicability; see CITATIONS.md."};
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  bethe::cli::option(app, "N", args.sites, "Number of particles or sites")->required();
  bethe::cli::option(app, "--ferromagnetic", args.ferromagnetic,
                     "use H=+sum (S.S)^2 and the exact ferro ground reference");
  bethe::cli::option(app, "--one-defect", args.one_defect, "complete analytic ferro ell=N-2 band; odd/even N>=2");
  bethe::cli::option(app, "--through-lines", args.through_lines,
                     "select TL module; standalone ferro minima require ELL=N or N-2");
  bethe::cli::option(app, "--sectors", args.sectors, "AF: lowest level in every TL module");
  bethe::cli::all_count_option(
      app, "--excitations", args.excitations,
      "lowest COUNT, or all, supported real-root levels; default ELL=2 (ferro: N-2); excludes complex roots");
  bethe::cli::option(app, "--max-candidates", args.max_candidates,
                     "real-family scan / analytic band row limit (default: 10000)");
  bethe::cli::option(app, "--quantum-numbers", args.numbers, "explicit integer labels; none for vacuum");
  bethe::cli::option(
      app, "--q-spectrum", args.q_spectrum,
      "budgeted Q-system search including complex roots; default ELL=0 (ferro: N-2); larger N is experimental");
  bethe::cli::option(app, "--q-seed", args.q_seed,
                     "selected Q-system branch: c0,c1,... for monic Q(x), x=cosh(2u); none for vacuum");
  bethe::cli::option(app, "--max-attempts", args.max_attempts, "Q-system seed budget (default: 4000)");
  bethe::cli::option(app, "--singlet-excitation", args.singlet_excitation,
                     "AF: target the low-lying complex-root singlet (one two-string); even N>=4; not a spectrum scan");
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
  if (args.one_defect && !args.ferromagnetic) throw std::invalid_argument("--one-defect requires --ferromagnetic");
  if (args.one_defect && (args.q_seed || args.q_spectrum || args.sectors || args.excitations || args.numbers ||
                          args.through_lines || args.singlet_excitation || args.max_attempts))
    throw std::invalid_argument("--one-defect cannot combine with other state selections");
  if (args.ferromagnetic)
  {
    if (args.sites < 2) throw std::invalid_argument("ferromagnetic free-end chain requires N>=2");
    if (args.sectors || args.singlet_excitation)
      throw std::invalid_argument(
          "ferromagnetic --sectors and --singlet-excitation are not implemented; use --q-spectrum --through-lines ELL");
    bool const analytic = !args.numbers && !args.excitations && !args.q_seed && !args.q_spectrum;
    if (analytic && args.through_lines && *args.through_lines != args.sites && *args.through_lines != args.sites - 2)
      throw std::invalid_argument(
          "ferromagnetic module minima are analytic only for ELL=N or N-2; use --q-spectrum --through-lines ELL");
    if (analytic && args.roots)
      throw std::invalid_argument(
          "analytic ferro modes have no root table; use --quantum-numbers or --excitations for XXZ roots");
  }
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
  if (args.max_candidates && !args.excitations && !args.one_defect)
    throw std::invalid_argument("--max-candidates requires --excitations or --one-defect");
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
int run_qsystem(Arguments const& args, bethe::SolverOptions<Real> const& options, cli::RunReport report)
{
  auto& context = report.context();
  if (args.sites > 8)
    std::cerr << "Warning: Q-system spectrum validation covers N<=8. Larger chains are experimental; "
                 "memory, search cost and string conditioning may prevent convergence or completeness.\n";
  std::vector<model::qsystem::State<Real>> states;
  auto const ell = args.through_lines.value_or(args.ferromagnetic ? args.sites - 2 : 0);
  bool complete = false;
  if (args.q_seed)
  {
    states.push_back(context.measure([&] {
      auto const seed = parse_q_seed<Real>(*args.q_seed);
      return args.ferromagnetic ? model::ferromagnetic::qsystem::solve<Real>(args.sites, seed, options)
                                : model::qsystem::solve<Real>(args.sites, seed, options);
    }));
    complete = states.front().reference.converged;
    report.field("calculation", "Calculation", "selected Q-system branch; not necessarily a lowest level");
  }
  else
  {
    auto scan = context.measure([&] {
      bethe::xxz::quantum_group::qsystem::SearchOptions search{.max_attempts = args.max_attempts.value_or(4000)};
      return args.ferromagnetic ? model::ferromagnetic::qsystem::spectrum<Real>(args.sites, ell, search, options)
                                : model::qsystem::spectrum<Real>(args.sites, ell, search, options);
    });
    complete = scan.complete();
    report.field("calculation", "Calculation", "Q-system spectrum search (real and complex roots)")
        .field("expected_module_dimension", "Expected module dimension", scan.expected_count,
               {.missing = "overflow (>size_t); completeness unavailable"})
        .field("discovered_levels", "Discovered levels", scan.states.size())
        .field("attempts", "Attempts", scan.attempts)
        .field("attempt_budget", "Attempt budget", args.max_attempts.value_or(4000))
        .field("failed_attempts", "Failed attempts", scan.failed_attempts)
        .field("ordering", "Ordering",
               complete ? "numerically complete TL module; count matched, not a rigorous certificate"
                        : "incomplete discoveries; NOT guaranteed lowest levels");
    states = std::move(scan.states);
  }
  auto const ground = context.measure([&] {
    return args.ferromagnetic ? model::ferromagnetic::solve_real<Real>(args.sites, {}, options)
                              : model::ground_state<Real>(args.sites, options);
  });
  complete = complete && ground.reference.converged;
  report.field("tl_through_lines", "TL through-lines", args.q_seed ? states.front().through_lines : ell)
      .field("tl_defects", "TL defects M", (args.sites - (args.q_seed ? states.front().through_lines : ell)) / 2)
      .field("q_system_validation", "Q-system validation",
             args.sites <= 8 ? "within small-chain regression range"
                             : "experimental beyond N=8; no completeness or convergence guarantee")
      .field("residual_convention", "Residual convention",
             "Q-system Wronskian coefficient backward error; not an energy-error bound")
      .field("root_coordinate", "Root coordinate", "x=cosh(2u)=cos(alpha); Bajnok u, alpha=-2iu")
      .field("multiplicity_meaning", "Multiplicity meaning", "physical states per TL eigenvector, not SU(2) multiplets")
      .field("gap_reference", "Gap reference",
             args.ferromagnetic           ? "E-(N-1); exact degenerate ferro ground space"
             : ground.reference.converged ? "E-E0; global singlet ground state"
                                          : "unavailable; ground solve failed")
      .result(complete, complete ? "converged" : "incomplete or unverified");
  std::vector<std::string> names{"states", "reference", "q_coefficients"};
  if (args.roots) names.push_back("roots");
  cli::ResultOutput output(report, args.output, names);
  output.table(
      "states", "Q-system levels",
      [&](auto& table) {
        for (std::size_t i = 0; i < states.size(); ++i)
        {
          auto const& s = states[i];
          auto const& r = s.reference;
          std::optional<Real> gap;
          if (r.converged && ground.reference.converged)
            gap =
                args.ferromagnetic ? s.tl_energy : std::optional<Real>{Real{2} * (*r.energy - ground.reference.energy)};
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
int run_singlet(Arguments const& args, bethe::SolverOptions<Real> const& options, cli::RunReport report)
{
  auto& context = report.context();
  auto computation = context.computation();
  auto const state = model::two_string::singlet<Real>(args.sites, options);
  auto const ground = model::ground_state<Real>(args.sites, options);
  computation.finish();
  auto const& r = state.reference;
  bool const complete = r.converged && ground.reference.converged;
  Real const eta = std::acosh(r.delta), d = std::exp(-r.log_deviation);
  std::optional<Real> gap;
  if (complete) gap = Real{2} * (r.energy - ground.reference.energy);
  report.field("calculation", "Calculation", "selected complex-root singlet excitation")
      .field("family", "Family", "one positive-deviation two-string; real I=1,...,N/2-2; string label 1")
      .field("ordering", "Ordering", "targeted branch, not an exhaustive search or a global first-excitation guarantee")
      .field("tl_through_lines", "TL through-lines", 0)
      .field("multiplicity_meaning", "Multiplicity meaning", "one physical singlet per TL eigenvector")
      .field("residual_convention", "Residual convention",
             "max phase/log-modulus equation residual divided by 2N; not an energy-error bound")
      .field("string_coordinate", "String coordinate", "u=(eta+d)/2 +/- i*a/2; d=exp(-L)>0; L is authoritative")
      .field("string_deviation", "String deviation",
             eta + d == eta ? "unresolved in rounded u; retained by L=-log(d)"
                            : "resolved in rounded u; L=-log(d) also retained")
      .field("gap_reference", "Gap reference",
             ground.reference.converged ? "E-E0; global singlet ground state" : "unavailable; ground solve failed")
      .result(complete, complete ? "converged" : "incomplete; unconverged estimate");
  std::vector<std::string> names{"states", "reference", "string"};
  if (args.roots) names.push_back("roots");
  cli::ResultOutput output(report, args.output, names);
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
void write_output(cli::RunReport& report, Arguments const& args, std::vector<model::State<Real> const*> const& states,
                  std::vector<std::optional<Real>> const& gaps, model::State<Real> const* reference = nullptr,
                  model::State<Real> const* failed = nullptr)
{
  std::vector<std::string> names{"states", "quantum_numbers"};
  if (reference) names.push_back("reference");
  if (failed) names.push_back("failed");
  if (args.roots) names.push_back("roots");
  cli::ResultOutput output(report, args.output, names);
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
               args.excitations ? (args.ferromagnetic ? "Real-root TL levels (complex levels excluded)"
                                                      : "Real-root TL levels (module minimum included)")
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

template <uni20::Real Real>
auto preamble(uni20::run_context& context, Arguments const& args, bethe::SolverOptions<Real> const& options)
{
  cli::RunReport report(context, "Spin-1 pure biquadratic chain (free ends)");
  report
      .field("hamiltonian", "Hamiltonian", args.ferromagnetic ? "H=+sum_i (S_i.S_(i+1))^2" : "H=-sum_i (S_i.S_(i+1))^2")
      .field("sites", "Sites", args.sites)
      .field("spin", "Spin", 1)
      .field("tl_loop_weight", "TL loop weight", 3)
      .field("xxz_delta", "XXZ Delta", Real{1.5})
      .field("xxz_reference", "XXZ reference", "spin-half exchange 1; +sqrt(5)/4*(sz_1-sz_N)")
      .field("precision", "Precision", args.precision)
      .field("residual_tolerance", "Residual tolerance", options.residual_tolerance);
  if (args.ferromagnetic)
    report.field("tl_convention", "TL convention", "tl_energy=+sum e_i=E-(N-1); auxiliary XXZ reference is unchanged")
        .field("defect_meaning", "Defect meaning", "M=(N-ell)/2; TL singlet defects, not physical spin flips")
        .field("spectral_gap", "Exact gap above ground space", model::ferromagnetic::spectral_gap<Real>(args.sites));
  return report;
}

template <uni20::Real Real> int run_ferro_analytic(Arguments const& args, cli::RunReport report)
{
  auto computation = report.context().computation();
  auto const ground = model::ferromagnetic::ground_space<Real>(args.sites);
  std::vector<model::ferromagnetic::OneDefectLevel<Real>> levels;
  if (args.one_defect)
  {
    if (args.sites - 1 > args.max_candidates.value_or(10000))
      throw std::length_error("one-defect band exceeds max_candidates; raise --max-candidates");
    levels.reserve(args.sites - 1);
    for (std::size_t j = args.sites - 1; j > 0; --j)
      levels.push_back(model::ferromagnetic::one_defect_level<Real>(args.sites, j));
  }
  else if (args.through_lines && *args.through_lines == args.sites - 2)
    levels.push_back(model::ferromagnetic::one_defect_level<Real>(args.sites, args.sites - 1));
  computation.finish();
  report
      .field("calculation", "Calculation",
             args.one_defect  ? "complete one-defect TL module"
             : levels.empty() ? "exact ferro ground space"
                              : "exact one-defect module minimum")
      .field("coverage", "Coverage",
             args.one_defect  ? "all N-1 eigenvalues of ell=N-2; NOT the full excited spectrum"
             : levels.empty() ? "ground space only; zero modes are not individually enumerated"
                              : "first positive level above the entire ground space")
      .field("wave_number", "Wave number", "k=pi*j/N is an OBC standing-wave coordinate, not lattice momentum")
      .field("ground_energy", "Ground energy", ground.energy)
      .field("ground_multiplicity", "Ground-space dimension", ground.multiplicity, {.missing = "overflow (>uint64)"})
      .field("gap_reference", "Gap reference", "E-(N-1); exact degenerate ferro ground space")
      .field("multiplicity_meaning", "Multiplicity meaning",
             "physical states per TL eigenvector; not SU(2) multiplets or accidental-degeneracy sums")
      .result(true, "exact spectral rules");
  cli::ResultOutput output(report, args.output, {"states", "reference"});
  output.table(
      "states", "Analytic ferromagnetic levels",
      [&](auto& table) {
        if (levels.empty())
          table.append(std::size_t{0}, ground.through_lines, std::size_t{0}, ground.multiplicity,
                       std::optional<std::size_t>{}, std::optional<Real>{}, ground.energy, Real{0});
        for (std::size_t i = 0; i < levels.size(); ++i)
        {
          auto const& level = levels[i];
          table.append(i, level.through_lines, std::size_t{1}, level.multiplicity,
                       std::optional<std::size_t>{level.mode}, std::optional<Real>{level.wave_number}, level.energy,
                       level.gap);
        }
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("through_lines"),
      cli::column<std::size_t>("defects"), cli::column<std::optional<std::uint64_t>>("multiplicity"),
      cli::column<std::optional<std::size_t>>("mode", "j"),
      cli::column<std::optional<Real>>("wave_number", "Standing-wave k"), cli::column<Real>("energy", "Energy"),
      cli::column<Real>("gap", "E-E0"));
  output.table(
      "reference", "Exact ferromagnetic ground space",
      [&](auto& table) {
        table.append(std::max(std::size_t{1}, levels.size()), ground.through_lines, ground.energy, ground.multiplicity);
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("through_lines"),
      cli::column<Real>("energy", "Energy"), cli::column<std::optional<std::uint64_t>>("multiplicity"));
  output.finish();
  return 0;
}

int finish(bool converged)
{
  if (!converged) std::cerr << "Biquadratic solve incomplete; consider a larger budget or higher precision.\n";
  return converged ? 0 : 2;
}

template <uni20::Real Real> int run(Arguments const& args, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  bethe::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");
  auto report = preamble(context, args, options);
  if (args.q_seed || args.q_spectrum) return run_qsystem(args, options, std::move(report));
  if (args.singlet_excitation) return run_singlet(args, options, std::move(report));
  if (args.ferromagnetic && !args.excitations && !args.numbers)
    return run_ferro_analytic<Real>(args, std::move(report));
  auto computation = context.computation();
  if (args.excitations)
  {
    auto const ell = args.through_lines.value_or(args.ferromagnetic ? args.sites - 2 : 2);
    bethe::RealExcitationOptions const enumeration{.count = *args.excitations,
                                                   .max_candidates = args.max_candidates.value_or(10000)};
    auto const scan = args.ferromagnetic
                          ? model::ferromagnetic::real_excitations<Real>(args.sites, ell, enumeration, options)
                          : model::real_excitations<Real>(args.sites, ell, enumeration, options);
    computation.finish();
    auto const& ground = scan.ground_state;
    report.field("calculation", "Calculation", "restricted real-root excitations")
        .field("family", "Family", "positive finite roots; complex-root levels excluded; NOT the complete TL spectrum")
        .field("tl_through_lines", "TL through-lines", ell)
        .field("tl_defects", "TL defects M", (args.sites - ell) / 2)
        .field("multiplicity_per_tl_eigenvector", "Multiplicity per TL eigenvector",
               bethe::temperley_lieb::spin_chain_multiplicity(3, ell), {.missing = "overflow (>uint64)"})
        .field("multiplicity_meaning", "Multiplicity meaning",
               "physical states, not SU(2) multiplets; not an accidental-degeneracy sum")
        .field("candidates", "Candidates", scan.candidate_count)
        .field("converged_candidates", "Converged candidates", scan.converged_count)
        .field("returned_levels", "Returned levels", scan.levels.size())
        .field("ordering", "Ordering",
               scan.family_converged() ? "complete within supported family" : "incomplete; failed candidates excluded")
        .field("ground_energy", "Ground energy", ground.energy)
        .field("ground_status", "Ground status", status(ground.reference.status))
        .field("ground_residual", "Ground residual", ground.reference.residual_norm)
        .field("ground_iterations", "Ground iterations", ground.reference.iterations)
        .field("gap_reference", "Gap reference",
               args.ferromagnetic           ? "E-(N-1); exact degenerate ferro ground space"
               : ground.reference.converged ? "E-E0; global singlet ground state"
                                            : "unavailable; ground solve failed")
        .result(scan.converged(), scan.converged() ? "converged" : "incomplete scan or ground reference");
    if (scan.first_unconverged)
      report
          .field("first_failed_i", "First failed I",
                 cli::quantum_number_text(scan.first_unconverged->reference.quantum_numbers))
          .field("first_failed_status", "First failed status", status(scan.first_unconverged->reference.status))
          .field("first_failed_residual", "First failed residual", scan.first_unconverged->reference.residual_norm);
    std::vector<model::State<Real> const*> states;
    std::vector<std::optional<Real>> gaps;
    for (auto const& level : scan.levels)
    {
      states.push_back(&level.state);
      gaps.push_back(level.gap);
    }
    write_output(report, args, states, gaps, &ground, scan.first_unconverged ? &*scan.first_unconverged : nullptr);
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
    computation.finish();
    report.field("calculation", "Calculation", "TL module minima (not physical-spin sectors)")
        .field("multiplicity_meaning", "Multiplicity meaning",
               "physical states per TL eigenvector, not SU(2) multiplets")
        .result(converged, converged ? "converged" : "incomplete; unconverged estimates");
    std::vector<model::State<Real> const*> rows;
    std::vector<std::optional<Real>> gaps;
    for (auto const& state : states)
    {
      rows.push_back(&state);
      gaps.push_back(state.reference.converged && states.front().reference.converged
                         ? std::optional<Real>{Real{2} * (state.reference.energy - states.front().reference.energy)}
                         : std::nullopt);
    }
    write_output(report, args, rows, gaps);
    return finish(converged);
  }
  auto const state =
      args.numbers ? (args.ferromagnetic ? model::ferromagnetic::solve_real<Real>(args.sites, *args.numbers, options)
                                         : model::solve_real<Real>(args.sites, *args.numbers, options))
                   : model::sector_ground_state<Real>(args.sites, args.through_lines.value_or(0), options);
  computation.finish();
  auto const& reference = state.reference;
  report
      .status(reference.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning,
              status(reference.status))
      .field("calculation", "Calculation",
             args.numbers               ? "specified real-root TL level"
             : state.through_lines == 0 ? "even-chain singlet ground state"
                                        : "TL module minimum")
      .field("tl_through_lines", "TL through-lines", state.through_lines)
      .field("tl_defects", "TL defects M", (args.sites - state.through_lines) / 2)
      .field("multiplicity",
             !args.ferromagnetic && state.through_lines == 0 ? "Ground-state multiplicity"
                                                             : "Multiplicity per TL eigenvector",
             state.multiplicity, {.missing = "overflow (>uint64)"})
      .result(reference.converged, status(reference.status))
      .field("total_energy", "Total energy", state.energy)
      .field("energy_per_site", "Energy per site", state.energy / Real(args.sites))
      .field("tl_energy_sum_e_i", args.ferromagnetic ? "TL energy (+sum e_i)" : "TL energy (-sum e_i)", state.tl_energy)
      .field("xxz_reference_energy", "XXZ reference energy", reference.energy)
      .field("residual_norm", "Residual norm", reference.residual_norm)
      .field("iterations", "Iterations", reference.iterations);
  if (state.through_lines == 0) report.field("total_spin", "Total spin", 0);
  if (state.through_lines != 0)
    report.field("multiplicity_meaning", "Multiplicity meaning", "physical states, not a physical-spin label");
  if (args.ferromagnetic)
    report.field("gap_reference", "Gap reference", "E-(N-1); exact degenerate ferro ground space");
  write_output<Real>(report, args, {&state},
                     {args.ferromagnetic && reference.converged ? std::optional<Real>{state.tl_energy} : std::nullopt});
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
