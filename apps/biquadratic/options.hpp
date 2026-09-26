// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "data-output-options.hpp"
#include "program-options.hpp"
#include <bethe/xxz_common.hpp>

namespace bethe::apps::biquadratic
{
namespace cli = bethe::cli;
struct Arguments
{
    std::size_t sites = 0, max_iterations = 10000;
    std::optional<std::size_t> through_lines, excitations, max_candidates, real_window;
    std::optional<std::size_t> bound_pairs, bound_triples, bound_quartets;
    std::optional<std::size_t> pair_defects, mixed_window, real_defects;
    std::optional<std::string> pair_defect;
    std::optional<std::string> two_pairs;
    std::optional<std::size_t> two_pair_states, pair_window;
    std::optional<std::string> triple_defect;
    std::optional<std::size_t> triple_defects;
    std::optional<bethe::xxz::QuantumNumbers> numbers;
    std::optional<std::string> tolerance;
    std::optional<std::string> q_seed;
    std::optional<std::size_t> max_attempts;
    std::string precision = "fp64";
    cli::DataOutputOptions output;
    bool roots = false, sectors = false, spin_content = false;
    bool q_spectrum = false, singlet_excitation = false;
    bool ferromagnetic = false, one_defect = false;
};
inline auto program_info()
{
  auto info = bethe::cli::program_info("bethe-biquadratic-obc",
                                       "Spin-1 pure biquadratic chain, free ends: H=-sum_i (S_i.S_(i+1))^2 by default.",
                                       bethe::citations::Tool::biquadratic_obc);
  info.notes = {
      "Default: unique singlet ground state; even N>=2, coefficient -1.",
      "--ferromagnetic reverses the sign: H=+sum (S.S)^2; gaps use the exact degenerate ground space.",
      "--ferromagnetic --one-defect gives the complete ell=N-2 band analytically, including odd N.",
      "--ferromagnetic --bound-pairs COUNT|all targets two-string modes in ell=N-4 on odd/even chains.",
      "--ferromagnetic --bound-triples COUNT|all targets three-string droplets in ell=N-6 on odd/even chains.",
      "--pair-defect I1,...,Ir,J selects one pair plus real roots; --pair-defects COUNT|all scans that family.",
      "Ferromagnetic real-root scans exclude complex-root levels; NOT general module minima.",
      "--ferromagnetic --excitations COUNT|all --real-window WIDTH scans only the highest WIDTH real labels.",
      "TL loop weight 3; reference XXZ Delta=3/2 with opposite end fields.",
      "Real-root scans are NOT complete spectra: complex-root levels are excluded.",
      "--q-spectrum searches real and complex levels of one TL module; validated through N=8.",
      "Q-system modes have no site cutoff; larger sizes are experimental and may be costly or unresolved.",
      "--q-seed selects a polynomial branch, not necessarily a low-lying state.",
      "--singlet-excitation targets one two-string singlet above a real sea, including on long chains.",
      "Multiplicity counts physical states per TL eigenvector, not SU(2) multiplets.",
      "TL through-lines are not physical spin; AF ground/sector helpers and Q-system modes require even N.",
      "Selected real roots and ferro real-root scans accept odd/even N; no lattice momentum.",
      "This is not the TB point, ULS point, or zero-boundary-field XXZ chain.",
      "See docs/biquadratic.md for the TL mapping and representation multiplicities.",
      "Use --references for literature and applicability; see CITATIONS.md."};
  return info;
}
inline void add_options(CLI::App& app, Arguments& args)
{
  bethe::cli::option(app, "N", args.sites, "Number of particles or sites")->required();
  bethe::cli::option(app, "--ferromagnetic", args.ferromagnetic,
                     "use H=+sum (S.S)^2 and the exact ferro ground reference");
  bethe::cli::option(app, "--one-defect", args.one_defect, "complete analytic ferro ell=N-2 band; odd/even N>=2");
  bethe::cli::all_count_option(
      app, "--bound-pairs", args.bound_pairs,
      "ferro: first COUNT, or all N-3, targeted two-string modes; NOT full two-defect spectrum");
  bethe::cli::all_count_option(
      app, "--bound-triples", args.bound_triples,
      "ferro: first COUNT, or all N-5, targeted three-string modes; NOT full three-defect spectrum");
  bethe::cli::all_count_option(
      app, "--bound-quartets", args.bound_quartets,
      "ferro: first COUNT, or all N-7, targeted four-string modes; NOT full four-defect spectrum");
  bethe::cli::option(app, "--pair-defect", args.pair_defect,
                     "ferro: selected pair plus real roots, integer labels I1,...,Ir,J; r>=1, ell=N-2(r+2)");
  bethe::cli::option(app, "--triple-defect", args.triple_defect,
                     "ferro: selected bound triple plus one real root, integer labels I,J; ell=N-8");
  bethe::cli::all_count_option(app, "--triple-defects", args.triple_defects,
                               "ferro: retain lowest COUNT or all converged triple-plus-real-root candidates");
  bethe::cli::option(app, "--two-pairs", args.two_pairs,
                     "ferro: two scattering pairs with ordered labels J1,J2; ell=N-8");
  bethe::cli::all_count_option(app, "--two-pair-states", args.two_pair_states,
                               "ferro: retain lowest COUNT or all converged two-pair candidates");
  bethe::cli::option(app, "--pair-window", args.pair_window,
                     "--two-pair-states: highest WIDTH pair labels; choose(WIDTH,2) candidates");
  bethe::cli::all_count_option(app, "--pair-defects", args.pair_defects,
                               "ferro: retain lowest COUNT or all converged one-pair-family candidates");
  bethe::cli::option(app, "--real-defects", args.real_defects,
                     "--pair-defects: number R>=1 of real roots alongside one pair (default: 1)");
  bethe::cli::option(
      app, "--mixed-window", args.mixed_window,
      "--pair-defects/--triple-defects: highest WIDTH real and string labels (triple: WIDTH^2 candidates)");
  bethe::cli::option(app, "--through-lines", args.through_lines,
                     "select TL module; standalone ferro minima require ELL=N or N-2");
  bethe::cli::option(app, "--sectors", args.sectors, "AF: lowest level in every TL module");
  bethe::cli::all_count_option(
      app, "--excitations", args.excitations,
      "lowest COUNT, or all, supported real-root levels; default ELL=2 (ferro: N-2); excludes complex roots");
  bethe::cli::option(app, "--max-candidates", args.max_candidates,
                     "real-family scan / analytic band / bound-cluster row limit (default: 10000)");
  bethe::cli::option(
      app, "--real-window", args.real_window,
      "ferro --excitations: highest WIDTH integer labels only; choose(WIDTH,M) candidates, not a full sector");
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
  bethe::cli::option(app, "--spin-content", args.spin_content,
                     "add physical SU(2) multiplets per TL eigenvector; null if total dimension exceeds uint64");
  bethe::cli::option(
      app, "--tolerance", args.tolerance,
      "equation tolerance (default: 32 epsilon); Q-system uses coefficient backward error, not phase error");
  bethe::cli::option(app, "--max-iterations", args.max_iterations, "accepted Newton updates (default: 10000)")
      ->capture_default_str();
  bethe::cli::precision_option(app, args.precision);
  cli::add_data_output_options(app, args.output, true);
}

inline void validate(Arguments const& args)
{
  args.output.validate();
  if (args.triple_defect || args.triple_defects)
  {
    if (!args.ferromagnetic || args.sites < 8)
      throw std::invalid_argument("triple plus defect requires --ferromagnetic and N>=8");
    if (args.one_defect || args.q_seed || args.q_spectrum || args.sectors || args.excitations || args.numbers ||
        args.through_lines || args.singlet_excitation || args.max_attempts || args.bound_pairs || args.bound_triples ||
        args.pair_defect || args.pair_defects || args.two_pairs || args.two_pair_states || args.real_defects ||
        (args.triple_defect && args.triple_defects))
      throw std::invalid_argument(
          "triple-defect selections cannot combine with other state selections or --real-defects");
    if (args.triple_defects && *args.triple_defects == 0)
      throw std::invalid_argument("--triple-defects requires a positive count or all");
    if (args.triple_defect && args.max_candidates)
      throw std::invalid_argument("--max-candidates is not used by --triple-defect I,J");
    if (args.mixed_window && (*args.mixed_window == 0 || *args.mixed_window > args.sites - 7))
      throw std::invalid_argument("triple --mixed-window requires 1<=WIDTH<=N-7");
  }
  if (args.pair_window && !args.two_pair_states)
    throw std::invalid_argument("--pair-window requires --two-pair-states COUNT|all");
  if (args.two_pairs || args.two_pair_states)
  {
    if (!args.ferromagnetic || args.sites < 8)
      throw std::invalid_argument("two-pair states require --ferromagnetic and N>=8");
    if (args.one_defect || args.q_seed || args.q_spectrum || args.sectors || args.excitations || args.numbers ||
        args.through_lines || args.singlet_excitation || args.max_attempts || args.bound_pairs || args.bound_triples ||
        args.pair_defect || args.pair_defects || (args.two_pairs && args.two_pair_states))
      throw std::invalid_argument("two-pair selections cannot combine with other state selections");
    if (args.two_pair_states && *args.two_pair_states == 0)
      throw std::invalid_argument("--two-pair-states requires a positive count or all");
    if (args.two_pairs && args.max_candidates)
      throw std::invalid_argument("--max-candidates is not used by --two-pairs J1,J2");
    if (args.pair_window && (*args.pair_window < 2 || *args.pair_window > args.sites - 6))
      throw std::invalid_argument("--pair-window requires 2<=WIDTH<=N-6");
  }
  if (args.mixed_window && !args.pair_defects && !args.triple_defects)
    throw std::invalid_argument("--mixed-window requires --pair-defects or --triple-defects COUNT|all");
  if (args.real_defects && !args.pair_defects)
    throw std::invalid_argument(
        "--real-defects requires --pair-defects COUNT|all; selected states infer R from labels");
  if (args.pair_defect || args.pair_defects)
  {
    if (!args.ferromagnetic || args.sites < 6)
      throw std::invalid_argument("pair plus defect requires --ferromagnetic and N>=6");
    if (args.one_defect || args.q_seed || args.q_spectrum || args.sectors || args.excitations || args.numbers ||
        args.through_lines || args.singlet_excitation || args.max_attempts || args.bound_pairs || args.bound_triples ||
        (args.pair_defect && args.pair_defects))
      throw std::invalid_argument(
          "pair-defect selections determine the TL module; cannot combine with other state selections");
    if (args.pair_defects && *args.pair_defects == 0)
      throw std::invalid_argument("--pair-defects requires a positive count or all");
    if (args.pair_defect && args.max_candidates)
      throw std::invalid_argument("--max-candidates is not used by --pair-defect I,J");
    auto const sea = args.real_defects.value_or(1);
    if (sea == 0 || sea > args.sites / 2 - 2) throw std::invalid_argument("--real-defects requires 1<=R<=N/2-2");
    if (args.mixed_window && (*args.mixed_window < sea || *args.mixed_window > args.sites - 2 * (sea + 2) + 1))
      throw std::invalid_argument("--mixed-window requires R<=WIDTH<=N-2(R+2)+1");
  }
  if (args.real_window && (!args.ferromagnetic || !args.excitations))
    throw std::invalid_argument("--real-window requires --ferromagnetic --excitations COUNT|all");
  if (args.bound_pairs || args.bound_triples || args.bound_quartets)
  {
    std::string const option = args.bound_pairs     ? "--bound-pairs"
                               : args.bound_triples ? "--bound-triples"
                                                    : "--bound-quartets";
    std::size_t const defects = args.bound_pairs ? 2 : args.bound_triples ? 3 : 4;
    if (!args.ferromagnetic) throw std::invalid_argument(option + " requires --ferromagnetic");
    if (args.sites < 2 * defects) throw std::invalid_argument(option + " requires N>=" + std::to_string(2 * defects));
    if ((args.bound_pairs ? *args.bound_pairs : args.bound_triples ? *args.bound_triples : *args.bound_quartets) == 0)
      throw std::invalid_argument(option + " requires a positive count or all");
    if (args.one_defect || args.q_seed || args.q_spectrum || args.sectors || args.excitations || args.numbers ||
        args.through_lines || args.singlet_excitation || args.max_attempts || args.pair_defect || args.pair_defects ||
        args.two_pairs || args.two_pair_states || args.triple_defect || args.triple_defects ||
        (int(args.bound_pairs.has_value()) + int(args.bound_triples.has_value()) +
             int(args.bound_quartets.has_value()) >
         1))
      throw std::invalid_argument(option + " determines ell=N-" + std::to_string(2 * defects) +
                                  "; cannot combine with other state selections");
  }
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
    bool const analytic = !args.numbers && !args.excitations && !args.q_seed && !args.q_spectrum && !args.bound_pairs &&
                          !args.bound_triples && !args.bound_quartets && !args.pair_defect && !args.pair_defects &&
                          !args.two_pairs && !args.two_pair_states && !args.triple_defect && !args.triple_defects;
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
  if (args.max_candidates && !args.excitations && !args.one_defect && !args.bound_pairs && !args.bound_triples &&
      !args.bound_quartets && !args.pair_defects && !args.two_pair_states && !args.triple_defects)
    throw std::invalid_argument("--max-candidates requires a state scan");
}
} // namespace bethe::apps::biquadratic
