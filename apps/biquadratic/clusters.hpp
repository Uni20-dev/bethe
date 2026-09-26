// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "report.hpp"

namespace bethe::apps::biquadratic
{
enum class ClusterFamily
{
  pair,
  triple,
  quartet,
  pair_real,
  two_pairs,
  triple_real
};

template <uni20::Real Real, ClusterFamily Family>
int run_bound_clusters(Arguments const& args, std::size_t requested, bethe::SolverOptions<Real> const& options,
                       cli::RunReport report)
{
  constexpr bool TwoPairs = Family == ClusterFamily::two_pairs;
  constexpr bool Quartet = Family == ClusterFamily::quartet;
  constexpr bool TripleMixed = Family == ClusterFamily::triple_real;
  constexpr bool PairMixed = Family == ClusterFamily::pair_real;
  constexpr bool Mixed = PairMixed || TripleMixed;
  constexpr std::size_t Defects = TwoPairs || TripleMixed || Quartet ? 4 : Family == ClusterFamily::pair ? 2 : 3;
  constexpr bool Scan = Mixed || TwoPairs;
  auto computation = report.context().computation();
  auto const& selected = TripleMixed ? args.triple_defect : args.pair_defect;
  std::size_t available = args.sites - (2 * Defects - 1), count = std::min(requested, available);
  std::size_t real_width = 0, pair_width = 0, selected_j = 0;
  std::size_t sea = args.real_defects.value_or(1), defects = Defects;
  std::size_t real_last = 0, pair_last = 0;
  std::vector<std::size_t> real_labels;
  std::array<std::size_t, 2> pair_labels{};
  if constexpr (TwoPairs)
  {
    (void)bethe::xxz::detail::checked_sites(args.sites);
    pair_last = args.sites - 6;
    pair_width = args.pair_window.value_or(pair_last);
    if (args.two_pairs)
    {
      std::string_view text = *args.two_pairs;
      auto const comma = text.find(',');
      if (comma == std::string_view::npos || text.find(',', comma + 1) != std::string_view::npos)
        throw std::invalid_argument("--two-pairs requires J1,J2");
      pair_labels = {cli::parse_size(text.substr(0, comma)), cli::parse_size(text.substr(comma + 1))};
      count = available = 1;
    }
    else
    {
      count = available = bethe::detail::bounded_binomial(pair_width, 2, args.max_candidates.value_or(10000));
      pair_labels = {pair_last - pair_width + 1, pair_last - pair_width + 2};
    }
  }
  if constexpr (Mixed)
  {
    (void)bethe::xxz::detail::checked_sites(args.sites);
    if (selected)
    {
      std::string_view text = *selected;
      sea = std::count(text.begin(), text.end(), ',');
      if (TripleMixed && sea != 1) throw std::invalid_argument("--triple-defect requires exactly I,J");
      if (sea == 0 || sea > args.sites / 2 - (TripleMixed ? 3 : 2))
        throw std::invalid_argument("mixed selection requires real-root labels followed by one string label");
      defects = sea + (TripleMixed ? 3 : 2);
      bethe::xxz::quantum_group::detail::check_matrix_size<Real>(defects);
      real_labels.reserve(sea);
      for (std::size_t i = 0; i < sea; ++i)
      {
        auto const comma = text.find(',');
        real_labels.push_back(cli::parse_size(text.substr(0, comma)));
        text.remove_prefix(comma + 1);
      }
      selected_j = cli::parse_size(text);
      count = available = 1;
    }
    else
    {
      defects = sea + (TripleMixed ? 3 : 2);
      bethe::xxz::quantum_group::detail::check_matrix_size<Real>(defects);
      real_last = args.sites - (TripleMixed ? 3 : defects);
      pair_last = args.sites - 2 * defects + 1;
      real_width = args.mixed_window.value_or(real_last);
      pair_width = args.mixed_window.value_or(pair_last);
      auto const budget = args.max_candidates.value_or(10000);
      auto const combinations = bethe::detail::bounded_binomial(real_width, sea, budget / pair_width);
      count = available = combinations * pair_width;
      real_labels.resize(sea);
      std::iota(real_labels.begin(), real_labels.end(), real_last - real_width + 1);
    }
  }
  if (count > args.max_candidates.value_or(10000))
    throw std::length_error("bound-cluster rows exceed max_candidates; raise --max-candidates");
  auto const ground = model::ferromagnetic::ground_space<Real>(args.sites);
  auto solve = [&](std::size_t mode) {
    if constexpr (TwoPairs)
      return model::ferromagnetic::two_bound_pairs<Real>(args.sites, pair_labels, options);
    else if constexpr (Quartet)
      return model::ferromagnetic::bound_quartet<Real>(args.sites, mode, options);
    else if constexpr (TripleMixed)
      return model::ferromagnetic::triple_defect<Real>(args.sites, real_labels.front(), selected_j, options);
    else if constexpr (Mixed)
      return model::ferromagnetic::pair_with_real_roots<Real>(args.sites, real_labels, selected_j, options);
    else if constexpr (Defects == 2)
      return model::ferromagnetic::bound_pair<Real>(args.sites, mode, options);
    else
      return model::ferromagnetic::bound_triple<Real>(args.sites, mode, options);
  };
  std::vector<decltype(solve(1))> states;
  states.reserve(count);
  bool complete = true;
  std::size_t converged_count = 0;
  auto append = [&](std::size_t mode) {
    states.push_back(solve(mode));
    complete = complete && states.back().reference.converged;
    converged_count += states.back().reference.converged;
  };
  if constexpr (TwoPairs)
  {
    do
    {
      append(1);
    }
    while (!args.two_pairs && bethe::detail::advance_combination<std::size_t>(pair_labels, pair_last));
  }
  else if constexpr (Mixed)
  {
    if (selected)
      append(1);
    else
    {
      do
      {
        for (std::size_t j = 0; j < pair_width; ++j)
        {
          selected_j = pair_last - j;
          append(1);
        }
      }
      while (bethe::detail::advance_combination<std::size_t>(real_labels, real_last));
    }
  }
  else
    for (std::size_t mode = 1; mode <= count; ++mode)
      append(mode);
  if constexpr (Scan)
  {
    std::stable_sort(states.begin(), states.end(), [](auto const& a, auto const& b) {
      if (a.reference.converged != b.reference.converged) return a.reference.converged;
      if (!a.reference.converged) return false;
      // Rank direct gaps, not extensive energies that can tie by rounding.
      if (a.tl_energy != b.tl_energy) return a.tl_energy < b.tl_energy;
      if constexpr (TwoPairs)
        return a.reference.string_labels < b.reference.string_labels;
      else if constexpr (TripleMixed)
        return std::pair{a.reference.real_label, a.reference.string_label} <
               std::pair{b.reference.real_label, b.reference.string_label};
      else
      {
        if (a.reference.real_labels != b.reference.real_labels)
          return a.reference.real_labels < b.reference.real_labels;
        return a.reference.string_label < b.reference.string_label;
      }
    });
    auto const retained = std::min(requested, converged_count);
    if (!complete)
    {
      auto failure = states[converged_count];
      states.resize(retained);
      states.push_back(std::move(failure)); // One explicitly unconverged diagnostic row.
    }
    else
      states.resize(retained);
    report.field("candidate_count", "Scanned candidates", count)
        .field("converged_candidates", "Converged candidates", converged_count)
        .field("failed_candidates", "Failed candidates", count - converged_count)
        .field("retained_levels", "Retained converged levels", retained)
        .field("failure_rows", "Failure rows", "at most one unconverged diagnostic follows the converged levels");
    if constexpr (TwoPairs)
      report.field("pair_label_first", "First allowed J", args.two_pairs ? pair_labels[0] : pair_last - pair_width + 1)
          .field("pair_label_last", "Last allowed J", args.two_pairs ? pair_labels[1] : pair_last);
    else
      report.field("real_defects", TripleMixed ? "Real roots alongside triple" : "Real roots alongside pair", sea)
          .field("real_label_first", "First allowed I", selected ? real_labels.front() : real_last - real_width + 1)
          .field("real_label_last", "Last allowed I", selected ? real_labels.back() : real_last)
          .field(TripleMixed ? "triple_label_first" : "pair_label_first", "First allowed J",
                 selected ? selected_j : pair_last - pair_width + 1)
          .field(TripleMixed ? "triple_label_last" : "pair_label_last", "Last allowed J",
                 selected ? selected_j : pair_last);
  }
  computation.finish();
  report
      .field("calculation", "Calculation",
             TwoPairs       ? "two scattering bound pairs"
             : Quartet      ? "targeted four-defect bound-quartet modes"
             : TripleMixed  ? "mixed bound triple plus one real defect"
             : Mixed        ? "mixed bound pair plus selected real defects"
             : Defects == 2 ? "targeted two-defect bound-pair modes"
                            : "targeted three-defect bound-triple modes")
      .field("family", "Family",
             TwoPairs ? "two two-strings, no real roots; ordered J1,J2; sign(d_i)=(-1)^(N-J_i-i), i=1,2"
             : Quartet
                 ? "one four-string, no real sea; J=N-6-mode, sign(d)=(-1)^(mode+1), complex z=exp(-L_outer+i*phi)"
             : TripleMixed  ? "one three-string plus one real root; I,J labels, complex deviation z=exp(-L+i*phi)"
             : Mixed        ? "one two-string plus real roots; I1,...,Ir,J labels, sign(d)=(-1)^(N-J-1)"
             : Defects == 2 ? "one two-string, no real sea; J=N-2-mode, sign(d)=(-1)^(mode+1)"
                            : "one three-string, no real sea; J=N-4-mode, complex deviation z=exp(-L+i*phi)")
      .field("coverage", "Coverage",
             TwoPairs      ? "selected two-pair family only; NOT the full four-defect or excited spectrum"
             : Quartet     ? "selected four-string family only; NOT the full four-defect or excited spectrum"
             : TripleMixed ? "selected triple-plus-real-root family only; NOT the full four-defect or excited spectrum"
             : Mixed
                 ? (sea == 1
                        ? "selected mixed-family label rectangle only; NOT the full three-defect or excited spectrum"
                        : "selected one-pair plus real-root family only; NOT the full module or excited spectrum")
             : Defects == 2 ? "selected two-string family only; NOT the full two-defect or excited spectrum"
                            : "selected three-string family only; NOT the full three-defect or excited spectrum")
      .field("ordering", "Ordering",
             Scan ? "direct gap order within scanned candidates, ties by Bethe labels; not a global excitation rank"
                  : "mode order from the low-energy branch edge; not a global excitation rank")
      .field("tl_through_lines", "TL through-lines", args.sites - 2 * defects)
      .field("tl_defects", "TL singlet insertions M", defects)
      .field("available_modes",
             Scan           ? "Selected label combinations"
             : Quartet      ? "Available four-string labels"
             : Defects == 2 ? "Available two-string labels"
                            : "Available three-string labels",
             available)
      .field("selected_modes", Scan ? "Returned rows" : "Selected modes", states.size())
      .field("residual_convention", "Residual convention",
             "max phase/log-modulus residual divided by 2N; not an energy-error bound")
      .field("string_coordinate", "String coordinate",
             Quartet ? "u=(eta+d+i*a)/2, w=u+eta+z, and conjugates; both logarithms, sign and phi are authoritative"
             : (Defects == 2 || PairMixed || TwoPairs)
                 ? "pair u=(eta+d)/2 +/- i*a/2; d=sign*exp(-L); real roots u=i*alpha/2; L is authoritative"
                 : "u0=i*a/2; u+/-=eta+Re(z) +/- i*(a/2+Im(z)); L and phi are authoritative")
      .field("wave_number", "Wave number", "mode and string center are branch coordinates, not lattice momentum")
      .field("multiplicity_meaning", "Multiplicity meaning", "physical states per TL eigenvector, not SU(2) multiplets")
      .field("gap_reference", "Gap reference", "E-(N-1); exact degenerate ferro ground space")
      .field(TwoPairs       ? "bulk_two_pair_threshold"
             : Quartet      ? "bulk_quartet_threshold"
             : TripleMixed  ? "bulk_triple_defect_threshold"
             : Mixed        ? "bulk_mixed_threshold"
             : Defects == 2 ? "bulk_pair_threshold"
                            : "bulk_triple_threshold",
             TwoPairs       ? "Bulk two-pair threshold"
             : Quartet      ? "Bulk bound-quartet threshold"
             : TripleMixed  ? "Bulk triple-plus-defect threshold"
             : Mixed        ? "Bulk pair-plus-defect threshold"
             : Defects == 2 ? "Bulk bound-pair threshold"
                            : "Bulk bound-triple threshold",
             TwoPairs       ? Real{10} / Real{3}
             : Quartet      ? Real{15} / Real{7}
             : TripleMixed  ? Real{3}
             : Mixed        ? Real{5} / Real{3} + Real(sea)
             : Defects == 2 ? Real{5} / Real{3}
                            : Real{2})
      .result(complete, complete ? "converged targeted branches" : "incomplete; unconverged estimates");
  std::vector<std::string> names{"states", "reference", "string"};
  if constexpr (Mixed) names.push_back("labels");
  if (args.roots) names.push_back("roots");
  cli::ResultOutput output(report, args.output, table_names(args, names));
  output.table(
      "states",
      TwoPairs       ? "Two scattering bound pairs"
      : Quartet      ? "Four-defect bound-quartet modes"
      : TripleMixed  ? "Mixed triple-plus-defect states"
      : Mixed        ? "Mixed pair-plus-defect states"
      : Defects == 2 ? "Two-defect bound-pair modes"
                     : "Three-defect bound-triple modes",
      [&](auto& table) {
        for (std::size_t i = 0; i < states.size(); ++i)
        {
          auto const& s = states[i];
          auto const& r = s.reference;
          table.append(i, Scan ? i + 1 : s.mode, s.through_lines, s.multiplicity, s.energy,
                       r.converged ? std::optional<Real>{s.tl_energy} : std::nullopt, s.tl_energy, r.energy,
                       r.residual_norm, r.phase_residual, r.modulus_residual, r.iterations, r.converged,
                       std::string(status(r.status)));
        }
      },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>(Scan ? "row" : "mode"),
      cli::column<std::size_t>("through_lines"), cli::column<std::optional<std::uint64_t>>("multiplicity"),
      cli::column<Real>("energy", "Energy"), cli::column<std::optional<Real>>("gap", "E-E0"),
      cli::column<Real>("tl_energy"), cli::column<Real>("reference_energy"), cli::column<Real>("residual"),
      cli::column<Real>("phase_residual"), cli::column<Real>("modulus_residual"),
      cli::column<std::size_t>("iterations"), cli::column<bool>("converged"), cli::column<std::string>("status"));
  output.table(
      "reference", "Exact ferromagnetic ground space",
      [&](auto& table) { table.append(states.size(), ground.through_lines, ground.energy, ground.multiplicity); },
      cli::column<std::size_t>("state_id"), cli::column<std::size_t>("through_lines"),
      cli::column<Real>("energy", "Energy"), cli::column<std::optional<std::uint64_t>>("multiplicity"));
  if constexpr (Mixed)
    output.table(
        "labels", "Mixed-family Bethe labels (not momenta or energy ranks)",
        [&](auto& table) {
          for (std::size_t i = 0; i < states.size(); ++i)
          {
            auto const& r = states[i].reference;
            if constexpr (TripleMixed)
              table.append(i, r.real_label, r.string_label, r.rapidity);
            else
              for (std::size_t j = 0; j < r.real_labels.size(); ++j)
                table.append(i, r.real_labels[j], r.string_label, r.rapidities[j]);
          }
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("real_label", "I"),
        cli::column<std::size_t>("string_label", "J"), cli::column<Real>("real_rapidity", "alpha"));
  if constexpr (Quartet)
    output.table(
        "string", "Four-string parameters (both logarithms, sign and phi are authoritative)",
        [&](auto& table) {
          for (std::size_t i = 0; i < states.size(); ++i)
          {
            auto const& r = states[i].reference;
            Real const d = Real(r.inner_deviation_sign) * std::exp(-r.inner_log_deviation);
            Real const z = std::exp(-r.outer_log_deviation);
            table.append(i, r.string_label, r.center, r.inner_deviation_sign, r.inner_log_deviation,
                         d != Real{0} ? std::optional<Real>{d} : std::nullopt, r.outer_log_deviation,
                         r.outer_deviation_phase,
                         z != Real{0} ? std::optional<Real>{z * std::cos(r.outer_deviation_phase)} : std::nullopt,
                         z != Real{0} ? std::optional<Real>{z * std::sin(r.outer_deviation_phase)} : std::nullopt);
          }
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("string_label", "J"),
        cli::column<Real>("center", "a"), cli::column<int>("inner_deviation_sign", "sign(d)"),
        cli::column<Real>("inner_log_deviation", "L_inner"),
        cli::column<std::optional<Real>>("inner_deviation", "d (null if underflow)"),
        cli::column<Real>("outer_log_deviation", "L_outer"), cli::column<Real>("outer_deviation_phase", "phi"),
        cli::column<std::optional<Real>>("outer_deviation_real", "Re(z) (null if underflow)"),
        cli::column<std::optional<Real>>("outer_deviation_imag", "Im(z) (null if underflow)"));
  else if constexpr (TwoPairs)
    output.table(
        "string", "Two signed strings (L and sign remain authoritative)",
        [&](auto& table) {
          for (std::size_t i = 0; i < states.size(); ++i)
          {
            auto const& r = states[i].reference;
            for (std::size_t p = 0; p < 2; ++p)
            {
              Real const d = Real(r.deviation_signs[p]) * std::exp(-r.log_deviations[p]);
              table.append(i, p, r.string_labels[p], r.centers[p], r.deviation_signs[p], r.log_deviations[p],
                           d != Real{0} ? std::optional<Real>{d} : std::nullopt);
            }
          }
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("pair_index"),
        cli::column<std::size_t>("string_label", "J"), cli::column<Real>("center", "a"),
        cli::column<int>("deviation_sign", "sign(d)"), cli::column<Real>("log_deviation", "L=-log|d|"),
        cli::column<std::optional<Real>>("deviation", "d (null if underflow)"));
  else if constexpr (Defects == 2 || PairMixed)
    output.table(
        "string", "Signed two-string parameters (L is authoritative even on underflow)",
        [&](auto& table) {
          for (std::size_t i = 0; i < states.size(); ++i)
          {
            auto const& r = states[i].reference;
            Real const d = Real(r.deviation_sign) * std::exp(-r.log_deviation);
            table.append(i, r.string_label, r.center, r.deviation_sign, r.log_deviation,
                         d != Real{0} ? std::optional<Real>{d} : std::nullopt);
          }
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("string_label", "J"),
        cli::column<Real>("center", "a"), cli::column<int>("deviation_sign", "sign(d)"),
        cli::column<Real>("log_deviation", "L=-log|d|"),
        cli::column<std::optional<Real>>("deviation", "d (null if underflow)"));
  else
    output.table(
        "string", "Complex three-string parameters (L and phi are authoritative even on underflow)",
        [&](auto& table) {
          for (std::size_t i = 0; i < states.size(); ++i)
          {
            auto const& r = states[i].reference;
            Real const magnitude = std::exp(-r.log_deviation);
            table.append(
                i, r.string_label, r.center, r.log_deviation, r.deviation_phase,
                magnitude != Real{0} ? std::optional<Real>{magnitude * std::cos(r.deviation_phase)} : std::nullopt,
                magnitude != Real{0} ? std::optional<Real>{magnitude * std::sin(r.deviation_phase)} : std::nullopt);
          }
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("string_label", "J"),
        cli::column<Real>("center", "a"), cli::column<Real>("log_deviation", "L=-log|z|"),
        cli::column<Real>("deviation_phase", "phi (radians)"),
        cli::column<std::optional<Real>>("deviation_real", "Re(z) (null if underflow)"),
        cli::column<std::optional<Real>>("deviation_imag", "Im(z) (null if underflow)"));
  if (args.roots)
    output.table(
        "roots", "Reference cluster roots (rounded u; string parameters retain the deviation)",
        [&](auto& table) {
          for (std::size_t i = 0; i < states.size(); ++i)
          {
            auto const& r = states[i].reference;
            if constexpr (Quartet)
            {
              Real const eta = std::acosh(r.delta), d = Real(r.inner_deviation_sign) * std::exp(-r.inner_log_deviation);
              Real const z = std::exp(-r.outer_log_deviation), inner = (eta + d) / Real{2};
              Real const outer = inner + eta + z * std::cos(r.outer_deviation_phase);
              Real const imag = r.center / Real{2} + z * std::sin(r.outer_deviation_phase);
              table.append(i, 0, inner, r.center / Real{2});
              table.append(i, 1, inner, -r.center / Real{2});
              table.append(i, 2, outer, imag);
              table.append(i, 3, outer, -imag);
            }
            else if constexpr (TwoPairs)
            {
              for (std::size_t p = 0; p < 2; ++p)
              {
                Real const width = std::acosh(r.delta) + Real(r.deviation_signs[p]) * std::exp(-r.log_deviations[p]);
                table.append(i, 2 * p, width / Real{2}, r.centers[p] / Real{2});
                table.append(i, 2 * p + 1, width / Real{2}, -r.centers[p] / Real{2});
              }
            }
            else if constexpr (Defects == 2 || PairMixed)
            {
              Real const width = std::acosh(r.delta) + Real(r.deviation_sign) * std::exp(-r.log_deviation);
              if constexpr (Mixed)
                for (std::size_t j = 0; j < r.rapidities.size(); ++j)
                  table.append(i, j, Real{0}, r.rapidities[j] / Real{2});
              for (std::size_t j = 0; j < 2; ++j)
                table.append(i, j + r.rapidities.size(), width / Real{2}, (j ? -r.center : r.center) / Real{2});
            }
            else
            {
              Real const magnitude = std::exp(-r.log_deviation);
              Real const real = std::acosh(r.delta) + magnitude * std::cos(r.deviation_phase);
              Real const imag = r.center / Real{2} + magnitude * std::sin(r.deviation_phase);
              table.append(i, 0, Real{0}, r.center / Real{2});
              table.append(i, 1, real, imag);
              table.append(i, 2, real, -imag);
              if constexpr (TripleMixed) table.append(i, 3, Real{0}, r.rapidity / Real{2});
            }
          }
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("index"), cli::column<Real>("u_real"),
        cli::column<Real>("u_imag"));
  write_spin_content(output, args, states.size() + 1,
                     [&](std::size_t i) { return i < states.size() ? states[i].through_lines : ground.through_lines; });
  output.finish();
  if (!complete) std::cerr << "Bound-cluster calculation incomplete; unconverged modes have no verified gap.\n";
  return complete ? 0 : 2;
}

} // namespace bethe::apps::biquadratic
