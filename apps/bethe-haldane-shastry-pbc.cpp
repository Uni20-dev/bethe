// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "data-output-options.hpp"
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/haldane_shastry.hpp>

namespace
{
namespace cli = bethe::cli;
namespace data = cli::data;
namespace model = bethe::haldane_shastry;
struct Arguments
{
    std::size_t sites = 0, max_motifs = 100000;
    std::size_t max_spin_updates = 1000000;
    bool spin_content = false;
    std::optional<std::string> motif, levels, sz;
    std::string precision = "fp64";
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-haldane-shastry-pbc", "Spin-1/2 inverse-chord-square ring, J=1, zero field.",
                                bethe::citations::Tool::haldane_shastry_pbc);
  info.notes = {"H=(pi/N)^2 sum_{i<j} S_i.S_j/sin^2(pi*(i-j)/N); 2<=N<=1000000. "
                "Default: ground multiplets, including both odd-ring momenta.",
                "Motifs obey 1<=m<=N-1; neighboring positions differ by at least 2. A row is a Yangian multiplet, "
                "not one state or one distinct energy. S_max is its largest SU(2) spin; degeneracy counts all spin "
                "projections even with --sz. Missing degeneracy means uint64 overflow.",
                "Gaps reference the global ground energy. No Newton solves or wavefunctions. Table: levels. "
                "Enumeration-budget refusal publishes no levels and exits 2.",
                "--spin-content adds SU(2) irreps for each selected Yangian multiplet, including all projections "
                "even with --sz. Decomposition failure leaves energies valid but omits spin counts (exit 2).",
                "See docs/haldane-shastry.md and docs/output.md. Use --references for literature and applicability."};
  info.examples = {{"bethe-haldane-shastry-pbc 8 --levels all --csv levels.csv", "All Yangian multiplets"}};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::count_option(app, "N", a.sites, "Number of sites")->required();
  auto* modes = app.add_option_group("State selection")->require_option(0, 1);
  cli::text_option(*modes, "--sz", a.sz, "Lowest energy in this exact integer/half-integer Sz sector");
  cli::text_option(*modes, "--motif", a.motif, "Explicit positions, e.g. 1,3,5; empty list is allowed");
  cli::text_option(*modes, "--levels", a.levels, "Lowest COUNT Yangian multiplets, or all")->type_name("COUNT|all");
  cli::count_option(app, "--max-motifs", a.max_motifs, "Enumeration budget; refusal exits 2")->capture_default_str();
  app.add_flag("--spin-content", a.spin_content, "Add the spin_content table of SU(2) irrep multiplicities");
  cli::count_option(app, "--max-spin-updates", a.max_spin_updates, "Clebsch-Gordan addition budget per motif")
      ->capture_default_str()
      ->needs("--spin-content");
  cli::precision_option(app, a.precision);
  cli::add_data_output_options(app, a.output, true);
}
std::vector<std::size_t> parse_motif(std::string_view text)
{
  std::vector<std::size_t> out;
  if (text.empty()) return out;
  for (;;)
  {
    auto const comma = text.find(',');
    out.push_back(cli::parse_size(text.substr(0, comma)));
    if (comma == std::string_view::npos) return out;
    text.remove_prefix(comma + 1);
  }
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  auto computation = context.computation();
  std::vector<model::Level<Real>> levels;
  std::optional<model::Spectrum<Real>> scan;
  if (a.motif)
    levels = {model::evaluate<Real>(a.sites, parse_motif(*a.motif))};
  else if (a.levels)
  {
    auto const count = *a.levels == "all" ? std::optional<std::size_t>{} : cli::parse_size(*a.levels);
    scan = model::spectrum<Real>(a.sites, count, a.max_motifs);
    levels = std::move(scan->levels);
  }
  else if (a.sz)
    levels = model::sector_ground_levels<Real>(a.sites, uni20::half_int::parse(*a.sz));
  else
    levels = model::ground_levels<Real>(a.sites);
  std::vector<bethe::spin::Decomposition> contents;
  bool complete = !scan || scan->complete;
  if (a.spin_content)
    for (auto const& level : levels)
    {
      contents.push_back(model::spin_decomposition(a.sites, level.motif, {.max_updates = a.max_spin_updates}));
      complete = complete && contents.back().complete;
    }
  computation.finish();
  cli::RunReport report(context, "Haldane-Shastry Yangian multiplets");
  report.field("hamiltonian", "Hamiltonian", "H=(pi/N)^2 sum_{i<j} S_i.S_j/sin^2(pi*(i-j)/N); J=1")
      .field("sites", "Sites", a.sites)
      .field("precision", "Precision", a.precision)
      .field("calculation", "Calculation",
             a.motif    ? "specified motif"
             : a.levels ? "energy-ordered motif spectrum"
             : a.sz     ? "Sz-sector minimum"
                        : "global ground multiplets")
      .field("momentum", "Momentum", "P=2*pi*momentum_index/N modulo 2*pi")
      .field("degeneracy", "Degeneracy", "whole Yangian multiplet; S_max is not a unique total spin")
      .result(complete, complete                  ? "exact spectral rules"
                        : scan && !scan->complete ? "motif budget exceeded; no levels published"
                                                  : "spin decomposition incomplete; energies remain valid");
  if (a.spin_content)
    report.field("spin_content", "Spin content", "whole Yangian multiplet; multiplicity counts SU(2) irreps")
        .field("max_spin_updates", "Max spin updates per motif", a.max_spin_updates);
  if (a.sz) report.field("selected_sz", "Selected Sz", uni20::half_int::parse(*a.sz));
  if (scan && scan->complete) report.field("motifs_enumerated", "Motifs enumerated", scan->total_motifs);
  std::vector<std::string> names{"levels"};
  if (a.spin_content) names.push_back("spin_content");
  cli::ResultOutput output(report, a.output, std::move(names), false);
  output.table(
      "levels", "Haldane-Shastry Yangian multiplets",
      [&](auto& table) {
        for (std::size_t i = 0; i < levels.size(); ++i)
        {
          auto const& s = levels[i];
          std::string motif;
          for (auto m : s.motif)
          {
            if (!motif.empty()) motif += ' ';
            motif += std::to_string(m);
          }
          if (motif.empty()) motif = "empty";
          table.append(i, motif, s.energy, s.gap, s.momentum_index, s.momentum, s.spinons, s.maximum_spin,
                       s.degeneracy);
        }
      },
      data::data_column<std::size_t>("state_id"), data::data_column<std::string>("motif"),
      data::data_column<Real>("energy").round_trip(), data::data_column<Real>("gap").round_trip(),
      data::data_column<std::size_t>("momentum_index"), data::data_column<Real>("p").unit("radians").round_trip(),
      data::data_column<std::size_t>("spinons"), data::data_column<uni20::half_int>("s_max"),
      data::data_column<std::optional<std::uint64_t>>("degeneracy"));
  if (a.spin_content)
    output.table(
        "spin_content", "SU(2) content of selected Yangian multiplets",
        [&](auto& table) {
          using Spin = std::optional<uni20::half_int>;
          using Count = std::optional<std::uint64_t>;
          for (std::size_t i = 0; i < contents.size(); ++i)
          {
            auto const& content = contents[i];
            if (content.complete)
              for (auto const& term : content.multiplets)
                table.append(i, Spin(term.spin), Count(term.multiplicity), content.updates, true,
                             std::string("complete"));
            else
              table.append(i, Spin{}, Count{}, content.updates, false,
                           std::string(content.status == bethe::spin::DecompositionStatus::count_overflow
                                           ? "count_overflow"
                                           : "work_limit"));
          }
        },
        cli::column<std::size_t>("state_id"), cli::column<std::optional<uni20::half_int>>("spin"),
        cli::column<std::optional<std::uint64_t>>("multiplicity"), cli::column<std::size_t>("updates"),
        cli::column<bool>("complete"), cli::column<std::string>("status"));
  output.finish();
  if (!complete)
    std::cerr << (scan && !scan->complete
                      ? "Motif budget exceeded; raise --max-motifs. No lowest levels claimed.\n"
                      : "Spin decomposition incomplete; see spin_content status. Motif energies remain valid.\n");
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
        return cli::dispatch_precision(a.precision, [&]<typename Real>() { return run<Real>(a, argc, argv); });
      });
}
