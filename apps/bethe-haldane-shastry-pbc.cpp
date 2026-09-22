// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "citation-report.hpp"
#include "report-common.hpp"
#include <array>
#include <bethe/haldane_shastry.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::haldane_shastry;
struct Arguments
{
    std::size_t sites = 0, max_motifs = 100000;
    std::optional<std::string_view> motif, levels, sz;
    std::string_view precision = "fp64", format = "auto";
};
void usage(std::ostream& out)
{
  out << "Usage: bethe-haldane-shastry-pbc N [options]\n"
      << "Spin-1/2 inverse-chord-square ring, J=1, zero field, 2<=N<=1000000.\n"
      << "H=(pi/N)^2 sum_{i<j} S_i.S_j / sin^2(pi*(i-j)/N).\n"
      << "Default: ground multiplet(s), including both odd-ring momenta.\n"
      << "  --sz VALUE             lowest energy in this Sz sector\n"
      << "  --motif LIST           explicit positions, e.g. 1,3,5; empty list is allowed\n"
      << "  --levels COUNT|all     lowest COUNT Yangian multiplets, or every motif\n"
      << "  --max-motifs COUNT     enumeration budget (100000); preflight refusal, exit 2\n"
      << "  --precision fp64|long-double|fp128 (fp64; fp128 requires MPLAPACK)\n"
      << "  --format auto|pretty|plain|csv|tsv (auto)\n"
      << "  --help                 show this help and references\n"
      << "Modes --sz, --motif and --levels are mutually exclusive.\n"
      << "Motifs obey 1<=m<=N-1 and neighboring positions differ by at least 2.\n"
      << "A row is a Yangian multiplet, not one state or one distinct energy.\n"
      << "S_max is its largest SU(2) spin; degeneracy counts all spin projections,\n"
      << "even with --sz. A missing degeneracy means uint64 count overflow.\n"
      << "Gaps reference the global ground energy. No Newton solves or wavefunctions.\n"
      << "See docs/haldane-shastry.md for normalization and motif counting.\n";
  cli::print_citations(out, bethe::citations::Tool::haldane_shastry_pbc);
}
Arguments parse(int argc, char** argv)
{
  Arguments a;
  a.sites = cli::parse_size(argv[1]);
  for (int i = 2; i < argc; ++i)
  {
    std::string_view const option = argv[i];
    if (option != "--sz" && option != "--motif" && option != "--levels" && option != "--max-motifs" &&
        option != "--precision" && option != "--format")
      throw std::invalid_argument("unknown option: " + std::string(option));
    if (++i == argc) throw std::invalid_argument("missing value for " + std::string(option));
    std::string_view const value = argv[i];
    if (option == "--sz")
      a.sz = value;
    else if (option == "--motif")
      a.motif = value;
    else if (option == "--levels")
      a.levels = value;
    else if (option == "--max-motifs")
      a.max_motifs = cli::parse_size(value);
    else if (option == "--precision")
      a.precision = value;
    else
      a.format = value;
  }
  if (int(a.sz.has_value()) + int(a.motif.has_value()) + int(a.levels.has_value()) > 1)
    throw std::invalid_argument("--sz, --motif and --levels are mutually exclusive");
  if (a.format != "auto" && a.format != "plain" && a.format != "pretty" && a.format != "csv" && a.format != "tsv")
    throw std::invalid_argument("unknown output format: " + std::string(a.format));
  return a;
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
template <uni20::Real Real> int run(Arguments const& a)
{
  cli::CpuTimer timer;
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
  auto const cpu = timer.elapsed_text();
  bool const complete = !scan || scan->complete;
  cli::report_builder report("Haldane-Shastry spin ring");
  report.field("Hamiltonian", "H=(pi/N)^2 sum_{i<j} S_i.S_j/sin^2(pi*(i-j)/N); J=1")
      .field("Sites", a.sites)
      .field("Precision", a.precision)
      .field("Calculation", a.motif    ? "specified motif"
                            : a.levels ? "energy-ordered motif spectrum"
                            : a.sz     ? "Sz-sector minimum"
                                       : "global ground multiplets")
      .field("Status", complete ? "exact spectral rules" : "motif budget exceeded; no levels published")
      .field("Momentum", "P=2*pi*momentum_index/N modulo 2*pi")
      .field("Degeneracy", "whole Yangian multiplet; S_max is not a unique total spin")
      .field("CPU time", cpu);
  if (a.sz) report.field("Selected Sz", *a.sz);
  if (scan && scan->complete) report.field("Motifs enumerated", scan->total_motifs);
  bool const separated = a.format == "csv" || a.format == "tsv";
  char const separator = a.format == "tsv" ? '\t' : ',';
  if (separated)
  {
    for (auto const& [key, value] : report.fields())
      std::cout << "# " << key << ": " << value << '\n';
    std::cout << "motif" << separator << "energy" << separator << "gap" << separator << "momentum_index" << separator
              << "p" << separator << "spinons" << separator << "s_max" << separator << "degeneracy\n";
  }
  auto& table = report.table("Yangian multiplets");
  table.header_separator()
      .column("Motif")
      .column("Energy")
      .column("Gap")
      .column("Momentum index")
      .column("P")
      .column("Spinons")
      .column("S_max")
      .column("Degeneracy");
  for (auto const& s : levels)
  {
    std::string motif;
    for (auto m : s.motif)
    {
      if (!motif.empty()) motif += ' ';
      motif += std::to_string(m);
    }
    if (motif.empty()) motif = "empty";
    std::array<std::string, 8> const cells{motif,
                                           uni20::format_real(s.energy),
                                           uni20::format_real(s.gap),
                                           std::to_string(s.momentum_index),
                                           uni20::format_real(s.momentum),
                                           std::to_string(s.spinons),
                                           uni20::to_string_fraction(s.maximum_spin),
                                           s.degeneracy ? std::to_string(*s.degeneracy) : ""};
    if (separated)
    {
      for (std::size_t i = 0; i < cells.size(); ++i)
      {
        if (i) std::cout << separator;
        std::cout << cells[i];
      }
      std::cout << '\n';
    }
    else
      table.row(cells[0], cells[1], cells[2], cells[3], cells[4], cells[5], cells[6], cells[7]);
  }
  if (!separated) cli::print_report(report, a.format);
  if (!complete) std::cerr << "Motif budget exceeded; raise --max-motifs. No lowest levels claimed.\n";
  return complete ? 0 : 2;
}
} // namespace
int main(int argc, char** argv)
{
  if (argc == 1)
  {
    usage(std::cerr);
    return 1;
  }
  for (int i = 1; i < argc; ++i)
    if (std::string_view(argv[i]) == "--help")
    {
      usage(std::cout);
      return 0;
    }
  try
  {
    auto const a = parse(argc, argv);
    return cli::dispatch_precision(a.precision, [&]<typename Real>() { return run<Real>(a); });
  }
  catch (std::exception const& e)
  {
    std::cerr << "Error: " << e.what() << '\n';
    return 1;
  }
}
