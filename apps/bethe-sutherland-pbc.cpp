// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "citation-report.hpp"
#include "report-common.hpp"
#include <bethe/sutherland.hpp>
#include <charconv>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::sutherland;
struct Arguments
{
    std::size_t particles = 0, max_states = 100000;
    std::optional<std::size_t> window;
    std::optional<std::string_view> length, lambda, labels, levels;
    std::string_view precision = "fp64", format = "auto";
    bool pseudomomenta = false;
};
void usage(std::ostream& out)
{
  out << "Usage: bethe-sutherland-pbc N --length L --lambda VALUE [options]\n"
      << "Scalar bosonic trigonometric inverse-square gas on a periodic circle.\n"
      << "H=-sum d_i^2+2*lambda*(lambda-1)*(pi/L)^2 sum_{i<j} csc^2(pi*(x_i-x_j)/L).\n"
      << "hbar=2m=1; 1<=N<=1000000, L>0, lambda>=0; default: ground state.\n"
      << "The collision branch is psi ~ |x_i-x_j|^lambda: lambda=0 is free bosons;\n"
      << "lambda=1 is hard-core bosons, despite the same zero potential coefficient.\n"
      << "  --labels LIST          exactly N nondecreasing integers, e.g. -1,0,0,2\n"
      << "  --levels COUNT|all     lowest COUNT states, or all, WITHIN the label window\n"
      << "  --window W             explicit scan window -W<=n_j<=W (0<=W<=1000000)\n"
      << "  --max-states COUNT     enumeration budget (100000); refusal exits 2\n"
      << "  --pseudomomenta        also display the exact-rule k_j values\n"
      << "  --precision fp64|long-double|fp128 (fp64; fp128 requires MPLAPACK)\n"
      << "  --format auto|pretty|plain|csv|tsv (auto)\n"
      << "  --help                 show this help and references\n"
      << "--levels requires --window and vice versa; --labels excludes both.\n"
      << "Labels may repeat and include uniform integer boosts. P is NOT modulo 2*pi.\n"
      << "A finite window is not the infinite spectrum or a global low-energy guarantee.\n"
      << "No root solving, wavefunctions, spin or alternative collision domains.\n"
      << "See docs/sutherland.md for conventions and resource limits.\n";
  cli::print_citations(out, bethe::citations::Tool::sutherland_pbc);
}
Arguments parse(int argc, char** argv)
{
  Arguments a;
  a.particles = cli::parse_size(argv[1]);
  for (int i = 2; i < argc; ++i)
  {
    std::string_view const option = argv[i];
    if (option == "--pseudomomenta")
    {
      a.pseudomomenta = true;
      continue;
    }
    if (option != "--length" && option != "--lambda" && option != "--labels" && option != "--levels" &&
        option != "--window" && option != "--max-states" && option != "--precision" && option != "--format")
      throw std::invalid_argument("unknown option: " + std::string(option));
    if (++i == argc) throw std::invalid_argument("missing value for " + std::string(option));
    std::string_view const value = argv[i];
    if (option == "--length")
      a.length = value;
    else if (option == "--lambda")
      a.lambda = value;
    else if (option == "--labels")
      a.labels = value;
    else if (option == "--levels")
      a.levels = value;
    else if (option == "--window")
      a.window = cli::parse_size(value);
    else if (option == "--max-states")
      a.max_states = cli::parse_size(value);
    else if (option == "--precision")
      a.precision = value;
    else
      a.format = value;
  }
  if (!a.length || !a.lambda) throw std::invalid_argument("--length and --lambda are required");
  if (a.levels.has_value() != a.window.has_value())
    throw std::invalid_argument("--levels and --window must be supplied together");
  if (a.labels && a.levels) throw std::invalid_argument("--labels and --levels are mutually exclusive");
  if (a.format != "auto" && a.format != "plain" && a.format != "pretty" && a.format != "csv" && a.format != "tsv")
    throw std::invalid_argument("unknown output format: " + std::string(a.format));
  return a;
}
std::vector<std::int64_t> parse_labels(std::string_view text)
{
  std::vector<std::int64_t> out;
  for (;;)
  {
    auto const comma = text.find(',');
    auto const token = text.substr(0, comma);
    std::int64_t label;
    auto const [end, ec] = std::from_chars(token.data(), token.data() + token.size(), label);
    if (ec != std::errc{} || end != token.data() + token.size())
      throw std::invalid_argument("labels must be comma-separated int64 integers");
    out.push_back(label);
    if (comma == std::string_view::npos) return out;
    text.remove_prefix(comma + 1);
  }
}
template <typename Range, typename Format> std::string join(Range const& values, Format format)
{
  std::string text;
  for (auto const& value : values)
  {
    if (!text.empty()) text += ' ';
    text += format(value);
  }
  return text;
}
template <uni20::Real Real> int run(Arguments const& a)
{
  cli::CpuTimer timer;
  Real const length = uni20::parse_real<Real>(*a.length), lambda = uni20::parse_real<Real>(*a.lambda);
  std::vector<model::State<Real>> states;
  std::optional<model::Spectrum<Real>> scan;
  if (a.labels)
    states = {model::evaluate(a.particles, length, lambda, parse_labels(*a.labels))};
  else if (a.levels)
  {
    auto const count = *a.levels == "all" ? std::optional<std::size_t>{} : cli::parse_size(*a.levels);
    scan = model::spectrum(a.particles, length, lambda, *a.window, count, a.max_states);
    states = std::move(scan->states);
  }
  else
    states = {model::ground_state(a.particles, length, lambda)};
  auto const cpu = timer.elapsed_text();
  bool const complete = !scan || scan->complete;
  cli::report_builder report("Sutherland gas on a circle");
  report.field("Hamiltonian", "H=-sum d_i^2+2*lambda*(lambda-1)*(pi/L)^2 sum_{i<j} csc^2(pi*(x_i-x_j)/L)")
      .field("Units", "hbar=2m=1")
      .field("Statistics/domain", "periodic bosons; collision branch |x_i-x_j|^lambda")
      .field("Particles", a.particles)
      .field("Length", uni20::format_real(length))
      .field("Lambda", uni20::format_real(lambda))
      .field("Precision", a.precision)
      .field("Calculation", a.labels   ? "specified integer labels"
                            : a.levels ? "label-window spectrum"
                                       : "ground state")
      .field("Status", complete ? "exact spectral rules" : "state budget exceeded; no levels published")
      .field("Momentum", "P=2*pi*momentum_index/L; not reduced modulo 2*pi")
      .field("CPU time", cpu);
  if (!states.empty()) report.field("Ground energy", uni20::format_real(states.front().ground_energy));
  if (scan)
  {
    report.field("Label window", "[-" + std::to_string(*a.window) + "," + std::to_string(*a.window) + "]")
        .field("Coverage", "window only; no global low-energy completeness claimed");
    if (complete) report.field("States enumerated", scan->total_states);
  }
  bool const separated = a.format == "csv" || a.format == "tsv";
  char const separator = a.format == "tsv" ? '\t' : ',';
  if (separated)
  {
    for (auto const& [key, value] : report.fields())
      std::cout << "# " << key << ": " << value << '\n';
    std::cout << "labels" << separator << "energy" << separator << "gap" << separator << "momentum_index" << separator
              << "p";
    if (a.pseudomomenta) std::cout << separator << "pseudomomenta";
    std::cout << '\n';
  }
  auto& table = report.table("States");
  table.header_separator().column("Labels").column("Energy").column("Gap").column("Momentum index").column("P");
  if (a.pseudomomenta) table.column("Pseudomomenta");
  for (auto const& s : states)
  {
    std::vector<std::string> cells{join(s.labels, [](auto n) { return std::to_string(n); }),
                                   uni20::format_real(s.energy), uni20::format_real(s.gap),
                                   std::to_string(s.momentum_index), uni20::format_real(s.momentum)};
    if (a.pseudomomenta) cells.push_back(join(s.pseudomomenta, [](Real k) { return uni20::format_real(k); }));
    if (separated)
    {
      for (std::size_t i = 0; i < cells.size(); ++i)
      {
        if (i) std::cout << separator;
        std::cout << cells[i];
      }
      std::cout << '\n';
    }
    else if (a.pseudomomenta)
      table.row(cells[0], cells[1], cells[2], cells[3], cells[4], cells[5]);
    else
      table.row(cells[0], cells[1], cells[2], cells[3], cells[4]);
  }
  if (!separated) cli::print_report(report, a.format);
  if (!complete)
    std::cerr << "State budget exceeded; raise --max-states or narrow --window. No lowest levels claimed.\n";
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
