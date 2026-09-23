// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "data-output-options.hpp"
#include "program-options.hpp"
#include <bethe/sutherland.hpp>
#include <charconv>

namespace
{
namespace cli = bethe::cli;
namespace data = cli::data;
namespace model = bethe::sutherland;
struct Arguments
{
    std::size_t particles = 0, max_states = 100000;
    std::optional<std::size_t> window;
    std::optional<std::string> labels, levels;
    std::string length, lambda, precision = "fp64";
    bool pseudomomenta = false;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-sutherland-pbc", "Scalar bosonic trigonometric inverse-square gas on a circle.",
                                bethe::citations::Tool::sutherland_pbc);
  info.notes = {"H=-sum d_i^2+2*lambda*(lambda-1)*(pi/L)^2 sum_{i<j} csc^2(pi*(x_i-x_j)/L). "
                "hbar=2m=1; 1<=N<=1000000, L>0, lambda>=0; default: ground state.",
                "Collision branch psi ~ |x_i-x_j|^lambda: lambda=0 is free bosons; lambda=1 is hard-core bosons, "
                "despite the same zero potential coefficient. Labels may repeat and include uniform integer boosts. "
                "P is NOT modulo 2*pi.",
                "A finite label window is not the infinite spectrum or a global low-energy guarantee. "
                "Enumeration-budget refusal publishes no states and exits 2. No root solving, wavefunctions or spin.",
                "Tables: states; pseudomomenta with --pseudomomenta (joined by state_id). "
                "See docs/sutherland.md and docs/output.md. Use --references for literature and applicability."};
  info.examples = {{"bethe-sutherland-pbc 3 --length 4 --lambda 2 --levels all --window 2 --json spectrum.json",
                    "All states within a finite label window"}};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::count_option(app, "N", a.particles, "Number of particles")->required();
  app.add_option("--length", a.length, "Circle length L>0")->required()->type_name("REAL");
  app.add_option("--lambda", a.lambda, "Collision exponent lambda>=0")->required()->type_name("REAL");
  auto* labels = cli::text_option(app, "--labels", a.labels, "Exactly N nondecreasing integer labels, e.g. -1,0,2");
  auto* levels = cli::text_option(app, "--levels", a.levels, "Lowest COUNT states, or all, within the label window")
                     ->type_name("COUNT|all")
                     ->excludes(labels);
  auto* window = cli::count_option(app, "--window", a.window, "Scan window -W<=n_j<=W; W<=1000000")
                     ->needs(levels)
                     ->excludes(labels);
  levels->needs(window);
  cli::count_option(app, "--max-states", a.max_states, "Enumeration budget; refusal exits 2")->capture_default_str();
  app.add_flag("--pseudomomenta", a.pseudomomenta, "Include exact-rule k_j values in a separate table");
  cli::precision_option(app, a.precision);
  cli::add_data_output_options(app, a.output, true);
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
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  cli::CpuTimer timer;
  Real const length = uni20::parse_real<Real>(a.length), lambda = uni20::parse_real<Real>(a.lambda);
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
  auto metadata = cli::provenance("bethe-sutherland-pbc", argc, argv);
  metadata.insert({{"Hamiltonian", "H=-sum d_i^2+2*lambda*(lambda-1)*(pi/L)^2 sum_{i<j} csc^2(pi*(x_i-x_j)/L)"},
                   {"Units", "hbar=2m=1"},
                   {"Statistics/domain", "periodic bosons; collision branch |x_i-x_j|^lambda"},
                   {"Particles", std::to_string(a.particles)},
                   {"Length", uni20::format_real(length)},
                   {"Lambda", uni20::format_real(lambda)},
                   {"Precision", a.precision},
                   {"Calculation", a.labels   ? "specified integer labels"
                                   : a.levels ? "label-window spectrum"
                                              : "ground state"},
                   {"Momentum", "P=2*pi*momentum_index/L; not reduced modulo 2*pi"}});
  if (!states.empty()) metadata.emplace("Ground energy", uni20::format_real(states.front().ground_energy));
  if (scan)
  {
    metadata.emplace("Label window", "[-" + std::to_string(*a.window) + "," + std::to_string(*a.window) + "]");
    metadata.emplace("Coverage", "window only; no global low-energy completeness claimed");
    if (complete) metadata.emplace("States enumerated", std::to_string(scan->total_states));
  }
  data::data_table_options table_options{.retain = a.output.retain ? data::retention::all : data::retention::none,
                                         .metadata = metadata};
  auto table = data::make_data_table(
      "Sutherland states", table_options, data::data_column<std::size_t>("state_id"),
      data::data_column<std::string>("labels"), data::data_column<Real>("energy").round_trip(),
      data::data_column<Real>("gap").round_trip(), data::data_column<std::int64_t>("momentum_index"),
      data::data_column<Real>("p").round_trip());
  std::vector<std::string> names{"states"};
  if (a.pseudomomenta) names.push_back("pseudomomenta");
  cli::DataOutput output(a.output, names);
  data::table_metadata summary{
      {"Status", complete ? "exact spectral rules" : "state budget exceeded; no levels published"}, {"CPU time", cpu}};
  output.write_table(
      "states", table,
      [&](auto& table) {
        for (std::size_t i = 0; i < states.size(); ++i)
        {
          auto const& s = states[i];
          std::string labels;
          for (auto n : s.labels)
          {
            if (!labels.empty()) labels += ' ';
            labels += std::to_string(n);
          }
          table.append(i, labels, s.energy, s.gap, s.momentum_index, s.momentum);
        }
      },
      summary);
  if (a.pseudomomenta)
  {
    auto roots =
        data::make_data_table("Sutherland pseudomomenta", table_options, data::data_column<std::size_t>("state_id"),
                              data::data_column<std::size_t>("index"), data::data_column<std::int64_t>("label"),
                              data::data_column<Real>("k").round_trip());
    output.write_table(
        "pseudomomenta", roots,
        [&](auto& table) {
          for (std::size_t i = 0; i < states.size(); ++i)
            for (std::size_t j = 0; j < states[i].pseudomomenta.size(); ++j)
              table.append(i, j, states[i].labels[j], states[i].pseudomomenta[j]);
        },
        summary);
  }
  output.finish_document();
  if (!complete) std::cerr << "State budget exceeded; raise --max-states. No lowest levels claimed.\n";
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
