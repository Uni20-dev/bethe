// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/sine_gordon_particles.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::sine_gordon;
struct Arguments
{
    std::string mass = "1", precision = "fp64", branch = "all";
    std::optional<std::string> coupling, momentum, maximum;
    std::size_t points = 65;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info =
      cli::program_info("bethe-sine-gordon-dispersion", "Sine-Gordon particle lines and two-particle thresholds.",
                        bethe::citations::Tool::sine_gordon_dispersion);
  info.examples = {
      {"bethe-sine-gordon-dispersion --p 0.4 --csv sg.csv", "Attractive particle spectrum and pair onsets"},
      {"bethe-sine-gordon-dispersion --p 1 --momentum 0 --precision long-double", "Free Dirac masses"}};
  info.notes = {"Physical soliton mass M; p=beta^2/(8*pi-beta^2), kinetic term (1/2)(partial phi)^2; velocity=hbar=1.",
                "Soliton/antisoliton charges +/-1; neutral B_n exists only for n*p<1, mass=2*M*sin(n*pi*p/2).",
                "Energy=sqrt(mass^2+k^2). Pair rows are lower thresholds, not upper bounds or spectral intensities.",
                "Default grid: k=0..5*M. --momentum accepts either sign. No lattice Brillouin zone or folding.",
                "At p=1/3, s, anti-s and B1 have equal mass; B2/M=sqrt(3). This is the leading two-flavour "
                "Schwinger light-sector benchmark at theta=0, not the full massive Schwinger model.",
                "All stable species and unordered pairs are enumerated; limits: 256 breathers, 1000000 output rows.",
                "See docs/sine-gordon-excitations.md; --references prints literature."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::text_option(app, "--p", a.coupling, "Positive sine-Gordon coupling p")->required()->type_name("REAL");
  app.add_option("--mass", a.mass, "Soliton mass M>0")->capture_default_str()->type_name("REAL");
  app.add_option("--branch", a.branch, "Rows to emit")
      ->check(CLI::IsMember({"particles", "thresholds", "all"}))
      ->capture_default_str();
  auto* count = cli::count_option(app, "--points", a.points, "Grid points, 2..1000000")->capture_default_str();
  auto* maximum =
      cli::text_option(app, "--max-momentum", a.maximum, "Positive grid endpoint; default 5*M")->type_name("REAL");
  cli::text_option(app, "--momentum", a.momentum, "Single physical momentum k, either sign")
      ->type_name("REAL")
      ->excludes(count)
      ->excludes(maximum);
  cli::precision_option(app, a.precision);
  cli::add_data_output_options(app, a.output, true);
}
std::string name(model::Particle particle)
{
  switch (particle.kind)
  {
    case model::ParticleKind::soliton:
      return "soliton";
    case model::ParticleKind::antisoliton:
      return "antisoliton";
    case model::ParticleKind::breather:
      return "B" + std::to_string(particle.index);
  }
  throw std::invalid_argument("unknown particle");
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const mass = uni20::parse_real<Real>(a.mass), p = uni20::parse_real<Real>(*a.coupling);
  model::ParticleSpectrum<Real> const spectrum(mass, p);
  auto const particles = spectrum.particles(256);
  std::optional<Real> momentum;
  if (a.momentum) momentum = uni20::parse_real<Real>(*a.momentum);
  Real const maximum = momentum ? Real{0} : a.maximum ? uni20::parse_real<Real>(*a.maximum) : Real{5} * mass;
  if ((momentum && !uni20::isfinite(*momentum)) || (!momentum && (!uni20::isfinite(maximum) || maximum <= Real{0})))
    throw std::invalid_argument("require finite momentum or finite positive grid endpoint");
  struct Branch
  {
      std::string name, kind;
      int charge;
      std::vector<model::Particle> content;
      Real rest;
  };
  std::vector<Branch> branches;
  for (std::size_t i = 0; i < particles.size(); ++i)
  {
    auto const first = particles[i];
    if (a.branch != "thresholds")
      branches.push_back({name(first), "particle", first.charge(), {first}, spectrum.mass(first)});
    if (a.branch != "particles")
      for (std::size_t j = i; j < particles.size(); ++j)
      {
        std::vector<model::Particle> pair{first, particles[j]};
        branches.push_back({name(first) + "+" + name(particles[j]), "threshold", first.charge() + particles[j].charge(),
                            pair, spectrum.threshold(pair, Real{0})});
      }
  }
  std::size_t const count = momentum ? 1 : a.points;
  if (branches.size() > 1000000 / count) throw std::length_error("output exceeds 1000000 rows");
  cli::RunReport report(context, "Sine-Gordon infinite-volume excitations");
  report.field("mass", "Soliton mass M", mass)
      .field("p", "Coupling p", p)
      .field("convention", "Convention", "p=beta^2/(8*pi-beta^2); canonical kinetic term; velocity=hbar=1")
      .field("energy_reference", "Energy reference", "excitation energy above the infinite-volume vacuum")
      .field("scope", "Scope",
             "particle dispersions and pair lower thresholds; no spectral weights or finite-volume levels")
      .field("species", "Stable species", particles.size())
      .field("precision", "Precision", a.precision);
  if (momentum)
    report.field("momentum", "Requested momentum", *momentum);
  else
    report.field("maximum", "Maximum momentum", maximum).field("points", "Points per branch", count);
  using Optional = std::optional<Real>;
  auto table = cli::data::make_data_table(
      "Sine-Gordon dispersion and thresholds",
      {.retain = a.output.retain ? cli::data::retention::all : cli::data::retention::none,
       .metadata = report.metadata()},
      cli::column<std::string>("branch"), cli::column<std::string>("kind"), cli::column<int>("charge"),
      cli::column<Real>("rest_energy"), cli::column<Real>("momentum"), cli::column<Optional>("energy"),
      cli::column<std::string>("status"));
  bool const complete = cli::stream_result_table(report, context, a.output, "dispersion", table, [&](auto& table) {
    bool success = true;
    for (auto const& branch : branches)
      for (std::size_t j = 0; j < count; ++j)
      {
        Real const k = momentum ? *momentum : j == count - 1 ? maximum : maximum * (Real(j) / Real(count - 1));
        Optional energy;
        std::string status = "converged";
        try
        {
          context.measure([&] {
            energy = branch.content.size() == 1 ? spectrum.energy(branch.content[0], k)
                                                : spectrum.threshold(branch.content, k);
          });
        }
        catch (std::overflow_error const&)
        {
          status = "precision_limit";
          success = false;
        }
        catch (std::underflow_error const&)
        {
          status = "precision_limit";
          success = false;
        }
        table.append(branch.name, branch.kind, branch.charge, branch.rest, k, energy, status);
      }
    return success;
  });
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
        if (!a.momentum && (a.points < 2 || a.points > 1000000))
          throw std::invalid_argument("require 2<=points<=1000000");
        return cli::dispatch_precision(a.precision, [&]<typename Real>() { return run<Real>(a, argc, argv); });
      });
}
