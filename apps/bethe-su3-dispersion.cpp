// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/su3_dispersion.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::su3;
struct Arguments
{
    std::string exchange = "1", branch = "all", precision = "fp64";
    std::optional<std::string> momentum;
    std::size_t points = 65;
    bool folded = false;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info =
      cli::program_info("bethe-su3-dispersion", "SU(3)/ULS elementary lines and multiparticle continuum bounds.",
                        bethe::citations::Tool::su3_dispersion);
  info.examples = {
      {"bethe-su3-dispersion --points 129 --csv uls.csv", "One-site momentum reference curves"},
      {"bethe-su3-dispersion --branch four-soliton --folded --json edges.json", "Three-site iMPS envelope"}};
  info.notes = {
      "Infinite zero-field H=J sum P, J>0. Spin-1 ULS H=J sum [S.S+(S.S)^2] differs only by +J per bond.",
      "3: p in [0,4*pi/3]; bar3: p in [0,2*pi/3]. An isolated soliton is not a balanced periodic-chain state.",
      "Two-soliton content: 3 x bar3 = 1 + 8. Four-soliton content: (3 x bar3)^2. Local spin and quadrupolar "
      "responses from a singlet transform in 8. Bounds do not predict spectral weights.",
      "Continuum Q in [0,2*pi]. --folded takes the union of three momentum images spaced by 2*pi/3; "
      "cell_momentum is the phase 3*p modulo 2*pi. Elementary lines are never folded into a single band.",
      "The two-soliton lower edge is not the full low-energy threshold between 2*pi/3 and 4*pi/3. "
      "Four-soliton upper bounds are not upper bounds on the full spectrum.",
      "See docs/su3-dispersion.md. Use --references for literature; no finite-size levels or spectral weights."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  app.add_option("--exchange", a.exchange, "Permutation exchange J>0")->capture_default_str()->type_name("REAL");
  app.add_option("--branch", a.branch, "Branch or specified particle content")
      ->check(CLI::IsMember({"3", "bar3", "two-soliton", "four-soliton", "all"}))
      ->capture_default_str();
  app.add_flag("--folded", a.folded, "Three-site folding of continuum envelopes only");
  auto* points = cli::count_option(app, "--points", a.points, "Samples per branch, 2..1000000")->capture_default_str();
  cli::text_option(app, "--momentum", a.momentum, "Single momentum (radians/site), valid for every selected branch")
      ->type_name("REAL")
      ->excludes(points);
  cli::precision_option(app, a.precision);
  cli::add_data_output_options(app, a.output, true);
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  model::ExcitationDispersion<Real> const band(uni20::parse_real<Real>(a.exchange));
  Real const pi = bethe::detail::pi<Real>(), soft = Real{2} * pi / Real{3};
  std::optional<Real> momentum;
  if (a.momentum)
  {
    momentum = uni20::parse_real<Real>(*a.momentum);
    Real const end = a.branch == "all" || a.branch == "bar3" ? soft : a.branch == "3" ? Real{2} * soft : Real{2} * pi;
    if (!uni20::isfinite(*momentum) || *momentum < Real{0} || *momentum > end)
      throw std::invalid_argument("momentum outside a selected branch's range");
  }
  cli::RunReport report(context, "SU(3)/ULS thermodynamic excitations");
  report.field("exchange", "Exchange J", band.exchange())
      .field("hamiltonian", "Hamiltonian", "H=J sum P; spin-1 ULS energy is E+J*L")
      .field("energy_reference", "Energy reference", "excitation energy above the zero-field thermodynamic vacuum")
      .field("velocity", "Soft-mode velocity", band.velocity())
      .field("precision", "Precision", a.precision)
      .field("branches", "Branches", a.branch)
      .field("momentum_convention", "Continuum momentum", a.folded ? "three-site folded envelope" : "one-site unfolded")
      .field("cell_momentum", "Three-site translation phase", "3*p modulo 2*pi; not radians/site")
      .field("local_response", "Local response sectors", "spin and quadrupole: SU(3) adjoint 8; weights not calculated")
      .field("spectrum_scope", "Scope",
             "elementary lines; 3 x bar3 and (3 x bar3)^2 envelopes, not full spectral support");
  if (momentum)
    report.field("momentum", "Requested momentum", *momentum);
  else
    report.field("points", "Points per branch", a.points);
  namespace data = cli::data;
  using Optional = std::optional<Real>;
  auto table = data::make_data_table(
      "SU(3) elementary lines and continuum bounds",
      {.retain = a.output.retain ? data::retention::all : data::retention::none, .metadata = report.metadata()},
      cli::column<std::string>("branch"), cli::column<std::string>("representations"),
      cli::column<Real>("p").unit("radians/site"), cli::column<Real>("p_over_pi"),
      cli::column<Real>("cell_momentum").unit("radians/cell"), cli::column<Optional>("energy"),
      cli::column<Optional>("lower"), cli::column<Optional>("upper"), cli::column<std::string>("status"));
  cli::DataOutput output(a.output, {"dispersion"});
  bool complete = true;
  try
  {
    output.attach(table, "dispersion");
    for (std::string const branch : {"3", "bar3", "two-soliton", "four-soliton"})
    {
      if (a.branch != "all" && a.branch != branch) continue;
      bool const elementary = branch == "3" || branch == "bar3";
      auto const particle = branch == "3" ? model::Particle::fundamental : model::Particle::antifundamental;
      Real const end = elementary ? band.momentum_max(particle) : Real{2} * pi;
      auto const count = momentum ? 1 : a.points;
      for (std::size_t i = 0; i < count; ++i)
      {
        Real const p = momentum ? *momentum : i == count - 1 ? end : end * (Real(i) / Real(count - 1));
        Optional energy, lower, upper;
        std::string status = "converged";
        try
        {
          context.measure([&] {
            if (elementary)
              energy = band.energy(particle, p);
            else
            {
              auto const bounds = band.continuum(
                  p, branch == "two-soliton" ? model::Continuum::two_soliton : model::Continuum::four_soliton,
                  a.folded ? model::Momentum::three_site_folded : model::Momentum::unfolded);
              lower = bounds.lower;
              upper = bounds.upper;
            }
          });
        }
        catch (std::overflow_error const&)
        {
          status = "precision_limit";
          complete = false;
        }
        catch (std::underflow_error const&)
        {
          status = "precision_limit";
          complete = false;
        }
        Real cell = p < soft ? p : p < Real{2} * soft ? p - soft : p < Real{2} * pi ? p - Real{2} * soft : Real{0};
        cell *= Real{3};
        table.append(branch,
                     elementary                ? branch
                     : branch == "two-soliton" ? "1+8"
                                               : "(1+8)x(1+8)",
                     p, p / pi, cell, energy, lower, upper, status);
      }
    }
    report.result(complete, complete ? "converged" : "precision_limit; unavailable quantities omitted");
    auto const summary = report.finish();
    output.overview(report.overview(summary));
    output.finish(table, summary);
    output.finish_document();
  }
  catch (...)
  {
    data::table_metadata aborted{{"Status", "aborted"}};
    if (!context.finished()) try
      {
        aborted = cli::run_summary(context.finish(uni20::run_outcome::failed));
      }
      catch (...)
      {}
    output.abort(table, std::move(aborted));
    throw;
  }
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
        if (a.folded && (a.branch == "3" || a.branch == "bar3"))
          throw std::invalid_argument("--folded changes continuum envelopes only");
        return cli::dispatch_precision(a.precision, [&]<typename Real>() { return run<Real>(a, argc, argv); });
      });
}
