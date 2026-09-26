// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/tb_dispersion.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::takhtajan_babujian;
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
  auto info = cli::program_info("bethe-tb-dispersion", "Spin-1 Takhtajan-Babujian spinons and continuum bounds.",
                                bethe::citations::Tool::tb_dispersion);
  info.examples = {
      {"bethe-tb-dispersion --points 129 --csv tb.csv", "Spinon and two-/four-spinon reference curves"},
      {"bethe-tb-dispersion --branch four-spinon --folded --precision long-double", "Two-site iMPS envelope"}};
  info.notes = {"Infinite zero-field spin-1 H=J sum [S.S-(S.S)^2], J>0, no additive constant. Our J=1 is J_paper=4 "
                "in Vlijm-Caux. For cos(theta),sin(theta) normalization at theta=-pi/4, set J=1/sqrt(2).",
                "Spinon: spin 1/2, p in [0,pi], energy=2*pi*J*sin(p). Integer-spin periodic chains require even "
                "spinon number. Two-spinon sectors S=0,1; four-spinon sectors include S=0,1,2.",
                "Local spin response from a singlet selects S=1; local quadrupole selects S=2, requiring at least "
                "four spinons. Shared sine-band kinematics do not give XXX state counting: TB has SU(2)_2 content.",
                "Continuum Q in [0,2*pi]. --folded unions the Q and Q+pi images; cell_momentum=2*p modulo 2*pi. "
                "Folding is an output convention, not a claim of a gapped dimerized vacuum at the critical TB point.",
                "Bounds are for specified particle contents, not spectral intensities or finite-ring excited levels. "
                "The four-spinon upper edge does not bound higher-particle continua.",
                "See docs/tb-dispersion.md; --references shows literature and applicability."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  app.add_option("--exchange", a.exchange, "Bilinear exchange J>0 (biquadratic coefficient -J)")
      ->capture_default_str()
      ->type_name("REAL");
  app.add_option("--branch", a.branch, "Elementary or multiparticle branch")
      ->check(CLI::IsMember({"spinon", "two-spinon", "four-spinon", "all"}))
      ->capture_default_str();
  app.add_flag("--folded", a.folded, "Two-site folding of continuum envelopes only");
  auto* points = cli::count_option(app, "--points", a.points, "Samples per branch, 2..1000000")->capture_default_str();
  cli::text_option(app, "--momentum", a.momentum, "Single momentum, in radians/site")
      ->type_name("REAL")
      ->excludes(points);
  cli::precision_option(app, a.precision);
  cli::add_data_output_options(app, a.output, true);
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const exchange = uni20::parse_real<Real>(a.exchange), pi = bethe::detail::pi<Real>();
  model::SpinonDispersion<Real> const band(exchange);
  std::optional<Real> momentum;
  if (a.momentum)
  {
    momentum = uni20::parse_real<Real>(*a.momentum);
    Real const end = a.branch == "spinon" || a.branch == "all" ? pi : Real{2} * pi;
    if (!uni20::isfinite(*momentum) || *momentum < Real{0} || *momentum > end)
      throw std::invalid_argument("momentum outside a selected branch's range");
  }
  cli::RunReport report(context, "Takhtajan-Babujian thermodynamic spinons");
  report.field("exchange", "Exchange J", exchange)
      .field("hamiltonian", "Hamiltonian", "H=J sum [S.S-(S.S)^2]; spin 1, no constant")
      .field("energy_reference", "Energy reference",
             "excitation energy above the zero-field thermodynamic ground state")
      .field("velocity", "Spinon velocity", band.velocity())
      .field("spinon_spin", "Spinon spin", band.spin())
      .field("precision", "Precision", a.precision)
      .field("branches", "Branches", a.branch)
      .field("momentum_convention", "Continuum momentum", a.folded ? "two-site folded envelope" : "one-site unfolded")
      .field("cell_momentum", "Two-site translation phase", "2*p modulo 2*pi; not radians/site")
      .field("local_response", "Local response sectors", "spin: S=1; quadrupole: S=2; spectral weights not calculated")
      .field("spectrum_scope", "Scope",
             "spinon lines and two-/four-spinon bounds, not finite-size levels or full spectral support");
  if (momentum)
    report.field("momentum", "Requested momentum", *momentum);
  else
    report.field("points", "Points per branch", a.points);
  namespace data = cli::data;
  using Optional = std::optional<Real>;
  using Spin = std::optional<uni20::half_int>;
  auto table = data::make_data_table(
      "TB spinon lines and continuum bounds",
      {.retain = a.output.retain ? data::retention::all : data::retention::none, .metadata = report.metadata()},
      cli::column<std::string>("branch"), cli::column<Spin>("spin"), cli::column<std::string>("sectors"),
      cli::column<Real>("p").unit("radians/site"), cli::column<Real>("p_over_pi"),
      cli::column<Real>("cell_momentum").unit("radians/cell"), cli::column<Optional>("energy"),
      cli::column<Optional>("lower"), cli::column<Optional>("upper"), cli::column<std::string>("status"));
  bool const complete = cli::stream_result_table(report, context, a.output, "dispersion", table, [&](auto& table) {
    bool success = true;
    for (std::string const branch : {"spinon", "two-spinon", "four-spinon"})
    {
      if (a.branch != "all" && a.branch != branch) continue;
      bool const elementary = branch == "spinon";
      Real const end = elementary ? pi : Real{2} * pi;
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
              energy = band.energy(p);
            else
            {
              auto const convention = a.folded ? bethe::SpinonMomentum::folded : bethe::SpinonMomentum::unfolded;
              auto const bounds =
                  branch == "two-spinon" ? band.continuum(p, convention) : band.four_spinon_continuum(p, convention);
              lower = bounds.lower;
              upper = bounds.upper;
            }
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
        Real const cell = Real{2} * (p < pi ? p : p < Real{2} * pi ? p - pi : Real{0});
        table.append(branch, elementary ? Spin{band.spin()} : Spin{},
                     elementary               ? "S=1/2"
                     : branch == "two-spinon" ? "S=0,1"
                                              : "S=0,1,2",
                     p, p / pi, cell, energy, lower, upper, status);
      }
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
        if (a.folded && a.branch == "spinon") throw std::invalid_argument("--folded changes continuum envelopes only");
        return cli::dispatch_precision(a.precision, [&]<typename Real>() { return run<Real>(a, argc, argv); });
      });
}
