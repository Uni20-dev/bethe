// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "execution-options.hpp"
#include "result-output.hpp"
#include <bethe/xxx_thermodynamic_structure_factor.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::heisenberg;
struct Arguments
{
    std::string exchange = "1", channel = "zz", precision = "fp64", omega_min = "0";
    std::optional<std::string> momentum, omega, omega_max, tolerance;
    std::size_t points = 201, momentum_points = 65, evaluations = 200000;
    int threads = 1;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-xxx-structure-factor-thermo", "Infinite-chain XXX two-spinon spectral density.",
                                bethe::citations::Tool::xxx_structure_factor_thermo);
  info.examples = {
      {"bethe-xxx-structure-factor-thermo --momentum 1.5 --points 501", "Unbroadened frequency cut"},
      {"bethe-xxx-structure-factor-thermo --momentum-points 129 --points 501 --threads 4 --csv spectrum.csv",
       "Momentum-frequency grid for a heat map"}};
  info.notes = {
      "Infinite spin-1/2 XXX chain, H=J sum S_j.S_(j+1), J>0, zero field and temperature.",
      "Exact two-spinon contribution only, NOT the full DSF. No broadening or sum-rule rescaling.",
      "Density is S(q,w) itself, not S/(2*pi) and not a finite-ring line weight; raising (S-+) = 2*zz.",
      "Momenta in [0,2*pi], radians/site. Grid includes endpoints; do not count the periodic seam twice.",
      "At the divergent lower threshold density/error are absent with status lower_threshold (not failure).",
      "At q=0,2*pi the density is zero. Upper-edge limit is zero. Negative frequencies are allowed and zero.",
      "Error estimates cover kernel quadrature/roundoff, not input/threshold conditioning or missing higher spinons.",
      "See docs/xxx-structure-factor-thermo.md; --references for literature."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  auto* nq = cli::count_option(app, "--momentum-points", a.momentum_points, "Uniform q grid, 2..1000000")
                 ->capture_default_str();
  cli::text_option(app, "--momentum", a.momentum, "Single momentum, radians/site")->excludes(nq)->type_name("REAL");
  auto* nw =
      cli::count_option(app, "--points", a.points, "Uniform frequency samples, 2..1000000")->capture_default_str();
  auto* lo = cli::option(app, "--omega-min", a.omega_min, "Frequency grid lower endpoint")->capture_default_str();
  auto* hi = cli::text_option(app, "--omega-max", a.omega_max, "Frequency grid upper endpoint (default pi*J)");
  cli::text_option(app, "--omega", a.omega, "Single frequency instead of a grid")
      ->excludes(nw)
      ->excludes(lo)
      ->excludes(hi);
  cli::option(app, "--exchange", a.exchange, "Antiferromagnetic exchange J>0")->capture_default_str();
  cli::option(app, "--channel", a.channel, "zz or raising (S-+)")
      ->check(CLI::IsMember({"zz", "raising"}))
      ->capture_default_str();
  cli::text_option(app, "--tolerance", a.tolerance, "Absolute I(rho) error target (default 4096 epsilon)");
  cli::count_option(app, "--max-evaluations", a.evaluations, "Kernel evaluation budget per point")
      ->capture_default_str();
  cli::precision_option(app, a.precision);
  cli::threads_option(app, a.threads)->description("Scheduler concurrency limit for independent spectral points");
  cli::add_data_output_options(
      app, a.output,
      {{.name = "spectrum", .description = "Unbroadened two-spinon density and numerical status", .primary = true},
       {.name = "continuum", .description = "Two-spinon lower/upper edges at each sampled q"}});
}
char const* status(model::SpectralDensityStatus s)
{
  switch (s)
  {
    case model::SpectralDensityStatus::converged:
      return "converged";
    case model::SpectralDensityStatus::outside_continuum:
      return "outside_continuum";
    case model::SpectralDensityStatus::lower_threshold:
      return "lower_threshold";
    case model::SpectralDensityStatus::numerical_failure:
      return "numerical_failure";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const pi = bethe::detail::pi<Real>(), exchange = uni20::parse_real<Real>(a.exchange);
  model::StructureFactorOptions<Real> options;
  options.max_evaluations = a.evaluations;
  if (a.tolerance) options.tolerance = uni20::parse_real<Real>(*a.tolerance);
  model::ThermodynamicTwoSpinonStructureFactor<Real> sf(exchange, options);
  auto grid = [](std::optional<std::string> const& one, Real lo, Real hi, std::size_t n) {
    if (one) return std::vector<Real>{uni20::parse_real<Real>(*one)};
    if (!uni20::isfinite(lo) || !uni20::isfinite(hi) || !(hi > lo) || n < 2 || n > 1000000)
      throw std::invalid_argument("invalid grid: finite increasing endpoints and 2..1000000 points required");
    std::vector<Real> values(n);
    for (std::size_t i = 0; i < n; ++i)
      values[i] = (Real{1} - Real(i) / Real(n - 1)) * lo + Real(i) / Real(n - 1) * hi;
    return values;
  };
  auto const q = grid(a.momentum, Real{0}, Real{2} * pi, a.momentum_points);
  auto const w = grid(a.omega, uni20::parse_real<Real>(a.omega_min),
                      a.omega_max ? uni20::parse_real<Real>(*a.omega_max) : pi * exchange, a.points);
  for (Real p : q)
    (void)sf.boundaries(p);
  for (Real f : w)
    if (!uni20::isfinite(f)) throw std::invalid_argument("frequency must be finite");
  if (q.size() > 1000000 / w.size()) throw std::length_error("spectral grid exceeds 1000000 points");
  std::vector<model::SpectralDensity<Real>> values(q.size() * w.size());
  Real const scale = a.channel == "zz" ? Real{1} : Real{2};
  uni20::async::TbbScheduler scheduler(a.threads);
  context.measure([&] {
    for (std::size_t offset = 0; offset < values.size(); offset += 128)
      scheduler.execute_batch(std::min(std::size_t{128}, values.size() - offset), [&](std::size_t j) {
        auto const i = offset + j;
        values[i] = sf(q[i / w.size()], w[i % w.size()]);
        auto& value = values[i];
        if (value.value && (!uni20::isfinite(scale * *value.value) || !uni20::isfinite(scale * *value.error)))
        {
          value.value.reset();
          value.error.reset();
          value.status = model::SpectralDensityStatus::numerical_failure;
        }
      });
  });
  std::size_t failures = 0, thresholds = 0;
  for (auto const& v : values)
  {
    failures += !v.converged();
    thresholds += v.status == model::SpectralDensityStatus::lower_threshold;
  }
  cli::RunReport report(context, "XXX thermodynamic two-spinon structure factor");
  report.field("hamiltonian", "Hamiltonian", "H=J sum S_j.S_(j+1); infinite chain; zero field and temperature")
      .field("exchange", "Exchange J", exchange)
      .field("channel", "Channel", a.channel)
      .field("family", "Family", "Thermodynamic two-spinon contribution; partial DSF; no rescaling")
      .field("density_convention", "Density convention", "S(q,w) itself; integrals use dw/(2*pi); raising=2*zz")
      .field("precision", "Precision", a.precision)
      .field("tolerance", "Kernel tolerance", options.tolerance)
      .field("max_evaluations", "Max kernel evaluations per point", options.max_evaluations)
      .field("threads", "Scheduler concurrency limit", a.threads)
      .field("points", "Spectral points", values.size())
      .field("failures", "Failed points", failures)
      .field("thresholds", "Divergent threshold points", thresholds)
      .result(failures == 0, failures == 0 ? "converged" : "numerical_failure");
  cli::ResultOutput output(report, a.output, {"spectrum", "continuum"});
  using Optional = std::optional<Real>;
  using cli::column;
  output.table(
      "spectrum", "Unbroadened spectral density",
      [&](auto& table) {
        for (std::size_t i = 0; i < values.size(); ++i)
        {
          auto const& v = values[i];
          table.append(q[i / w.size()], w[i % w.size()], v.value ? Optional(scale * *v.value) : Optional{},
                       v.error ? Optional(scale * *v.error) : Optional{}, v.evaluations, std::string(status(v.status)));
        }
      },
      column<Real>("q"), column<Real>("omega"), column<Optional>("density"), column<Optional>("estimated_error"),
      column<std::size_t>("evaluations"), column<std::string>("status"));
  output.table(
      "continuum", "Two-spinon continuum",
      [&](auto& table) {
        for (Real p : q)
        {
          auto const [lo, hi] = sf.boundaries(p);
          table.append(p, lo, hi);
        }
      },
      column<Real>("q"), column<Real>("lower"), column<Real>("upper"));
  output.finish();
  return failures == 0 ? 0 : 2;
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
