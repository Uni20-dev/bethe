// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "execution-options.hpp"
#include "result-output.hpp"
#include "solver-precision.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/xxz_thermodynamic_structure_factor.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::xxz;
struct Arguments
{
    std::string delta, exchange = "1", precision = "fp64", omega_min = "0";
    std::optional<std::string> momentum, omega, omega_max, tolerance;
    std::size_t points = 201, momentum_points = 65, terms = 200000;
    int threads = 1;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-xxz-structure-factor-thermo",
                                "Massive XXZ thermodynamic transverse two-spinon spectral density.",
                                bethe::citations::Tool::xxz_structure_factor_thermo);
  info.examples = {
      {"bethe-xxz-structure-factor-thermo --delta 1.5 --momentum 1.5 --points 501", "Unbroadened frequency cut"},
      {"bethe-xxz-structure-factor-thermo --delta 1.5 --momentum-points 129 --points 501 --threads 4 --csv "
       "spectrum.csv",
       "Momentum-frequency grid"},
      {"bethe-xxz-structure-factor-thermo --delta 1.001 --momentum 1.5 --omega 1.8 --precision bits-160",
       "Approach the XXX raising-channel limit with MPFR enabled"}};
  info.notes = {
      "Infinite spin-1/2 chain, H=J sum(SxSx+SySy+Delta SzSz), J>0, Delta>1; zero field and temperature.",
      "Transverse S-+ two-spinon contribution only; no broadening or sum-rule rescaling. Szz is a different response.",
      "Density is S(k,w) itself, not S/(2*pi); frequency integrals use dw/(2*pi).",
      "Physical momentum k in [0,2*pi], radians/site; the paper uses Q=pi-k. The response is not pi-periodic.",
      "Both spectral sheets contribute with unequal weights. Exact sheet endpoints have absent values and status "
      "threshold.",
      "Error estimates cover the kernel tail and arithmetic conditioning, not missing sectors or input uncertainty.",
      "Raise precision near thresholds or the isotropic limit; unresolved arithmetic is not reported as zero.",
      "Numerical failures retain their grid rows and exit 2; threshold and outside_continuum are valid statuses.",
      "Streaming output is ordered and bounded when used with --no-retain. Do not count the periodic seam twice.",
      "See docs/xxz-structure-factor-thermo.md; --references for literature."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::option(app, "--delta", a.delta, "Massive antiferromagnetic anisotropy Delta>1")->required();
  cli::option(app, "--exchange", a.exchange, "Antiferromagnetic exchange J>0")->capture_default_str();
  auto* nk = cli::count_option(app, "--momentum-points", a.momentum_points, "Uniform momentum grid, 2..1000000")
                 ->capture_default_str();
  cli::text_option(app, "--momentum", a.momentum, "Single physical momentum, radians/site")->excludes(nk);
  auto* nw =
      cli::count_option(app, "--points", a.points, "Uniform frequency samples, 2..1000000")->capture_default_str();
  auto* lo = cli::option(app, "--omega-min", a.omega_min, "Frequency grid lower endpoint")->capture_default_str();
  auto* hi =
      cli::text_option(app, "--omega-max", a.omega_max, "Frequency grid upper endpoint (default twice spinon maximum)");
  cli::text_option(app, "--omega", a.omega, "Single frequency instead of a grid")
      ->excludes(nw)
      ->excludes(lo)
      ->excludes(hi);
  cli::text_option(app, "--tolerance", a.tolerance, "Absolute log-kernel tail target (default 4096 epsilon)");
  cli::count_option(app, "--max-terms", a.terms, "Series budget per spectral sheet")->capture_default_str();
  cli::precision_option(app, a.precision, cli::ArbitraryPrecision::real);
  cli::threads_option(app, a.threads)->description("Shared scheduler concurrency limit");
  cli::add_data_output_options(
      app, a.output,
      {{.name = "spectrum", .description = "Transverse two-spinon density, error estimate and status", .primary = true},
       {.name = "continuum", .description = "Separate C_plus and C_minus spectral-sheet bounds"}});
}
char const* status(model::SpectralDensityStatus value)
{
  switch (value)
  {
    case model::SpectralDensityStatus::converged:
      return "converged";
    case model::SpectralDensityStatus::outside_continuum:
      return "outside_continuum";
    case model::SpectralDensityStatus::threshold:
      return "threshold";
    case model::SpectralDensityStatus::numerical_failure:
      return "numerical_failure";
  }
  return "unknown";
}
template <uni20::Real R>
int run(Arguments const& a, model::TwoSpinonStructureFactorOptions<R> options, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  R const zero = uni20::scalar_at_precision<R>(0, options.precision()), one = uni20::scalar_like(zero, 1);
  R const pi = bethe::detail::spinon_pi(zero);
  auto parse = [&](std::string const& s) { return cli::parse_solver_real<R>(s, options); };
  R const delta = parse(a.delta), exchange = parse(a.exchange);
  options.max_terms = a.terms;
  if (a.tolerance) options.tolerance = parse(*a.tolerance);
  model::ThermodynamicTwoSpinonStructureFactor<R> const sf(delta, exchange, options);
  auto grid = [&](std::optional<std::string> const& point, R lo, R hi, std::size_t n) {
    if (point) return std::vector<R>{parse(*point)};
    if (!uni20::isfinite(lo) || !uni20::isfinite(hi) || !(hi > lo) || n < 2 || n > 1000000)
      throw std::invalid_argument("grid requires finite increasing endpoints and 2..1000000 points");
    std::vector<R> result;
    result.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
    {
      R const fraction = uni20::scalar_like(zero, i) / R(n - 1);
      result.push_back((one - fraction) * lo + fraction * hi);
    }
    return result;
  };
  auto const k = grid(a.momentum, zero, 2 * pi, a.momentum_points);
  auto const w =
      grid(a.omega, parse(a.omega_min), a.omega_max ? parse(*a.omega_max) : 2 * sf.maximum_energy(), a.points);
  for (R p : k)
    (void)sf.sheets(p);
  for (R f : w)
    if (!uni20::isfinite(f)) throw std::invalid_argument("frequency must be finite");
  if (k.size() > 1000000 / w.size()) throw std::length_error("spectral grid exceeds 1000000 points");
  std::size_t const count = k.size() * w.size();
  uni20::async::TbbScheduler scheduler(a.threads);
  std::size_t const block_size = uni20::numeric_limits<R>::digits(options.precision()) > 64 ? 1 : 16;
  bethe::ExecutionOptions const execution{
      .scheduler = &scheduler, .max_in_flight = std::max(std::size_t{128 / block_size}, std::size_t(a.threads) * 2)};
  std::size_t failures = 0, thresholds = 0;
  cli::RunReport report(context, "XXZ thermodynamic transverse two-spinon structure factor");
  report.field("hamiltonian", "Hamiltonian", "H=J sum(SxSx+SySy+Delta SzSz); infinite chain; zero field/temperature")
      .field("delta", "Delta", delta)
      .field("exchange", "Exchange J", exchange)
      .field("channel", "Channel", "raising (S-+)")
      .field("family", "Family", "Thermodynamic two-spinon contribution; partial DSF; no rescaling")
      .field("density_convention", "Density convention", "S(k,w); integrals use dw/(2*pi)")
      .field("momentum_convention", "Momentum convention", "Physical k; paper Q=pi-k")
      .field("spinon_gap", "Single-spinon gap", sf.gap())
      .precision(a.precision)
      .field("tolerance", "Log-kernel tail tolerance", options.tolerance)
      .field("max_terms", "Max series terms per sheet", options.max_terms)
      .field("threads", "Scheduler concurrency limit", a.threads)
      .field("max_in_flight", "Maximum points in flight", execution.max_in_flight * block_size)
      .field("block_size", "Points per task", block_size)
      .field("cpu_timing_scope", "CPU timing scope", "spectral grid including ordered row output")
      .field("points", "Spectral points", count);
  namespace data = cli::data;
  using Optional = std::optional<R>;
  using cli::column;
  data::data_table_options const table_options{.retain = a.output.retain ? data::retention::all : data::retention::none,
                                               .metadata = report.metadata()};
  auto spectrum =
      data::make_data_table("Unbroadened transverse spectral density", table_options, column<R>("k"),
                            column<R>("omega"), column<Optional>("density"), column<Optional>("estimated_error"),
                            column<std::size_t>("evaluations"), column<std::string>("status"));
  cli::DataOutput output(a.output, {"spectrum", "continuum"});
  try
  {
    output.attach(spectrum, "spectrum");
    context.measure([&] {
      bethe::evaluate_observable_grid((count + block_size - 1) / block_size,
                                      [&](std::size_t block) {
                                        std::size_t const first = block * block_size,
                                                          size = std::min(block_size, count - first);
                                        std::vector<model::SpectralDensity<R>> values;
                                        values.reserve(size);
                                        for (std::size_t j = 0; j < size; ++j)
                                          values.push_back(sf(k[(first + j) / w.size()], w[(first + j) % w.size()]));
                                        return values;
                                      },
                                      [&](std::size_t block, auto const& values) {
                                        for (std::size_t j = 0; j < values.size(); ++j)
                                        {
                                          auto const i = block * block_size + j;
                                          auto const& value = values[j];
                                          failures += !value.converged();
                                          thresholds += value.status == model::SpectralDensityStatus::threshold;
                                          spectrum.append(k[i / w.size()], w[i % w.size()], value.value, value.error,
                                                          value.evaluations, std::string(status(value.status)));
                                        }
                                      },
                                      execution);
    });
    report.result(failures == 0, failures == 0 ? "converged" : "numerical_failure");
    auto summary = report.finish();
    summary.emplace("Failed points", std::to_string(failures));
    summary.emplace("Sheet endpoint points", std::to_string(thresholds));
    output.overview(report.overview(summary));
    output.finish(spectrum, summary);
    auto continuum =
        data::make_data_table("Two-spinon spectral sheets", table_options, column<R>("k"), column<std::string>("sheet"),
                              column<bool>("present"), column<Optional>("lower"), column<Optional>("upper"));
    output.write_table(
        "continuum", continuum,
        [&](auto& table) {
          for (R p : k)
            for (auto const& sheet : sf.sheets(p))
              table.append(p, std::string(sheet.sign == 1 ? "C_plus" : "C_minus"), sheet.present,
                           sheet.present ? Optional(sheet.lower) : Optional{},
                           sheet.present ? Optional(sheet.upper) : Optional{});
        },
        summary);
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
    output.abort(spectrum, std::move(aborted));
    throw;
  }
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
        return cli::dispatch_options_precision<cli::ArbitraryPrecision::real, model::TwoSpinonStructureFactorOptions>(
            a.precision, [&](auto options) { return run(a, std::move(options), argc, argv); });
      });
}
