// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "result-output.hpp"
#include "solver-precision.hpp"
#include <bethe/xxx_four_spinon_thermodynamics.hpp>
#include <bethe/xxx_six_spinon_structure_factor.hpp>
#include <bethe/xxx_thermodynamic_structure_factor.hpp>

namespace bethe::cli::spinon
{
namespace model = bethe::heisenberg;
struct Controls
{
    std::string sectors = "2", four_cutoff = "16", six_cutoff = "8", relative_tolerance = "1e-6";
    std::string absolute_tolerance = "0";
    std::string six_soft = "0.5", background_weight = "0.5", one_pair_weight = "0.25", two_pairs_weight = "0.25";
    std::string one_pair_scale = "4", one_pair_gap = "0.5", two_pairs_scale = "8", two_pairs_gap = "1";
    std::size_t order = 8, nodes = 100000, terms = 64, tail_order = 16, detail_rows = 100000;
    std::optional<std::size_t> fine_order;
    std::optional<std::string> extended_cutoff;
    bethe::MonteCarloOptions sampling{8, 32, 0, 0, 52};
    double recovery_seconds = 60;
    bool diagnostics = false, replicates = false;
};
inline std::vector<int> sectors(std::string_view text)
{
  std::vector<int> result;
  while (!text.empty())
  {
    auto const comma = text.find(',');
    auto token = text.substr(0, comma);
    int n = token == "2" ? 2 : token == "4" ? 4 : token == "6" ? 6 : 0;
    if (!n || std::ranges::find(result, n) != result.end())
      throw std::invalid_argument("--sectors requires distinct entries from 2,4,6");
    result.push_back(n);
    if (comma == std::string_view::npos) break;
    text.remove_prefix(comma + 1);
    if (text.empty()) throw std::invalid_argument("empty trailing spinon sector");
  }
  if (result.empty()) throw std::invalid_argument("empty spinon sector list");
  std::ranges::sort(result);
  return result;
}
inline void add_options(CLI::App& app, Controls& c)
{
  option(app, "--sectors", c.sectors, "Comma-separated spinon sectors from 2,4,6; partial response")
      ->capture_default_str();
  count_option(app, "--four-order", c.order, "Four-spinon Gauss order per panel, 2..256")->capture_default_str();
  option(app, "--four-cutoff", c.four_cutoff, "Four-spinon log-radius cutoff (not an absolute rapidity bound)")
      ->capture_default_str();
  count_option(app, "--four-max-nodes", c.nodes, "Four-spinon geometry node budget per control")->capture_default_str();
  auto* fine = count_option(app, "--four-refine-order", c.fine_order, "Fine order for four order/cutoff controls");
  auto* extended =
      text_option(app, "--four-extended-cutoff", c.extended_cutoff, "Extended log-radius cutoff for refinement");
  fine->needs(extended);
  extended->needs(fine);
  option(app, "--six-cutoff", c.six_cutoff, "Six-spinon absolute rapidity cutoff |beta|<R")->capture_default_str();
  option(app, "--six-soft-fraction", c.six_soft, "Soft component fraction of the native rapidity marginal")
      ->capture_default_str();
  option(app, "--six-background-weight", c.background_weight, "Background joint-proposal weight; weights sum to 1")
      ->capture_default_str();
  option(app, "--six-one-pair-weight", c.one_pair_weight, "One-soft-pair joint-proposal weight")->capture_default_str();
  option(app, "--six-two-pairs-weight", c.two_pairs_weight, "Two-soft-pair joint-proposal weight")
      ->capture_default_str();
  option(app, "--six-one-pair-scale", c.one_pair_scale, "One-pair center scale")->capture_default_str();
  option(app, "--six-one-pair-gap", c.one_pair_gap, "One-pair rapidity-gap scale")->capture_default_str();
  option(app, "--six-two-pairs-scale", c.two_pairs_scale, "Two-pair center scale")->capture_default_str();
  option(app, "--six-two-pairs-gap", c.two_pairs_gap, "Two-pair rapidity-gap scale")->capture_default_str();
  app.add_option("--uniform-bits", c.sampling.uniform_bits,
                 "Exact midpoint uniform grid bits; separate from arithmetic precision")
      ->capture_default_str();
  count_option(app, "--replicates", c.sampling.replicates, "Independent six-spinon replicates")->capture_default_str();
  count_option(app, "--samples", c.sampling.samples_per_replicate, "Configurations per six-spinon replicate")
      ->capture_default_str();
  option(app, "--seed", c.sampling.seed, "Six-spinon reproducible stream seed")->capture_default_str();
  option(app, "--density-relative-tolerance", c.relative_tolerance, "W4/W6 recovery relative tolerance")
      ->capture_default_str();
  option(app, "--density-absolute-tolerance", c.absolute_tolerance, "W4/W6 recovery absolute tolerance")
      ->capture_default_str();
  count_option(app, "--density-terms", c.terms, "First contour residue budget; guard uses 16 more")
      ->capture_default_str();
  count_option(app, "--density-tail-order", c.tail_order, "First contour tail order; guard uses 4 more")
      ->capture_default_str();
  app.add_option("--density-seconds", c.recovery_seconds, "Wall-time budget per density recovery")
      ->capture_default_str();
  count_option(app, "--detail-row-limit", c.detail_rows, "Maximum buffered diagnostic/replicate rows when requested")
      ->capture_default_str();
  app.add_flag("--diagnostics", c.diagnostics, "Show separate numerical diagnostics on screen");
  app.add_flag("--replicate-output", c.replicates, "Show six-spinon replicate records on screen");
}

template <class R>
using DensityLow = std::conditional_t<uni20::has_runtime_precision_v<R>, R,
                                      std::conditional_t<std::is_same_v<R, double>, float, double>>;
template <uni20::Real R> model::SpinonRecoveryPolicy recovery(Controls const& c, R const& one)
{
  int const guard = uni20::numeric_limits<R>::digits(one);
  int low;
  if constexpr (uni20::has_runtime_precision_v<R>)
    low = std::max(24, guard - 64);
  else
    low = uni20::numeric_limits<DensityLow<R>>::digits(DensityLow<R>{1});
  if (guard <= low)
    throw std::invalid_argument("selected precision has no distinct lower density comparison; choose a wider format");
  if (c.terms > std::size_t(std::numeric_limits<int>::max() - 22) ||
      c.tail_order > std::size_t(std::numeric_limits<int>::max() - 8))
    throw std::invalid_argument("density contour budgets exceed implementation range");
  return {c.relative_tolerance,
          c.absolute_tolerance,
          {{low, guard, {{{int(c.terms), int(c.tail_order)}, {int(c.terms + 16), int(c.tail_order + 4)}}}}},
          c.recovery_seconds};
}

template <uni20::Real R> struct Row
{
    std::string sector, status, reason;
    std::optional<R> value, kernel_error, quadrature, cutoff, contour, arithmetic, sampling;
    std::size_t evaluations = 0;
};
template <uni20::Real R> struct Point
{
    std::vector<Row<R>> rows;
    std::optional<model::SixSpinonSpectrumEstimate<R>> six;
};
template <uni20::Real R> struct Plan
{
    using Four = model::ThermodynamicFourSpinon<R, DensityLow<R>, R>;
    using Six = model::ThermodynamicSixSpinonStructureFactor<R, DensityLow<R>, R>;
    Controls controls;
    std::vector<int> selected;
    model::SpinonRecoveryPolicy policy;
    model::ThermodynamicTwoSpinonStructureFactor<R> two;
    std::optional<Four> four;
    std::optional<Six> six;
    R one, four_cutoff, six_cutoff;
    std::optional<R> extended_cutoff;
    Plan(Controls c, R exchange, model::StructureFactorOptions<R> options)
        : controls(std::move(c)), selected(sectors(controls.sectors)), two(exchange, options),
          one(uni20::scalar_like(exchange, 1)), four_cutoff(parse_solver_real<R>(controls.four_cutoff, options)),
          six_cutoff(parse_solver_real<R>(controls.six_cutoff, options))
    {
      policy = recovery(controls, one);
      auto has = [&](int n) { return std::ranges::find(selected, n) != selected.end(); };
      if (has(4))
      {
        four.emplace(policy, exchange);
        if (controls.extended_cutoff) extended_cutoff = parse_solver_real<R>(*controls.extended_cutoff, options);
        // Validate controls before any destination is opened, without constructing a grid.
        if (!controls.nodes || controls.order < 2 || controls.order > 256 || !uni20::isfinite(four_cutoff) ||
            four_cutoff <= R{0})
          throw std::invalid_argument("invalid four-spinon quadrature controls");
        if (bool(controls.fine_order) != bool(extended_cutoff) ||
            (extended_cutoff && (*controls.fine_order <= controls.order || *controls.fine_order > 256 ||
                                 !uni20::isfinite(*extended_cutoff) || *extended_cutoff <= four_cutoff)))
          throw std::invalid_argument("four-spinon refinement requires increasing orders and cutoffs");
      }
      if (has(6))
      {
        using Proposal = typename Six::Proposal;
        using Kind = model::detail::SpinonProposalKind;
        auto parse = [&](std::string const& text) { return parse_solver_real<R>(text, options); };
        Proposal proposal(six_cutoff, parse(controls.six_soft),
                          {{Kind::background, parse(controls.background_weight), one, one},
                           {Kind::one_pair, parse(controls.one_pair_weight), parse(controls.one_pair_scale),
                            parse(controls.one_pair_gap)},
                           {Kind::two_pairs, parse(controls.two_pairs_weight), parse(controls.two_pairs_scale),
                            parse(controls.two_pairs_gap)}});
        six.emplace(std::move(proposal), policy, exchange);
        bethe::detail::validate_monte_carlo(controls.sampling, 3, one);
        if (controls.sampling.samples_per_replicate >
            std::numeric_limits<std::size_t>::max() / controls.sampling.replicates / 6 / (2 * policy.stages.size()))
          throw std::invalid_argument("six-spinon work counters overflow");
      }
    }
};

template <uni20::Real R> Row<R> four_row(model::FourSpinonIntegralResult<R, 1> const& result)
{
  Row<R> row;
  row.sector = "4";
  row.evaluations = result.density_attempts;
  row.reason = result.reason;
  switch (result.status)
  {
    case model::FourSpinonIntegralStatus::estimated:
      row.status = "estimated";
      break;
    case model::FourSpinonIntegralStatus::physical_zero:
      row.status = "physical_zero";
      break;
    case model::FourSpinonIntegralStatus::threshold:
      row.status = "threshold";
      break;
    case model::FourSpinonIntegralStatus::cancelled:
      row.status = "cancelled";
      break;
    default:
      row.status = "numerical_failure";
  }
  if (result.estimates)
  {
    auto const& v = (*result.estimates)[0];
    row.value = v.value;
    row.contour = v.contour_diagnostic;
    row.arithmetic = v.arithmetic_diagnostic;
  }
  return row;
}

template <uni20::Real R>
uni20::async::AsyncTask point_task(std::shared_ptr<Plan<R> const> plan, R q, R omega, std::size_t index,
                                   bethe::ExecutionOptions execution, uni20::async::WriteBuffer<Point<R>> output)
{
  auto destination = co_await output;
  using F = model::FourSpinonIntegralResult<R, 1>;
  using Refined = model::FourSpinonRefinementResult<R, 1>;
  using S = model::SixSpinonSpectrumEstimate<R>;
  std::optional<uni20::async::Async<F>> four;
  std::optional<uni20::async::Async<Refined>> refined;
  std::optional<uni20::async::Async<S>> six;
  auto const& c = plan->controls;
  Point<R> result;
  bethe::detail::ExecutionControl control(execution.stop);
  execution.stop = control.token();
  std::exception_ptr error;
  auto fail = [&] {
    if (!error) error = std::current_exception();
    control.fail(error);
  };
  try
  {
    if (plan->four)
    {
      if (c.fine_order)
        refined.emplace(
            plan->four
                ->refined_spectrum(q, omega, c.order, *c.fine_order, plan->four_cutoff, *plan->extended_cutoff, c.nodes)
                .evaluate_async(execution));
      else
        four.emplace(plan->four->spectrum(q, omega, c.order, plan->four_cutoff, c.nodes).evaluate_async(execution));
    }
    if (plan->six)
    {
      auto sampling = c.sampling;
      sampling.first_replicate = std::uint64_t(index) * sampling.replicates;
      six.emplace(plan->six->at(q, omega).evaluate_async(sampling, execution));
    }
    if (std::ranges::find(plan->selected, 2) != plan->selected.end())
    {
      auto v = plan->two(q, omega);
      Row<R> row;
      row.sector = "2";
      row.value = v.value;
      row.kernel_error = v.error;
      row.evaluations = v.evaluations;
      row.status = v.status == model::SpectralDensityStatus::converged           ? "converged"
                   : v.status == model::SpectralDensityStatus::outside_continuum ? "physical_zero"
                   : v.status == model::SpectralDensityStatus::lower_threshold   ? "lower_threshold"
                                                                                 : "numerical_failure";
      result.rows.push_back(std::move(row));
    }
  }
  catch (...)
  {
    fail();
  }
  try
  {
    if (four)
    {
      auto const& v = co_await four->read();
      result.rows.push_back(four_row(v.get()));
    }
  }
  catch (...)
  {
    fail();
  }
  try
  {
    if (refined)
    {
      auto const proxy = co_await refined->read();
      auto const& v = proxy.get();
      Row<R> row;
      row.sector = "4";
      if (v.estimates)
      {
        row = four_row(*v.controls[3]);
        row.quadrature = (*v.quadrature_difference)[0];
        row.cutoff = (*v.cutoff_change)[0];
        row.evaluations = 0;
        for (auto const& control : v.controls)
          if (control) row.evaluations += control->density_attempts;
      }
      else
      {
        row.status = v.status == model::FourSpinonIntegralStatus::threshold   ? "threshold"
                     : v.status == model::FourSpinonIntegralStatus::cancelled ? "cancelled"
                                                                              : "numerical_failure";
      }
      row.reason = v.reason;
      result.rows.push_back(std::move(row));
    }
  }
  catch (...)
  {
    fail();
  }
  try
  {
    if (six)
    {
      auto proxy = co_await six->read();
      result.six = proxy.get();
      auto const& v = *result.six;
      Row<R> row;
      row.sector = "6";
      row.value = v.value;
      row.contour = v.contour_diagnostic;
      row.arithmetic = v.arithmetic_diagnostic;
      row.sampling = v.sampling_error;
      row.evaluations = v.density_attempts;
      row.reason = v.reason;
      switch (v.status)
      {
        case model::SixSpinonSpectrumStatus::estimated:
          row.status = "estimated";
          break;
        case model::SixSpinonSpectrumStatus::outside_continuum:
          row.status = "physical_zero";
          break;
        case model::SixSpinonSpectrumStatus::lower_threshold:
          row.status = "lower_threshold";
          break;
        case model::SixSpinonSpectrumStatus::upper_edge:
          row.status = "physical_zero";
          break;
        case model::SixSpinonSpectrumStatus::no_support_sampled:
          row.status = "no_support_sampled";
          break;
        case model::SixSpinonSpectrumStatus::cancelled:
          row.status = "cancelled";
          break;
        default:
          row.status = "numerical_failure";
      }
      result.rows.push_back(std::move(row));
    }
  }
  catch (...)
  {
    fail();
  }
  if (error) std::rethrow_exception(error);
  destination = std::move(result);
}

template <uni20::Real R>
auto point_async(std::shared_ptr<Plan<R> const> plan, R q, R omega, std::size_t index,
                 bethe::ExecutionOptions execution)
{
  uni20::async::Async<Point<R>> output;
  auto task = point_task(std::move(plan), std::move(q), std::move(omega), index, execution, output.write());
  bethe::detail::schedule_owned(execution.continuous_scheduler(), std::move(task));
  return output;
}
} // namespace bethe::cli::spinon

namespace bethe::cli::spinon
{
template <uni20::Real R, class Arguments, class Info>
int run_spectrum(Arguments const& a, model::StructureFactorOptions<R> options, Info info, int argc, char** argv)
{
  R const one = uni20::scalar_at_precision<R>(1, options.precision()), zero = one - one,
          pi = R{4} * uni20::math::atan(one);
  R const exchange = parse_solver_real<R>(a.exchange, options), scale = a.channel == "zz" ? one : R{2} * one;
  options.max_evaluations = a.evaluations;
  if (a.max_order) options.max_order = *a.max_order;
  if (a.tolerance) options.tolerance = parse_solver_real<R>(*a.tolerance, options);
  std::shared_ptr<Plan<R> const> plan = std::make_shared<Plan<R>>(a.spinons, exchange, options);
  auto grid = [&](std::optional<std::string> const& single, R lo, R hi, std::size_t n) {
    if (single) return std::vector<R>{parse_solver_real<R>(*single, options)};
    if (!uni20::isfinite(lo) || !uni20::isfinite(hi) || !(hi > lo) || n < 2 || n > 1000000)
      throw std::invalid_argument("invalid grid: finite increasing endpoints and 2..1000000 points required");
    std::vector<R> values;
    values.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
    {
      R f = uni20::scalar_like(one, i) / uni20::scalar_like(one, n - 1);
      values.push_back((one - f) * lo + f * hi);
    }
    return values;
  };
  auto const q = grid(a.momentum, zero, R{2} * pi, a.momentum_points);
  R const upper_default = uni20::scalar_like(one, plan->selected.back() / 2) * pi * exchange;
  auto const w = grid(a.omega, parse_solver_real<R>(a.omega_min, options),
                      a.omega_max ? parse_solver_real<R>(*a.omega_max, options) : upper_default, a.points);
  for (auto const& x : q)
    (void)plan->two.boundaries(x);
  for (auto const& x : w)
    if (!uni20::isfinite(x)) throw std::invalid_argument("frequency must be finite");
  if (q.size() > 1000000 / w.size()) throw std::length_error("spectral grid exceeds 1000000 points");
  std::size_t const count = q.size() * w.size();
  if (plan->six && a.spinons.sampling.replicates > std::numeric_limits<std::uint64_t>::max() / count)
    throw std::invalid_argument("six-spinon replicate identities overflow on this grid");
  bool const want_diagnostics = a.output.requested("diagnostics") || a.output.screen_table("diagnostics");
  bool const want_replicates = plan->six && (a.output.requested("replicates") || a.output.screen_table("replicates"));
  std::size_t const per_point =
      (want_diagnostics ? plan->selected.size() : 0) + (want_replicates ? a.spinons.sampling.replicates : 0);
  if (per_point && count > a.spinons.detail_rows / per_point)
    throw std::invalid_argument(
        "requested diagnostic/replicate tables exceed --detail-row-limit; split the grid or raise the limit");
  std::vector<std::tuple<R, R, Row<R>>> diagnostic_rows;
  std::vector<std::tuple<R, R, model::SpinonSpectrumReplicate<R>>> replicate_rows;
  uni20::async::TbbScheduler scheduler(a.threads);
  bethe::ExecutionOptions execution{.scheduler = &scheduler,
                                    .max_in_flight = std::max(std::size_t{2}, std::size_t(a.threads) * 2)};
  uni20::run_context context(info, {.invocation = std::vector<std::string>(argv, argv + argc)});
  RunReport report(context, "XXX thermodynamic selected spinon contributions");
  report.field("hamiltonian", "Hamiltonian", "H=J sum S_j.S_(j+1); infinite chain; zero field and temperature")
      .field("exchange", "Exchange J", exchange)
      .field("channel", "Channel", a.channel)
      .field("sectors", "Selected spinon sectors", a.spinons.sectors)
      .field("family", "Family", "Partial DSF; physical normalization; no fitted rescaling")
      .field("density_convention", "Density convention", "S(q,w); moments use dw/(2*pi); raising=2*zz")
      .precision(a.precision)
      .field("geometry_bits", "Geometry and output bits", uni20::numeric_limits<R>::digits(one))
      .field("density_bits", "Density comparison bits", plan->policy.stages[0].bits)
      .field("density_guard_bits", "Density guard bits", plan->policy.stages[0].guard_bits)
      .field("two_tolerance", "Two-spinon transition-integral tolerance", options.tolerance)
      .field("two_max_order", "Two-spinon maximum quadrature order", options.max_order)
      .field("two_max_evaluations", "Two-spinon kernel budget", options.max_evaluations)
      .field("density_relative_tolerance", "Density relative tolerance", a.spinons.relative_tolerance)
      .field("density_absolute_tolerance", "Density absolute tolerance", a.spinons.absolute_tolerance)
      .field("density_terms", "First density residue budget", a.spinons.terms)
      .field("density_guard_terms", "Second density residue budget", a.spinons.terms + 16)
      .field("density_tail_order", "First density tail order", a.spinons.tail_order)
      .field("density_guard_tail_order", "Second density tail order", a.spinons.tail_order + 4)
      .field("density_seconds", "Recovery seconds per configuration", a.spinons.recovery_seconds)
      .field("threads", "Scheduler concurrency limit", a.threads)
      .field("points", "Spectral points", count)
      .field("detail_row_limit", "Maximum buffered detail rows", a.spinons.detail_rows)
      .field("cpu_timing_scope", "CPU timing scope", "spectral grid including geometry and ordered output");
  if (plan->four)
  {
    report.field("four_order", "Four-spinon order", a.spinons.order)
        .field("four_log_radius_cutoff", "Four-spinon log-radius cutoff", plan->four_cutoff)
        .field("four_max_nodes", "Four-spinon nodes per control", a.spinons.nodes);
    if (plan->extended_cutoff)
      report.field("four_extended_cutoff", "Four-spinon extended log-radius cutoff", *plan->extended_cutoff)
          .field("four_fine_order", "Four-spinon fine order", *a.spinons.fine_order);
  }
  if (plan->six)
    report.field("six_rapidity_cutoff", "Six-spinon absolute rapidity cutoff", plan->six_cutoff)
        .field("six_replicates", "Six-spinon replicates", a.spinons.sampling.replicates)
        .field("six_samples", "Samples per replicate", a.spinons.sampling.samples_per_replicate)
        .field("six_seed", "Six-spinon seed", a.spinons.sampling.seed)
        .field("six_streams", "Six-spinon stream identities",
               "first_replicate=grid_index*replicates; distinct streams at each point")
        .field("six_uniform_bits", "Six-spinon uniform grid bits", a.spinons.sampling.uniform_bits)
        .field("six_soft_fraction", "Six-spinon native soft fraction", a.spinons.six_soft)
        .field("six_background_weight", "Six-spinon background weight", a.spinons.background_weight)
        .field("six_one_pair_weight", "Six-spinon one-pair weight", a.spinons.one_pair_weight)
        .field("six_two_pairs_weight", "Six-spinon two-pair weight", a.spinons.two_pairs_weight)
        .field("six_one_pair_scale", "Six-spinon one-pair center scale", a.spinons.one_pair_scale)
        .field("six_one_pair_gap", "Six-spinon one-pair gap scale", a.spinons.one_pair_gap)
        .field("six_two_pairs_scale", "Six-spinon two-pair center scale", a.spinons.two_pairs_scale)
        .field("six_two_pairs_gap", "Six-spinon two-pair gap scale", a.spinons.two_pairs_gap);
  using Optional = std::optional<R>;
  data::data_table_options table_options{.retain = a.output.retain ? data::retention::all : data::retention::none,
                                         .metadata = report.metadata()};
  auto spectrum = data::make_data_table("Selected spinon spectral densities", table_options, column<R>("q"),
                                        column<R>("omega"), column<std::string>("sector"), column<Optional>("density"),
                                        column<Optional>("estimated_error"), column<std::size_t>("evaluations"),
                                        column<std::string>("status"));
  auto diagnostics = data::make_data_table(
      "Separate numerical diagnostics", table_options, column<R>("q"), column<R>("omega"),
      column<std::string>("sector"), column<Optional>("two_kernel_error"), column<Optional>("quadrature_difference"),
      column<Optional>("cutoff_change"), column<Optional>("contour_diagnostic"),
      column<Optional>("arithmetic_diagnostic"), column<Optional>("sampling_error"), column<std::size_t>("evaluations"),
      column<std::string>("reason"));
  auto replicates =
      data::make_data_table("Six-spinon replicate records", table_options, column<R>("q"), column<R>("omega"),
                            column<std::uint64_t>("replicate"), column<std::size_t>("completed_samples"),
                            column<std::size_t>("supported_roots"), column<std::size_t>("density_attempts"),
                            column<Optional>("maximum_absolute_density"), column<Optional>("density"),
                            column<Optional>("contour_diagnostic"), column<Optional>("arithmetic_diagnostic"),
                            column<std::optional<std::size_t>>("failed_sample"), column<std::string>("reason"));
  std::vector<std::string> tables{"spectrum", "continuum", "diagnostics"};
  if (plan->six) tables.push_back("replicates");
  DataOutput output(a.output, tables);
  std::size_t failures = 0, thresholds = 0;
  try
  {
    output.attach(spectrum, "spectrum");
    context.measure([&] {
      for (std::size_t i = 0; i < count; ++i)
      {
        R const& momentum = q[i / w.size()];
        R const& frequency = w[i % w.size()];
        auto pending = point_async(plan, momentum, frequency, i, execution);
        auto value = pending.get_wait(scheduler);
        bool failed = false, threshold = false;
        bethe::detail::CompensatedSum<R> total(zero);
        std::string label = "sum:";
        for (auto& row : value.rows)
        {
          for (auto* field : {&row.value, &row.kernel_error, &row.quadrature, &row.cutoff, &row.contour,
                              &row.arithmetic, &row.sampling})
            if (*field)
            {
              **field *= scale;
              if (!uni20::isfinite(**field))
              {
                row.value.reset();
                row.status = "numerical_failure";
                row.reason = "channel conversion exceeds exponent range";
                break;
              }
            }
          if (row.status == "numerical_failure" && !row.value)
          {
            row.kernel_error.reset();
            row.quadrature.reset();
            row.cutoff.reset();
            row.contour.reset();
            row.arithmetic.reset();
            row.sampling.reset();
          }
          bool const edge = row.status == "threshold" || row.status == "lower_threshold";
          threshold |= edge;
          failed |= !row.value && !edge;
          if (row.value) total.add(*row.value);
          label += (label == "sum:" ? "" : "+") + row.sector;
          spectrum.append(momentum, frequency, row.sector, row.value, row.kernel_error, row.evaluations, row.status);
          if (want_diagnostics) diagnostic_rows.emplace_back(momentum, frequency, row);
        }
        if (value.rows.size() > 1)
        {
          Optional sum;
          if (!failed && !threshold)
          {
            sum = total.value();
            if (!uni20::isfinite(*sum))
            {
              sum.reset();
              failed = true;
            }
          }
          spectrum.append(momentum, frequency, label, sum, Optional{}, std::size_t{0},
                          std::string(failed      ? "incomplete"
                                      : threshold ? "threshold"
                                                  : "partial_sum"));
        }
        failures += failed;
        thresholds += threshold;
        if (want_replicates && value.six)
          for (auto const& row : value.six->replicates)
            if (row) replicate_rows.emplace_back(momentum, frequency, *row);
      }
    });
    report.result(failures == 0, failures ? "incomplete" : "finite_domain_estimate");
    auto summary = report.finish();
    summary.emplace("Failed points", std::to_string(failures));
    summary.emplace("Threshold points", std::to_string(thresholds));
    output.overview(report.overview(summary));
    output.finish(spectrum, summary);
    output.write_table(
        "diagnostics", diagnostics,
        [&](auto& table) {
          for (auto const& [momentum, frequency, row] : diagnostic_rows)
            table.append(momentum, frequency, row.sector, row.kernel_error, row.quadrature, row.cutoff, row.contour,
                         row.arithmetic, row.sampling, row.evaluations, row.reason);
        },
        summary);
    diagnostic_rows.clear();
    if (plan->six)
      output.write_table(
          "replicates", replicates,
          [&](auto& table) {
            for (auto const& [momentum, frequency, row] : replicate_rows)
              table.append(momentum, frequency, row.id, row.completed_samples, row.supported_roots,
                           row.density_attempts,
                           row.maximum_absolute_observation ? Optional(scale * (*row.maximum_absolute_observation)[0])
                                                            : Optional{},
                           row.means ? Optional(scale * (*row.means)[0]) : Optional{},
                           row.means ? Optional(scale * (*row.means)[1]) : Optional{},
                           row.means ? Optional(scale * (*row.means)[2]) : Optional{}, row.failed_sample, row.reason);
          },
          summary);
    replicate_rows.clear();
    auto continuum = data::make_data_table("Spinon kinematic continua", table_options, column<R>("q"),
                                           column<int>("sector"), column<R>("lower"), column<R>("upper"));
    output.write_table(
        "continuum", continuum,
        [&](auto& table) {
          for (auto const& momentum : q)
            for (int sector : plan->selected)
            {
              auto bounds = plan->two.boundaries(momentum);
              R folded = std::min(momentum, R{2} * pi - momentum);
              if (sector == 4) bounds.second = R{2} * pi * exchange * uni20::math::cos(folded / R{4});
              if (sector == 6) bounds.second = R{3} * pi * exchange * uni20::math::cos((pi - folded) / R{6});
              table.append(momentum, sector, bounds.first, bounds.second);
            }
        },
        summary);
    output.finish_document();
  }
  catch (...)
  {
    data::table_metadata aborted{{"Status", "aborted"}};
    if (!context.finished()) try
      {
        aborted = run_summary(context.finish(uni20::run_outcome::failed));
      }
      catch (...)
      {}
    output.abort(spectrum, aborted);
    output.abort(diagnostics, aborted);
    if (plan->six) output.abort(replicates, aborted);
    throw;
  }
  return failures ? 2 : 0;
}
} // namespace bethe::cli::spinon
