// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "execution-options.hpp"
#include "spinon-thermo.hpp"
#include <bethe/xxx_six_spinon_moments.hpp>
#include <bethe/xxx_two_spinon_moments.hpp>

namespace
{
namespace cli = bethe::cli;
namespace spinon = cli::spinon;
namespace model = bethe::heisenberg;
struct Arguments
{
    bool global = false, covariance = false;
    std::optional<std::string> momentum, tolerance;
    std::string exchange = "1", channel = "zz", precision = "fp64", two_cutoff = "48";
    std::size_t two_order = 16, evaluations = 200000;
    std::optional<std::size_t> max_order;
    std::vector<std::string> gaps{"16", "8", "16"}, extended_gaps;
    int threads = 1;
    cli::DataOutputOptions output;
    spinon::Controls controls;
};
auto program_info()
{
  auto info = cli::program_info("bethe-xxx-moments-thermo", "Direct thermodynamic XXX spinon zeroth and first moments.",
                                bethe::citations::Tool::xxx_moments_thermo);
  info.examples = {{"bethe-xxx-moments-thermo --global", "Two-spinon global moments"},
                   {"bethe-xxx-moments-thermo --momentum 1.5 --sectors 2,4 --precision bits-192 --threads 2",
                    "Selected fixed-momentum contributions"}};
  info.notes = {
      "Infinite spin-1/2 XXX chain, J>0, zero field and temperature; choose --global or --momentum.",
      "Moments include dw/(2*pi); global moments also include dk/(2*pi). No display-grid integration or rescaling.",
      "Default sector 2. --sectors 2,4,6 requests a partial response with physically normalized contributions.",
      "Two-spinons use a relative-rapidity gap cutoff; four-spinons use log radius (fixed k) or three gaps (global).",
      "Six-spinons use an absolute-rapidity box. Equal cutoff numbers do not define equal domains.",
      "Finite-cutoff estimates exclude omitted tails. Vary quadrature, sampling, precision and cutoff separately.",
      "At k=pi the unrestricted static response diverges. At k=0 the first-moment fraction is undefined.",
      "The raising channel is twice zz; values, diagnostics and exact reference denominators are converted together.",
      "For sectors 4/6 prefer --precision bits-192/bits-384 respectively. See --references and "
      "docs/xxx-spinon-moments.md."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  auto* global = app.add_flag("--global", a.global, "Integrate over frequency and momentum");
  auto* q = cli::text_option(app, "--momentum", a.momentum, "Fixed momentum in radians/site, [0,2*pi]");
  global->excludes(q);
  cli::option(app, "--exchange", a.exchange, "Positive exchange J")->capture_default_str();
  cli::option(app, "--channel", a.channel, "zz or raising (S-+)")
      ->check(CLI::IsMember({"zz", "raising"}))
      ->capture_default_str();
  cli::precision_option(app, a.precision, cli::ArbitraryPrecision::real);
  cli::threads_option(app, a.threads)->description("Shared scheduler concurrency limit");
  cli::count_option(app, "--two-order", a.two_order, "Two-spinon Gauss order per gap panel, 2..256")
      ->capture_default_str();
  cli::option(app, "--two-gap-cutoff", a.two_cutoff, "Two-spinon relative-rapidity gap cutoff")->capture_default_str();
  cli::text_option(app, "--tolerance", a.tolerance, "Two-spinon transition-kernel tolerance (default 4096 epsilon)");
  cli::count_option(app, "--max-evaluations", a.evaluations, "Two-spinon transition-kernel evaluation budget")
      ->capture_default_str();
  cli::count_option(app, "--max-quadrature-order", a.max_order,
                    "Two-spinon kernel Gauss order limit, power of two 32..2048");
  spinon::add_options(app, a.controls, true);
  app.add_option("--four-gap-cutoffs", a.gaps, "Three global four-spinon gap cutoffs, A,B,C")
      ->delimiter(',')
      ->expected(3)
      ->capture_default_str();
  app.add_option("--four-extended-gap-cutoffs", a.extended_gaps,
                 "Extended global gap cutoffs paired with --four-refine-order")
      ->delimiter(',')
      ->expected(3);
  app.add_flag("--covariance", a.covariance, "Show covariance of the six-spinon moment means");
  cli::add_data_output_options(
      app, a.output,
      {{.name = "moments", .description = "Requested sector moments and their partial sum", .primary = true},
       {.name = "diagnostics", .description = "Separate numerical diagnostics", .screen_option = "--diagnostics"},
       {.name = "replicates",
        .description = "Six-spinon replicate means and maxima",
        .availability = "sector 6 selected",
        .screen_option = "--replicate-output"},
       {.name = "covariance",
        .description = "Covariance of the six-spinon moment means",
        .availability = "sector 6 selected",
        .screen_option = "--covariance"}});
}

template <class R> struct Row
{
    std::string sector, status, reason;
    std::array<std::optional<R>, 2> values, denominators, fractions, kernel, quadrature, cutoff, contour, arithmetic,
        shift_quadrature, shift_tail, sampling;
    std::size_t evaluations = 0;
};
template <class R> struct Result
{
    std::vector<Row<R>> rows;
    std::optional<bethe::MonteCarloEstimate<R, model::SpinonSampleMetadata>> six;
};
template <class R> struct Prepared
{
    using Four = model::ThermodynamicFourSpinon<R, spinon::DensityLow<R>, R>;
    using Six = model::ThermodynamicSixSpinonMoments<R, spinon::DensityLow<R>, R>;
    R one, exchange, two_cutoff, four_cutoff, six_cutoff;
    std::optional<R> momentum;
    std::array<R, 3> gaps;
    std::optional<std::array<R, 3>> extended_gaps;
    std::optional<R> extended_cutoff;
    std::vector<int> selected;
    std::optional<model::SpinonRecoveryPolicy> policy;
    std::optional<typename model::ThermodynamicTwoSpinonMoments<R>::Query> two;
    std::optional<typename Four::template Query<2>> four;
    std::optional<typename Four::template RefinedQuery<2>> refined;
    std::optional<Six> six;
    std::optional<typename Six::FixedQuery> fixed_six;
    bethe::MonteCarloOptions sampling;
    Prepared(Arguments const& a, model::StructureFactorOptions<R> options)
        : one(uni20::scalar_at_precision<R>(1, options.precision())),
          exchange(cli::parse_solver_real<R>(a.exchange, options)),
          two_cutoff(cli::parse_solver_real<R>(a.two_cutoff, options)),
          four_cutoff(cli::parse_solver_real<R>(a.controls.four_cutoff, options)),
          six_cutoff(cli::parse_solver_real<R>(a.controls.six_cutoff, options)),
          selected(spinon::sectors(a.controls.sectors)), sampling(a.controls.sampling)
    {
      auto parse = [&](std::string const& text) { return cli::parse_solver_real<R>(text, options); };
      if (a.momentum) momentum = parse(*a.momentum);
      options.max_evaluations = a.evaluations;
      if (a.max_order) options.max_order = *a.max_order;
      if (a.tolerance) options.tolerance = parse(*a.tolerance);
      model::ThermodynamicTwoSpinonMoments<R> reference(exchange, options);
      auto has = [&](int n) { return std::ranges::find(selected, n) != selected.end(); };
      if (has(2))
        two.emplace(momentum ? reference.fixed_moments(*momentum, a.two_order, two_cutoff)
                             : reference.global_moments(a.two_order, two_cutoff));
      if (has(4) || has(6)) policy = spinon::recovery(a.controls, one);
      if (has(4))
      {
        Four model(*policy, exchange);
        if (a.gaps.size() != 3) throw std::invalid_argument("require three global four-spinon gap cutoffs");
        for (int i = 0; i < 3; ++i)
          gaps[i] = parse(a.gaps[i]);
        if (!a.extended_gaps.empty())
        {
          std::array<R, 3> parsed;
          for (int i = 0; i < 3; ++i)
            parsed[i] = parse(a.extended_gaps.at(i));
          extended_gaps = parsed;
        }
        if (a.controls.extended_cutoff) extended_cutoff = parse(*a.controls.extended_cutoff);
        if (a.controls.fine_order)
        {
          if (momentum)
            refined.emplace(model.refined_fixed_moments(*momentum, a.controls.order, *a.controls.fine_order,
                                                        four_cutoff, *extended_cutoff, a.controls.nodes));
          else
            refined.emplace(model.refined_global_moments(a.controls.order, *a.controls.fine_order, gaps, *extended_gaps,
                                                         a.controls.nodes));
        }
        else
        {
          if (momentum)
            four.emplace(model.fixed_moments(*momentum, a.controls.order, four_cutoff, a.controls.nodes));
          else
            four.emplace(model.global_moments(a.controls.order, gaps, a.controls.nodes));
        }
      }
      if (has(6))
      {
        six.emplace(spinon::joint_proposal(a.controls, one, options), *policy, exchange);
        bethe::detail::validate_monte_carlo(sampling, 6, one);
        if (sampling.samples_per_replicate >
            std::numeric_limits<std::size_t>::max() / sampling.replicates / (2 * policy->stages.size()))
          throw std::invalid_argument("six-spinon work counters overflow");
        if (momentum) fixed_six.emplace(six->fixed_moments(*momentum));
      }
    }
};
template <class R> Row<R> four_row(model::FourSpinonIntegralResult<R, 2> const& r)
{
  Row<R> row;
  row.sector = "4";
  row.reason = r.reason;
  row.evaluations = r.density_attempts;
  row.status = r.status == model::FourSpinonIntegralStatus::estimated       ? "estimated"
               : r.status == model::FourSpinonIntegralStatus::physical_zero ? "physical_zero"
               : r.status == model::FourSpinonIntegralStatus::cancelled     ? "cancelled"
                                                                            : "numerical_failure";
  row.denominators = r.sum_rule_denominators;
  row.fractions = r.sum_rule_fractions;
  if (r.estimates)
    for (int i = 0; i < 2; ++i)
    {
      auto const& v = (*r.estimates)[i];
      row.values[i] = v.value;
      row.contour[i] = v.contour_diagnostic;
      row.arithmetic[i] = v.arithmetic_diagnostic;
      row.shift_quadrature[i] = v.shift_quadrature_diagnostic;
      row.shift_tail[i] = v.shift_tail_bound;
    }
  return row;
}
template <class R>
void six_row(Result<R>& out, bethe::MonteCarloEstimate<R, model::SpinonSampleMetadata> integral,
             std::array<std::optional<R>, 2> denominators, std::string status, std::string reason)
{
  Row<R> row;
  row.sector = "6";
  row.status = std::move(status);
  row.reason = std::move(reason);
  row.denominators = std::move(denominators);
  row.evaluations = integral.metadata.density_attempts;
  if (integral.mean)
    for (int i = 0; i < 2; ++i)
    {
      row.values[i] = (*integral.mean)[i];
      row.contour[i] = (*integral.mean)[i + 2];
      row.arithmetic[i] = (*integral.mean)[i + 4];
      if (row.denominators[i])
        row.fractions[i] =
            model::detail::spinon_weighted_value(*row.values[i], -uni20::math::log(*row.denominators[i]));
      if (integral.covariance) row.sampling[i] = uni20::math::sqrt((*integral.covariance)[7 * i]);
    }
  out.rows.push_back(std::move(row));
  out.six = std::move(integral);
}
template <class R>
uni20::async::AsyncTask calculate(std::shared_ptr<Prepared<R> const> p, bethe::ExecutionOptions execution,
                                  uni20::async::WriteBuffer<Result<R>> output)
{
  auto destination = co_await output;
  std::optional<uni20::async::Async<model::TwoSpinonMoments<R>>> two;
  std::optional<uni20::async::Async<model::FourSpinonIntegralResult<R, 2>>> four;
  std::optional<uni20::async::Async<model::FourSpinonRefinementResult<R, 2>>> refined;
  std::optional<uni20::async::Async<model::SixSpinonFixedMoments<R>>> fixed_six;
  std::optional<uni20::async::Async<model::SixSpinonGlobalMoments<R>>> global_six;
  Result<R> result;
  bethe::detail::ExecutionControl control(execution.stop);
  execution.stop = control.token();
  std::exception_ptr error;
  auto fail = [&] {
    if (!error) error = std::current_exception();
    control.fail(error);
  };
  try
  {
    if (p->two) two.emplace(p->two->evaluate_async(execution));
    if (p->four) four.emplace(p->four->evaluate_async(execution));
    if (p->refined) refined.emplace(p->refined->evaluate_async(execution));
    if (p->fixed_six)
      fixed_six.emplace(p->fixed_six->evaluate_async(p->sampling, execution));
    else if (p->six)
      global_six.emplace(p->six->global_async(p->sampling, execution));
  }
  catch (...)
  {
    fail();
  }
  try
  {
    if (two)
    {
      auto proxy = co_await two->read();
      auto const& r = proxy.get();
      Row<R> row;
      row.sector = "2";
      row.reason = r.integral.reason;
      row.evaluations = r.integral.work_units;
      row.status = r.status == model::TwoSpinonMomentStatus::estimated       ? "estimated"
                   : r.status == model::TwoSpinonMomentStatus::physical_zero ? "physical_zero"
                   : r.status == model::TwoSpinonMomentStatus::cancelled     ? "cancelled"
                                                                             : "numerical_failure";
      row.denominators = r.sum_rule_denominators;
      row.fractions = r.sum_rule_fractions;
      if (r.integral.sums)
        for (int i = 0; i < 2; ++i)
        {
          row.values[i] = (*r.integral.sums)[i];
          row.kernel[i] = (*r.integral.sums)[i + 2];
        }
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
      auto proxy = co_await four->read();
      result.rows.push_back(four_row(proxy.get()));
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
      auto proxy = co_await refined->read();
      auto const& r = proxy.get();
      Row<R> row;
      row.sector = "4";
      row.status = "numerical_failure";
      row.reason = r.reason;
      if (r.estimates)
      {
        row = four_row(*r.controls[3]);
        row.quadrature = {(*r.quadrature_difference)[0], (*r.quadrature_difference)[1]};
        row.cutoff = {(*r.cutoff_change)[0], (*r.cutoff_change)[1]};
        row.evaluations = 0;
        for (auto const& c : r.controls)
          if (c) row.evaluations += c->density_attempts;
      }
      else if (r.status == model::FourSpinonIntegralStatus::cancelled)
        row.status = "cancelled";
      result.rows.push_back(std::move(row));
    }
  }
  catch (...)
  {
    fail();
  }
  try
  {
    if (fixed_six)
    {
      auto proxy = co_await fixed_six->read();
      auto const& r = proxy.get();
      std::string status = r.status == model::SixSpinonMomentStatus::estimated            ? "estimated"
                           : r.status == model::SixSpinonMomentStatus::physical_zero      ? "physical_zero"
                           : r.status == model::SixSpinonMomentStatus::no_support_sampled ? "no_support_sampled"
                           : r.status == model::SixSpinonMomentStatus::cancelled          ? "cancelled"
                                                                                          : "numerical_failure";
      six_row(result, r.integral, {std::optional<R>{}, r.first_moment_denominator}, std::move(status), r.reason);
    }
  }
  catch (...)
  {
    fail();
  }
  try
  {
    if (global_six)
    {
      auto proxy = co_await global_six->read();
      auto const& r = proxy.get();
      R const exact = R{2} * p->exchange / R{3} * (uni20::math::log(R{2} * p->one) - p->one / R{4});
      std::string status = r.integral.status == bethe::MonteCarloStatus::complete    ? "estimated"
                           : r.integral.status == bethe::MonteCarloStatus::cancelled ? "cancelled"
                                                                                     : "numerical_failure";
      six_row(result, r.integral, {std::optional<R>(p->one / R{4}), std::optional<R>(exact)}, std::move(status),
              r.integral.reason);
    }
  }
  catch (...)
  {
    fail();
  }
  if (error) std::rethrow_exception(error);
  destination = std::move(result);
}

template <class R> int run(Arguments const& a, model::StructureFactorOptions<R> options, int argc, char** argv)
{
  std::shared_ptr<Prepared<R> const> p = std::make_shared<Prepared<R>>(a, options);
  R const scale = a.channel == "zz" ? p->one : R{2} * p->one;
  bool const want_replicates = a.output.requested("replicates") || a.output.screen_table("replicates");
  if (want_replicates && p->sampling.replicates > a.controls.detail_rows)
    throw std::invalid_argument("replicate records exceed --detail-row-limit");
  std::vector<std::string> tables{"moments", "diagnostics"};
  if (p->six)
  {
    tables.push_back("replicates");
    tables.push_back("covariance");
  }
  cli::DataOutput output(a.output, tables);
  uni20::async::TbbScheduler scheduler(a.threads);
  bethe::ExecutionOptions execution{.scheduler = &scheduler,
                                    .max_in_flight = std::max(std::size_t{2}, std::size_t(a.threads) * 2)};
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  cli::RunReport report(context, "XXX thermodynamic spinon moments");
  report.field("hamiltonian", "Hamiltonian", "H=J sum S_j.S_(j+1); infinite chain; zero field and temperature")
      .field("exchange", "Exchange J", p->exchange)
      .field("channel", "Channel", a.channel)
      .field("sectors", "Selected sectors", a.controls.sectors)
      .field("domain", "Moment domain", a.global ? "global" : "fixed momentum")
      .field("momentum", "Momentum", p->momentum)
      .field("normalization", "Normalization",
             "Physical Szz; dw/(2*pi), and dk/(2*pi) for global moments; raising=2*zz")
      .field("family", "Family", "Finite-domain partial response; no fitted rescaling")
      .precision(a.precision)
      .field("threads", "Scheduler concurrency limit", a.threads)
      .field("two_order", "Two-spinon gap order", a.two_order)
      .field("two_gap_cutoff", "Two-spinon gap cutoff", p->two_cutoff)
      .field("two_kernel_tolerance", "Two-spinon kernel tolerance",
             a.tolerance ? cli::parse_solver_real<R>(*a.tolerance, options) : options.tolerance)
      .field("two_kernel_max_order", "Two-spinon kernel order limit", a.max_order.value_or(options.max_order))
      .field("two_kernel_max_evaluations", "Two-spinon kernel budget", a.evaluations)
      .field("four_order", "Four-spinon quadrature order", a.controls.order)
      .field("four_max_nodes", "Four-spinon node budget", a.controls.nodes);
  if (p->four || p->refined)
  {
    if (a.global)
      for (int i = 0; i < 3; ++i)
        report.field("four_gap_" + std::to_string(i), "Four-spinon gap cutoff " + std::to_string(i), p->gaps[i]);
    else
      report.field("four_log_radius_cutoff", "Four-spinon log-radius cutoff", p->four_cutoff);
    if (a.controls.fine_order)
    {
      report.field("four_fine_order", "Four-spinon fine order", *a.controls.fine_order);
      if (p->extended_gaps)
        for (int i = 0; i < 3; ++i)
          report.field("four_extended_gap_" + std::to_string(i), "Extended gap cutoff " + std::to_string(i),
                       (*p->extended_gaps)[i]);
      if (p->extended_cutoff) report.field("four_extended_cutoff", "Extended log-radius cutoff", *p->extended_cutoff);
    }
  }
  if (p->policy)
    report.field("density_bits", "Density comparison bits", p->policy->stages[0].bits)
        .field("density_guard_bits", "Density guard bits", p->policy->stages[0].guard_bits)
        .field("density_relative_tolerance", "Density relative tolerance", a.controls.relative_tolerance)
        .field("density_absolute_tolerance", "Density absolute tolerance", a.controls.absolute_tolerance)
        .field("density_terms", "First density terms", a.controls.terms)
        .field("density_guard_terms", "Second density terms", a.controls.terms + 16)
        .field("density_tail_order", "First density tail order", a.controls.tail_order)
        .field("density_guard_tail_order", "Second density tail order", a.controls.tail_order + 4)
        .field("density_seconds", "Recovery seconds per configuration", a.controls.recovery_seconds);
  if (p->six)
    report.field("six_rapidity_cutoff", "Six-spinon absolute rapidity cutoff", p->six_cutoff)
        .field("six_replicates", "Six-spinon replicates", p->sampling.replicates)
        .field("six_samples", "Samples per replicate", p->sampling.samples_per_replicate)
        .field("six_seed", "Six-spinon seed", p->sampling.seed)
        .field("six_uniform_bits", "Uniform grid bits", p->sampling.uniform_bits)
        .field("six_soft_fraction", "Native soft fraction", a.controls.six_soft)
        .field("six_background_weight", "Background weight", a.controls.background_weight)
        .field("six_one_pair_weight", "One-pair weight", a.controls.one_pair_weight)
        .field("six_two_pairs_weight", "Two-pair weight", a.controls.two_pairs_weight)
        .field("six_one_pair_scale", "One-pair center scale", a.controls.one_pair_scale)
        .field("six_one_pair_gap", "One-pair gap scale", a.controls.one_pair_gap)
        .field("six_two_pairs_scale", "Two-pair center scale", a.controls.two_pairs_scale)
        .field("six_two_pairs_gap", "Two-pair gap scale", a.controls.two_pairs_gap);
  using O = std::optional<R>;
  namespace data = cli::data;
  data::data_table_options table_options{.retain = a.output.retain ? data::retention::all : data::retention::none,
                                         .metadata = report.metadata()};
  auto moments = data::make_data_table(
      "Selected spinon moments", table_options, cli::column<O>("q"), cli::column<std::string>("sector"),
      cli::column<int>("moment"), cli::column<O>("value"), cli::column<O>("sum_rule"), cli::column<O>("fraction"),
      cli::column<std::string>("status"), cli::column<bool>("unrestricted_static_divergence"));
  auto diagnostics = data::make_data_table(
      "Separate moment diagnostics", table_options, cli::column<std::string>("sector"), cli::column<int>("moment"),
      cli::column<O>("two_kernel_error"), cli::column<O>("quadrature_difference"), cli::column<O>("cutoff_change"),
      cli::column<O>("contour_diagnostic"), cli::column<O>("arithmetic_diagnostic"),
      cli::column<O>("shift_quadrature_diagnostic"), cli::column<O>("shift_tail_bound"),
      cli::column<O>("sampling_error"), cli::column<std::size_t>("evaluations"), cli::column<std::string>("reason"));
  bool failed = false;
  Result<R> result;
  try
  {
    output.attach(moments, "moments");
    context.measure([&] {
      uni20::async::Async<Result<R>> pending;
      auto task = calculate(p, execution, pending.write());
      bethe::detail::schedule_owned(scheduler, std::move(task));
      result = pending.get_wait(scheduler);
    });
    bool const divergence = p->momentum && *p->momentum == R{4} * uni20::math::atan(p->one);
    std::array<O, 2> denominators;
    std::array<bethe::detail::CompensatedSum<R>, 2> sum{bethe::detail::CompensatedSum<R>(p->one - p->one),
                                                        bethe::detail::CompensatedSum<R>(p->one - p->one)};
    std::string label = "sum:";
    for (auto& row : result.rows)
    {
      for (auto* fields : {&row.values, &row.denominators, &row.kernel, &row.quadrature, &row.cutoff, &row.contour,
                           &row.arithmetic, &row.shift_quadrature, &row.shift_tail, &row.sampling})
        for (auto& value : *fields)
          if (value)
          {
            *value *= scale;
            if (!uni20::isfinite(*value)) throw std::overflow_error("channel conversion exceeds exponent range");
          }
      bool const complete = row.values[0] && row.values[1];
      failed |= !complete;
      label += (label == "sum:" ? "" : "+") + row.sector;
      for (int i = 0; i < 2; ++i)
      {
        if (row.denominators[i]) denominators[i] = row.denominators[i];
        if (row.values[i]) sum[i].add(*row.values[i]);
        moments.append(p->momentum, row.sector, i, row.values[i], row.denominators[i], row.fractions[i], row.status,
                       divergence && i == 0);
      }
    }
    if (result.rows.size() > 1)
      for (int i = 0; i < 2; ++i)
      {
        O value, fraction;
        if (!failed)
        {
          value = sum[i].value();
          if (!uni20::isfinite(*value)) throw std::overflow_error("partial moment sum exceeds exponent range");
          if (denominators[i])
            fraction = model::detail::spinon_weighted_value(*value, -uni20::math::log(*denominators[i]));
        }
        moments.append(p->momentum, label, i, value, denominators[i], fraction,
                       std::string(failed ? "incomplete" : "partial_sum"), divergence && i == 0);
      }
    report.result(!failed, failed ? "incomplete" : "finite_domain_estimate");
    auto summary = report.finish();
    output.overview(report.overview(summary));
    output.finish(moments, summary);
    output.write_table(
        "diagnostics", diagnostics,
        [&](auto& t) {
          for (auto const& row : result.rows)
            for (int i = 0; i < 2; ++i)
              t.append(row.sector, i, row.kernel[i], row.quadrature[i], row.cutoff[i], row.contour[i],
                       row.arithmetic[i], row.shift_quadrature[i], row.shift_tail[i], row.sampling[i], row.evaluations,
                       row.reason);
        },
        summary);
    auto converted = [&](O value, int power = 1) -> O {
      if (value)
        for (int i = 0; i < power; ++i)
        {
          *value *= scale;
          if (!uni20::isfinite(*value))
            throw std::overflow_error("replicate/covariance conversion exceeds exponent range");
        }
      return value;
    };
    if (result.six)
    {
      auto reps = data::make_data_table(
          "Six-spinon moment replicates", table_options, cli::column<std::uint64_t>("replicate"),
          cli::column<std::size_t>("completed_samples"), cli::column<std::size_t>("supported_configurations"),
          cli::column<O>("m0"), cli::column<O>("m1"), cli::column<O>("maximum_m0"), cli::column<O>("maximum_m1"),
          cli::column<std::optional<std::size_t>>("failed_sample"), cli::column<std::string>("reason"));
      output.write_table(
          "replicates", reps,
          [&](auto& t) {
            for (auto const& r : result.six->replicates)
              if (r)
                t.append(r->id, r->completed_samples, r->metadata.supported_configurations,
                         converted(r->means ? O((*r->means)[0]) : O{}), converted(r->means ? O((*r->means)[1]) : O{}),
                         converted(r->maximum_absolute_observation.size() >= 2 ? O(r->maximum_absolute_observation[0])
                                                                               : O{}),
                         converted(r->maximum_absolute_observation.size() >= 2 ? O(r->maximum_absolute_observation[1])
                                                                               : O{}),
                         r->failed_sample, r->reason);
          },
          summary);
      auto cov =
          data::make_data_table("Covariance of six-spinon moment means", table_options, cli::column<int>("moment_i"),
                                cli::column<int>("moment_j"), cli::column<O>("covariance"));
      output.write_table(
          "covariance", cov,
          [&](auto& t) {
            for (int i = 0; i < 2; ++i)
              for (int j = 0; j < 2; ++j)
                t.append(i, j, converted(result.six->covariance ? O((*result.six->covariance)[6 * i + j]) : O{}, 2));
          },
          summary);
    }
    output.finish_document();
  }
  catch (...)
  {
    cli::data::table_metadata aborted{{"Status", "aborted"}};
    if (!context.finished()) try
      {
        aborted = cli::run_summary(context.finish(uni20::run_outcome::failed));
      }
      catch (...)
      {}
    output.abort(moments, aborted);
    output.abort(diagnostics, aborted);
    throw;
  }
  return failed ? 2 : 0;
}
} // namespace
int main(int argc, char** argv)
{
  Arguments a;
  return cli::program_main(
      argc, argv, program_info(), [&](auto& app) { add_options(app, a); },
      [&](auto&) {
        if (a.global == bool(a.momentum)) throw std::invalid_argument("choose exactly one of --global and --momentum");
        auto selected = spinon::sectors(a.controls.sectors);
        bool const six = std::ranges::find(selected, 6) != selected.end();
        if ((a.controls.replicates || a.covariance) && !six)
          throw std::invalid_argument("replicates/covariance require sector 6");
        if (a.global && a.controls.extended_cutoff)
          throw std::invalid_argument("global refinement uses --four-extended-gap-cutoffs");
        if (!a.global && !a.extended_gaps.empty())
          throw std::invalid_argument("fixed-momentum refinement uses --four-extended-cutoff");
        if (bool(a.controls.fine_order) != (a.global ? !a.extended_gaps.empty() : bool(a.controls.extended_cutoff)))
          throw std::invalid_argument("four-spinon refinement requires both a fine order and extended domain");
        if (a.controls.diagnostics) a.output.screen_options.push_back("--diagnostics");
        if (a.controls.replicates) a.output.screen_options.push_back("--replicate-output");
        if (a.covariance) a.output.screen_options.push_back("--covariance");
        a.output.validate();
        return cli::dispatch_options_precision<cli::ArbitraryPrecision::real, model::StructureFactorOptions>(
            a.precision, [&](auto options) { return run(a, std::move(options), argc, argv); });
      });
}
