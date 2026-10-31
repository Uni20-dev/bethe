// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/detail/spinon_numerics.hpp>
#include <bethe/observable_grid.hpp>
#include <bethe/sampling.hpp>
#include <chrono>

namespace bethe
{
struct MonteCarloOptions
{
    std::size_t replicates = 8, samples_per_replicate = 256;
    std::uint64_t seed = 0, first_replicate = 0;
    unsigned uniform_bits = 52;
};

struct EmptySampleMetadata
{
    EmptySampleMetadata& operator+=(EmptySampleMetadata const&) { return *this; }
};

/// Values include the integrand's proposal/Jacobian weights. The engine does
/// not infer a physical measure or adjust an importance-sampling distribution.
/// Metadata supplies operator+= for work counters and model failure details.
template <uni20::Real R, class Metadata = EmptySampleMetadata> struct MonteCarloObservation
{
    using metadata_type = Metadata;
    std::optional<std::vector<R>> values;
    Metadata metadata;
    bool cancelled = false;
    std::string reason;
};

template <uni20::Real R, class Metadata = EmptySampleMetadata> struct MonteCarloReplicate
{
    std::uint64_t id = 0;
    std::size_t completed_samples = 0;
    std::optional<std::vector<R>> means;
    std::vector<R> maximum_absolute_observation;
    std::optional<std::size_t> failed_sample;
    Metadata metadata;
    bool cancelled = false;
    std::string reason;
};

enum class MonteCarloStatus
{
  complete,
  numerical_failure,
  cancelled
};

template <uni20::Real R, class Metadata = EmptySampleMetadata> struct MonteCarloEstimate
{
    MonteCarloStatus status = MonteCarloStatus::numerical_failure;
    MonteCarloOptions sampling;
    std::optional<std::vector<R>> mean, covariance; // covariance of the mean, row-major
    std::vector<std::optional<MonteCarloReplicate<R, Metadata>>> replicates;
    Metadata metadata;
    std::size_t completed_samples = 0;
    double wall_seconds = 0;
    std::string reason;
};

namespace detail
{
template <uni20::Real R>
void validate_monte_carlo(MonteCarloOptions const& options, std::size_t dimension, R const& exemplar)
{
  R const zero = precision_zero(exemplar);
  if (!options.replicates || !options.samples_per_replicate || !dimension || !options.uniform_bits ||
      options.uniform_bits >= unsigned(std::numeric_limits<int>::max()) ||
      uni20::numeric_limits<R>::digits(zero) < options.uniform_bits + 1 ||
      options.replicates - 1 > std::numeric_limits<std::uint64_t>::max() - options.first_replicate ||
      options.samples_per_replicate > std::numeric_limits<std::size_t>::max() / options.replicates ||
      dimension > std::numeric_limits<std::size_t>::max() / dimension)
    throw std::invalid_argument("invalid Monte Carlo budget, dimensions or input precision");
}

template <std::size_t D, class Observe>
using SampleMetadata =
    typename std::invoke_result_t<Observe const&, UniformSample<D> const&, std::stop_token>::metadata_type;
} // namespace detail

/// Evaluate one independent replicate. Observe must be const-callable and own
/// any state needed for replay. Numerical failures retain the exact sample
/// index; invalid interface contracts and unexpected exceptions propagate.
template <std::size_t D, uni20::Real R, class Observe>
auto monte_carlo_replicate(MonteCarloOptions options, std::size_t index, std::size_t dimension, R const& exemplar,
                           Observe const& observe, std::stop_token stop = {})
{
  static_assert(D > 0);
  detail::validate_monte_carlo(options, dimension, exemplar);
  if (index >= options.replicates) throw std::out_of_range("replicate outside requested work");
  using Metadata = detail::SampleMetadata<D, Observe>;
  MonteCarloReplicate<R, Metadata> result;
  result.id = options.first_replicate + index;
  MonteCarloStream<D> stream(options.seed, result.id, options.uniform_bits);
  R const zero = detail::precision_zero(exemplar);
  result.maximum_absolute_observation.assign(dimension, zero);
  std::vector<detail::SpinonSum<R>> sums(dimension, detail::SpinonSum<R>(zero));
  for (std::size_t i = 0; i < options.samples_per_replicate; ++i)
  {
    if (stop.stop_requested())
    {
      result.cancelled = true;
      result.failed_sample = i;
      result.reason = "replicate cancelled";
      return result;
    }
    auto value = std::invoke(observe, stream.next(), stop);
    result.metadata += value.metadata;
    if (!value.values || value.cancelled)
    {
      result.failed_sample = i;
      result.cancelled = value.cancelled;
      result.reason = value.reason.empty() ? "observation failed without a value" : value.reason;
      return result;
    }
    if (value.values->size() != dimension) throw std::invalid_argument("observation has wrong dimension");
    for (auto const& x : *value.values)
    {
      detail::require_same_precision(x, zero);
      if (!uni20::isfinite(x))
      {
        result.failed_sample = i;
        result.reason = "observation contains a nonfinite value";
        return result;
      }
    }
    for (std::size_t field = 0; field < dimension; ++field)
    {
      R const& x = (*value.values)[field];
      sums[field].add(x);
      result.maximum_absolute_observation[field] =
          std::max(result.maximum_absolute_observation[field], uni20::math::abs(x));
    }
    ++result.completed_samples;
  }
  R const count = uni20::scalar_like(zero, options.samples_per_replicate);
  std::vector<R> means;
  means.reserve(dimension);
  for (auto const& sum : sums)
  {
    R const total = sum.value(), mean = total / count;
    if (!uni20::isfinite(mean) || (total != zero && mean == zero))
    {
      result.reason = "replicate mean exceeded exponent range";
      return result;
    }
    means.push_back(mean);
  }
  result.means = std::move(means);
  return result;
}

namespace detail
{
template <std::size_t D, uni20::Real R, class Observe> struct MonteCarloRun
{
    using Metadata = SampleMetadata<D, Observe>;
    using Result = MonteCarloEstimate<R, Metadata>;
    MonteCarloOptions options;
    std::size_t dimension;
    R exemplar;
    Observe observe;
    Result result;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    MonteCarloRun(MonteCarloOptions o, std::size_t d, R r, Observe f)
        : options(o), dimension(d), exemplar(std::move(r)), observe(std::move(f))
    {
      result.sampling = options;
      result.replicates.resize(options.replicates);
    }

    static uni20::async::AsyncTask finish(uni20::async::ReadBuffer<std::size_t> work,
                                          std::shared_ptr<MonteCarloRun> run, uni20::async::WriteBuffer<Result> output)
    {
      auto destination = co_await output;
      bool cancelled = false;
      try
      {
        (void)co_await work;
      }
      catch (ExecutionCancelled const&)
      {
        cancelled = true;
      }
      auto& r = run->result;
      ReplicateStatistics<R> statistics(run->options.replicates, run->options.samples_per_replicate, run->dimension,
                                        run->exemplar);
      for (std::size_t i = 0; i < r.replicates.size(); ++i)
        if (auto const& row = r.replicates[i])
        {
          r.completed_samples += row->completed_samples;
          r.metadata += row->metadata;
          cancelled = cancelled || row->cancelled;
          if (row->means)
            statistics.accept(i, row->completed_samples, *row->means);
          else
            statistics.reject(i, row->reason);
        }
      if (cancelled)
      {
        r.status = MonteCarloStatus::cancelled;
        r.reason = "cancelled after draining admitted work";
      }
      else
      {
        try
        {
          auto summary = statistics.summary();
          if (summary.mean)
          {
            r.status = MonteCarloStatus::complete;
            r.mean = std::move(summary.mean);
            r.covariance = std::move(summary.covariance);
          }
          else
            r.reason = "incomplete or failed replicate; no survivor-only estimate";
        }
        catch (std::runtime_error const& error)
        {
          r.reason = error.what();
        }
      }
      r.wall_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - run->start).count();
      destination = std::move(r);
    }
};
} // namespace detail

/// Bounded Async execution of independent replicates, with ordered reduction
/// and retained terminal records. Observe receives exact uniform identities;
/// proposal selection and adaptation policies belong to its caller.
template <std::size_t D, uni20::Real R, class Observe>
[[nodiscard]] auto monte_carlo_async(MonteCarloOptions options, std::size_t dimension, R exemplar, Observe observe,
                                     ExecutionOptions execution = {})
{
  static_assert(D > 0);
  detail::validate_monte_carlo(options, dimension, exemplar);
  auto& scheduler = execution.continuous_scheduler();
  using Run = detail::MonteCarloRun<D, R, Observe>;
  auto run = std::make_shared<Run>(options, dimension, std::move(exemplar), std::move(observe));
  auto work = evaluate_observable_grid_async(
      options.replicates,
      [run, stop = execution.stop](std::size_t i) {
        run->result.replicates[i] =
            monte_carlo_replicate<D>(run->options, i, run->dimension, run->exemplar, std::as_const(run->observe), stop);
        return i;
      },
      [](std::size_t, std::size_t) {}, execution);
  uni20::async::Async<typename Run::Result> output;
  auto task = Run::finish(work.read(), std::move(run), output.write());
  detail::schedule_owned(scheduler, std::move(task));
  return output;
}
} // namespace bethe
