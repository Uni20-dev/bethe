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

namespace detail
{
template <uni20::Real R, class Metadata> struct MonteCarloAccumulator
{
    MonteCarloReplicate<R, Metadata> row;
    R zero;
    std::vector<SpinonSum<R>> sums;

    MonteCarloAccumulator(std::uint64_t id, std::size_t dimension, R const& exemplar)
        : zero(precision_zero(exemplar)), sums(dimension, SpinonSum<R>(zero))
    {
      row.id = id;
      row.maximum_absolute_observation.assign(dimension, zero);
    }

    void fail(std::size_t index, std::string reason, bool cancelled = false)
    {
      if (!row.failed_sample || index < *row.failed_sample)
      {
        row.failed_sample = index;
        row.reason = std::move(reason);
      }
      row.cancelled = row.cancelled || cancelled;
    }

    bool accept(std::size_t index, MonteCarloObservation<R, Metadata> value)
    {
      // Numerical failures retain the valid prefix, independent of speculative
      // work in later blocks. Cancellation retains completed observations but
      // can never produce a successful estimate.
      if (row.failed_sample && !row.cancelled) return false;
      row.metadata += value.metadata;
      if (!value.values || value.cancelled)
      {
        fail(index, value.reason.empty() ? "observation failed without a value" : value.reason, value.cancelled);
        return false;
      }
      if (value.values->size() != sums.size()) throw std::invalid_argument("observation has wrong dimension");
      for (auto const& x : *value.values)
      {
        require_same_precision(x, zero);
        if (!uni20::isfinite(x))
        {
          fail(index, "observation contains a nonfinite value");
          return false;
        }
      }
      for (std::size_t field = 0; field < sums.size(); ++field)
      {
        R const& x = (*value.values)[field];
        sums[field].add(x);
        row.maximum_absolute_observation[field] =
            std::max(row.maximum_absolute_observation[field], uni20::math::abs(x));
      }
      ++row.completed_samples;
      return true;
    }

    void finish(std::size_t samples)
    {
      if (row.failed_sample || row.cancelled) return;
      if (row.completed_samples != samples)
      {
        row.reason = "incomplete replicate";
        return;
      }
      R const count = uni20::scalar_like(zero, samples);
      std::vector<R> means;
      means.reserve(sums.size());
      for (auto const& sum : sums)
      {
        R const total = sum.value(), mean = total / count;
        if (!uni20::isfinite(mean) || (total != zero && mean == zero))
        {
          row.reason = "replicate mean exceeded exponent range";
          return;
        }
        means.push_back(mean);
      }
      row.means = std::move(means);
    }
};
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
  detail::MonteCarloAccumulator<R, Metadata> accumulator(options.first_replicate + index, dimension, exemplar);
  MonteCarloStream<D> stream(options.seed, accumulator.row.id, options.uniform_bits);
  for (std::size_t i = 0; i < options.samples_per_replicate; ++i)
  {
    if (stop.stop_requested())
    {
      accumulator.fail(i, "replicate cancelled", true);
      break;
    }
    if (!accumulator.accept(i, std::invoke(observe, stream.next(), stop))) break;
  }
  accumulator.finish(options.samples_per_replicate);
  return std::move(accumulator.row);
}

namespace detail
{
template <std::size_t D, uni20::Real R, class Observe> struct MonteCarloRun
{
    using Metadata = SampleMetadata<D, Observe>;
    using Result = MonteCarloEstimate<R, Metadata>;
    using Observation = MonteCarloObservation<R, Metadata>;
    using Accumulator = MonteCarloAccumulator<R, Metadata>;
    struct Block
    {
        std::size_t index, replicate, first, count;
        bool entered = false, cancelled = false;
        std::vector<Observation> observations;
    };
    struct Input
    {
        std::shared_ptr<Block> block;
        MonteCarloStream<D> stream;
    };
    MonteCarloOptions options;
    std::size_t dimension, grain, blocks_per_replicate;
    R exemplar;
    Observe observe;
    Result result;
    std::vector<std::unique_ptr<Accumulator>> accumulators;
    // Only the coordinator accesses this ring. Workers own their individual
    // Block through Input; completion/drain precedes collection or destruction.
    std::vector<std::shared_ptr<Block>> pending;
    std::optional<MonteCarloStream<D>> stream;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    MonteCarloRun(MonteCarloOptions o, std::size_t d, R r, Observe f, ExecutionOptions execution)
        : options(o), dimension(d), grain(std::min(o.samples_per_replicate, execution.batch_size)),
          blocks_per_replicate(1 + (o.samples_per_replicate - 1) / grain), exemplar(std::move(r)),
          observe(std::move(f)), accumulators(o.replicates),
          pending(std::min(o.replicates * blocks_per_replicate, execution.max_in_flight))
    {
      result.sampling = options;
      result.replicates.resize(options.replicates);
    }

    Input prepare(std::size_t index)
    {
      auto const replicate = index / blocks_per_replicate;
      auto const first = (index % blocks_per_replicate) * grain;
      auto const count = std::min(grain, options.samples_per_replicate - first);
      if (first == 0) stream.emplace(options.seed, options.first_replicate + replicate, options.uniform_bits);
      auto block = std::make_shared<Block>(Block{index, replicate, first, count});
      auto& slot = pending[index % pending.size()];
      if (slot) throw std::logic_error("Monte Carlo block credit was not released");
      Input input{block, *stream};
      stream->discard(count);
      slot = std::move(block);
      return input;
    }

    std::shared_ptr<Block> evaluate(Input input, std::stop_token stop) const
    {
      auto& block = *input.block;
      block.entered = true;
      block.observations.reserve(block.count);
      for (std::size_t i = 0; i < block.count; ++i)
      {
        if (stop.stop_requested())
        {
          block.cancelled = true;
          break;
        }
        block.observations.push_back(std::invoke(std::as_const(observe), input.stream.next(), stop));
        auto const& observation = block.observations.back();
        if (!observation.values || observation.cancelled) break;
      }
      return std::move(input.block);
    }

    void collect(std::size_t index, std::shared_ptr<Block> block)
    {
      auto& slot = pending[index % pending.size()];
      if (slot != block || index != block->index) throw std::logic_error("Monte Carlo block identity changed");
      auto& accumulator = accumulators[block->replicate];
      if (!accumulator)
        accumulator = std::make_unique<Accumulator>(options.first_replicate + block->replicate, dimension, exemplar);
      for (std::size_t j = 0; j < block->observations.size(); ++j)
        if (!accumulator->accept(block->first + j, std::move(block->observations[j]))) break;
      if (block->cancelled || !block->entered)
        accumulator->fail(block->first + block->observations.size(), "replicate cancelled", true);
      slot.reset();
    }

    void collect_drained()
    {
      // Cancellation stops ordered grid collection, but completed observations
      // remain owned here. Recover them only after every admitted task drains.
      std::vector<std::shared_ptr<Block>> remaining;
      for (auto const& block : pending)
        if (block) remaining.push_back(block);
      std::ranges::sort(remaining, {}, [](auto const& block) { return block->index; });
      for (auto& block : remaining)
      {
        auto index = block->index;
        collect(index, std::move(block));
      }
      for (std::size_t i = 0; i < accumulators.size(); ++i)
        if (auto& accumulator = accumulators[i])
        {
          accumulator->finish(options.samples_per_replicate);
          result.replicates[i] = std::move(accumulator->row);
        }
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
      catch (...)
      {
        // The grid has drained. Release retained observations before exposing
        // its exception, even if the scheduler keeps this coroutine frame alive.
        run->pending.clear();
        run->accumulators.clear();
        throw;
      }
      run->collect_drained();
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

/// Bounded Async execution of sample blocks within independent replicates.
/// execution.batch_size is the maximum samples per task; max_in_flight bounds
/// admitted blocks and their retained observations. Sample identities and
/// accumulation order do not depend on task size or completion order.
/// Observe receives exact uniform identities; physical proposals remain its responsibility.
template <std::size_t D, uni20::Real R, class Observe>
[[nodiscard]] auto monte_carlo_async(MonteCarloOptions options, std::size_t dimension, R exemplar, Observe observe,
                                     ExecutionOptions execution = {})
{
  static_assert(D > 0);
  detail::validate_monte_carlo(options, dimension, exemplar);
  auto& scheduler = execution.continuous_scheduler();
  using Run = detail::MonteCarloRun<D, R, Observe>;
  auto run = std::make_shared<Run>(options, dimension, std::move(exemplar), std::move(observe), execution);
  auto work = evaluate_prepared_observable_grid_async(
      options.replicates * run->blocks_per_replicate, [run](std::size_t i) { return run->prepare(i); },
      [run, stop = execution.stop](typename Run::Input input) { return run->evaluate(std::move(input), stop); },
      [run](std::size_t i, std::shared_ptr<typename Run::Block> block) { run->collect(i, std::move(block)); },
      execution);
  uni20::async::Async<typename Run::Result> output;
  auto task = Run::finish(work.read(), std::move(run), output.write());
  detail::schedule_owned(scheduler, std::move(task));
  return output;
}
} // namespace bethe
