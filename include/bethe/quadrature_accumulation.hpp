// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <array>
#include <bethe/observable_grid.hpp>
#include <bethe/solver.hpp>
#include <chrono>
#include <variant>

namespace bethe
{
enum class QuadratureStatus
{
  complete,
  numerical_failure,
  cancelled
};

template <uni20::Real R, std::size_t N, class Metadata = std::monostate> struct QuadratureObservation
{
    std::optional<std::array<R, N>> values;
    // Optional adapter-defined work count and failure details. The numerical
    // engine does not interpret kernels, geometry, units or uncertainty fields.
    std::size_t work_units = 0;
    Metadata metadata{};
    bool cancelled = false;
    std::string reason;
};

template <uni20::Real R, std::size_t N, class Metadata = std::monostate> struct QuadratureReduction
{
    QuadratureStatus status = QuadratureStatus::numerical_failure;
    std::optional<std::array<R, N>> sums;
    std::size_t total_nodes = 0, completed_nodes = 0, work_units = 0;
    std::optional<std::size_t> failed_node;
    std::optional<Metadata> failure_metadata;
    std::string reason;
    double wall_seconds = 0;
};

namespace detail
{
struct QuadratureStop
{};

template <uni20::Real R, std::size_t N, class Metadata, class Evaluate> struct QuadratureRun
{
    using Result = QuadratureReduction<R, N, Metadata>;
    using Observation = QuadratureObservation<R, N, Metadata>;
    Evaluate evaluate;
    R zero;
    ExecutionControl control;
    Result result;
    std::vector<CompensatedSum<R>> sums;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    QuadratureRun(std::size_t count, Evaluate f, R const& exemplar, std::stop_token stop)
        : evaluate(std::move(f)), zero(precision_zero(exemplar)), control(stop)
    {
      static_assert(N > 0);
      result.total_nodes = count;
      for (std::size_t i = 0; i < N; ++i)
        sums.emplace_back(zero);
    }
    void collect(std::size_t index, Observation observation)
    {
      if (observation.work_units > std::numeric_limits<std::size_t>::max() - result.work_units)
        throw std::overflow_error("quadrature work counter overflow");
      result.work_units += observation.work_units;
      // Cancellation is terminal even when an adapter supplied a last value.
      if (observation.cancelled) observation.values.reset();
      if (observation.values)
      {
        bool finite = true;
        for (auto const& value : *observation.values)
        {
          require_same_precision(value, zero);
          finite &= uni20::isfinite(value);
        }
        if (!finite)
        {
          observation.values.reset();
          observation.reason = "nonfinite quadrature contribution";
        }
      }
      if (!observation.values)
      {
        result.failed_node = index;
        result.failure_metadata = std::move(observation.metadata);
        result.status = observation.cancelled ? QuadratureStatus::cancelled : QuadratureStatus::numerical_failure;
        result.reason = observation.reason.empty()
                            ? (observation.cancelled ? "observation cancelled" : "observation has no numerical value")
                            : std::move(observation.reason);
        control.fail(std::make_exception_ptr(QuadratureStop{}));
        throw QuadratureStop{};
      }
      for (std::size_t i = 0; i < N; ++i)
        sums[i].add((*observation.values)[i]);
      ++result.completed_nodes;
    }

    static uni20::async::AsyncTask execute(std::shared_ptr<QuadratureRun> run, ExecutionOptions execution,
                                           uni20::async::WriteBuffer<Result> output)
    {
      auto destination = co_await output;
      bool complete = false;
      try
      {
        execution.stop = run->control.token();
        auto evaluate = [run](std::size_t index) {
          try
          {
            return std::invoke(run->evaluate, index, run->control.token());
          }
          catch (...)
          {
            run->control.fail();
            throw;
          }
        };
        auto collect = [run](std::size_t index, auto observation) {
          try
          {
            run->collect(index, std::move(observation));
          }
          catch (...)
          {
            run->control.fail();
            throw;
          }
        };
        using Window = ObservableWindow<decltype(evaluate)>;
        auto window = std::make_shared<Window>(execution.continuous_scheduler(), std::move(evaluate), execution,
                                               run->result.total_nodes);
        uni20::async::Async<std::size_t> work;
        // Construct the grid inside its sole admitted owner. The directly
        // awaited child drains all accepted observations before propagating.
        co_await observable_grid_task(std::move(window), std::move(collect), work.write());
        complete = true;
      }
      catch (QuadratureStop const&)
      {}
      catch (ExecutionCancelled const&)
      {
        if (!run->result.failed_node)
        {
          run->result.status = QuadratureStatus::cancelled;
          run->result.reason = "cancelled after draining admitted work";
        }
      }
      auto& result = run->result;
      if (complete)
      {
        std::array<R, N> values;
        bool finite = true;
        for (std::size_t i = 0; i < N; ++i)
        {
          values[i] = run->sums[i].value();
          finite &= uni20::isfinite(values[i]);
        }
        if (finite)
        {
          result.status = QuadratureStatus::complete;
          result.sums = std::move(values);
        }
        else
        {
          result.status = QuadratureStatus::numerical_failure;
          result.reason = "quadrature sum exceeded exponent range";
        }
      }
      result.wall_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - run->start).count();
      destination = std::move(result);
    }
};

// Prepare owned state without constructing or admitting a task, so model
// adapters can directly await the reduction inside their own final-result task.
template <uni20::Real R, std::size_t N, class Metadata = std::monostate, class Evaluate>
auto prepare_quadrature_run(std::size_t count, Evaluate evaluate, R const& exemplar, ExecutionOptions execution)
{
  (void)execution.continuous_scheduler();
  using Run = QuadratureRun<R, N, Metadata, Evaluate>;
  return std::make_shared<Run>(count, std::move(evaluate), exemplar, execution.stop);
}
} // namespace detail

/// Bounded ordered reduction of already weighted quadrature observations.
/// evaluate(index, stop) owns its per-call scratch and may run concurrently.
/// Domain maps, weights, physical kernels, diagnostics and refinement policies
/// belong to the caller. No failed node is replaced or dropped; incomplete
/// prefixes never become sums. Unexpected exceptions propagate after draining.
/// The scheduler remains borrowed through completion; callables are owned.
/// The public call admits only its final-result owner; child admission failures
/// are reported through that result after accepted observations have drained.
template <uni20::Real R, std::size_t N, class Metadata = std::monostate, class Evaluate>
[[nodiscard]] uni20::async::Async<QuadratureReduction<R, N, Metadata>>
accumulate_quadrature_async(std::size_t count, Evaluate evaluate, R const& exemplar, ExecutionOptions execution = {})
{
  using Run = detail::QuadratureRun<R, N, Metadata, Evaluate>;
  auto& scheduler = execution.continuous_scheduler();
  execution.scheduler = &scheduler;
  auto run = detail::prepare_quadrature_run<R, N, Metadata>(count, std::move(evaluate), exemplar, execution);
  uni20::async::Async<typename Run::Result> output;
  auto task = Run::execute(std::move(run), execution, output.write());
  detail::schedule_owned(scheduler, std::move(task));
  return output;
}
} // namespace bethe
