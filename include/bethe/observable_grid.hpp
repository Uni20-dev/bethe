// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <algorithm>
#include <bethe/detail/async_execution.hpp>
#include <functional>
#include <memory>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

namespace bethe
{
namespace detail
{
struct UnobservedGridWait
{
    int operator()(std::size_t, std::size_t) const { return 0; }
};

template <typename Value> struct ResetObservableValue
{
    Value& value;
    ~ResetObservableValue() { value.reset(); }
};

// Slots are private to the coordinator. Once its read wait has returned, no
// other reader or writer can prevent this exclusive access. Acquiring the
// writer explicitly also preserves Uni20's exception and debug-proxy checks.
template <typename Value>
auto writable_observable_result(uni20::async::Async<Value>& result, uni20::async::IAsyncScheduler& scheduler)
{
  (void)result.get_wait(scheduler);
  auto writer = result.write();
  if (!writer.await_ready()) throw std::logic_error("observable-grid result has an unexpected dependency");
  return writer;
}

template <typename Value>
Value take_observable_result(uni20::async::Async<Value>& result, uni20::async::IAsyncScheduler& scheduler)
{
  auto writer = writable_observable_result(result, scheduler);
  auto source = writer.await_resume();
  ResetObservableValue reset{source.get()};
  return std::move(source.get());
}

template <typename Value>
void discard_observable_result(uni20::async::Async<Value>& result, uni20::async::IAsyncScheduler& scheduler)
{
  auto writer = writable_observable_result(result, scheduler);
  writer.await_resume()->reset();
}

/// Ordered output requires a different admission policy from unordered/top-K
/// scanners: uncollected results retain their credits. This is a bounded ring,
/// not a second thread pool or an unbounded reorder queue.
template <typename Evaluate> struct ObservableWindow : std::enable_shared_from_this<ObservableWindow<Evaluate>>
{
    using Result = std::remove_cv_t<std::invoke_result_t<Evaluate&, std::size_t>>;
    static_assert(std::is_object_v<Result> && !std::is_array_v<Result>,
                  "grid evaluation must return an owning value, not void or a reference");
    using Value = std::optional<Result>;
    using Pending = uni20::async::Async<Value>;
    uni20::async::IAsyncScheduler& scheduler;
    Evaluate evaluate;
    ExecutionControl control;
    std::vector<std::unique_ptr<Pending>> slots;
    std::size_t count;

    ObservableWindow(uni20::async::IAsyncScheduler& s, Evaluate e, ExecutionOptions options, std::size_t n)
        : scheduler(s), evaluate(std::move(e)), control(options.stop), slots(std::min(n, options.max_in_flight)),
          count(n)
    {}

    static uni20::async::AsyncTask point(std::shared_ptr<ObservableWindow> self, std::size_t i,
                                         uni20::async::WriteBuffer<Value> output)
    {
      try
      {
        auto destination = co_await output;
        destination.emplace(); // explicit empty result when cancelled before entry
        if (!self->control.stopped()) destination->emplace(std::invoke(self->evaluate, i));
      }
      catch (...)
      {
        self->control.fail();
        throw;
      }
    }

    void start(std::size_t i)
    {
      if (control.stopped()) return;
      auto pending = std::make_unique<Pending>();
      auto task = point(this->shared_from_this(), i, pending->write());
      schedule_owned(scheduler, std::move(task));
      slots[i % slots.size()] = std::move(pending); // noexcept, after successful admission
    }

    void prime()
    {
      for (std::size_t i = 0; i < slots.size() && !control.stopped(); ++i)
        start(i);
    }

    void refill(std::size_t i)
    {
      // Subtraction avoids wrapping when a caller supplies a very large count.
      if (slots.size() < count - i) start(i + slots.size());
    }

    void finish()
    {
      control.rethrow_error(); // an execution failure takes precedence over cancellation
      if (control.stopped()) throw ExecutionCancelled{};
    }
};

template <typename Evaluate, typename Collect>
uni20::async::AsyncTask observable_grid_task(std::shared_ptr<ObservableWindow<Evaluate>> window, Collect collect,
                                             uni20::async::WriteBuffer<std::size_t> output)
{
  auto destination = co_await output;
  try
  {
    window->prime();
    for (std::size_t i = 0; i < window->count && !window->control.stopped(); ++i)
    {
      auto& slot = window->slots[i % window->slots.size()];
      {
        auto value = co_await slot->write();
        ResetObservableValue reset{value.get()};
        if (!window->control.stopped())
        {
          if (!value->has_value()) throw std::logic_error("missing observable-grid value");
          std::invoke(collect, i, std::move(value->value()));
        }
      }
      slot.reset();
      window->refill(i);
    }
  }
  catch (...)
  {
    window->control.fail();
  }
  // A collector, evaluator or admission failure must not orphan earlier work.
  // Await each accepted task and destroy discarded values before completion.
  // Readiness can precede release of the producer's last storage reference, so
  // dropping our slot alone does not guarantee prompt payload destruction.
  for (auto& slot : window->slots)
    if (slot)
    {
      try
      {
        auto value = co_await slot->write();
        value->reset();
      }
      catch (...)
      {
        window->control.fail();
      }
      slot.reset();
    }
  window->finish();
  destination = window->count;
}
} // namespace detail

/// Continuously evaluate independent indices with bounded ordered collection.
/// Owns decay-copied callables; references captured inside them remain borrowed.
/// evaluate(index) may run concurrently. collect(index,value) is serialized in
/// increasing index order ON THE SCHEDULER, not necessarily the initiating thread.
/// Returns the complete point count; execution failure/cancellation throws after
/// all admitted tasks finish. Previously collected ordered prefixes are retained.
/// Numerical failure values are ordinary results, never silently filtered out.
///
/// max_in_flight bounds active AND completed-but-uncollected points. A slow early
/// point/collector applies backpressure. No whole-batch join or count-sized task
/// graph is created. Collector-retained storage is outside this bound. Scheduler
/// lifetime must cover completion, even when the returned handle is discarded.
template <typename Evaluate, typename Collect>
[[nodiscard]] uni20::async::Async<std::size_t>
evaluate_observable_grid_async(std::size_t count, Evaluate evaluate, Collect collect, ExecutionOptions execution = {})
{
  auto& scheduler = execution.continuous_scheduler();
  auto window = std::make_shared<detail::ObservableWindow<Evaluate>>(scheduler, std::move(evaluate), execution, count);
  uni20::async::Async<std::size_t> result;
  auto task = detail::observable_grid_task(std::move(window), std::move(collect), result.write());
  detail::schedule_owned(scheduler, std::move(task));
  return result;
}

/// Synchronous adapter with the SAME bounded window, but collection remains on
/// the calling thread (safe for CLI tables). Evaluators may overlap collection.
/// Callables are borrowed only until all owned work has joined before return.
/// Cooperative waits may execute a kernel on the caller, delaying collection;
/// use the Async entry point when the coordinator must progress on other threads.
///
/// Optional observe_wait(index,1) runs on the calling thread and returns a scope
/// guard around waiting/moving that point, excluding collection. Unlike the old
/// public batched adapter this is not a compute-only batch timer: workers can
/// compute during collection too. Time the whole call for end-to-end costs.
/// batch_size is retained for older batched clients; this uses max_in_flight.
template <typename Evaluate, typename Collect, typename ObserveWait = detail::UnobservedGridWait>
void evaluate_observable_grid(std::size_t count, Evaluate&& evaluate, Collect&& collect,
                              ExecutionOptions execution = {}, ObserveWait&& observe_wait = {})
{
  auto& scheduler = execution.continuous_scheduler();
  using Window = detail::ObservableWindow<std::reference_wrapper<std::remove_reference_t<Evaluate>>>;
  auto window = std::make_shared<Window>(scheduler, std::ref(evaluate), execution, count);
  try
  {
    window->prime();
    for (std::size_t i = 0; i < count && !window->control.stopped(); ++i)
    {
      auto& slot = window->slots[i % window->slots.size()];
      {
        auto value = [&] {
          [[maybe_unused]] auto scope = std::invoke(observe_wait, i, std::size_t{1});
          return detail::take_observable_result(*slot, scheduler);
        }();
        slot.reset();
        if (!window->control.stopped())
        {
          if (!value) throw std::logic_error("missing observable-grid value");
          std::invoke(collect, i, std::move(*value));
        }
      }
      window->refill(i);
    }
  }
  catch (...)
  {
    window->control.fail();
  }
  for (auto& slot : window->slots)
    if (slot)
    {
      try
      {
        detail::discard_observable_result(*slot, scheduler);
      }
      catch (...)
      {
        window->control.fail();
      }
      slot.reset();
    }
  window->finish();
}
} // namespace bethe
