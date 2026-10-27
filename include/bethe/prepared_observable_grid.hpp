// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/observable_grid.hpp>
#include <functional>
#include <stop_token>
#include <type_traits>
#include <utility>

namespace bethe
{
namespace detail
{
template <typename Prepare, typename Background>
uni20::async::AsyncTask prepare_background_task(Prepare prepare, std::stop_token stop,
                                                uni20::async::WriteBuffer<Background> output)
{
  auto destination = co_await output;
  if (stop.stop_requested()) throw ExecutionCancelled{};
  destination.emplace(std::invoke(prepare));
  if (stop.stop_requested()) throw ExecutionCancelled{};
}

template <typename Background, typename Evaluate, typename Collect>
uni20::async::AsyncTask prepared_observable_grid_task(uni20::async::ReadBuffer<Background> input, std::size_t count,
                                                      Evaluate evaluate, Collect collect, ExecutionOptions execution,
                                                      uni20::async::WriteBuffer<std::size_t> output)
{
  auto destination = co_await output;
  // Await an lvalue reader by reference: `auto background` would copy the
  // prepared meshes. Keep this reader enrolled until the entire grid drains.
  auto const& background = co_await input;
  if (execution.stop.stop_requested()) throw ExecutionCancelled{};
  auto query = [&background, evaluate = std::move(evaluate)](std::size_t i) mutable {
    return std::invoke(evaluate, background, i);
  };
  using Window = ObservableWindow<decltype(query)>;
  auto window = std::make_shared<Window>(*execution.scheduler, std::move(query), execution, count);
  uni20::async::Async<std::size_t> pending;
  auto writer = pending.write();
  // Enroll completion before admitting children. A reader-allocation failure
  // must not release this background while already-admitted queries use it.
  auto complete = pending.read();
  auto task = observable_grid_task(std::move(window), std::move(collect), std::move(writer));
  schedule_owned(*execution.scheduler, std::move(task));
  destination = co_await complete;
}
} // namespace detail

/// Schedule one synchronous preparation, owning the decay-copied callable.
/// Its result must own its data and support concurrent const queries. Numerical
/// status remains in that result; only execution errors/cancellation throw.
/// Cancellation is checked before and after preparation, not inside its kernel.
/// The borrowed scheduler and any references captured by prepare must outlive
/// completion, even if the returned handle is dropped.
template <typename Prepare>
[[nodiscard]] auto prepare_background_async(Prepare prepare, ExecutionOptions execution = {})
{
  using Background = std::remove_cv_t<std::invoke_result_t<Prepare&>>;
  static_assert(std::is_object_v<Background> && !std::is_array_v<Background>,
                "background preparation must return an owning value, not void or a reference");
  auto& scheduler = execution.continuous_scheduler();
  uni20::async::Async<Background> result;
  auto task = detail::prepare_background_task(std::move(prepare), execution.stop, result.write());
  detail::schedule_owned(scheduler, std::move(task));
  return result;
}

/// Submit a bounded ordered grid before its shared background is ready.
/// evaluate(background const&, index) must preserve model preparation failures
/// in its numerical results (or explicitly reject them); it may run concurrently.
/// The input reader retains the background without copying it, including after
/// the original Async handle is destroyed. Different grids enroll separate readers.
/// collect(index, value) runs serially ON THE SCHEDULER, as in the ordinary grid.
///
/// Cancellation belongs to this grid, not to the shared preparation or sibling
/// grids. Dependency waits finish before cancellation is observed; kernels are
/// not interrupted. Preparation exceptions propagate even to an empty grid.
/// Completion, including failure, covers all admitted queries before releasing
/// the background reader. Each grid has its own max_in_flight retention bound;
/// the scheduler controls shared CPU capacity. Captured references and schedulers
/// remain borrowed, including when the result handle is dropped.
template <typename Background, typename Evaluate, typename Collect>
[[nodiscard]] uni20::async::Async<std::size_t>
evaluate_observable_grid_async(uni20::async::ReadBuffer<Background> background, std::size_t count, Evaluate evaluate,
                               Collect collect, ExecutionOptions execution = {})
{
  auto& scheduler = execution.continuous_scheduler();
  // Resolve once: a resumed coroutine must not consult another thread's active
  // scheduler when it submits the dependent grid.
  execution.scheduler = &scheduler;
  uni20::async::Async<std::size_t> result;
  auto task = detail::prepared_observable_grid_task(std::move(background), count, std::move(evaluate),
                                                    std::move(collect), execution, result.write());
  detail::schedule_owned(scheduler, std::move(task));
  return result;
}
} // namespace bethe
