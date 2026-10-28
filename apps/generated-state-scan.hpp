// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/detail/continuous_work.hpp>
#include <functional>
#include <mutex>
#include <optional>
#include <type_traits>
#include <uni20/async/tbb_scheduler.hpp>
#include <utility>
#include <vector>

namespace bethe::cli
{
namespace detail
{
template <typename Next, typename Solve>
using GeneratedState = std::invoke_result_t<Solve&, std::invoke_result_t<Next&>>;

template <typename Next, typename Solve>
uni20::async::AsyncTask
generated_states_task(uni20::async::IAsyncScheduler& scheduler, std::size_t workers, std::size_t count, Next next,
                      Solve solve, uni20::async::WriteBuffer<std::vector<GeneratedState<Next, Solve>>> output)
{
  auto destination = co_await output;
  using State = GeneratedState<Next, Solve>;
  std::vector<std::optional<State>> states(count);
  std::mutex generation;
  std::size_t next_id = 0;
  bethe::detail::continuous_work(scheduler, workers, count, [&](std::size_t, std::stop_token) {
    auto request = [&] {
      std::lock_guard lock(generation);
      auto value = std::invoke(next);
      return std::pair{next_id++, std::move(value)};
    }();
    // IDs follow generation order, never runner arrival or completion order.
    // Each runner owns its request and writes a different non-packed slot.
    states[request.first].emplace(std::invoke(solve, std::move(request.second)));
  });
  std::vector<State> result;
  result.reserve(count);
  for (auto& state : states)
    result.push_back(std::move(*state));
  destination = std::move(result);
}
} // namespace detail

/// CLI adapter for scans whose labels must be generated in sequence. Next
/// returns an owning request; Solve must permit simultaneous calls. Only label
/// generation is serialized. Callbacks are owned until all admitted solves
/// finish, including after an exception. Final states are retained for ranking
/// and calling-thread output; transient requests are bounded by concurrency.
template <typename Next, typename Solve>
auto scan_generated_states(std::size_t count, Next next, Solve solve, int threads)
{
  using State = detail::GeneratedState<Next, Solve>;
  static_assert(std::is_object_v<State>);
  static_assert(!std::is_reference_v<std::invoke_result_t<Next&>>, "requests must own their labels");
  if (threads < 1) throw std::invalid_argument("scan threads must be positive");
  if (threads == 1 || count < 2)
  {
    std::vector<State> states;
    states.reserve(count);
    for (std::size_t i = 0; i < count; ++i)
      states.push_back(std::invoke(solve, std::invoke(next)));
    return states;
  }
  auto const workers = std::min(count, static_cast<std::size_t>(threads));
  uni20::async::TbbScheduler scheduler(workers);
  uni20::async::Async<std::vector<State>> result;
  auto task =
      detail::generated_states_task(scheduler, workers, count, std::move(next), std::move(solve), result.write());
  bethe::detail::schedule_owned(scheduler, std::move(task));
  return bethe::detail::take_result(result, scheduler);
}
} // namespace bethe::cli
