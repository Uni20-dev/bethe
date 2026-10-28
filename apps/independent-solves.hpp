// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/observable_grid.hpp>
#include <optional>
#include <uni20/async/tbb_scheduler.hpp>
#include <variant>

namespace bethe::cli
{
// Two independent solves, possibly with different result types. Use the same
// bounded grid as longer scans: on failure it waits for every admitted solve
// before returning. No worker touches run_context or output tables.
template <typename First, typename Second> auto independent_solves(First first, Second second, int threads)
{
  using A = std::invoke_result_t<First&>;
  using B = std::invoke_result_t<Second&>;
  if (threads < 1) throw std::invalid_argument("solver concurrency must be positive");
  if (threads == 1)
  {
    auto a = std::invoke(first);
    auto b = std::invoke(second);
    return std::pair<A, B>{std::move(a), std::move(b)};
  }
  uni20::async::TbbScheduler scheduler(2);
  std::optional<A> a;
  std::optional<B> b;
  using Value = std::variant<A, B>;
  bethe::evaluate_observable_grid(
      2,
      [&](std::size_t i) -> Value {
        if (i == 0) return Value(std::in_place_index<0>, std::invoke(first));
        return Value(std::in_place_index<1>, std::invoke(second));
      },
      [&](std::size_t i, Value value) {
        if (i == 0)
          a.emplace(std::get<0>(std::move(value)));
        else
          b.emplace(std::get<1>(std::move(value)));
      },
      {.scheduler = &scheduler, .max_in_flight = 2});
  return std::pair<A, B>{std::move(*a), std::move(*b)};
}
} // namespace bethe::cli
