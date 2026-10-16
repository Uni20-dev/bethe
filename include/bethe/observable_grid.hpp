// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <algorithm>
#include <bethe/execution.hpp>
#include <functional>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

namespace bethe
{
namespace detail
{
struct UnobservedGridBatch
{
    int operator()(std::size_t, std::size_t) const { return 0; }
};
} // namespace detail

/// Evaluate independent logical points [0,count) in bounded synchronous batches.
/// evaluate(index) returns an owning value and must permit concurrent calls;
/// collect(index, result) runs on the calling thread in increasing index order.
/// Numerical failure/status values are collected normally, never filtered here.
/// No tables, floating-point reductions or process-global scheduler changes.
///
/// At most execution.batch_size result slots are retained by this adapter.
/// Callers must validate the full request before opening output, share only
/// immutable prepared data, and keep worker scratch private. Continuation chains
/// are not independent points. Any storage retained by collect is caller-owned.
///
/// Optional observe_batch(offset,size) runs on the calling thread and returns a
/// scope guard covering only the synchronous scheduler call, not collection.
/// This permits computation timing without a run-context dependency or worker
/// writes to shared timing state. Exceptions propagate after workers have joined;
/// a failed evaluation batch emits no rows, leaving earlier batches collected.
template <typename Evaluate, typename Collect, typename ObserveBatch = detail::UnobservedGridBatch>
void evaluate_observable_grid(std::size_t count, Evaluate&& evaluate, Collect&& collect,
                              ExecutionOptions execution = {}, ObserveBatch&& observe_batch = {})
{
  using Result = std::remove_cv_t<std::invoke_result_t<Evaluate&, std::size_t>>;
  static_assert(std::is_object_v<Result>, "grid evaluation must return a value, not void or a reference");
  auto& scheduler = execution.resolved_scheduler();
  for (std::size_t offset = 0; offset < count;)
  {
    auto const size = std::min(execution.batch_size, count - offset);
    std::vector<std::optional<Result>> results(size);
    {
      [[maybe_unused]] auto scope = std::invoke(observe_batch, offset, size);
      scheduler.execute_batch(size, [&](std::size_t j) { results[j].emplace(std::invoke(evaluate, offset + j)); });
    }
    for (std::size_t j = 0; j < size; ++j)
    {
      std::invoke(collect, offset + j, std::move(*results[j]));
      results[j].reset();
    }
    offset += size;
  }
}
} // namespace bethe
