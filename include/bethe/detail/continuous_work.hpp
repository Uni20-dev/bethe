// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <algorithm>
#include <atomic>
#include <bethe/detail/async_execution.hpp>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <stop_token>

namespace bethe::detail
{
/// A fixed number of scheduler-owned runners pull indices continuously. There
/// are no intermediate joins: a slow index cannot hold up the next logical
/// batch. The callback owns loading, computing and publishing ONE result; no
/// completed-result queue or per-input scheduler task is retained here.
///
/// Exceptions stop admission, allow other active callbacks to checkpoint using
/// the shared stop token, then propagate after every runner has joined. The
/// callback must cooperate with cancellation; it is never forcibly interrupted.
template <typename Function>
std::size_t continuous_work(uni20::async::IAsyncScheduler& scheduler, std::size_t workers, std::size_t count,
                            Function&& function, std::stop_token external_stop = {})
{
  if (!workers) throw std::invalid_argument("continuous work requires a positive runner count");
  ExecutionControl control(external_stop);
  std::atomic<std::size_t> next{0}, completed{0};
  scheduler.execute_batch(std::min(workers, count), [&](std::size_t) {
    try
    {
      while (!control.stopped())
      {
        auto i = next.fetch_add(1, std::memory_order_relaxed);
        if (i >= count || control.stopped()) break;
        function(i, control.token());
        completed.fetch_add(1, std::memory_order_relaxed);
      }
    }
    catch (...)
    {
      control.fail();
    }
  });
  control.rethrow_error();
  return completed.load(std::memory_order_relaxed);
}
} // namespace bethe::detail
