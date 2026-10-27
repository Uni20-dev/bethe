// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <cstddef>
#include <stdexcept>
#include <stop_token>
#include <uni20/async/debug_scheduler.hpp>

namespace bethe
{
/// Execution resources. The borrowed scheduler must outlive all submitted work.
/// A null scheduler uses Uni20's active scheduler (serial by default). Bethe
/// never changes the global scheduler. Callbacks must permit concurrent calls.
struct ExecutionOptions
{
    uni20::async::IAsyncScheduler* scheduler = nullptr;
    std::size_t batch_size = 128;
    /// Independent bound for continuous work; not a scheduler thread count.
    /// Ordered grids include completed-but-not-collected results in this bound.
    /// batch_size remains the granularity of clients that still use batches.
    std::size_t max_in_flight = 128;
    /// Stop admitting work and finish active kernels before reporting cancellation.
    /// Dropping an Async handle does not request cancellation.
    std::stop_token stop = {};

    uni20::async::IAsyncScheduler& resolved_scheduler() const
    {
      if (batch_size == 0) throw std::invalid_argument("execution batch_size must be positive");
      return *(scheduler ? scheduler : uni20::async::get_global_scheduler());
    }

    uni20::async::IAsyncScheduler& continuous_scheduler() const
    {
      if (!max_in_flight) throw std::invalid_argument("execution max_in_flight must be positive");
      return resolved_scheduler();
    }
};

class ExecutionCancelled : public std::runtime_error {
  public:
    explicit ExecutionCancelled(char const* message = "Bethe operation cancelled") : std::runtime_error(message) {}
};

/// No partial scan is returned as a successful/completed family on cancellation.
class ScanCancelled : public ExecutionCancelled {
  public:
    ScanCancelled() : ExecutionCancelled("real-root scan cancelled") {}
};
} // namespace bethe
