// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <cstddef>
#include <stdexcept>
#include <uni20/async/debug_scheduler.hpp>

namespace bethe
{
/// Synchronous, bounded batches. The borrowed scheduler must outlive the call.
/// A null scheduler uses Uni20's active scheduler (serial by default). Bethe
/// never changes the global scheduler. Callbacks must permit concurrent calls.
struct ExecutionOptions
{
    uni20::async::IAsyncScheduler* scheduler = nullptr;
    std::size_t batch_size = 128;

    uni20::async::IAsyncScheduler& resolved_scheduler() const
    {
      if (batch_size == 0) throw std::invalid_argument("execution batch_size must be positive");
      return *(scheduler ? scheduler : uni20::async::get_global_scheduler());
    }
};
} // namespace bethe
