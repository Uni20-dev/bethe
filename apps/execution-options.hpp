// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include "program-options.hpp"
#include <bethe/execution.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace bethe::cli
{
inline CLI::Option* threads_option(CLI::App& app, int& threads)
{
  return app.add_option("--threads", threads, "Scheduler concurrency limit for independent states (1 = serial)")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();
}
} // namespace bethe::cli
