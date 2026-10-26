// SPDX-License-Identifier: GPL-3.0-or-later
// Minimal dependency-only sanitizer reproducers; no Bethe headers or kernels.
#include <gtest/gtest.h>
#include <oneapi/tbb/task.h>
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/task_group.h>
#include <uni20/async/async.hpp>

TEST(DependencyOnly, Uni20DebugNodeLifetime) { uni20::async::Async<int> ready(7); }

TEST(DependencyOnly, TbbSuspendResume)
{
  oneapi::tbb::task_arena arena(2);
  for (int i = 0; i < 10; ++i)
    arena.execute([] {
      oneapi::tbb::task_group group;
      group.run(
          [] { oneapi::tbb::task::suspend([](oneapi::tbb::task::suspend_point p) { oneapi::tbb::task::resume(p); }); });
      group.wait();
    });
}
