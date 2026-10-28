// SPDX-License-Identifier: GPL-3.0-or-later
#include "../apps/generated-state-scan.hpp"
#include <chrono>
#include <future>
#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>

TEST(GeneratedStateScan, OwnsGeneratedRequestsAndPreservesEnumeration)
{
  for (int threads : {1, 2, 4})
  {
    // The generator reuses mutable storage; each returned request is independent.
    auto states = bethe::cli::scan_generated_states(
        257,
        [buffer = std::vector<int>{-1, 7}]() mutable {
          ++buffer[0];
          return std::make_unique<std::vector<int>>(buffer);
        },
        [scale = std::make_unique<int>(3)](auto request) {
          return std::pair{request->front() * *scale, request->back()};
        },
        threads);
    ASSERT_EQ(states.size(), 257u);
    for (int i = 0; i != 257; ++i)
      EXPECT_EQ(states[i], (std::pair{i * 3, 7}));
  }
}

TEST(GeneratedStateScan, ExceptionCannotReturnWhileAnAdmittedSolveIsHeldOpen)
{
  using namespace std::chrono_literals;
  std::promise<void> entered, throwing, release;
  auto entered_future = entered.get_future().share();
  auto throwing_future = throwing.get_future();
  auto release_future = release.get_future().share();
  auto operation = std::async(std::launch::async, [&] {
    return bethe::cli::scan_generated_states(
        20, [i = 0]() mutable { return i++; },
        [&](int i) {
          if (i == 0)
          {
            entered.set_value();
            release_future.wait();
          }
          else
          {
            if (entered_future.wait_for(10s) != std::future_status::ready)
              throw std::runtime_error("first solve did not start");
            throwing.set_value();
            throw std::runtime_error("solve failed");
          }
          return i;
        },
        2);
  });
  // Release on every path, including a failed expectation; no global drain can
  // conceal early completion of the operation under test.
  EXPECT_EQ(throwing_future.wait_for(10s), std::future_status::ready);
  EXPECT_EQ(operation.wait_for(50ms), std::future_status::timeout);
  release.set_value();
  EXPECT_THROW(operation.get(), std::runtime_error);
}

TEST(GeneratedStateScan, OutOfOrderCompletionKeepsTheOriginalIds)
{
  using namespace std::chrono_literals;
  std::promise<void> second;
  auto second_ready = second.get_future().share();
  auto states = bethe::cli::scan_generated_states(
      64, [i = 0]() mutable { return i++; },
      [&](int i) {
        if (i == 0) EXPECT_EQ(second_ready.wait_for(10s), std::future_status::ready);
        if (i == 1) second.set_value();
        return i;
      },
      2);
  for (int i = 0; i != 64; ++i)
    EXPECT_EQ(states[i], i);
}

TEST(GeneratedStateScan, EmptyAndSingletonScansAndGeneratorFailures)
{
  EXPECT_TRUE(bethe::cli::scan_generated_states(
                  0,
                  [] {
                    throw std::logic_error("unused");
                    return 0;
                  },
                  [](int i) { return i; }, 4)
                  .empty());
  EXPECT_EQ(bethe::cli::scan_generated_states(1, [] { return 7; }, [](int i) { return i; }, 64), std::vector{7});
  for (int threads : {1, 4})
    EXPECT_THROW(bethe::cli::scan_generated_states(
                     20,
                     [i = 0]() mutable {
                       if (i == 3) throw std::runtime_error("labels failed");
                       return i++;
                     },
                     [](int i) { return i; }, threads),
                 std::runtime_error);
  EXPECT_THROW(
      bethe::cli::scan_generated_states(2, [] { return 0; }, [](int i) { return i; }, 0), std::invalid_argument);
}
