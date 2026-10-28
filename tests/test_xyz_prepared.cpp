// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include <atomic>
#include <bethe/prepared_observable_grid.hpp>
#include <bethe/xyz_dispersion.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
template <typename Real> class XYZPrepared : public ::testing::Test {};
TYPED_TEST_SUITE(XYZPrepared, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XYZPrepared, OwningBranchAndUnchangedQueries)
{
  using R = TypeParam;
  using Band = bethe::xyz::SpinonDispersion<R>;
  R const pi = bethe::detail::pi<R>();
  for (R eta : {R{0.75}, R{0.875}})
    for (R t : {R{0.25}, R{1}, R{4}})
    {
      Band original(eta, t);
      for (unsigned s = 1; original.bound_exists(s); ++s)
      {
        // Returned branch owns all needed data after its temporary parent dies.
        auto prepared = Band(eta, t).prepare_bound(s);
        EXPECT_EQ(prepared.index(), s);
        EXPECT_EQ(prepared.x_parity(), original.bound_x_parity(s));
        EXPECT_EQ(prepared.gap(), original.bound_gap(s));
        auto copy = prepared;
        for (unsigned i = 0; i <= 32; ++i)
        {
          R const q = 2 * pi * R(i) / R{32};
          EXPECT_EQ(prepared.energy(q), original.bound_energy(s, q));
          EXPECT_EQ(copy.energy(q), prepared.energy(q));
        }
        for (R invalid : {-R{1}, 3 * pi, uni20::numeric_limits<R>::quiet_NaN(), uni20::numeric_limits<R>::infinity()})
          EXPECT_THROW(prepared.energy(invalid), std::invalid_argument);
      }
    }
  Band band(R{0.75}, R{1});
  EXPECT_THROW(band.prepare_bound(0), std::invalid_argument);
  EXPECT_THROW(band.prepare_bound(3), std::invalid_argument); // strict merger
}

TYPED_TEST(XYZPrepared, SharedAsyncPreparationAndIndependentFailures)
{
  using R = TypeParam;
  using Band = bethe::xyz::SpinonDispersion<R>;
  using Branch = typename Band::BoundBranch;
  R const pi = bethe::detail::pi<R>();
  constexpr std::size_t count = 129;
  auto query = [pi](Branch const& branch, std::size_t i, unsigned copy) {
    R q = 2 * pi * R(i) / R(count - 1);
    if (copy) q = q < pi ? q + pi : q - pi;
    return branch.energy(q);
  };
  std::array<std::vector<R>, 2> reference;
  auto bound = Band(R{0.75}, R{1}).prepare_bound(1);
  for (unsigned b = 0; b < 2; ++b)
    for (std::size_t i = 0; i < count; ++i)
      reference[b].push_back(query(bound, i, b));
  auto check = [&](auto& scheduler, bool reject_preparation) {
    std::atomic<unsigned> prepared_count{0};
    std::array<std::vector<R>, 2> rows;
    std::vector<R> failed;
    auto results = [&] {
      auto background = bethe::prepare_background_async(
          [&] {
            ++prepared_count;
            return Band(R{0.75}, R{1}).prepare_bound(reject_preparation ? 3 : 1);
          },
          {&scheduler});
      auto submit = [&](unsigned b) {
        return bethe::evaluate_observable_grid_async(
            background.read(), count, [=](Branch const& branch, std::size_t i) { return query(branch, i, b); },
            [&, b](std::size_t i, R value) {
              EXPECT_EQ(i, rows[b].size());
              rows[b].push_back(value);
            },
            {&scheduler, 128, 3});
      };
      auto invalid = bethe::evaluate_observable_grid_async(
          background.read(), count,
          [=](Branch const& branch, std::size_t i) {
            return i == count - 1 ? branch.energy(-R{1}) : query(branch, i, 0);
          },
          [&](std::size_t i, R value) {
            EXPECT_EQ(i, failed.size());
            failed.push_back(value);
          },
          {&scheduler, 128, 2});
      return std::tuple{submit(0), submit(1), std::move(invalid)};
    }(); // No preparation handle or parent band remains in the submitting scope.
    EXPECT_THROW((void)std::get<2>(results).get_wait(scheduler), std::invalid_argument);
    if (reject_preparation)
    {
      EXPECT_THROW((void)std::get<0>(results).get_wait(scheduler), std::invalid_argument);
      EXPECT_THROW((void)std::get<1>(results).get_wait(scheduler), std::invalid_argument);
      EXPECT_TRUE(rows[0].empty() && rows[1].empty() && failed.empty());
    }
    else
    {
      EXPECT_EQ(std::get<0>(results).get_wait(scheduler), count);
      EXPECT_EQ(std::get<1>(results).get_wait(scheduler), count);
      EXPECT_EQ(rows, reference);
    }
    EXPECT_EQ(prepared_count, 1u);
  };
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  check(deferred, false);
  check(deferred, true);
  for (int workers : {1, 4})
  {
    uni20::async::TbbScheduler scheduler(workers);
    check(scheduler, false);
    check(scheduler, true);
  }
}
} // namespace
