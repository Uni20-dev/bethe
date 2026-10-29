// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include "xyz_dispersion_mpfr_references.hpp"
#include <atomic>
#include <bethe/prepared_observable_grid.hpp>
#include <bethe/xyz.hpp>
#include <bethe/xyz_dispersion.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using P = uni20::Precision;
using Band = bethe::xyz::SpinonDispersion<R>;
namespace math = uni20::math;

TEST(XYZDispersionMpfr, IndependentJacobiAccuracyBeyondBinary128)
{
  auto const hp = P::bits(768);
  R previous(1, hp);
  auto const ambient = mpfr_get_default_prec(), emin = mpfr_get_emin(), emax = mpfr_get_emax();
  for (int bits : {80, 160, 256, 333})
  {
    auto p = P::bits(bits);
    R worst(0, hp);
    auto check = [&](R const& actual, char const* expected) {
      EXPECT_EQ(actual.precision(), p);
      R ref(expected, hp);
      worst = std::max(worst, math::abs(actual.at(hp) - ref) / (1 + math::abs(ref)));
    };
    for (auto const& f : xyz_mpfr_reference::fixtures)
    {
      SCOPED_TRACE(f.eta);
      Band band(R(f.eta, p), R(f.t, p));
      auto j = bethe::xyz::couplings(R(f.eta, p), R(f.t, p));
      check(j.x, f.x);
      check(j.y, f.y);
      check(j.z, f.z);
      check(band.maximum_energy(), f.maximum);
      check(band.gap(), f.gap);
      for (auto const& b : f.bounds)
      {
        auto branch = band.prepare_bound(b.s);
        check(branch.gap(), b.gap);
        check(branch.energy(R("0.7", p)), b.energy);
      }
    }
    EXPECT_LT(worst, 1024 * uni20::numeric_limits<R>::epsilon(p).at(hp));
    EXPECT_LT(worst, previous * R("1e-15", hp));
    if (bits >= 256)
    {
      EXPECT_LT(worst, R("1e-70", hp));
    }
    previous = worst;
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
  EXPECT_EQ(mpfr_get_emin(), emin);
  EXPECT_EQ(mpfr_get_emax(), emax);
}

TEST(XYZDispersionMpfr, EndpointsThresholdsAndExponentRange)
{
  for (int bits : {80, 160, 256})
  {
    auto p = P::bits(bits);
    R zero(0, p), one(1, p), pi = bethe::detail::spinon_pi(zero);
    Band band(R("0.75", p), one), scaled(R("0.75", p), one, R(3, p));
    EXPECT_EQ(band.energy(zero), band.gap());
    EXPECT_EQ(band.energy(pi), band.gap());
    EXPECT_EQ(band.energy(pi / 2), band.maximum_energy());
    auto edges = band.continuum(zero, bethe::SpinonMomentum::folded);
    EXPECT_EQ(edges.lower, 2 * band.gap());
    EXPECT_EQ(edges.upper, 2 * band.maximum_energy());
    EXPECT_EQ(edges.lower.precision(), p);
    EXPECT_EQ(edges.upper.precision(), p);
    EXPECT_TRUE(band.bound_exists(2));
    EXPECT_FALSE(band.bound_exists(3));
    for (unsigned s : {1u, 2u})
    {
      auto branch = band.prepare_bound(s);
      EXPECT_EQ(branch.energy(zero), branch.gap());
      EXPECT_EQ(branch.energy(2 * pi), branch.gap());
      EXPECT_EQ(branch.x_parity(), s % 2 ? -1 : 1);
      for (int i = 0; i <= 8; ++i)
      {
        R q = pi * R(i, p) / 4;
        auto value = branch.energy(q);
        EXPECT_EQ(value.precision(), p);
        EXPECT_LT(math::abs(scaled.bound_energy(s, q) - 3 * value), 4096 * uni20::numeric_limits<R>::epsilon(p));
      }
    }
    // The finite precision of eta decides a strict merger; do not substitute a fuzzy threshold.
    Band above(R("0.5", p) + uni20::numeric_limits<R>::epsilon(p), one);
    EXPECT_TRUE(above.bound_exists(1));
    EXPECT_FALSE(Band(R("0.5", p), one).bound_exists(1));
    Band tiny(R("0.25", p), R(10000, p));
    EXPECT_GT(tiny.gap(), zero);
    EXPECT_LT(tiny.gap(), R("1e-10000", p));
    EXPECT_EQ(tiny.energy(pi), tiny.gap());
    EXPECT_EQ(tiny.gap().precision(), p);
  }
}

TEST(XYZDispersionMpfr, InvalidPrecisionAndQueries)
{
  auto p = P::bits(160);
  R eta("0.75", p), one(1, p);
  Band band(eta, one);
  auto branch = band.prepare_bound(1);
  for (R const& bad : {R{}, R(1), R(1, P::bits(80)), R("nan", p)})
  {
    EXPECT_THROW((void)Band(bad, one), std::invalid_argument);
    EXPECT_THROW((void)Band(eta, bad), std::invalid_argument);
    EXPECT_THROW(band.energy(bad), std::invalid_argument);
    EXPECT_THROW(band.continuum(bad), std::invalid_argument);
    EXPECT_THROW(branch.energy(bad), std::invalid_argument);
    EXPECT_THROW(bethe::xyz::couplings(bad, one), std::invalid_argument);
  }
  EXPECT_THROW((void)Band(eta, one, R(1, P::bits(80))), std::invalid_argument);
  EXPECT_THROW(band.prepare_bound(0), std::invalid_argument);
  EXPECT_THROW(band.prepare_bound(3), std::invalid_argument);
  EXPECT_THROW(branch.energy(R(-1, p)), std::invalid_argument);
  EXPECT_THROW(branch.energy(R(10, p)), std::invalid_argument);
}

TEST(XYZDispersionMpfr, OwningMixedPrecisionAsyncReadersAndFailures)
{
  constexpr std::size_t count = 65;
  auto ambient = mpfr_get_default_prec();
  auto check = [&](auto& scheduler) {
    using Branch = Band::BoundBranch;
    std::array<std::array<std::vector<R>, 2>, 3> rows;
    std::vector<uni20::async::Async<std::size_t>> results, failures;
    std::atomic<unsigned> preparations{0};
    for (std::size_t b = 0; b < 3; ++b)
    {
      auto p = P::bits(80 + 80 * b);
      auto background = bethe::prepare_background_async(
          [p, &preparations] {
            ++preparations;
            return Band(R("0.875", p), R("0.75", p)).prepare_bound(2);
          },
          {&scheduler});
      R const pi = bethe::detail::spinon_pi(R(0, p));
      for (unsigned c = 0; c < 2; ++c)
        results.push_back(bethe::evaluate_observable_grid_async(
            background.read(), count,
            [=](Branch const& branch, std::size_t i) {
              return branch.energy(2 * pi * R(c ? count - 1 - i : i, p) / R(count - 1, p));
            },
            [&, b, c](std::size_t i, R value) {
              EXPECT_EQ(i, rows[b][c].size());
              rows[b][c].push_back(std::move(value));
            },
            {&scheduler, 8, 3}));
      failures.push_back(bethe::evaluate_observable_grid_async(
          background.read(), 1,
          [=](Branch const& branch, std::size_t) { return branch.energy(R(1, P::bits(p.bit_count() + 1))); },
          [](std::size_t, R) { FAIL() << "invalid precision must not publish a value"; }, {&scheduler, 1, 1}));
    }
    // The submitting background handles have gone; wait on each owning operation.
    for (auto& failure : failures)
      EXPECT_THROW((void)failure.get_wait(scheduler), std::invalid_argument);
    for (auto& result : results)
      EXPECT_EQ(result.get_wait(scheduler), count);
    EXPECT_EQ(preparations, 3u);
    for (std::size_t b = 0; b < 3; ++b)
    {
      auto p = P::bits(80 + 80 * b);
      auto branch = Band(R("0.875", p), R("0.75", p)).prepare_bound(2);
      R pi = bethe::detail::spinon_pi(R(0, p));
      for (unsigned c = 0; c < 2; ++c)
        for (std::size_t i = 0; i < count; ++i)
        {
          EXPECT_EQ(rows[b][c][i].precision(), p);
          EXPECT_EQ(rows[b][c][i], branch.energy(2 * pi * R(c ? count - 1 - i : i, p) / R(count - 1, p)));
        }
    }
  };
  uni20::async::DebugScheduler debug({.order = uni20::async::DebugSchedulerOrder::reverse});
  check(debug);
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    check(scheduler);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
} // namespace
