// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/tb_dispersion.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using Band = bethe::takhtajan_babujian::SpinonDispersion<R>;
using bethe::SpinonMomentum;
class TBDispersionMpfr : public ::testing::TestWithParam<int> {
  protected:
    uni20::Precision p() const { return uni20::Precision::bits(GetParam()); }
    R v(int n) const { return R(n, p()); }
};
INSTANTIATE_TEST_SUITE_P(WorkingPrecisions, TBDispersionMpfr, ::testing::Values(128, 256, 333));
void finite(R const& x, uni20::Precision p)
{
  EXPECT_FALSE(x.is_exact());
  EXPECT_EQ(x.precision(), p);
  EXPECT_TRUE(uni20::isfinite(x));
}
std::pair<R, R> sea_reference(R lambda)
{
  R pi = uni20::pi<R>.at(lambda.precision());
  // Dressed two-string hole momentum/energy in rapidity space.
  return {2 * uni20::math::atan(uni20::math::exp(-pi * lambda)), 2 * pi / uni20::math::cosh(pi * lambda)};
}

TEST_P(TBDispersionMpfr, IndependentSeaAndTypedEndpoints)
{
  Band band(v(1));
  auto const eps = uni20::numeric_limits<R>::epsilon(p());
  R pi = uni20::pi<R>.at(p());
  EXPECT_EQ(band.spin(), uni20::from_twice(1));
  finite(band.gap(), p());
  EXPECT_EQ(band.gap(), 0);
  for (R lambda : {-v(2), -v(1) / 2, v(0), v(1) / 2, v(2)})
  {
    auto [momentum, energy] = sea_reference(lambda);
    auto actual = band.energy(momentum);
    finite(actual, p());
    EXPECT_REAL_NEAR(actual, energy, 64 * eps);
  }
  for (R q : {v(0), pi})
  {
    auto energy = band.energy(q);
    finite(energy, p());
    EXPECT_EQ(energy, 0);
  }
  EXPECT_EQ(band.energy(pi / 2), band.maximum_energy());
  R tiny("1e-7000", p());
  for (auto convention : {SpinonMomentum::unfolded, SpinonMomentum::folded})
  {
    auto two = band.continuum(tiny, convention), four = band.four_spinon_continuum(tiny, convention);
    finite(two.lower, p());
    finite(two.upper, p());
    finite(four.lower, p());
    finite(four.upper, p());
    EXPECT_GT(two.lower, 0);
    EXPECT_EQ(two.lower, four.lower);
    EXPECT_REAL_NEAR(two.lower / tiny, band.velocity(), 32 * eps);
    if (convention == SpinonMomentum::unfolded) EXPECT_REAL_NEAR(two.upper / tiny, band.velocity(), 32 * eps);
  }
  Band huge(R("1e7000", p()));
  EXPECT_REAL_NEAR(huge.energy(tiny), band.velocity(), 64 * eps);
  finite(huge.four_spinon_continuum(v(0)).upper, p());
}

TEST_P(TBDispersionMpfr, PreserveMomentumWhenHalvingUnderflows)
{
  // Include nonzero underflow: between min and 2*min, q/2 can round UP to
  // min, so merely checking q/2==0 would corrupt the upper edge by O(1).
  // Construct values at the current exponent range without changing it.
  Band band(v(1));
  for (unsigned numerator : {4u, 5u, 6u, 7u, 8u, 12u})
  {
    SCOPED_TRACE(numerator);
    mpfr_t raw;
    mpfr_init2(raw, p().bit_count());
    mpfr_set_ui_2exp(raw, numerator, mpfr_get_emin() - 3, MPFR_RNDN);
    R q(raw, p());
    mpfr_clear(raw);
    ASSERT_GT(q, 0);
    if (numerator == 4)
      EXPECT_EQ(q / 2, 0);
    else
      EXPECT_GT(q / 2, 0);
    EXPECT_EQ((q / 2) * 2 != q, numerator < 8);
    R expected = band.velocity() * q;
    ASSERT_GT(expected, 0);
    EXPECT_EQ(band.energy(q), expected);
    auto two = band.continuum(q);
    EXPECT_EQ(two.lower, expected);
    EXPECT_EQ(two.upper, expected);
    EXPECT_EQ(band.continuum(q, SpinonMomentum::folded).lower, expected);
    EXPECT_EQ(band.four_spinon_continuum(q, SpinonMomentum::folded).lower, expected);
  }
}

TEST_P(TBDispersionMpfr, ContinuumExtremaAndFoldedUnions)
{
  Band band(v(1)), scaled(v(7));
  R pi = uni20::pi<R>.at(p()), eps = uni20::numeric_limits<R>::epsilon(p());
  for (unsigned j = 0; j <= 12; ++j)
  {
    R q = 2 * pi * (R(j) / R(12)), reduced = std::min(q, 2 * pi - q);
    auto two = band.continuum(q), four = band.four_spinon_continuum(q);
    // Equal sharing attains each upper bound. Two soft spectators attain
    // the four-particle minimum; the maximum uses the alias nearest 2*pi.
    EXPECT_REAL_NEAR(two.upper, 2 * band.energy(q / 2), 64 * eps);
    EXPECT_REAL_NEAR(four.upper, 4 * band.energy((2 * pi + reduced) / 4), 128 * eps);
    R lo = std::max(v(0), q - pi), hi = std::min(pi, q);
    for (unsigned i = 0; i <= 32; ++i)
    {
      R a = lo + (hi - lo) * (R(i) / R(32));
      R b = std::clamp(q - a, v(0), pi);
      R sample = band.energy(a) + band.energy(b);
      EXPECT_GE(sample + 128 * eps, two.lower);
      EXPECT_LE(sample - 128 * eps, two.upper);
      if (i == 0 || i == 32) EXPECT_REAL_NEAR(sample, two.lower, 128 * eps);
    }
    EXPECT_EQ(four.lower, two.lower);
    R shifted = q < pi ? q + pi : q - pi;
    for (bool four_particles : {false, true})
    {
      auto get = [&](Band const& x, R k, SpinonMomentum c) {
        return four_particles ? x.four_spinon_continuum(k, c) : x.continuum(k, c);
      };
      auto a = get(band, q, SpinonMomentum::unfolded), b = get(band, shifted, SpinonMomentum::unfolded);
      auto folded = get(band, q, SpinonMomentum::folded), large = get(scaled, q, SpinonMomentum::folded);
      finite(folded.lower, p());
      finite(folded.upper, p());
      EXPECT_REAL_NEAR(folded.lower, std::min(a.lower, b.lower), 128 * eps);
      EXPECT_REAL_NEAR(folded.upper, std::max(a.upper, b.upper), 128 * eps);
      EXPECT_REAL_NEAR(large.lower, 7 * folded.lower, 512 * eps);
      EXPECT_REAL_NEAR(large.upper, 7 * folded.upper, 512 * eps);
    }
  }
  for (R q : {v(0), pi, 2 * pi})
  {
    auto two = band.continuum(q), four = band.four_spinon_continuum(q);
    EXPECT_EQ(two.lower, 0);
    EXPECT_EQ(four.lower, 0);
    finite(two.lower, p());
    finite(four.lower, p());
  }
}

TEST_P(TBDispersionMpfr, SharedMassiveBandAndPrecisionContracts)
{
  R pi = uni20::pi<R>.at(p()), eps = uni20::numeric_limits<R>::epsilon(p());
  using Shared = bethe::detail::SpinonBand<R>;
  for (R gap : {v(0), R("1e-80", p()), v(1) / 5, v(1)})
  {
    Shared band(v(1), gap);
    EXPECT_EQ(band.energy(v(0)), gap);
    EXPECT_EQ(band.energy(pi), gap);
    EXPECT_EQ(band.energy(pi / 2), 1);
    for (R q : {pi / 8, pi / 3, 3 * pi / 4})
    {
      auto bounds = band.continuum(q);
      finite(bounds.lower, p());
      finite(bounds.upper, p());
      for (unsigned i = 0; i <= 64; ++i)
      {
        R a = q * (R(i) / R(64)), b = q - a;
        // Independent squared dispersion, avoiding the production hypot.
        auto energy = [&](R k) {
          R c = uni20::math::cos(k), s = uni20::math::sin(k);
          return uni20::math::sqrt(gap * gap * c * c + s * s);
        };
        R sample = energy(a) + energy(b);
        EXPECT_GE(sample + 64 * eps, bounds.lower);
        EXPECT_LE(sample - 64 * eps, bounds.upper);
      }
    }
  }
  Band band(v(1));
  static_assert(!std::is_default_constructible_v<Band>);
  for (R bad : {R(uni20::uninitialized), R(1), R(1, uni20::Precision::bits(GetParam() + 1))})
  {
    EXPECT_THROW(band.energy(bad), std::invalid_argument);
    EXPECT_THROW(band.continuum(bad), std::invalid_argument);
    EXPECT_THROW(band.four_spinon_continuum(bad), std::invalid_argument);
    EXPECT_THROW((Shared{v(1), bad}), std::invalid_argument);
  }
  EXPECT_THROW((Band{R(1)}), std::invalid_argument);
  EXPECT_THROW((Band{R(uni20::uninitialized)}), std::invalid_argument);
  for (R bad : {v(0), -v(1), R("nan", p()), R("inf", p())})
    EXPECT_THROW((Band{bad}), std::invalid_argument);
  EXPECT_THROW((Shared{v(1), v(2)}), std::invalid_argument);
  EXPECT_THROW(band.energy(-v(1)), std::invalid_argument);
  EXPECT_THROW(band.energy(v(4)), std::invalid_argument);
  EXPECT_THROW(band.continuum(v(7)), std::invalid_argument);
  EXPECT_THROW(band.four_spinon_continuum(R("nan", p())), std::invalid_argument);
  EXPECT_THROW(band.continuum(v(0), static_cast<SpinonMomentum>(7)), std::invalid_argument);
}

TEST(TBDispersionMpfrShared, AccuracyGrowsBeyondBinary128)
{
  auto reference_precision = uni20::Precision::bits(512);
  auto [momentum, reference] = sea_reference(R(1, reference_precision) / 7);
  R previous(1, reference_precision);
  for (int bits : {128, 256, 333})
  {
    auto p = uni20::Precision::bits(bits);
    Band band(R(1, p));
    R error = uni20::math::abs(band.energy(momentum.at(p)).at(reference_precision) - reference);
    EXPECT_LT(error, previous / R("1e15", reference_precision));
    EXPECT_LT(error, (128 * uni20::numeric_limits<R>::epsilon(p)).at(reference_precision));
    previous = error;
  }
}

TEST(TBDispersionMpfrShared, OwningMixedPrecisionQueriesAndFailureIsolation)
{
  struct Point
  {
      R energy;
      bethe::SpinonContinuum<R> two, four;
  };
  auto const ambient = mpfr_get_default_prec();
  auto submit = [](uni20::async::IAsyncScheduler& scheduler, std::vector<Point>& rows, bool fail) {
    struct Request
    {
        Band band;
        R momentum;
    };
    std::vector<Request> requests;
    for (int bits : {128, 256, 333})
    {
      auto p = uni20::Precision::bits(bits);
      requests.push_back({Band(R(7, p)), R("0.3", p)});
    }
    if (fail) requests.back().momentum = R(1);
    return bethe::evaluate_observable_grid_async(
        requests.size(),
        [requests](std::size_t i) {
          auto const& r = requests[i];
          return Point{r.band.energy(r.momentum), r.band.continuum(r.momentum),
                       r.band.four_spinon_continuum(r.momentum, SpinonMomentum::folded)};
        },
        [&rows](std::size_t i, Point value) {
          EXPECT_EQ(i, rows.size());
          rows.push_back(std::move(value));
        },
        {&scheduler, 128, 2});
  };
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  std::vector<Point> reference;
  auto baseline = submit(deferred, reference, false);
  EXPECT_EQ(baseline.get_wait(deferred), 3u);
  for (int workers : {1, 4})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<Point> actual, failed;
    auto invalid = submit(scheduler, failed, true);
    auto complete = submit(scheduler, actual, false);
    EXPECT_THROW((void)invalid.get_wait(scheduler), std::invalid_argument);
    EXPECT_EQ(complete.get_wait(scheduler), 3u);
    ASSERT_EQ(actual.size(), reference.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      auto const& a = actual[i];
      auto const& b = reference[i];
      auto p = b.energy.precision();
      for (auto x : {a.energy, a.two.lower, a.two.upper, a.four.lower, a.four.upper})
        finite(x, p);
      EXPECT_EQ(a.energy, b.energy);
      EXPECT_EQ(a.two.lower, b.two.lower);
      EXPECT_EQ(a.two.upper, b.two.upper);
      EXPECT_EQ(a.four.lower, b.four.lower);
      EXPECT_EQ(a.four.upper, b.four.upper);
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
} // namespace
