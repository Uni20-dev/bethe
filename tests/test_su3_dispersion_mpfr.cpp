// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/su3_dispersion.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using Band = bethe::su3::ExcitationDispersion<R>;
using bethe::su3::Continuum;
using bethe::su3::Momentum;
using bethe::su3::Particle;
class SU3DispersionMpfr : public ::testing::TestWithParam<int> {
  protected:
    uni20::Precision p() const { return uni20::Precision::bits(GetParam()); }
    R v(int n) const { return R(n, p()); }
};
INSTANTIATE_TEST_SUITE_P(WorkingPrecisions, SU3DispersionMpfr, ::testing::Values(128, 256, 333));

void finite(R const& x, uni20::Precision p)
{
  EXPECT_FALSE(x.is_exact());
  EXPECT_EQ(x.precision(), p);
  EXPECT_TRUE(uni20::isfinite(x));
}

// Independent nested-sea density and its integrated hole momentum, rather
// than the implementation's trigonometric dispersion in momentum space.
std::pair<R, R> sea_reference(unsigned sea, R lambda)
{
  auto p = lambda.precision();
  R pi = uni20::pi<R>.at(p), angle = sea * pi / 3;
  R density = uni20::math::sin(angle) / (3 * (uni20::math::cosh(2 * pi * lambda / 3) - uni20::math::cos(angle)));
  R momentum = pi - angle - 2 * uni20::math::atan(uni20::math::tanh(pi * lambda / 3) / uni20::math::tan(angle / 2));
  return {momentum, 2 * pi * density};
}

TEST_P(SU3DispersionMpfr, IndependentSeasEndpointsAndTinyEnergies)
{
  Band band(v(1));
  R const eps = uni20::numeric_limits<R>::epsilon(p());
  for (unsigned sea : {1u, 2u})
  {
    auto particle = sea == 1 ? Particle::fundamental : Particle::antifundamental;
    for (R lambda : {-v(2), -v(1) / 4, v(0), v(1) / 4, v(2)})
    {
      auto [momentum, expected] = sea_reference(sea, lambda);
      auto energy = band.energy(particle, momentum);
      finite(energy, p());
      EXPECT_REAL_NEAR(energy, expected, 128 * eps);
    }
    R end = band.momentum_max(particle);
    for (R q : {v(0), end})
    {
      R zero = band.energy(particle, q);
      finite(zero, p());
      EXPECT_EQ(zero, 0);
    }
    R tiny("1e-7000", p());
    auto energy = band.energy(particle, tiny);
    finite(energy, p());
    EXPECT_GT(energy, 0);
    EXPECT_REAL_NEAR(energy / tiny, band.velocity(), 32 * eps);
    EXPECT_REAL_NEAR(band.energy(particle, end / 2), band.maximum_energy(particle), 32 * eps);
  }
  // MPFR's exponent range must not be restricted to a native float's range.
  Band huge(R("1e7000", p()));
  finite(huge.velocity(), p());
  EXPECT_REAL_NEAR(huge.energy(Particle::fundamental, R("1e-7000", p())), band.velocity(), 64 * eps);
}

TEST_P(SU3DispersionMpfr, NonzeroUnderflowWhenHalvingMomentum)
{
  Band band(v(1));
  R const eps = uni20::numeric_limits<R>::epsilon(p());
  for (unsigned numerator : {4u, 5u, 6u, 7u, 8u, 12u})
  {
    SCOPED_TRACE(numerator);
    mpfr_t raw;
    mpfr_init2(raw, p().bit_count());
    mpfr_set_ui_2exp(raw, numerator, mpfr_get_emin() - 3, MPFR_RNDN);
    R q(raw, p());
    mpfr_clear(raw);
    ASSERT_GT(q, 0);
    EXPECT_EQ((q / 2) * 2 != q, numerator < 8);
    if (numerator > 4) EXPECT_GT(q / 2, 0);
    for (auto particle : {Particle::fundamental, Particle::antifundamental})
    {
      R energy = band.energy(particle, q);
      finite(energy, p());
      EXPECT_GT(energy, 0);
      EXPECT_REAL_NEAR(energy / q, band.velocity(), 32 * eps);
    }
    for (auto content : {Continuum::two_soliton, Continuum::four_soliton})
      for (auto convention : {Momentum::unfolded, Momentum::three_site_folded})
      {
        auto edge = band.continuum(q, content, convention);
        finite(edge.lower, p());
        EXPECT_GT(edge.lower, 0);
        EXPECT_REAL_NEAR(edge.lower / q, band.velocity(), 32 * eps);
      }
  }
}

TEST_P(SU3DispersionMpfr, ContinuumFoldingAndExchangeScaling)
{
  Band band(v(1)), scaled(v(7));
  R const pi = uni20::pi<R>.at(p()), eps = uni20::numeric_limits<R>::epsilon(p());
  for (auto content : {Continuum::two_soliton, Continuum::four_soliton})
    for (unsigned i = 0; i <= 19; ++i)
    {
      R q = 2 * pi * (R(i) / R(19));
      auto edge = band.continuum(q, content), large = scaled.continuum(q, content);
      finite(edge.lower, p());
      finite(edge.upper, p());
      EXPECT_GE(edge.lower, 0);
      EXPECT_LE(edge.lower, edge.upper);
      EXPECT_REAL_NEAR(large.lower, 7 * edge.lower, 512 * eps);
      EXPECT_REAL_NEAR(large.upper, 7 * edge.upper, 512 * eps);
      R low(100, p()), high(0, p());
      for (unsigned image = 0; image < 3; ++image)
      {
        R r = q + 2 * image * pi / 3;
        if (r > 2 * pi) r -= 2 * pi;
        auto e = band.continuum(r, content);
        low = std::min(low, e.lower);
        high = std::max(high, e.upper);
      }
      auto folded = band.continuum(q, content, Momentum::three_site_folded);
      finite(folded.lower, p());
      finite(folded.upper, p());
      EXPECT_REAL_NEAR(folded.lower, low, 128 * eps);
      EXPECT_REAL_NEAR(folded.upper, high, 128 * eps);
    }
  auto two = band.continuum(pi), four = band.continuum(pi, Continuum::four_soliton);
  EXPECT_REAL_NEAR(two.lower, 2 * four.lower, 32 * eps);
  for (R q : {v(0), 2 * pi / 3, 2 * pi})
  {
    auto edge = band.continuum(q);
    finite(edge.lower, p());
    EXPECT_EQ(edge.lower, 0);
  }
}

TEST_P(SU3DispersionMpfr, ExplicitPrecisionAndDomainChecks)
{
  static_assert(!std::is_default_constructible_v<Band>);
  Band band(v(1));
  for (R bad : {R(uni20::uninitialized), R(1), R(1, uni20::Precision::bits(GetParam() + 1))})
  {
    EXPECT_THROW(band.energy(Particle::fundamental, bad), std::invalid_argument);
    EXPECT_THROW(band.continuum(bad), std::invalid_argument);
  }
  EXPECT_THROW((Band{R(1)}), std::invalid_argument);
  EXPECT_THROW((Band{R(uni20::uninitialized)}), std::invalid_argument);
  EXPECT_THROW((Band{v(0)}), std::invalid_argument);
  EXPECT_THROW((Band{-v(1)}), std::invalid_argument);
  EXPECT_THROW((Band{R("nan", p())}), std::invalid_argument);
  EXPECT_THROW((Band{R("inf", p())}), std::invalid_argument);
  EXPECT_THROW(band.energy(Particle::antifundamental, v(3)), std::invalid_argument);
  EXPECT_THROW(band.continuum(v(7)), std::invalid_argument);
  EXPECT_THROW(band.continuum(R("nan", p())), std::invalid_argument);
}

TEST(SU3DispersionMpfrShared, AccuracyGrowsBeyondBinary128)
{
  auto p = uni20::Precision::bits(512);
  for (unsigned sea : {1u, 2u})
  {
    auto [momentum, reference] = sea_reference(sea, R(1, p) / 7);
    R previous(1, p);
    for (int bits : {128, 256, 333})
    {
      auto q = uni20::Precision::bits(bits);
      Band band(R(1, q));
      R value = band.energy(sea == 1 ? Particle::fundamental : Particle::antifundamental, momentum.at(q));
      R error = uni20::math::abs(value.at(p) - reference);
      EXPECT_LT(error, previous / R("1e15", p));
      EXPECT_LT(error, (128 * uni20::numeric_limits<R>::epsilon(q)).at(p));
      previous = error;
    }
  }
}

TEST(SU3DispersionMpfrShared, OwningMixedPrecisionQueriesAndFailureIsolation)
{
  auto const ambient = mpfr_get_default_prec();
  auto submit = [](uni20::async::IAsyncScheduler& scheduler, std::vector<R>& rows, bool fail) {
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
    if (fail) requests.back().momentum = R(1); // exact input is not a query precision
    return bethe::evaluate_observable_grid_async(
        requests.size(),
        [requests](std::size_t i) { return requests[i].band.energy(Particle::fundamental, requests[i].momentum); },
        [&rows](std::size_t i, R value) {
          EXPECT_EQ(i, rows.size());
          rows.push_back(std::move(value));
        },
        {&scheduler, 128, 2});
  };
  std::vector<R> reference;
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  auto serial = submit(deferred, reference, false);
  EXPECT_EQ(serial.get_wait(deferred), 3u);
  for (int workers : {1, 4})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<R> actual, failed;
    auto invalid = submit(scheduler, failed, true);
    auto complete = submit(scheduler, actual, false);
    EXPECT_THROW((void)invalid.get_wait(scheduler), std::invalid_argument);
    EXPECT_EQ(complete.get_wait(scheduler), 3u);
    ASSERT_EQ(actual.size(), reference.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      finite(actual[i], reference[i].precision());
      EXPECT_EQ(actual[i], reference[i]);
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
} // namespace
