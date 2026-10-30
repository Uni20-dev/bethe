// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/detail/xxx_four_spinon_geometry.hpp>
#include <bethe/xxx_spinon_density.hpp>

namespace
{
namespace math = uni20::math;
template <class R> using Patch = bethe::heisenberg::detail::FourSpinonPatch<R>;
template <class R> using Geometry = bethe::heisenberg::detail::FourSpinonGeometry<R>;

template <class R> void point_checks(R const& one)
{
  R const zero = one - one, pi = R{4} * math::atan(one), eps = uni20::numeric_limits<R>::epsilon(one);
  for (int quarter : {1, 2, 4, 5, 7})
    for (bool mixed : {false, true})
    {
      if (mixed && quarter <= 4) continue;
      for (int fraction : (mixed ? std::vector{1, 3, 7} : std::vector{1, 3}))
      {
        R const total = R(quarter) * pi / R{4};
        Patch<R> patch(total, R(fraction) * one / R{8}, mixed);
        R const r = patch.radius() / R{3};
        auto const bounds = patch.bounds(r);
        R const omega = (bounds.lower_energy + bounds.upper_energy) / R{2};
        auto const support = patch.support(omega);
        ASSERT_TRUE(support);
        EXPECT_LT(support->first, r);
        EXPECT_GT(support->second, r);
        auto const point = patch.spectral_point(math::log(r), omega);
        R momentum = zero, energy = zero;
        std::array<R, 4> p;
        for (int j = 0; j < 4; ++j)
        {
          p[j] = R{2} * math::atan(math::exp(point.beta[j]));
          momentum += p[j];
          energy += pi / (R{2} * math::cosh(point.beta[j]));
        }
        EXPECT_TRUE(std::is_sorted(p.begin(), p.end()));
        EXPECT_REAL_NEAR(momentum, total, R{32768} * eps);
        EXPECT_REAL_NEAR(energy, omega, R{32768} * eps);
        int const first = mixed ? 1 : 2;
        R const direct = pi / R{2} * (math::cos(p[first]) - math::cos(p[first + 1]));
        EXPECT_REAL_NEAR(math::exp(-point.log_jacobian), direct, R{32768} * eps);
        auto const moment = patch.moment_point(math::log(r), one / R{3});
        R direct_energy = zero;
        for (auto const& b : moment.beta)
          direct_energy += pi / (R{2} * math::cosh(b));
        EXPECT_REAL_NEAR(moment.energy, direct_energy, R{32768} * eps);
        if constexpr (uni20::has_runtime_precision_v<R>)
        {
          EXPECT_EQ(point.beta[0].precision(), one.precision());
          EXPECT_EQ(point.log_jacobian.precision(), one.precision());
        }
      }
    }
}

template <class R> class FourSpinonGeometry : public ::testing::Test {};
TYPED_TEST_SUITE(FourSpinonGeometry, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(FourSpinonGeometry, ConstraintsAndActualSupport) { point_checks(TypeParam{1}); }

TYPED_TEST(FourSpinonGeometry, SoftLogDistancesOutliveExponentRange)
{
  using R = TypeParam;
  R const pi = R{4} * math::atan(R{1}), log_r{-20000};
  for (bool mixed : {false, true})
  {
    Patch<R> patch(mixed ? R{3} * pi / R{2} : pi, R{1} / R{4}, mixed);
    auto const point = patch.spectral_point(log_r, R{2});
    EXPECT_REAL_NEAR(point.beta[0], log_r - math::log(R{8}),
                     R{8} * uni20::numeric_limits<R>::epsilon() * math::abs(log_r));
    R const second = log_r + math::log(R{3} / R{8});
    EXPECT_REAL_NEAR(point.beta[mixed ? 3 : 1], mixed ? -second : second,
                     R{8} * uni20::numeric_limits<R>::epsilon() * math::abs(log_r));
    for (auto const& b : point.beta)
      EXPECT_TRUE(uni20::isfinite(b));
  }
}

TYPED_TEST(FourSpinonGeometry, EnergyDeltaRecoversHardInterval)
{
  using R = TypeParam;
  R const pi = R{4} * math::atan(R{1});
  auto const rule = bethe::detail::gauss_legendre<R>(32);
  for (bool mixed : {false, true})
  {
    Patch<R> patch(R{7} * pi / R{4}, R{3} / R{8}, mixed);
    R const r = patch.radius() / R{3};
    auto const bounds = patch.bounds(r);
    R const width = bounds.upper_energy - bounds.lower_energy;
    R integral{0};
    for (std::size_t i = 0; i < rule.x.size(); ++i)
    {
      R const angle = (R{1} + rule.x[i]) * pi / R{4}, sine = math::sin(angle);
      R const omega = bounds.upper_energy - width * sine * sine;
      auto const point = patch.spectral_point(math::log(r), omega);
      integral += width * math::sin(R{2} * angle) * pi / R{4} * rule.w[i] * math::exp(point.log_jacobian);
    }
    EXPECT_REAL_NEAR(integral, R(bounds.hard_sum / R{2} - bounds.hard_lower), uni20::parse_real<R>("2e-11"));
  }
}

TEST(FourSpinonMeasure, BothSectionsAndOrderedFactorMatchBoxSplines)
{
  double const pi = 4 * std::atan(1.);
  Geometry<double> geometry(24, 24.);
  auto volume = [pi](double q) {
    double excess = std::max(0., q - pi);
    return (q * q * q - 4 * excess * excess * excess) / 6;
  };
  auto energy_volume = [pi](double q) {
    auto f = [pi](double y) {
      if (y <= 0) return 0.;
      return y <= pi ? y * y + 2 * std::cos(y) - 2 : 2 * y * y - 2 * pi * y + pi * pi - 4;
    };
    return pi * (f(q) - 3 * f(q - pi));
  };
  for (double q : {pi / 4, pi / 2, pi, 7 * pi / 4})
  {
    bethe::detail::CompensatedSum<double> sum, energy;
    std::array<bool, 2> sections{};
    auto count = geometry.moments(q, [&](auto const& point) {
      double const weight = std::exp(point.log_measure);
      sum.add(weight);
      energy.add(weight * point.energy);
      if (!sections[0] || !sections[1])
      {
        double total = 0;
        for (double b : point.beta)
          total += 2 * std::atan(std::exp(b));
        int const section = int(std::round((total - q) / (2 * pi)));
        ASSERT_TRUE(section == 0 || section == 1);
        EXPECT_NEAR(total, q + 2 * pi * section, 1e-13);
        sections[section] = true;
      }
    });
    double folded = std::min(q, 2 * pi - q);
    EXPECT_GT(count, 0u);
    EXPECT_TRUE(sections[0] && sections[1]);
    EXPECT_NEAR(sum.value(), volume(folded) + volume(2 * pi - folded), 2e-9);
    EXPECT_NEAR(energy.value(), energy_volume(folded) + energy_volume(2 * pi - folded), 2e-8);
  }
}

TEST(FourSpinonMeasure, PhysicalDensityReconstructsIndependentSlice)
{
  using R = long double;
  R const pi = R{4} * math::atan(R{1});
  Geometry<R> geometry(8, R{16});
  bethe::heisenberg::XXXSpinonDensity<R> density;
  bethe::detail::CompensatedSum<R> sum, tail, arithmetic;
  auto const count = geometry.spectrum(pi, R{2}, [&](auto const& point) {
    auto const value = density.evaluate(point.beta, {64, 16});
    sum.add(math::exp(math::log(value.value) + point.log_measure));
    tail.add(value.contour_tail_error * math::exp(point.log_measure));
    arithmetic.add(value.arithmetic_indicator * math::exp(point.log_measure));
  });
  // Independently parameterized radial quadrature and independent 8F7 kernel,
  // with the same finite order/cutoff, in the longitudinal spin convention.
  R const reference = uni20::parse_real<R>("0.0981022159286235520475");
  EXPECT_EQ(count, 512u);
  EXPECT_REAL_NEAR(sum.value(), reference, R(tail.value() + arithmetic.value() + R{2e-10L}));
  EXPECT_LT(tail.value() + arithmetic.value(), R{1e-7L});
}

TEST(FourSpinonControls, InvalidInputsFailBeforeSampling)
{
  EXPECT_THROW((Geometry<double>(0, 16)), std::invalid_argument);
  EXPECT_THROW((Geometry<double>(std::numeric_limits<std::size_t>::max(), 16)), std::invalid_argument);
  EXPECT_THROW((Geometry<double>(8, -1)), std::invalid_argument);
  EXPECT_THROW((Patch<double>(1, .5, false)), std::invalid_argument);
  Patch<double> patch(1., .25, false);
  EXPECT_FALSE(patch.support(20.));
  EXPECT_THROW(patch.support(patch.bounds(patch.radius()).upper_energy), std::runtime_error);
  EXPECT_THROW(patch.support(patch.bounds(0.).lower_energy), std::runtime_error);
  Geometry<double> geometry(4, 16);
  EXPECT_THROW(geometry.moments(0., [](auto const&) {}), std::invalid_argument);
  EXPECT_THROW(geometry.spectrum(1., -1., [](auto const&) {}), std::invalid_argument);
  EXPECT_EQ(geometry.spectrum(1., 20., [](auto const&) { ADD_FAILURE(); }), 0u);
}

#if UNI20_HAS_MPREAL
TEST(FourSpinonGeometryMpfr, ExplicitPrecisionBeyondBinary128)
{
  using R = uni20::mpreal;
  for (int bits : {160, 256})
  {
    auto const p = uni20::Precision::bits(bits);
    R const one(1, p);
    point_checks(one);
    Geometry<R> geometry(4, R(16, p));
    R const q = R{2} * math::atan(one);
    EXPECT_GT(geometry.spectrum(q, R(2, p),
                                [&](auto const& point) {
                                  EXPECT_EQ(point.log_measure.precision(), p);
                                  for (auto const& beta : point.beta)
                                    EXPECT_EQ(beta.precision(), p);
                                }),
              0u);
    EXPECT_THROW(geometry.moments(q.at(uni20::Precision::bits(bits + 16)), [](auto const&) {}), std::invalid_argument);
  }
  EXPECT_THROW((Geometry<R>(4, R{16})), std::invalid_argument);
}
#endif
} // namespace
