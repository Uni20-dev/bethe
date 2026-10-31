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

TEST(FourSpinonMeasure, NarrowUpperEdgeSupportIsNotMissedByFractionNodes)
{
  double const pi = 4 * std::atan(1.);
  Geometry<double> geometry(8, 16.);
  double weight = 0;
  auto const count =
      geometry.spectrum(pi, 4.4428825, [&](auto const& point) { weight += std::exp(point.log_measure); });
  EXPECT_GT(count, 0u);
  EXPECT_GT(weight, 0.);
  EXPECT_TRUE(std::isfinite(weight));
}

TEST(FourSpinonMeasure, HigherOrdersResolvePositiveSpectralSupport)
{
  double const pi = 4 * std::atan(1.);
  for (auto const& point : std::array{std::pair{.5, 3.117081}, std::pair{pi, 4.4428825}, std::pair{pi, 4.220739}})
    for (std::size_t order : {16u, 48u, 128u, 256u})
    {
      SCOPED_TRACE(::testing::Message() << "k=" << point.first << " omega=" << point.second << " order=" << order);
      Geometry<double> geometry(order, 16.);
      double measure = 0;
      std::size_t nodes = 0;
      EXPECT_NO_THROW(nodes = geometry.spectrum(point.first, point.second, [&](auto const& node) {
        EXPECT_TRUE(uni20::isfinite(node.log_measure));
        double energy = 0;
        for (double beta : node.beta)
          energy += pi / (2 * std::cosh(beta));
        EXPECT_NEAR(energy, point.second, 3e-13);
        measure += std::exp(node.log_measure);
      }));
      EXPECT_GT(nodes, 0u);
      EXPECT_GT(measure, 0.);
      EXPECT_TRUE(std::isfinite(measure));
    }
}

TEST(FourSpinonMeasure, AngularAndRadialCoareaMeasuresAgree)
{
  using R = long double;
  R const pi = R{4} * math::atan(R{1});
  auto const rule = bethe::detail::gauss_legendre<R>(64);
  for (bool mixed : {false, true})
  {
    Patch<R> patch(R{7} * pi / R{4}, R{3} / R{8}, mixed);
    R const omega = (patch.bounds(R{0}).upper_energy + patch.bounds(patch.radius()).upper_energy) / R{2};
    auto const radial = patch.support(omega);
    ASSERT_TRUE(radial);
    ASSERT_GT(radial->first, R{0});
    auto const angular = patch.angular_support(omega);
    std::array<R, 2> radial_sum{}, angular_sum{};
    auto accumulate = [](auto const& point, R weight, auto& sum) {
      R observable{0};
      for (R beta : point.beta)
      {
        R const momentum = R{2} * math::atan(math::exp(beta));
        observable += momentum * momentum;
      }
      sum[0] += weight;
      sum[1] += weight * observable;
    };
    for (std::size_t i = 0; i < rule.x.size(); ++i)
    {
      // Independent coarea coordinates: eliminate the hard pair at fixed r,
      // or eliminate r at fixed hard-pair angle. Their Jacobians differ.
      R const angle = pi * (R{1} + rule.x[i]) / R{4};
      R const sine = math::sin(angle), width = radial->second - radial->first;
      R const r = radial->first + width * sine * sine;
      auto const old_point = patch.spectral_point(math::log(r), omega);
      accumulate(old_point,
                 r * width * math::sin(R{2} * angle) * pi / R{4} * rule.w[i] * math::exp(old_point.log_jacobian),
                 radial_sum);
      auto const point = patch.spectral_angle_point(omega, angular, (R{1} + rule.x[i]) / R{2});
      accumulate(point, angular.angle * rule.w[i] / R{2} * math::exp(point.log_jacobian), angular_sum);
    }
    for (std::size_t j = 0; j < radial_sum.size(); ++j)
      EXPECT_REAL_NEAR(angular_sum[j], radial_sum[j], R{2e-10L} * math::abs(radial_sum[j]));
  }
}

TEST(FourSpinonMeasure, NarrowMixedPatchSupportIncludesItsInteriorMaximum)
{
  long double const pi = 4 * std::atan(1.L);
  Geometry<long double> geometry(8, 16.);
  for (long double a : {pi / 2, 3 * pi / 4})
  {
    // At the mixed-patch radius the momenta are (a*v,b,b,pi-b),
    // b=a*(1-v). Their energy has a stationary maximum at this v.
    long double const peak = std::atan2(1 - 3 * std::cos(a), 3 * std::sin(a)) / a;
    long double const maximum = pi / 2 * (std::sin(a * peak) + 3 * std::sin(a * (1 - peak)));
    long double const omega = maximum - 1e-3L;
    long double weight = 0;
    std::size_t mixed_count = 0;
    geometry.spectrum(pi - a, omega, [&](auto const& point) {
      // Only the reflected section reaches these energies. Reflect back to
      // its ordered Q=pi+a patch and identify its two closest endpoint distances.
      std::array<long double, 4> p;
      for (int j = 0; j < 4; ++j)
        p[j] = 2 * std::atan(std::exp(-point.beta[3 - j]));
      if (pi - p[3] < p[1])
      {
        long double const v = p[0] / (p[0] + pi - p[3]);
        EXPECT_NEAR(v, peak, .03);
        EXPECT_NEAR(pi / 2 * (std::sin(p[0]) + std::sin(p[1]) + std::sin(p[2]) + std::sin(p[3])), omega, 1e-15);
        weight += std::exp(point.log_measure);
        ++mixed_count;
      }
    });
    EXPECT_GT(mixed_count, 0u);
    EXPECT_GT(weight, 0.);
    EXPECT_TRUE(std::isfinite(weight));
  }
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
  double const pi = 4 * std::atan(1.);
  double const upper = 2 * pi * std::sin(pi / 4);
  EXPECT_THROW(geometry.spectrum(pi, upper, [](auto const&) { ADD_FAILURE(); }), std::runtime_error);
  EXPECT_EQ(geometry.spectrum(pi, upper + 1e-8, [](auto const&) { ADD_FAILURE(); }), 0u);
}

#if UNI20_HAS_MPREAL
TEST(FourSpinonGeometryMpfr, NativeAngularMeasureMatchesWideArithmetic)
{
  using R = uni20::mpreal;
  auto const precision = uni20::Precision::bits(160);
  double const pi = 4 * std::atan(1.);
  for (auto const& point : std::array{std::pair{.5, 3.117081}, std::pair{pi, 4.4428825}})
    for (std::size_t order : {16u, 48u})
    {
      SCOPED_TRACE(::testing::Message() << "k=" << point.first << " omega=" << point.second << " order=" << order);
      Geometry<double> native(order, 16.);
      Geometry<R> wide(order, R(16, precision));
      R sum(0, precision);
      double native_sum = 0;
      auto const native_count = native.spectrum(point.first, point.second,
                                                [&](auto const& node) { native_sum += math::exp(node.log_measure); });
      auto const wide_count = wide.spectrum(R(point.first, precision), R(point.second, precision),
                                            [&](auto const& node) { sum += math::exp(node.log_measure); });
      EXPECT_EQ(native_count, wide_count);
      EXPECT_GT(sum, R(0, precision));
      // At the extremely narrow edge, input/constant rounding is amplified
      // by the inverse distance to the edge. Compare identical finite rules,
      // not an asserted continuum quadrature accuracy.
      EXPECT_REAL_NEAR(R(native_sum, precision), sum, R(2e-8, precision) * sum);
    }
}

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
    R const pi = R{2} * q;
    R const upper = R{2} * pi * math::sin(pi / R{4});
    // The allowed fraction interval is far narrower than binary64 can resolve.
    R const omega = upper - uni20::parse_real<R>("1e-30", p);
    R weight(0, p);
    EXPECT_GT(geometry.spectrum(pi, omega, [&](auto const& point) { weight += math::exp(point.log_measure); }), 0u);
    EXPECT_GT(weight, R(0, p));
    EXPECT_TRUE(uni20::isfinite(weight));
    EXPECT_THROW(geometry.moments(q.at(uni20::Precision::bits(bits + 16)), [](auto const&) {}), std::invalid_argument);
  }
  EXPECT_THROW((Geometry<R>(4, R{16})), std::invalid_argument);
}
#endif
} // namespace
