// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/detail/xxx_four_spinon_common_shift.hpp>
#include <bethe/xxx_spinon_density.hpp>

namespace
{
namespace math = uni20::math;
template <class R> using Shift = bethe::heisenberg::detail::FourSpinonCommonShift<R>;

template <class R> R momentum(std::array<R, 4> const& beta, R const& one)
{
  R const pi = R{4} * math::atan(one);
  R total = one - one;
  for (auto const& b : beta)
    total += b <= R{0} ? R{2} * math::atan(math::exp(b)) : pi - R{2} * math::atan(math::exp(-b));
  return total;
}

template <class R> void geometry_checks(R const& one)
{
  Shift<R> shift(one);
  R const pi = R{4} * math::atan(one), eps = uni20::numeric_limits<R>::epsilon(one);
  for (auto gaps :
       {std::array{one / R{4}, R{2} * one, R{3} * one}, std::array{R{1000} * one, R{800} * one, R{1200} * one}})
  {
    auto const representative = shift.representative(gaps);
    EXPECT_REAL_NEAR(momentum(representative.beta, one), pi, R{32768} * eps);
    for (int j = 0; j < 3; ++j)
      EXPECT_REAL_NEAR(representative.beta[j + 1] - representative.beta[j], gaps[j], R{32768} * eps);
    EXPECT_TRUE(uni20::isfinite(representative.log_rapidity_measure));
    EXPECT_LT(representative.log_rapidity_measure, one - one);
    if constexpr (uni20::has_runtime_precision_v<R>)
    {
      EXPECT_EQ(representative.log_rapidity_measure.precision(), one.precision());
      for (auto const& b : representative.beta)
        EXPECT_EQ(b.precision(), one.precision());
    }
  }
}

template <class R> class CommonShift : public ::testing::Test {};
TYPED_TEST_SUITE(CommonShift, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(CommonShift, RepresentativePreservesGapsAndPiSection) { geometry_checks(TypeParam{1}); }

TYPED_TEST(CommonShift, ClosedExpressionMatchesPositiveQuadrature)
{
  using R = TypeParam;
  R const one{1}, zero{0};
  Shift<R> shift(one);
  for (auto const& gaps :
       {std::array{zero, zero, zero}, std::array{one / R{4}, R{2}, R{3}}, std::array{R{2}, R{32}, one}})
  {
    auto const analytic = shift.zeroth(gaps), direct = shift.quadrature(gaps);
    EXPECT_FALSE(analytic.used_quadrature);
    EXPECT_TRUE(direct.used_quadrature);
    EXPECT_REAL_NEAR(analytic.value, direct.value,
                     R(analytic.arithmetic_indicator + direct.arithmetic_indicator + direct.quadrature_difference +
                       direct.omitted_tail_bound));
  }
  auto const close = shift.zeroth({uni20::parse_real<R>("1e-30"), zero, uni20::parse_real<R>("1e-30")});
  EXPECT_TRUE(close.used_quadrature);
  EXPECT_REAL_NEAR(close.value, R{8} / R{3},
                   R(close.arithmetic_indicator + close.quadrature_difference + close.omitted_tail_bound));
}

TEST(CommonShiftIntegral, FirstMomentFromIndependentShiftQuadrature)
{
  using R = long double;
  R const one{1}, pi = R{4} * math::atan(one);
  Shift<R> shift(one);
  auto const rule = bethe::detail::gauss_legendre<R>(24);
  for (std::array<R, 4> delta : {std::array<R, 4>{0, 0, 0, 0}, std::array<R, 4>{0, 1, 4, 8}})
  {
    bethe::detail::SpinonSum<R> sum(R{0});
    for (int panel = -40; panel < 40; ++panel)
      for (std::size_t j = 0; j < rule.x.size(); ++j)
      {
        R const eta = R(panel) + (one + rule.x[j]) / R{2};
        std::array<R, 4> beta;
        R energy{0};
        for (int k = 0; k < 4; ++k)
        {
          beta[k] = eta + delta[k];
          energy += pi / (R{2} * math::cosh(beta[k]));
        }
        R const sine = math::sin(momentum(beta, one) / R{2});
        sum.add(rule.w[j] / R{2} * sine * sine * energy);
      }
    EXPECT_REAL_NEAR(sum.value(), shift.first(), R{1e-15L});
  }
}

TEST(CommonShiftDensity, PhysicalWardFactorIsIndependentOfCommonShift)
{
  using R = long double;
  R const one{1};
  Shift<R> shift(one);
  bethe::heisenberg::XXXSpinonDensity<R> density;
  auto const representative = shift.representative({R{1}, R{2}, R{3}});
  auto const reference = density.evaluate(representative.beta, {64, 16});
  R const reference_scale = math::exp(representative.log_rapidity_measure);
  R const invariant = reference.value * reference_scale;
  for (R eta : {R{-1}, R{.5}, R{1}})
  {
    auto beta = representative.beta;
    R jacobian{1};
    for (auto& b : beta)
    {
      b += eta;
      jacobian /= math::cosh(b);
    }
    R const sine = math::sin(momentum(beta, one) / R{2});
    R const scale = jacobian / (sine * sine);
    auto const value = density.evaluate(beta, {64, 16});
    R const diagnostic = (reference.contour_tail_error + reference.arithmetic_indicator) * reference_scale +
                         (value.contour_tail_error + value.arithmetic_indicator) * scale;
    EXPECT_REAL_NEAR(value.value * scale, invariant, R{8} * diagnostic);
    EXPECT_LT(diagnostic, R{1e-6L} * invariant);
  }
}

TEST(CommonShiftControls, InvalidAndUnboundedInputsFail)
{
  Shift<double> shift(1.);
  EXPECT_THROW(shift.representative({-1., 2., 3.}), std::invalid_argument);
  EXPECT_THROW(shift.zeroth({1., std::numeric_limits<double>::infinity(), 3.}), std::invalid_argument);
  EXPECT_THROW(shift.quadrature({1., 9000., 3.}), std::invalid_argument);
}

TEST(CommonShiftGlobal, GapMeasureAndFirstMomentHaveNoAdditionalFactor)
{
  using R = long double;
  using Geometry = bethe::heisenberg::detail::FourSpinonGapGeometry<R>;
  Geometry geometry(16, {R{3}, R{2}, R{5}});
  R const pi = R{4} * math::atan(R{1});
  bethe::detail::SpinonSum<R> sum(R{0});
  std::size_t count = 0;
  geometry.integrate([&](auto const& point) {
    // Test density W=exp(-sum gaps)*sin^2(P/2)/prod(sech beta).
    // Its full first-moment shift integral is pi^2; after both observable
    // 1/(2*pi) measures and ordering, the gap integral has coefficient 6.
    R log_density = -(point.beta.back() - point.beta.front());
    for (auto const& b : point.beta)
      log_density += math::log(math::cosh(b));
    sum.add(math::exp(log_density + point.log_measure) * point.first_shift / (R{4} * pi * pi));
    ++count;
  });
  R const exact = R{6} * (-math::expm1(R{-3})) * (-math::expm1(R{-2})) * (-math::expm1(R{-5}));
  EXPECT_EQ(count, geometry.size());
  EXPECT_REAL_NEAR(sum.value(), exact, R{1e-14L});
  EXPECT_THROW((Geometry(256, {R{8192}, R{1}, R{1}}, 100)), std::invalid_argument);
  EXPECT_THROW((Geometry(4, {R{0}, R{1}, R{1}})), std::invalid_argument);
}

#if UNI20_HAS_MPREAL
TEST(CommonShiftMpfr, ExplicitPrecisionAndCancellationFallback)
{
  using R = uni20::mpreal;
  for (int bits : {160, 256})
  {
    auto const p = uni20::Precision::bits(bits);
    R const one(1, p), zero(0, p);
    geometry_checks(one);
    Shift<R> shift(one);
    R const small = uni20::parse_real<R>("1e-60", p);
    auto const value = shift.zeroth({small, zero, small});
    EXPECT_TRUE(value.used_quadrature);
    EXPECT_EQ(value.value.precision(), p);
    EXPECT_REAL_NEAR(value.value, R{8} * one / R{3},
                     R(value.arithmetic_indicator + value.quadrature_difference + value.omitted_tail_bound));
    EXPECT_THROW(shift.zeroth({one, one.at(uni20::Precision::bits(bits + 16)), one}), std::invalid_argument);
  }
  EXPECT_THROW((Shift<R>(R{1})), std::invalid_argument);
}
#endif
} // namespace
