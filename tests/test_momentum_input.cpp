// SPDX-License-Identifier: GPL-3.0-or-later
#include <bethe/momentum_input.hpp>
#include <cmath>
#include <gtest/gtest.h>
#include <limits>

using bethe::MomentumCoordinate;
using bethe::MomentumEndpoint;
using bethe::MomentumInput;
using Kind = bethe::MomentumInputKind;

TEST(MomentumInput, LogCoordinateSurvivesBeyondNativeExponentRange)
{
  for (double b : {-1000., 1000.})
  {
    auto c = MomentumCoordinate<double>::from_beta(b);
    ASSERT_TRUE(c.beta());
    EXPECT_EQ(*c.beta(), b);
    EXPECT_NEAR(*c.log_distance(), std::log(2.) - 1000., 1e-12);
    EXPECT_THROW(c.radians(), std::range_error);
  }
}

TEST(MomentumInput, ExplicitDistancePreservesBothSidesWithoutSubtraction)
{
  auto a = MomentumInput{Kind::radians, "1e-100"}.resolve(1.);
  auto b = MomentumInput{Kind::distance_from_pi, "1e-100"}.resolve(1.);
  EXPECT_NEAR(*a.coordinate.beta(), -*b.coordinate.beta(), 1e-13);
  EXPECT_NEAR(*b.coordinate.log_distance(), -100. * std::log(10.), 1e-12);
  EXPECT_NEAR(a.require_radians() / 1e-100, 1., 1e-13);
  EXPECT_FALSE(b.scalar_radians);
  EXPECT_THROW(b.require_radians(), std::range_error);
  EXPECT_FALSE(b.beta_accuracy_warning);
}

TEST(MomentumInput, ScalarDistanceDoesNotUnderflowThroughTheHalfAngle)
{
  double const smallest = std::numeric_limits<double>::denorm_min();
  for (auto text : {"3e-324", "5e-324", "7e-324"})
  {
    auto r = MomentumInput{Kind::radians, text}.resolve(1.);
    ASSERT_TRUE(r.scalar_radians) << text;
    EXPECT_EQ(*r.scalar_radians, smallest) << text;
    EXPECT_EQ(r.coordinate.radians(), smallest) << text;
  }
  auto two = MomentumInput{Kind::radians, "8e-324"}.resolve(1.);
  EXPECT_EQ(two.require_radians(), 2 * smallest);
  EXPECT_EQ(two.coordinate.radians(), 2 * smallest);
  auto too_small = MomentumInput{Kind::radians, "2e-324"}.resolve(1.);
  EXPECT_TRUE(too_small.coordinate.beta());
  EXPECT_FALSE(too_small.scalar_radians);
}

TEST(MomentumInput, OrdinaryAndReflectedAnglesRemainSupported)
{
  for (auto text : {"0.3", "1.5", "1.8", "2", "3.0", "3.2", "5.0"})
  {
    auto r = MomentumInput{Kind::radians, text}.resolve(1.);
    EXPECT_EQ(r.require_radians(), std::stod(text));
    EXPECT_EQ((MomentumInput{Kind::radians, text}.resolve(1.L).require_radians()), std::stold(text));
#if UNI20_HAS_FLOAT128
    EXPECT_EQ((MomentumInput{Kind::radians, text}.resolve(uni20::float128{1}).require_radians()),
              uni20::parse_real<uni20::float128>(text));
#endif
    EXPECT_EQ(r.coordinate.reflected(), std::stod(text) > std::acos(-1.));
    EXPECT_FALSE(r.beta_accuracy_warning);
    EXPECT_FALSE(r.scalar_accuracy_warning);
  }
}

TEST(MomentumInput, OrdinaryDecimalThatRoundsOntoPiIsRejected)
{
  auto r = MomentumInput{Kind::radians, "3.141592653589793"}.resolve(1.);
  ASSERT_TRUE(r.coordinate.beta());
  EXPECT_FALSE(r.scalar_radians);
  EXPECT_THROW(r.require_radians(), std::range_error);
}

TEST(MomentumInput, ExactEndpointsAreDistinctFromFiniteBeta)
{
  for (auto text : {"0", "pi", "2*pi"})
  {
    auto r = MomentumInput{Kind::radians, text}.resolve(1.);
    EXPECT_FALSE(r.coordinate.beta());
    EXPECT_FALSE(r.coordinate.log_distance());
    EXPECT_TRUE(r.scalar_radians);
  }
  auto r = MomentumInput{Kind::distance_from_pi, "0"}.resolve(1.);
  EXPECT_EQ(r.coordinate.endpoint(), MomentumEndpoint::pi);
  EXPECT_EQ(r.require_radians(), std::acos(-1.));
  auto zero = MomentumInput{Kind::distance_from_pi, "pi"}.resolve(1.);
  EXPECT_EQ(zero.coordinate.endpoint(), MomentumEndpoint::zero);
  EXPECT_EQ(zero.require_radians(), 0.);
  auto center = MomentumInput{Kind::beta, "0"}.resolve(1.);
  EXPECT_EQ(center.coordinate.endpoint(), MomentumEndpoint::interior);
  EXPECT_EQ(center.require_radians(), std::acos(-1.) / 2);
}

TEST(MomentumInput, RejectsInvalidAndNonfiniteInputs)
{
  for (auto kind : {Kind::radians, Kind::distance_from_pi, Kind::beta})
    for (auto text : {"nan", "inf", "-inf", "junk", "1.2junk"})
      EXPECT_THROW((MomentumInput{kind, text}.resolve(1.)), std::exception);
  for (auto text : {"-1", "7"})
    EXPECT_THROW((MomentumInput{Kind::radians, text}.resolve(1.)), std::invalid_argument);
  for (auto text : {"-1", "4"})
    EXPECT_THROW((MomentumInput{Kind::distance_from_pi, text}.resolve(1.)), std::invalid_argument);
  for (auto tol : {"0", "-1", "nan"})
    EXPECT_THROW((MomentumInput{Kind::beta, "1"}.resolve(1., tol)), std::invalid_argument);
}

TEST(MomentumInput, NonzeroDecimalCannotUnderflowToAnExactEndpoint)
{
  for (auto kind : {Kind::radians, Kind::distance_from_pi, Kind::beta})
    EXPECT_THROW((MomentumInput{kind, "1e-999999999999999"}.resolve(1.)), std::exception);
  auto zero = MomentumInput{Kind::radians, "0.000e-123"}.resolve(1.);
  EXPECT_EQ(zero.coordinate.endpoint(), MomentumEndpoint::zero);
}

TEST(MomentumInput, DistinguishesBetaAndScalarGeometryAccuracy)
{
  auto r = MomentumInput{Kind::distance_from_pi, "1e-10"}.resolve(1., "1e-12");
  EXPECT_FALSE(r.beta_accuracy_warning);
  EXPECT_TRUE(r.scalar_radians);
  EXPECT_TRUE(r.scalar_accuracy_warning);
}

TEST(MomentumInput, ScalarRoundtripIncludesItsOwnReferenceConditioning)
{
  // Independent high-precision evaluation of this scalar angle against pi
  // gives |delta beta| = 8.85411894716521739146e-5. A native reference can
  // accidentally reproduce the input beta instead; it must flag uncertainty.
  auto r = MomentumInput{Kind::beta, "35.80001810975594200500693453648182185134"}.resolve(1.);
  ASSERT_TRUE(r.scalar_radians);
  EXPECT_EQ(*r.scalar_radians, 0x1.921fb54442d17p+1);
  EXPECT_FALSE(r.beta_accuracy_warning);
  EXPECT_TRUE(r.scalar_accuracy_warning);
  EXPECT_GT(std::stod(r.scalar_beta_absolute_error), 8e-5);
  if constexpr (std::numeric_limits<long double>::digits == 64)
  {
    // The corresponding x86 extended-precision example has true error
    // 0.26311679938959533976, despite matching the native reference closely.
    auto extended = MomentumInput{Kind::beta, "43.6682723752765544929"}.resolve(1.L);
    ASSERT_TRUE(extended.scalar_radians);
    EXPECT_EQ(*extended.scalar_radians, 0xc.90fdaa22168c234p-2L);
    EXPECT_TRUE(extended.scalar_accuracy_warning);
    EXPECT_GT(std::stod(extended.scalar_beta_absolute_error), 0.26);
  }
}

#if UNI20_HAS_FLOAT128
TEST(MomentumInput, Binary128CannotCertifyItsOwnLastEndpointInterval)
{
  using R = uni20::float128;
  R const one = 1, pi = R{4} * uni20::math::atan(one);
  R const q = pi - R{2} * uni20::numeric_limits<R>::epsilon(one);
  // Deliberately derive beta from a same-precision approximation of pi.
  // A nearly perfect same-precision roundtrip is not an accurate oracle.
  R const beta = -uni20::math::log((pi - q) / R{2});
  auto r = MomentumInput{Kind::beta, uni20::format_scalar(beta)}.resolve(one);
  ASSERT_TRUE(r.scalar_radians);
  EXPECT_EQ(*r.scalar_radians, q);
  EXPECT_TRUE(r.scalar_accuracy_warning);
}
#endif

#if !UNI20_HAS_MPREAL
TEST(MomentumInput, RoundedEndpointDiagnosticNamesTheRequestedEndpoint)
{
#if UNI20_HAS_FLOAT128
  uni20::float128 const one = 1;
#else
  long double const one = 1;
#endif
  auto const pi = 4 * uni20::math::atan(one);
  for (auto const& [value, endpoint] : {std::pair{pi, "pi"}, std::pair{2 * pi, "2*pi"}})
  {
    try
    {
      (void)MomentumInput{Kind::radians, uni20::format_scalar(value)}.resolve(1.);
      FAIL() << "rounded endpoint was accepted";
    }
    catch (std::range_error const& error)
    {
      EXPECT_NE(std::string(error.what()).find(std::string("specify ") + endpoint), std::string::npos);
    }
  }
}
#endif

TEST(MomentumInput, NativeExtendedTypesKeepTheirWorkingType)
{
  auto a = MomentumInput{Kind::beta, "10.5"}.resolve(1.L);
  static_assert(std::is_same_v<decltype(a.require_radians()), long double>);
  EXPECT_FALSE(a.beta_accuracy_warning);
#if UNI20_HAS_FLOAT128
  auto b = MomentumInput{Kind::beta, "10.5"}.resolve(uni20::float128{1});
  static_assert(std::is_same_v<decltype(b.require_radians()), uni20::float128>);
  EXPECT_FALSE(b.beta_accuracy_warning);
#endif
}

#if UNI20_HAS_MPREAL
TEST(MomentumInputMpfr, DecimalMomentumIsConvertedBeforeWorkingPrecisionRounding)
{
  using uni20::mpreal;
  using uni20::Precision;
  auto before = mpfr_get_default_prec();
  mpreal const one(1, Precision::bits(1024)), pi = 4 * uni20::math::atan(one);
  for (auto angle : {"1.8", "2", "3.2"})
    EXPECT_EQ((MomentumInput{Kind::radians, angle}.resolve(one).require_radians()), mpreal(angle, one.precision()));
  mpreal const distance("1e-100", one.precision());
  std::string const text = (pi - distance).to_string(180);
  auto native = MomentumInput{Kind::radians, text}.resolve(1.);
  ASSERT_TRUE(native.coordinate.beta());
  EXPECT_NEAR(*native.coordinate.beta(), 100. * std::log(10.) + std::log(2.), 1e-12);
  EXPECT_FALSE(native.scalar_radians);
  EXPECT_TRUE(native.mpfr_input);
  EXPECT_GE(native.input_bits, 700);
  auto wide = MomentumInput{Kind::radians, text}.resolve(mpreal(1, Precision::bits(512)), "1e-40");
  ASSERT_TRUE(wide.scalar_radians);
  EXPECT_FALSE(wide.scalar_accuracy_warning);
  EXPECT_FALSE(wide.beta_accuracy_warning);
  EXPECT_EQ(wide.require_radians().precision(), Precision::bits(512));
  // Input-to-beta subtraction also has a precision cost, even before the
  // result is rounded to the working scalar. Do not report only beta ULPs.
  auto tight = MomentumInput{Kind::radians, text}.resolve(one, "1e-250");
  EXPECT_TRUE(tight.beta_accuracy_warning);
  EXPECT_EQ(mpfr_get_default_prec(), before);
}

TEST(MomentumInputMpfr, ArbitraryExponentAndExplicitConversionWarning)
{
  auto small = MomentumInput{Kind::radians, "1e-10000"}.resolve(1.);
  ASSERT_TRUE(small.coordinate.beta());
  EXPECT_TRUE(std::isfinite(*small.coordinate.beta()));
  EXPECT_FALSE(small.scalar_radians);
  auto rounded = MomentumInput{Kind::beta, "1000.00000000000000000000000001"}.resolve(1., "1e-30");
  EXPECT_TRUE(rounded.beta_accuracy_warning);
  EXPECT_EQ(*rounded.coordinate.beta(), 1000.);
}
#endif
