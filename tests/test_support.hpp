// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <uni20/common/gtest.hpp>

#include <cmath>
#include <string>
#include <string_view>
#include <type_traits>

namespace test_support
{
using RealTypes = ::testing::Types<double, long double
#if UNI20_HAS_FLOAT128
                                   ,
                                   uni20::float128
#endif
                                   >;

struct PrecisionNames
{
    template <typename T> static std::string GetName(int)
    {
      if constexpr (std::is_same_v<T, double>) return "fp64";
      if constexpr (std::is_same_v<T, long double>) return "long_double";
      return "fp128";
    }
};

// Solver residuals and independently rounded analytic/ED references use an
// absolute error bound, not a ULP count. Unlike GoogleTest's EXPECT_NEAR, this
// predicate never converts the operands or tolerance to double. All three
// arguments must have the same scalar type. The strict bound preserves the
// original regression contracts, including rejecting NaNs and zero tolerance.
template <uni20::Real Real>
::testing::AssertionResult real_near(char const* actual_expression, char const* expected_expression,
                                     char const* tolerance_expression, Real actual, Real expected, Real tolerance)
{
  using std::abs;
  Real const error = abs(actual - expected);
  if (uni20::isfinite(actual) && uni20::isfinite(expected) && uni20::isfinite(tolerance) && tolerance > Real{0} &&
      error < tolerance)
    return ::testing::AssertionSuccess();
  return ::testing::AssertionFailure() << actual_expression << " = " << uni20::format_scalar(actual) << '\n'
                                       << expected_expression << " = " << uni20::format_scalar(expected) << '\n'
                                       << "absolute error = " << uni20::format_scalar(error) << '\n'
                                       << "required error < " << tolerance_expression << " = "
                                       << uni20::format_scalar(tolerance)
                                       << " (finite values and a positive finite tolerance required)";
}

// Round-trip formatting promises exact recovery, so a zero-ULP comparison is
// appropriate. The fallback retains native equality on platforms whose long
// double is not one of Uni20's supported ULP layouts (e.g. IBM double-double).
template <uni20::Real Real> void expect_exact(Real actual, Real expected, std::string_view context)
{
  SCOPED_TRACE(std::string(context));
  if constexpr (uni20::check::UlpComparable<Real>)
  {
    EXPECT_FLOATING_EQ(actual, expected, 0);
  }
  else
  {
    EXPECT_TRUE(actual == expected) << uni20::format_scalar(actual) << " != " << uni20::format_scalar(expected);
  }
}
} // namespace test_support

#define EXPECT_REAL_NEAR(actual, expected, tolerance)                                                                  \
  EXPECT_PRED_FORMAT3(::test_support::real_near, actual, expected, tolerance)
#define ASSERT_REAL_NEAR(actual, expected, tolerance)                                                                  \
  ASSERT_PRED_FORMAT3(::test_support::real_near, actual, expected, tolerance)
