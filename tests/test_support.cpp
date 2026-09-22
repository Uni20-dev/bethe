// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <gtest/gtest-spi.h>

namespace
{
template <typename Real> class PrecisionAssertions : public ::testing::Test {};
TYPED_TEST_SUITE(PrecisionAssertions, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(PrecisionAssertions, NativeArithmeticAndStrictBound)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const adjacent = Real{1} + eps;
  EXPECT_REAL_NEAR(adjacent, Real{1}, Real{2} * eps);
  ASSERT_REAL_NEAR(adjacent, Real{1}, Real{2} * eps);
  EXPECT_FALSE(test_support::real_near("adjacent", "one", "epsilon", adjacent, Real{1}, eps));
  if constexpr (uni20::numeric_limits<Real>::digits > uni20::numeric_limits<double>::digits)
  {
    ASSERT_EQ(static_cast<double>(adjacent), 1.0);
    EXPECT_NONFATAL_FAILURE(EXPECT_REAL_NEAR(adjacent, Real{1}, eps), "absolute error");
  }
}

TYPED_TEST(PrecisionAssertions, InvalidBoundsAndNonfiniteValues)
{
  using Real = TypeParam;
  Real const nan = uni20::numeric_limits<Real>::quiet_NaN();
  Real const inf = uni20::numeric_limits<Real>::infinity();
  for (Real tolerance : {Real{0}, Real{-1}, nan, inf})
    EXPECT_FALSE(test_support::real_near("a", "b", "tolerance", Real{1}, Real{1}, tolerance));
  for (Real value : {nan, inf, -inf})
  {
    EXPECT_FALSE(test_support::real_near("a", "b", "tolerance", value, Real{1}, Real{1}));
    EXPECT_FALSE(test_support::real_near("a", "b", "tolerance", Real{1}, value, Real{1}));
    EXPECT_FALSE(test_support::real_near("a", "b", "tolerance", value, value, Real{1}));
  }
}

TYPED_TEST(PrecisionAssertions, DiagnosticPreservesPrecision)
{
  using Real = TypeParam;
  Real const adjacent = Real{1} + uni20::numeric_limits<Real>::epsilon();
  auto const result = test_support::real_near("actual", "expected", "bound", adjacent, Real{1}, Real{0});
  ASSERT_FALSE(result);
  std::string const message = result.message();
  EXPECT_NE(message.find(uni20::format_scalar(adjacent)), std::string::npos);
  EXPECT_NE(message.find("actual"), std::string::npos);
  EXPECT_NE(message.find("bound"), std::string::npos);
}

TYPED_TEST(PrecisionAssertions, ExactComparisonRejectsOneUlp)
{
  using Real = TypeParam;
  Real const adjacent = Real{1} + uni20::numeric_limits<Real>::epsilon();
  test_support::expect_exact(adjacent, adjacent, "identical values");
  EXPECT_NONFATAL_FAILURE(test_support::expect_exact(adjacent, Real{1}, "distinct values"), "");
}
} // namespace
