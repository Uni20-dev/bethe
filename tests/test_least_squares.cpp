// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/detail/newton.hpp>

namespace
{
template <typename Real> class NewtonLeastSquares : public ::testing::Test {};
TYPED_TEST_SUITE(NewtonLeastSquares, test_support::RealTypes, test_support::PrecisionNames);

template <typename Real> void check_rectangular_fit(uni20::Precision p)
{
  auto number = [p](int x) { return uni20::scalar_at_precision<Real>(x, p); };
  std::vector<Real> const a{number(1), number(0), number(0), number(1), number(1), number(1)};
  auto const original = a;
  std::vector<Real> b{number(1), number(2), number(4)};
  ASSERT_TRUE(bethe::detail::least_squares_step(a, b, 2, p));
  ASSERT_EQ(b.size(), 2u);
  Real const eps = uni20::numeric_limits<Real>::epsilon(p);
  EXPECT_REAL_NEAR(b[0], number(4) / number(3), number(64) * eps);
  EXPECT_REAL_NEAR(b[1], number(7) / number(3), number(64) * eps);
  EXPECT_EQ(a, original);
  EXPECT_EQ(uni20::precision_of(b[0]), p);
  EXPECT_EQ(uni20::precision_of(b[1]), p);
  // An inconsistent fit is successful. Physical residual acceptance belongs
  // to the caller, not the linear-algebra adapter.
  Real const r0 = b[0] - number(1), r1 = b[1] - number(2), r2 = b[0] + b[1] - number(4);
  EXPECT_REAL_NEAR(r0 * r0 + r1 * r1 + r2 * r2, number(1) / number(3), number(128) * eps);
}

template <typename Real> void check_precision(uni20::Precision p)
{
  auto number = [p](int x) { return uni20::scalar_at_precision<Real>(x, p); };
  Real const eps = uni20::numeric_limits<Real>::epsilon(p);
  Real const offset = number(256) * eps;
  std::vector<Real> const a{number(1), number(0), number(0), number(1), number(1), number(1)};
  std::vector<Real> b{number(1) + offset, number(2) - offset, number(3)};
  ASSERT_TRUE(bethe::detail::least_squares_step(a, b, 2, p));
  EXPECT_REAL_NEAR(b[0], number(1) + offset, number(32) * eps);
  EXPECT_REAL_NEAR(b[1], number(2) - offset, number(32) * eps);
  EXPECT_EQ(uni20::precision_of(b[0]), p);
  EXPECT_EQ(uni20::precision_of(b[1]), p);
}

TYPED_TEST(NewtonLeastSquares, InconsistentRectangularFitPreservesCoefficients)
{
  check_rectangular_fit<TypeParam>(uni20::native_precision<TypeParam>());
}

TYPED_TEST(NewtonLeastSquares, NativePrecisionIsNotNarrowed)
{
  check_precision<TypeParam>(uni20::native_precision<TypeParam>());
}

TYPED_TEST(NewtonLeastSquares, ResolvesColumnsLostByNormalEquations)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const delta = uni20::math::sqrt(eps) / Real{4};
  // A^T A would have identical rounded entries, but QR can resolve delta.
  EXPECT_EQ(Real{1} + delta * delta, Real{1});
  std::vector<Real> b{Real{2}, delta, Real{1}};
  ASSERT_TRUE(
      bethe::detail::least_squares_step(std::vector<Real>{Real{1}, Real{1}, Real{0}, delta, Real{0}, Real{0}}, b, 2));
  EXPECT_REAL_NEAR(b[0], Real{1}, Real{64} * eps / delta);
  EXPECT_REAL_NEAR(b[1], Real{1}, Real{64} * eps / delta);
}

TYPED_TEST(NewtonLeastSquares, RejectsBadShapesAndSingularityWithoutChangingRhs)
{
  using Real = TypeParam;
  std::vector<Real> empty;
  EXPECT_TRUE(bethe::detail::least_squares_step(empty, empty, 0));
  std::vector<Real> b{Real{1}, Real{2}, Real{4}};
  auto const original = b;
  EXPECT_FALSE(bethe::detail::least_squares_step(std::vector<Real>{Real{1}}, b, 2));
  EXPECT_EQ(b, original);
  EXPECT_FALSE(bethe::detail::least_squares_step(empty, b, 4));
  EXPECT_EQ(b, original);
  EXPECT_FALSE(
      bethe::detail::least_squares_step(std::vector<Real>{Real{1}, Real{2}, Real{2}, Real{4}, Real{3}, Real{6}}, b, 2));
  EXPECT_EQ(b, original);
  EXPECT_FALSE(
      bethe::detail::least_squares_step(std::vector<Real>{Real{1}, Real{0}, Real{2}, Real{0}, Real{3}, Real{0}}, b, 2));
  EXPECT_EQ(b, original);
  EXPECT_TRUE(bethe::detail::least_squares_step(empty, b, 0));
  EXPECT_TRUE(b.empty());
}

TYPED_TEST(NewtonLeastSquares, RetainsBetheDiagonalThresholdAtExtremeScales)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real scale : {Real{1}, uni20::numeric_limits<Real>::max() / Real{4}, uni20::numeric_limits<Real>::min() / eps})
    for (int multiple : {32, 64, 128})
    {
      Real const d = scale * (Real(multiple) * eps);
      std::vector<Real> b{scale, d, Real{0}};
      auto const original = b;
      bool const solved =
          bethe::detail::least_squares_step(std::vector<Real>{scale, Real{0}, Real{0}, d, Real{0}, Real{0}}, b, 2);
      EXPECT_EQ(solved, multiple > 64);
      if (solved)
      {
        EXPECT_REAL_NEAR(b[0], Real{1}, Real{16} * eps);
        EXPECT_REAL_NEAR(b[1], Real{1}, Real{16} * eps);
      }
      else
        EXPECT_EQ(b, original);
    }
}

TYPED_TEST(NewtonLeastSquares, NonfiniteInputAndOverflowPreserveRhs)
{
  using Real = TypeParam;
  for (Real invalid : {uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    std::vector<Real> b{Real{3}, Real{4}};
    auto const original = b;
    EXPECT_FALSE(bethe::detail::least_squares_step(std::vector<Real>{invalid, Real{1}}, b, 1));
    EXPECT_EQ(b, original);
    b[0] = invalid;
    EXPECT_FALSE(bethe::detail::least_squares_step(std::vector<Real>{Real{1}, Real{1}}, b, 1));
    EXPECT_FALSE(uni20::isfinite(b[0]));
    EXPECT_EQ(b.size(), 2u);
    EXPECT_EQ(b[1], Real{4});
  }
  std::vector<Real> b{uni20::numeric_limits<Real>::max() / Real{2}, Real{0}};
  auto const original = b;
  EXPECT_FALSE(bethe::detail::least_squares_step(std::vector<Real>{Real{0.25}, Real{0}}, b, 1));
  EXPECT_EQ(b, original);
}

#if UNI20_HAS_MPREAL
template <typename Real>
concept ImplicitLeastSquaresPrecision =
    requires(std::vector<Real> const& a, std::vector<Real>& b) { bethe::detail::least_squares_step(a, b, 1); };
static_assert(!ImplicitLeastSquaresPrecision<uni20::mpreal>);

TEST(NewtonLeastSquaresMpfr, ExplicitPrecisionAndNoAmbientDefaultChanges)
{
  auto const before = mpfr_get_default_prec();
  for (int bits : {128, 256})
  {
    auto const p = uni20::Precision::bits(bits);
    check_rectangular_fit<uni20::mpreal>(p);
    check_precision<uni20::mpreal>(p);
    EXPECT_EQ(mpfr_get_default_prec(), before);
  }
}

TEST(NewtonLeastSquaresMpfr, ExactInputsNeedFiniteRequestAndFailuresPreserveRhs)
{
  using Real = uni20::mpreal;
  auto const p = uni20::Precision::bits(128);
  std::vector<Real> const a{Real{1}, Real{0}, Real{0}, Real{1}, Real{1}, Real{1}};
  std::vector<Real> b{Real{1}, Real{2}, Real{3}};
  auto const original = b;
  // Uni20 rejects asking an exact precision for a finite bit count.
  EXPECT_THROW(bethe::detail::least_squares_step(a, b, 2, uni20::Precision::exact()), std::logic_error);
  EXPECT_EQ(b, original);
  ASSERT_TRUE(bethe::detail::least_squares_step(a, b, 2, p));
  EXPECT_EQ(b[0].precision(), p);
  EXPECT_EQ(b[1].precision(), p);
  EXPECT_TRUE(a[0].is_exact());
  Real const eps = uni20::numeric_limits<Real>::epsilon(p);
  for (int multiple : {32, 128})
  {
    Real const tiny = multiple * eps;
    b = {Real(1, p), tiny, Real(0, p)};
    auto const saved = b;
    bool const success = bethe::detail::least_squares_step(
        std::vector<Real>{Real{1}, Real{0}, Real{0}, tiny, Real{0}, Real{0}}, b, 2, p);
    EXPECT_EQ(success, multiple > 64);
    if (!success) EXPECT_EQ(b, saved);
  }
}
#endif
} // namespace
