// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "xxx_spinon_surface_references.hpp"
#include <bethe/detail/xxx_spinon_surface.hpp>

namespace
{
namespace math = uni20::math;
template <class R> using Surface = bethe::heisenberg::detail::SpinonEnergySurface<R>;

template <class R> R parse(char const* text, R const& one)
{
  if constexpr (uni20::has_runtime_precision_v<R>)
    return R(text, one.precision());
  else
    return uni20::parse_real<R>(text);
}

template <class R> void check_reference(spinon_surface_reference::Point const& row, R const& one)
{
  SCOPED_TRACE(row.n);
  SCOPED_TRACE(row.section);
  SCOPED_TRACE(row.cutoff);
  R const pi = R{4} * math::atan(one), eps = uni20::numeric_limits<R>::epsilon(one);
  std::vector<R> beta;
  for (int i = 0; i < row.n; ++i)
    beta.push_back(parse(row.beta[i], one));
  Surface<R> const surface(row.n, R(row.cutoff) * one, pi * parse(row.q_over_pi, one), row.section,
                           parse(row.omega, one));
  auto const result = surface.solve(std::span<R const>(beta).first(row.n - 2));
  ASSERT_TRUE(result.inside);
  EXPECT_REAL_NEAR(result.log_coarea_jacobian, parse(row.log_jacobian, one), R{65536} * eps);
  for (int root = 0; root < 2; ++root)
  {
    ASSERT_EQ(result.beta[root].size(), beta.size());
    for (int i = 0; i < row.n - 2; ++i)
      EXPECT_EQ(result.beta[root][i], beta[i]);
    // All reference hard pairs are supplied in ascending rapidity order.
    EXPECT_REAL_NEAR(result.beta[root][row.n - 2], beta[row.n - 1 - root], R{65536} * eps);
    EXPECT_REAL_NEAR(result.beta[root][row.n - 1], beta[row.n - 2 + root], R{65536} * eps);
    if constexpr (uni20::has_runtime_precision_v<R>)
      for (auto const& b : result.beta[root])
        EXPECT_EQ(b.precision(), one.precision());
  }
  if constexpr (uni20::has_runtime_precision_v<R>) EXPECT_EQ(result.log_coarea_jacobian.precision(), one.precision());
}

template <class R> class SpinonSurface : public ::testing::Test {};
TYPED_TEST_SUITE(SpinonSurface, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(SpinonSurface, IndependentForwardKinematicsAndBothRoots)
{
  for (auto const& row : spinon_surface_reference::points)
    if (row.cutoff == 8) check_reference(row, TypeParam{1});
}

TYPED_TEST(SpinonSurface, SupportZerosAndUnresolvedBoundariesDiffer)
{
  using R = TypeParam;
  R const pi = R{4} * math::atan(R{1}), q = pi / R{2};
  EXPECT_FALSE((Surface<R>(2, R{8}, q, 0, R{10}).solve({}).inside));
  EXPECT_FALSE((Surface<R>(2, R{8}, q, 0, R{1}).solve({}).inside));
  EXPECT_THROW((Surface<R>(2, R{8}, q, 0, pi * math::sin(q / R{2})).solve({})), std::runtime_error);
  EXPECT_THROW((Surface<R>(2, R{1040}, q, 0, R{2})), std::runtime_error);
  EXPECT_THROW((Surface<R>(2, R{8}, q, 1, R{2})), std::invalid_argument);
  EXPECT_THROW((Surface<R>(3, R{8}, q, 0, R{2})), std::invalid_argument);
  EXPECT_THROW((Surface<R>(6, R{8}, q, 0, R{2}).solve({})), std::invalid_argument);
  std::array<R, 4> bad{R{-8}, R{-2}, R{-1}, R{0}};
  EXPECT_THROW((Surface<R>(6, R{8}, q, 0, R{20}).solve(bad)), std::invalid_argument);
}

#if UNI20_HAS_MPREAL
TEST(SpinonSurfaceMpfr, ExplicitPrecisionAndDeepSoftGeometry)
{
  using R = uni20::mpreal;
  for (int bits : {128, 256, 1800, 2048})
  {
    SCOPED_TRACE(bits);
    auto const p = uni20::Precision::bits(bits);
    R const one(1, p), pi = R{4} * math::atan(one);
    for (auto const& row : spinon_surface_reference::points)
      if (row.cutoff == 8 || bits >= 1800) check_reference(row, one);
    EXPECT_THROW((Surface<R>(2, R(8, p), R(pi / R{2}, uni20::Precision::bits(bits + 32)), 0, R(2, p))),
                 std::invalid_argument);
    if (bits < 1800) EXPECT_THROW((Surface<R>(2, R(1040, p), pi / R{2}, 0, R(2, p))), std::runtime_error);
  }
  EXPECT_THROW((Surface<R>(2, R{8}, R{1}, 0, R{2})), std::invalid_argument);
}
#endif
} // namespace
