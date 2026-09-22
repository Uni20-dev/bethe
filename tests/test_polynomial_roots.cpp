// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/polynomial_roots.hpp>

#include "test_support.hpp"
#include <array>

namespace
{
namespace engine = bethe::detail;
using Status = engine::PolynomialRootStatus;
template <typename Real> class PolynomialRoots : public ::testing::Test {};
TYPED_TEST_SUITE(PolynomialRoots, test_support::RealTypes, test_support::PrecisionNames);

template <typename Real>
void expect_roots(engine::PolynomialRoots<Real> const& result, std::vector<std::complex<Real>> const& expected,
                  Real tolerance)
{
  ASSERT_EQ(result.status, Status::resolved) << "residual=" << uni20::format_scalar(result.residual_norm)
                                             << " reconstruction=" << uni20::format_scalar(result.reconstruction_error)
                                             << " separation=" << uni20::format_scalar(result.separation_ratio);
  ASSERT_EQ(result.roots.size(), expected.size());
  std::vector<bool> used(expected.size());
  for (auto x : result.roots)
  {
    std::size_t nearest = expected.size();
    Real distance = uni20::numeric_limits<Real>::infinity();
    for (std::size_t j = 0; j < expected.size(); ++j)
      if (!used[j] && std::abs(x - expected[j]) < distance)
      {
        nearest = j;
        distance = std::abs(x - expected[j]);
      }
    ASSERT_LT(nearest, expected.size());
    EXPECT_LT(distance, tolerance * std::max(Real{1}, std::abs(expected[nearest])));
    used[nearest] = true;
  }
}

TYPED_TEST(PolynomialRoots, VacuumLinearAndMixedExactRoots)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  auto const vacuum = engine::recover_polynomial_roots<Real>({});
  EXPECT_EQ(vacuum.status, Status::resolved);
  EXPECT_TRUE(vacuum.roots.empty());
  EXPECT_EQ(vacuum.iterations, 0U);
  for (Real root : {Real{0}, -Real{7} / Real{4}, Real{1024}})
  {
    auto const result = engine::recover_polynomial_roots<Real>(std::array<Real, 1>{-root});
    expect_roots(result, std::vector<C>{C{root}}, Real{8} * eps);
    EXPECT_EQ(result.iterations, 0U);
  }
  // x*(x-2)*(x+1)*(x^2+1) = x^5-x^4-x^3-x^2-2x.
  std::array<Real, 5> const c{Real{0}, -Real{2}, -Real{1}, -Real{1}, -Real{1}};
  auto const result = engine::recover_polynomial_roots<Real>(c);
  expect_roots(result, std::vector<C>{C{}, C{2}, C{-1}, C{0, 1}, C{0, -1}}, Real{1024} * eps);
  EXPECT_LE(result.residual_norm, engine::PolynomialRootOptions<Real>{}.tolerance);
  EXPECT_LE(result.reconstruction_error, engine::PolynomialRootOptions<Real>{}.tolerance);
  EXPECT_GT(result.separation_ratio, Real{8});
}

TYPED_TEST(PolynomialRoots, RootsOfUnityThroughDegreeTwentyFour)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (unsigned n : {2, 3, 8, 17, 24})
  {
    SCOPED_TRACE(::testing::Message() << "degree=" << n);
    std::vector<Real> c(n);
    c[0] = -Real{1};
    std::vector<C> expected;
    for (unsigned j = 0; j < n; ++j)
    {
      Real const angle = Real{2} * pi * Real(j) / Real(n);
      expected.emplace_back(std::cos(angle), std::sin(angle));
    }
    expect_roots(engine::recover_polynomial_roots<Real>(c), expected, Real{2048} * eps);
  }
}

TYPED_TEST(PolynomialRoots, ScalingAndCallerSeeds)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real scale : {Real{1} / Real{1024}, Real{1024}})
  {
    std::array<Real, 3> const c{-Real{6} * scale * scale * scale, Real{11} * scale * scale, -Real{6} * scale};
    expect_roots(engine::recover_polynomial_roots<Real>(c), std::vector<C>{C{scale}, C{2 * scale}, C{3 * scale}},
                 Real{8192} * eps);
  }
  std::array<Real, 2> const c{-Real{2}, Real{0}};
  std::array<C, 2> const seed{C{-1}, C{1}};
  auto const initial = engine::recover_polynomial_roots<Real>(c, {.max_iterations = 0}, seed);
  EXPECT_EQ(initial.status, Status::iteration_limit);
  EXPECT_EQ(initial.iterations, 0U);
  EXPECT_TRUE(initial.roots[0] == seed[0]);
  EXPECT_GT(initial.residual_norm, Real{1} / Real{10});
  auto const single = engine::recover_polynomial_roots<Real>(c, {.max_iterations = 1}, seed);
  EXPECT_EQ(single.status, Status::iteration_limit);
  EXPECT_EQ(single.iterations, 1U);
  expect_roots(engine::recover_polynomial_roots<Real>(c, {}, seed),
               std::vector<C>{C{-std::sqrt(Real{2})}, C{std::sqrt(Real{2})}}, Real{1024} * eps);
}

TYPED_TEST(PolynomialRoots, RepeatedRootsAndDuplicateSolutionsAreNotResolved)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  for (auto const& c :
       {std::vector<Real>{Real{0}, Real{0}}, {Real{1}, -Real{2}}, {Real{1}, -Real{4}, Real{6}, -Real{4}}})
  {
    auto const result = engine::recover_polynomial_roots<Real>(c);
    EXPECT_NE(result.status, Status::resolved);
  }
  std::array<Real, 2> const c{-Real{1}, Real{0}};
  std::array<C, 2> const duplicate{C{1}, C{1}};
  auto const result = engine::recover_polynomial_roots<Real>(c, {.max_iterations = 0}, duplicate);
  EXPECT_EQ(result.residual_norm, Real{0}); // Individually exact roots are insufficient.
  EXPECT_GT(result.reconstruction_error, Real{1});
  EXPECT_EQ(result.status, Status::iteration_limit);
  EXPECT_EQ(engine::recover_polynomial_roots<Real>(c, {}, duplicate).status, Status::unresolved_cluster);
  // A pair distinguishable at this precision, and an unresolved rounded pair.
  Real const gap = std::sqrt(std::sqrt(uni20::numeric_limits<Real>::epsilon()));
  std::array<Real, 2> const separated{Real{1} - gap * gap, -Real{2}};
  auto const resolved = engine::recover_polynomial_roots<Real>(separated);
  EXPECT_EQ(resolved.status, Status::resolved);
}

TYPED_TEST(PolynomialRoots, InputValidationAndUnrepresentableScaling)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  auto solve = [](auto const& c, engine::PolynomialRootOptions<Real> options = {}) {
    return engine::recover_polynomial_roots<Real>(c, options);
  };
  std::vector<Real> const c{Real{1}, Real{0}};
  EXPECT_THROW(solve(c, {.tolerance = Real{0}}), std::invalid_argument);
  EXPECT_THROW(solve(c, {.tolerance = uni20::numeric_limits<Real>::quiet_NaN()}), std::invalid_argument);
  EXPECT_THROW(solve(c, {.max_degree = 1}), std::length_error);
  EXPECT_THROW(solve(std::array<Real, 1>{uni20::numeric_limits<Real>::infinity()}), std::invalid_argument);
  EXPECT_THROW(engine::recover_polynomial_roots<Real>(c, {}, std::array<C, 1>{C{1}}), std::invalid_argument);
  EXPECT_THROW(engine::recover_polynomial_roots<Real>(
                   c, {}, std::array<C, 2>{C{1}, C{0, uni20::numeric_limits<Real>::infinity()}}),
               std::invalid_argument);
  std::array<Real, 2> const extreme{uni20::numeric_limits<Real>::min(), uni20::numeric_limits<Real>::max()};
  EXPECT_EQ(solve(extreme).status, Status::nonfinite_or_unrepresentable);
}
} // namespace
