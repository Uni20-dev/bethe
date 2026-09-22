// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxz_odd_continuation.hpp>
#include <bethe/xxz_wronskian.hpp>
#include <bit>

namespace
{
namespace engine = bethe::xxz::detail;
template <typename Real> class XXZWronskian : public ::testing::Test {};
TYPED_TEST_SUITE(XXZWronskian, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXZWronskian, ContinuedSmallSectors)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {3, 5, 7, 9})
    for (std::size_t m = 0; m <= n / 2; ++m)
      for (Real d : {-Real{1} / Real{5}, -Real{7} / Real{10}, -Real{9} / Real{10}})
      {
        SCOPED_TRACE(::testing::Message() << "N=" << n << " M=" << m << " Delta=" << uni20::format_scalar(d));
        auto const branch = engine::continue_odd_polynomial(n, d, uni20::from_twice(std::int64_t(n - 2 * m)));
        ASSERT_TRUE(branch.equations_converged);
        auto const result = engine::check_odd_wronskian<Real>(n, branch.coefficients, d, branch.center,
                                                              branch.coordinate_scale, Real{65536} * eps);
        EXPECT_EQ(result.status, engine::WronskianStatus::consistent)
            << "residual=" << uni20::format_scalar(result.residual_norm)
            << "rcond=" << uni20::format_scalar(result.reciprocal_condition);
      }
}

TYPED_TEST(XXZWronskian, SingularPairFailsOnOddRings)
{
  using Real = TypeParam;
  // u=+-eta/2 maps to z=+-i. The cleared TQ equation is polynomial,
  // but the extra singular-pair condition is 1=(-1)^N, impossible for odd N.
  std::vector<Real> const c{Real{1}, Real{0}};
  for (std::size_t n : {5, 7, 9})
  {
    auto const result = engine::check_odd_wronskian<Real>(n, c, -Real{3} / Real{5});
    EXPECT_EQ(result.status, engine::WronskianStatus::inconsistent);
    EXPECT_GT(result.residual_norm, Real{1} / Real{100});
  }
}

TYPED_TEST(XXZWronskian, LargerContinuedBranches)
{
  using Real = TypeParam;
  Real const d = -Real{999} / Real{1000};
  for (std::size_t n : {13, 17, 21})
  {
    SCOPED_TRACE(::testing::Message() << "N=" << n);
    auto const branch =
        engine::continue_odd_polynomial(n, d, uni20::from_twice(std::int64_t{1}), {.max_iterations = 500});
    ASSERT_TRUE(branch.equations_converged);
    auto const result =
        engine::check_odd_wronskian<Real>(n, branch.coefficients, d, branch.center, branch.coordinate_scale);
    EXPECT_EQ(result.status, engine::WronskianStatus::consistent)
        << "residual=" << uni20::format_scalar(result.residual_norm)
        << "rcond=" << uni20::format_scalar(result.reciprocal_condition);
  }
}

TYPED_TEST(XXZWronskian, NativeAnalyticPolynomialAndAffineInvariance)
{
  using Real = TypeParam;
  using Complex = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  Real const C = std::cos(pi / Real{5}), d = -Real{7} / Real{10};
  Real const e = (-Real{3} * d - C - std::sqrt((d + C) * (d + C) + Real{4} * C * C)) / Real{2};
  Real const u = Real{4} * C * C / (Real{2} * C * C - Real{2} * d - e);
  std::vector<Real> const c{Real{1} - u, std::tan(pi / Real{5}) * u};
  auto const direct = engine::check_odd_wronskian<Real>(5, c, d);
  ASSERT_EQ(direct.status, engine::WronskianStatus::consistent);
  Real const center = -Real{1} / Real{3}, scale = Real{3} / Real{7};
  std::vector<Real> const affine{(c[0] + c[1] * center + center * center) / (scale * scale),
                                 (c[1] + Real{2} * center) / scale};
  auto const shifted = engine::check_odd_wronskian<Real>(5, affine, d, center, scale);
  ASSERT_EQ(shifted.status, engine::WronskianStatus::consistent);
  for (std::size_t j = 0; j < direct.q.size(); ++j)
    EXPECT_REAL_NEAR(direct.q[j], shifted.q[j], Real{512} * eps);

  // Independently evaluate the Laurent identity at complex t, not by
  // reusing the diagnostic's coefficient matrix or selected rows.
  auto laurent = [](std::vector<Real> const& coefficients, Complex t) {
    Complex value{};
    for (std::size_t j = coefficients.size(); j-- > 0;)
      value = value * t * t + coefficients[j];
    for (std::size_t j = 1; j < coefficients.size(); ++j)
      value /= t;
    return value;
  };
  Real const gamma = std::acos(d);
  Complex const shift{std::cos(gamma / Real{2}), std::sin(gamma / Real{2})};
  for (Complex t : {Complex{Real{7} / Real{5}, Real{1} / Real{5}}, Complex{Real{4} / Real{5}, -Real{2} / Real{5}}})
  {
    Complex const w = (laurent(direct.p, t * shift) * laurent(direct.q, t / shift) -
                       laurent(direct.p, t / shift) * laurent(direct.q, t * shift)) /
                      Complex{Real{0}, Real{2}};
    Complex const a = t - Real{1} / t;
    EXPECT_LT(std::abs(w - a * a * a * a * a / Real{10}), Real{1024} * eps);
  }
}

TYPED_TEST(XXZWronskian, MomentumConstraintDoesNotReplaceWronskian)
{
  using Real = TypeParam;
  Real const d = -Real{7} / Real{10};
  engine::OddPolynomialCoordinates<Real> const coordinates(7, 3);
  auto branch = engine::continue_odd_polynomial(7, d, uni20::from_twice(std::int64_t{1}));
  ASSERT_TRUE(branch.equations_converged);
  branch.coefficients[0] += Real{1} / Real{100};
  coordinates.impose_momentum(branch.coefficients, coordinates.weights(branch.coordinate_scale));
  engine::PolynomialBetheSystem<Real> const system(7, 3, branch.center, branch.coordinate_scale);
  EXPECT_LT(system.momentum_defect(branch.coefficients), Real{256} * uni20::numeric_limits<Real>::epsilon());
  auto const result =
      engine::check_odd_wronskian<Real>(7, branch.coefficients, d, branch.center, branch.coordinate_scale);
  EXPECT_EQ(result.status, engine::WronskianStatus::inconsistent);
  EXPECT_GT(result.residual_norm, Real{1} / Real{10000});
}

TYPED_TEST(XXZWronskian, PhantomPointNeedsSeparateTreatment)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {5, 7, 9})
  {
    Real const d = -std::cos(pi / Real(n));
    Real const s = std::sqrt((Real{1} + d) / (Real{1} - d));
    // Q(x)=x^M with z=-s+x, so its Laurent transform is proportional to t^M.
    std::vector<Real> const c(n / 2, Real{0});
    auto const result = engine::check_odd_wronskian<Real>(n, c, d, -s);
    EXPECT_EQ(result.endpoint_ratio, Real{0});
    EXPECT_NE(result.status, engine::WronskianStatus::consistent);
    engine::PolynomialBetheSystem<Real> const system(n, c.size(), -s);
    EXPECT_REAL_NEAR(system.energy(c, d), Real(n) * d / Real{4}, Real{128} * Real(n) * eps);
  }
}

TYPED_TEST(XXZWronskian, OtherRootsOfUnityCanBeInconclusive)
{
  using Real = TypeParam;
  // At Delta=-1/2, q^3=1. Even the physical vacuum cannot have the
  // assumed Laurent P: the discrete derivative annihilates t^+-3 while
  // the RHS has nonzero coefficients there. This is NOT a bad state.
  std::vector<Real> const vacuum;
  for (std::size_t n : {3, 5, 7})
  {
    auto const result = engine::check_odd_wronskian<Real>(n, vacuum, -Real{1} / Real{2});
    EXPECT_EQ(result.status, engine::WronskianStatus::ill_conditioned);
  }
}

TYPED_TEST(XXZWronskian, ProjectedHelixIsNonzeroEigenvector)
{
  using Real = TypeParam;
  using Complex = std::complex<Real>;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  // Direct Hamiltonian action in the spin basis, independent of Bethe roots,
  // in NATIVE precision. Both chiralities and every magnetization are checked.
  for (unsigned n : {3, 5, 7, 9})
    for (int chirality : {-1, 1})
    {
      Real const pitch = Real(chirality) * (pi + pi / Real(n));
      Real const d = -std::cos(pi / Real(n)), energy = Real(n) * d / Real{4};
      for (unsigned m = 0; m <= n; ++m)
      {
        auto amplitude = [&](unsigned bits) {
          unsigned sum = 0;
          for (unsigned j = 0; j < n; ++j)
            if (bits & (1U << j)) sum += j;
          Real const angle = pitch * Real(sum);
          return Complex{std::cos(angle), std::sin(angle)};
        };
        std::size_t count = 0;
        for (unsigned bits = 0; bits < (1U << n); ++bits)
          if (std::popcount(bits) == int(m))
          {
            ++count;
            Complex const v = amplitude(bits);
            Complex action{};
            for (unsigned j = 0; j < n; ++j)
            {
              unsigned const next = (j + 1) % n;
              bool const opposite = ((bits >> j) & 1U) != ((bits >> next) & 1U);
              action += d * (opposite ? -Real{1} / Real{4} : Real{1} / Real{4}) * v;
              if (opposite) action += amplitude(bits ^ (1U << j) ^ (1U << next)) / Real{2};
            }
            EXPECT_LT(std::abs(action - energy * v), Real{512} * Real(n) * eps);
            EXPECT_REAL_NEAR(std::abs(v), Real{1}, Real{8} * eps);
          }
        EXPECT_GT(count, std::size_t{0});
      }
    }
}

TYPED_TEST(XXZWronskian, InvalidArguments)
{
  using Real = TypeParam;
  std::vector<Real> c{Real{1}, Real{0}};
  Real const d = -Real{3} / Real{5};
  EXPECT_THROW(engine::check_odd_wronskian<Real>(6, c, d), std::invalid_argument);
  EXPECT_THROW(engine::check_odd_wronskian<Real>(3, c, d), std::invalid_argument);
  for (Real bad :
       {-Real{1}, Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW(engine::check_odd_wronskian<Real>(5, c, bad), std::invalid_argument);
  EXPECT_THROW(engine::check_odd_wronskian<Real>(5, c, d, Real{0}, Real{0}), std::invalid_argument);
  EXPECT_THROW(engine::check_odd_wronskian<Real>(5, c, d, Real{0}, Real{1}, Real{0}), std::invalid_argument);
  c[0] = uni20::numeric_limits<Real>::quiet_NaN();
  EXPECT_THROW(engine::check_odd_wronskian<Real>(5, c, d), std::invalid_argument);
}
} // namespace
