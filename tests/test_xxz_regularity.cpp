// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxz_odd_continuation.hpp>
#include <bethe/xxz_regularity.hpp>

namespace
{
namespace engine = bethe::xxz::detail;
template <typename Real> class XXZRegularity : public ::testing::Test {};
TYPED_TEST_SUITE(XXZRegularity, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXZRegularity, RegularContinuedSectors)
{
  using Real = TypeParam;
  for (std::size_t n : {3, 5, 7, 9})
    for (std::size_t m = 0; m <= n / 2; ++m)
      for (Real d : {Real{0}, -Real{1} / Real{5}, -Real{7} / Real{10}, -Real{9} / Real{10}, -Real{999} / Real{1000}})
      {
        SCOPED_TRACE(::testing::Message() << "N=" << n << " M=" << m << " Delta=" << uni20::format_scalar(d));
        auto const branch = engine::continue_odd_polynomial(n, d, uni20::from_twice(std::int64_t(n - 2 * m)));
        ASSERT_TRUE(branch.equations_converged);
        engine::PolynomialBetheSystem<Real> const system(n, m, branch.center, branch.coordinate_scale);
        auto const result = engine::check_regular_polynomial<Real>(system, branch.coefficients, d);
        EXPECT_EQ(result.status, engine::RegularityStatus::regular_on_shell)
            << "endpoints=" << uni20::format_scalar(result.endpoint_margin)
            << "distinct=" << uni20::format_scalar(result.distinct_roots_margin)
            << "scattering=" << uni20::format_scalar(result.scattering_margin)
            << "jacobian=" << uni20::format_scalar(result.jacobian_margin);
      }
}

TYPED_TEST(XXZRegularity, ExceptionalFactorsAreNotCertified)
{
  using Real = TypeParam;
  Real const d = -Real{3} / Real{5};
  engine::PolynomialBetheSystem<Real> const system(7, 2);
  // Repeated finite root, singular driving pair z=+-i, and infinite
  // conventional rapidity z=-sqrt((1+d)/(1-d))=-1/2 respectively.
  for (std::vector<Real> const c :
       {std::vector<Real>{Real{0}, Real{0}}, {Real{1}, Real{0}}, {Real{0}, Real{1} / Real{2}}})
  {
    auto const result = engine::check_regular_polynomial<Real>(system, c, d);
    EXPECT_EQ(result.status, engine::RegularityStatus::exceptional_or_unresolved);
  }
  // A pair separated by i*gamma has a vanishing pair-scattering factor.
  // Conjugate roots z=+-i*r have 1+d-(1-d)*r^2+2*d*r=0.
  Real const r = (d + Real{1}) / (Real{1} - d);
  std::vector<Real> const exact_string{r * r, Real{0}};
  auto const result = engine::check_regular_polynomial<Real>(system, exact_string, d);
  EXPECT_EQ(result.status, engine::RegularityStatus::exceptional_or_unresolved);
  EXPECT_LE(result.scattering_margin, Real{64} * uni20::numeric_limits<Real>::epsilon());
}

TYPED_TEST(XXZRegularity, LargerRegularBranches)
{
  using Real = TypeParam;
  for (std::size_t n : {13, 17, 21})
    for (Real d : {-Real{97} / Real{100}, -Real{999} / Real{1000}})
    {
      SCOPED_TRACE(::testing::Message() << "N=" << n << " Delta=" << uni20::format_scalar(d));
      auto const branch =
          engine::continue_odd_polynomial(n, d, uni20::from_twice(std::int64_t{1}), {.max_iterations = 500});
      ASSERT_TRUE(branch.equations_converged);
      engine::PolynomialBetheSystem<Real> const system(n, n / 2, branch.center, branch.coordinate_scale);
      auto const result = engine::check_regular_polynomial<Real>(system, branch.coefficients, d);
      EXPECT_EQ(result.status, engine::RegularityStatus::regular_on_shell)
          << "endpoints=" << uni20::format_scalar(result.endpoint_margin)
          << "distinct=" << uni20::format_scalar(result.distinct_roots_margin)
          << "scattering=" << uni20::format_scalar(result.scattering_margin)
          << "jacobian=" << uni20::format_scalar(result.jacobian_margin);
    }
}

TYPED_TEST(XXZRegularity, InfiniteAndPhantomContinuedStatesNeedLimits)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1});
  for (std::size_t n : {5, 7, 9})
    for (Real d : {-Real{1} / Real{2}, -std::cos(pi / Real(n))})
    {
      auto const branch = engine::continue_odd_polynomial(n, d, uni20::from_twice(std::int64_t{1}));
      ASSERT_TRUE(branch.equations_converged);
      engine::PolynomialBetheSystem<Real> const system(n, n / 2, branch.center, branch.coordinate_scale);
      auto const result = engine::check_regular_polynomial<Real>(system, branch.coefficients, d);
      EXPECT_EQ(result.status, engine::RegularityStatus::exceptional_or_unresolved)
          << "N=" << n << " Delta=" << uni20::format_scalar(d)
          << "endpoints=" << uni20::format_scalar(result.endpoint_margin)
          << "distinct=" << uni20::format_scalar(result.distinct_roots_margin);
    }
  // Root-of-unity anisotropy alone is not disqualifying: the vacuum and
  // this finite one-magnon solution still satisfy the regular criterion.
  for (std::size_t m : {0, 1})
  {
    Real const d = -Real{1} / Real{2};
    auto const branch = engine::continue_odd_polynomial(7, d, uni20::from_twice(std::int64_t(7 - 2 * m)));
    engine::PolynomialBetheSystem<Real> const system(7, m, branch.center, branch.coordinate_scale);
    EXPECT_EQ(engine::check_regular_polynomial<Real>(system, branch.coefficients, d).status,
              engine::RegularityStatus::regular_on_shell);
  }
}

TYPED_TEST(XXZRegularity, OffShellAndRequestedTolerance)
{
  using Real = TypeParam;
  Real const d = -Real{7} / Real{10};
  auto branch = engine::continue_odd_polynomial(7, d, uni20::from_twice(std::int64_t{1}));
  ASSERT_TRUE(branch.equations_converged);
  engine::PolynomialBetheSystem<Real> const system(7, 3, branch.center, branch.coordinate_scale);
  branch.coefficients[0] += Real{1} / Real{1000};
  auto const result = engine::check_regular_polynomial<Real>(system, branch.coefficients, d);
  EXPECT_EQ(result.status, engine::RegularityStatus::off_shell);
  EXPECT_EQ(result.residual_norm, system.evaluate(branch.coefficients, d).norm);
  EXPECT_GT(result.residual_norm, bethe::SolverOptions<Real>{}.residual_tolerance);
  // The caller's tolerance is used literally, not replaced by the default.
  // This intentionally loose check is not a claim that the perturbed vector
  // is an exact eigenstate.
  auto const loose =
      engine::check_regular_polynomial<Real>(system, branch.coefficients, d, Real{2} * result.residual_norm);
  EXPECT_EQ(loose.status, engine::RegularityStatus::regular_on_shell);
  EXPECT_EQ(loose.residual_norm, result.residual_norm);
  auto const tight =
      engine::check_regular_polynomial<Real>(system, branch.coefficients, d, result.residual_norm / Real{2});
  EXPECT_EQ(tight.status, engine::RegularityStatus::off_shell);
}

TYPED_TEST(XXZRegularity, RegularFiveSiteVectorsAreNonzeroEigenstates)
{
  using Real = TypeParam;
  using Complex = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  Real const C = std::cos(pi / Real{5});
  auto power = [](Complex z, unsigned n) {
    Complex out{1};
    for (unsigned j = 0; j < n; ++j)
      out *= z;
    return out;
  };
  for (Real d : {-Real{1} / Real{5}, -Real{7} / Real{10}, -Real{9} / Real{10}, -Real{999} / Real{1000}})
  {
    Real const e = (-Real{3} * d - C - std::sqrt((d + C) * (d + C) + Real{4} * C * C)) / Real{2};
    Real const u = Real{4} * C * C / (Real{2} * C * C - Real{2} * d - e);
    std::vector<Real> const c{Real{1} - u, std::tan(pi / Real{5}) * u};
    engine::PolynomialBetheSystem<Real> const system(5, 2);
    auto const result = engine::check_regular_polynomial<Real>(system, c, d, Real{256} * eps);
    ASSERT_EQ(result.status, engine::RegularityStatus::regular_on_shell);
    Real const discriminant = c[1] * c[1] - Real{4} * c[0];
    Complex const root =
        discriminant >= Real{0} ? Complex{std::sqrt(discriminant), 0} : Complex{0, std::sqrt(-discriminant)};
    Complex const z1 = (-c[1] + root) / Real{2}, z2 = (-c[1] - root) / Real{2}, imaginary{0, 1};
    Complex const v1 = -(Real{1} - imaginary * z1) / (Real{1} + imaginary * z1);
    Complex const v2 = -(Real{1} - imaginary * z2) / (Real{1} + imaginary * z2);
    Complex const a = Real{1} + v1 * v2 - Real{2} * d * v1;
    Complex const b = -(Real{1} + v1 * v2 - Real{2} * d * v2);
    std::vector<Complex> wave(32);
    Real largest = Real{0};
    for (unsigned x = 0; x < 5; ++x)
      for (unsigned y = x + 1; y < 5; ++y)
      {
        auto const bits = (1U << x) | (1U << y);
        wave[bits] = a * power(v1, x) * power(v2, y) + b * power(v2, x) * power(v1, y);
        largest = std::max(largest, std::abs(wave[bits]));
      }
    ASSERT_GT(largest, Real{1} / Real{10000});
    for (unsigned x = 0; x < 5; ++x)
      for (unsigned y = x + 1; y < 5; ++y)
      {
        auto const bits = (1U << x) | (1U << y);
        Complex action{};
        for (unsigned j = 0; j < 5; ++j)
        {
          auto const k = (j + 1) % 5;
          bool const anti = ((bits >> j) & 1U) != ((bits >> k) & 1U);
          action += d * (anti ? -Real{1} / Real{4} : Real{1} / Real{4}) * wave[bits];
          if (anti) action += wave[bits ^ (1U << j) ^ (1U << k)] / Real{2};
        }
        EXPECT_LT(std::abs(action - (Real{5} * d / Real{4} + e) * wave[bits]) / largest, Real{8192} * eps);
      }
    // Check the self-removed scattering polynomial independently, including
    // a nontrivial affine coordinate change, at both analytic roots.
    Real const center = -Real{1} / Real{3}, scale = Real{3} / Real{7};
    std::vector<Real> const affine{(c[0] + center * c[1] + center * center) / (scale * scale),
                                   (c[1] + Real{2} * center) / scale};
    engine::PolynomialBetheSystem<Real> const shifted(5, 2, center, scale);
    auto const k = shifted.scattering_remainder(affine, d);
    for (auto const z : {z1, z2})
    {
      Complex const other = -c[1] - z;
      Complex const expected = Real{1} + d - (Real{1} - d) * z * other - imaginary * d * (z - other);
      EXPECT_LT(std::abs(k[0] + k[1] * ((z - center) / scale) - expected), Real{512} * eps);
    }
  }
}

TYPED_TEST(XXZRegularity, MultiplicationMatricesDetectCommonRoots)
{
  using Real = TypeParam;
  using Complex = std::complex<Real>;
  std::vector<Real> const c{-Real{1}, Real{0}}; // x^2-1
  EXPECT_LE(engine::quotient_multiplication_margin<Real>(c, {Complex{-1}, Complex{1}}),
            Real{64} * uni20::numeric_limits<Real>::epsilon());
  EXPECT_GT(engine::quotient_multiplication_margin<Real>(c, {Complex{0}, Complex{1}}), Real{1} / Real{8});
  EXPECT_GT(engine::quotient_multiplication_margin<Real>(c, {Complex{0, 1}, Complex{1}}), Real{1} / Real{10});
}

TYPED_TEST(XXZRegularity, InvalidInput)
{
  using Real = TypeParam;
  engine::PolynomialBetheSystem<Real> const system(7, 2);
  std::vector<Real> c{-Real{1}, Real{0}};
  for (Real d : {-Real{1}, Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW(engine::check_regular_polynomial<Real>(system, c, d), std::invalid_argument);
  for (Real t : {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW(engine::check_regular_polynomial<Real>(system, c, -Real{1} / Real{2}, t), std::invalid_argument);
  c.pop_back();
  EXPECT_THROW(engine::check_regular_polynomial<Real>(system, c, -Real{1} / Real{2}), std::invalid_argument);
}
} // namespace
