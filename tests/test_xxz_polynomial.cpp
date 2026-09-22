// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "exact_spectrum.hpp"
#include "test_support.hpp"
#include <bethe/xxz_polynomial.hpp>
#include <uni20/linalg/ops/linear_solve.hpp>

namespace
{
namespace engine = bethe::xxz::detail;
template <typename Real> class XXZPolynomial : public ::testing::Test {};
TYPED_TEST_SUITE(XXZPolynomial, test_support::RealTypes, test_support::PrecisionNames);

template <typename Scalar> std::vector<Scalar> coefficients(std::span<Scalar const> roots)
{
  std::vector<Scalar> c{Scalar{1}};
  for (Scalar z : roots)
  {
    std::vector<Scalar> next(c.size() + 1);
    for (std::size_t i = 0; i < c.size(); ++i)
    {
      next[i] -= z * c[i];
      next[i + 1] += c[i];
    }
    c = std::move(next);
  }
  c.pop_back();
  return c;
}

template <typename Real> std::vector<Real> free_coefficients(std::size_t n, std::size_t m)
{
  Real const pi = Real{4} * std::atan(Real{1});
  std::vector<Real> roots(m);
  for (std::size_t j = 0; j < m; ++j)
    roots[j] = std::tan(pi * (Real(j) - Real(m - 1 + n % 2) / Real{2}) / Real(n));
  return coefficients<Real>(roots);
}

TYPED_TEST(XXZPolynomial, RootwiseEquationsAndObservables)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  // Deliberately not solutions: both parity branches must agree with direct
  // products at every root, including a complex conjugate pair.
  std::vector<std::vector<C>> const sets{
      {C{-Real{1} / Real{4}}, C{Real{1} / Real{8}}},
      {C{-Real{1} / Real{4}}, C{Real{1} / Real{8}, Real{1} / Real{5}}, C{Real{1} / Real{8}, -Real{1} / Real{5}}}};
  for (auto const& roots : sets)
    for (std::size_t n : {8, 9})
      for (Real d : {-Real{1} / Real{10}, -Real{1} / Real{2}, -Real{9} / Real{10}})
      {
        auto const complex_c = coefficients<C>(roots);
        std::vector<Real> c;
        for (C v : complex_c)
        {
          EXPECT_REAL_NEAR(v.imag(), Real{0}, Real{32} * eps);
          c.push_back(v.real());
        }
        engine::PolynomialBetheSystem<Real> const system(n, roots.size());
        auto const f = system.evaluate(c, d);
        C energy{Real(n) * d / Real{4}}, momentum{1};
        int const sign = (n - roots.size() - 1) % 2 ? -1 : 1;
        for (std::size_t i = 0; i < roots.size(); ++i)
        {
          C const z = roots[i], one{1}, imag{0, 1};
          C plus{1}, minus{1}, driving_plus{1}, driving_minus{1};
          for (std::size_t j = 0; j < roots.size(); ++j)
            if (i != j)
            {
              C const w = roots[j], common = Real{1} + d - (Real{1} - d) * z * w;
              plus *= common + imag * d * (z - w);
              minus *= common - imag * d * (z - w);
            }
          for (std::size_t j = 0; j < n; ++j)
          {
            driving_plus *= one + imag * z;
            driving_minus *= one - imag * z;
          }
          C const original = driving_plus * minus - Real(sign) * driving_minus * plus;
          C reduced{};
          for (std::size_t j = c.size(); j-- > 0;)
            reduced = reduced * z + f.residual[j];
          reduced *= sign == 1 ? C{0, 2} : C{2};
          EXPECT_REAL_NEAR(std::abs(original - reduced), Real{0}, Real{8192} * eps);
          energy -= (Real{1} + d - (Real{1} - d) * z * z) / (one + z * z);
          momentum *= -(one - imag * z) / (one + imag * z);
        }
        EXPECT_REAL_NEAR(energy.imag(), Real{0}, Real{128} * eps);
        EXPECT_REAL_NEAR(system.energy(c, d), energy.real(), Real{128} * eps);
        EXPECT_REAL_NEAR(std::abs(system.momentum_phase(c) - momentum), Real{0}, Real{128} * eps);
      }
}

TYPED_TEST(XXZPolynomial, AnalyticJacobian)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), h = std::cbrt(eps);
  for (std::size_t n : {8, 9})
    for (std::size_t m : {1, 2, 3, 4})
      for (Real d : {Real{0}, -Real{1} / Real{2}, -Real{99} / Real{100}})
      {
        engine::PolynomialBetheSystem<Real> const system(n, m);
        auto const c = free_coefficients<Real>(n, m);
        uni20::DenseMatrix<Real> jac(m, m);
        auto const f = system.evaluate(c, d, &jac);
        EXPECT_EQ(f.residual, system.evaluate(c, d).residual);
        for (std::size_t j = 0; j < m; ++j)
        {
          auto plus = c, minus = c;
          plus[j] += h;
          minus[j] -= h;
          auto const fp = system.evaluate(plus, d), fm = system.evaluate(minus, d);
          for (std::size_t i = 0; i < m; ++i)
          {
            Real const difference = (fp.residual[i] - fm.residual[i]) / (Real{2} * h);
            EXPECT_REAL_NEAR((jac[i, j]), difference, Real{8192} * h * h * (Real{1} + std::abs(jac[i, j])));
          }
        }
      }
}

TYPED_TEST(XXZPolynomial, FiveSiteExactRootCollision)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  Real const C = std::cos(pi / Real{5}), tangent = std::tan(pi / Real{5});
  engine::PolynomialBetheSystem<Real> const system(5, 2);
  // In momentum 2*pi/5, the relative-separation basis r=1,2 has H-E_F =
  // [[-Delta,C],[C,-2*Delta-C]]. Its smaller eigenvalue fixes Q directly.
  for (Real d : {-Real{1} / Real{2}, -C, -Real{9} / Real{10}, -Real{1} + Real{128} * eps})
  {
    Real const relative = (-Real{3} * d - C - std::sqrt((d + C) * (d + C) + Real{4} * C * C)) / Real{2};
    Real const u = Real{4} * C * C / (Real{2} * C * C - Real{2} * d - relative);
    std::vector<Real> const c{Real{1} - u, tangent * u};
    auto const f = system.evaluate(c, d);
    EXPECT_LT(f.norm, Real{1024} * eps);
    EXPECT_REAL_NEAR(system.energy(c, d), Real{5} * d / Real{4} + relative, Real{128} * eps);
    EXPECT_LT(system.momentum_defect(c), Real{128} * eps);
    Real const discriminant = c[1] * c[1] - Real{4} * c[0];
    if (d > -C)
      EXPECT_GT(discriminant, Real{0});
    else if (d < -C)
      EXPECT_LT(discriminant, Real{0});
    else
      EXPECT_REAL_NEAR(discriminant, Real{0}, Real{128} * eps);
    // Analytic coefficient derivatives remain nonsingular at the collision.
    uni20::DenseMatrix<Real> jac(2, 2);
    (void)system.evaluate(c, d, &jac);
    EXPECT_GT(std::abs(jac[0, 0] * jac[1, 1] - jac[0, 1] * jac[1, 0]), Real{1} / Real{100});
  }
}

TYPED_TEST(XXZPolynomial, SmallRingContinuationAgainstED)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  // Test-only continuation is deliberately not a production solver: no
  // adaptive steps, budget API, or general physical-state certification.
  for (unsigned n : {3, 5, 7, 9})
    for (unsigned m = 1; m <= n / 2; ++m)
    {
      engine::PolynomialBetheSystem<Real> const system(n, m);
      auto c = free_coefficients<Real>(n, m);
      for (unsigned coupling = 0; coupling <= 40; ++coupling)
      {
        Real const d = -Real(coupling) * Real{999} / Real{40000};
        SCOPED_TRACE(::testing::Message() << "N=" << n << " M=" << m << " Delta=" << uni20::format_scalar(d));
        bool converged = false;
        for (unsigned iteration = 0; iteration < 30; ++iteration)
        {
          uni20::DenseMatrix<Real> jac(m, m), step(m, 1);
          auto const f = system.evaluate(c, d, &jac);
          if (f.norm < Real{2048} * eps)
          {
            converged = true;
            break;
          }
          for (unsigned i = 0; i < m; ++i)
            step[i, 0] = -f.residual[i];
          uni20::linalg::solve_inplace(jac, step);
          for (unsigned i = 0; i < m; ++i)
            c[i] += step[i, 0];
        }
        ASSERT_TRUE(converged);
        // Crucially, a coefficient residual alone is not the test oracle.
        EXPECT_LT(system.momentum_defect(c), Real{65536} * eps);
        Real const energy = system.energy(c, d);
        auto const ed = test_support::exact_spectrum(n, m, 0, true, static_cast<double>(d));
        EXPECT_REAL_NEAR(static_cast<double>(energy), ed.front(), 3e-10);
        if (coupling % 10 == 0)
        {
          auto const resolved = test_support::exact_spectrum(n, m, .371, true, static_cast<double>(d));
          double const target = static_cast<double>(energy + Real{371} / Real{1000} * system.momentum_phase(c).real());
          EXPECT_TRUE(
              std::any_of(resolved.begin(), resolved.end(), [&](double e) { return std::abs(e - target) < 3e-10; }));
        }
      }
    }
}

TYPED_TEST(XXZPolynomial, EmptyAndInvalidInputs)
{
  using Real = TypeParam;
  engine::PolynomialBetheSystem<Real> const empty(5, 0), system(5, 2);
  EXPECT_EQ(empty.evaluate({}, -Real{1} / Real{2}).norm, Real{0});
  EXPECT_EQ(empty.energy({}, -Real{1} / Real{2}), -Real{5} / Real{8});
  EXPECT_EQ(empty.momentum_defect({}), Real{0});
  std::vector<Real> const c{Real{1} / Real{10}, Real{1} / Real{2}};
  EXPECT_THROW((void)engine::PolynomialBetheSystem<Real>(1, 0), std::invalid_argument);
  EXPECT_THROW((void)engine::PolynomialBetheSystem<Real>(5, 3), std::invalid_argument);
  EXPECT_THROW((void)system.evaluate({}, -Real{1} / Real{2}), std::invalid_argument);
  for (Real bad : {-Real{1}, Real{1}, uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW((void)system.evaluate(c, bad), std::invalid_argument);
  uni20::DenseMatrix<Real> wrong(1, 1);
  EXPECT_THROW((void)system.evaluate(c, -Real{1} / Real{2}, &wrong), std::invalid_argument);
  std::vector<Real> const pole{Real{1}, Real{0}}; // Q(z)=z^2+1
  EXPECT_THROW((void)system.energy(pole, -Real{1} / Real{2}), std::runtime_error);
  EXPECT_THROW((void)system.momentum_phase(pole), std::runtime_error);
  // Generic std::complex<fp128> division must not square these huge values.
  Real const large = uni20::numeric_limits<Real>::max() / Real{16};
  std::vector<Real> const huge{large, large / Real{2}};
  auto const phase = system.momentum_phase(huge);
  EXPECT_REAL_NEAR(phase.real(), Real{3} / Real{5}, Real{32} * uni20::numeric_limits<Real>::epsilon());
  EXPECT_REAL_NEAR(phase.imag(), -Real{4} / Real{5}, Real{32} * uni20::numeric_limits<Real>::epsilon());
  EXPECT_TRUE(uni20::isfinite(system.energy(huge, -Real{1} / Real{2})));
}

TYPED_TEST(XXZPolynomial, SelfScatteringFactorMustBeRemoved)
{
  using Real = TypeParam;
  // z=-1/2 at Delta=-3/5 makes 1+Delta-(1-Delta)*z^2 vanish.
  // Including this self factor would make the cleared Bethe equation zero
  // although the one-magnon momentum is not allowed on a five-site ring.
  engine::PolynomialBetheSystem<Real> const system(5, 1);
  std::vector<Real> const c{Real{1} / Real{2}};
  auto const f = system.evaluate(c, -Real{3} / Real{5});
  EXPECT_GT(f.norm, Real{1} / Real{2});
  EXPECT_GT(system.momentum_defect(c), Real{1} / Real{2});
}

TEST(XXZPolynomialAudit, SmallRawResidualDoesNotCertifyMomentum)
{
  // Recorded from an unguarded coefficient-Newton continuation. This is not
  // a reference ground state: the necessary translation condition fails.
  engine::PolynomialBetheSystem<double> const system(21, 10);
  std::vector<double> const c{-3.5966213246302427e-10, -1.1365731703277832e-8, 2.002100484436825e-7,
                              1.5977595890131712e-5,   .00036578804595449957,  .004660972853483087,
                              .03736415744276909,      .19379360494403983,     .6348494784002282,
                              1.1996157338248057};
  auto const f = system.evaluate(c, -.97);
  for (double value : f.residual)
    EXPECT_LT(std::abs(value), 1e-7);
  EXPECT_GT(system.momentum_defect(c), .1);
}
} // namespace
