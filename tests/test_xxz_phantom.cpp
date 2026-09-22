// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/xxz_phantom.hpp>

#include "test_support.hpp"
#include <bethe/xxz_odd_continuation.hpp>
#include <bit>

namespace
{
namespace engine = bethe::xxz::detail;
template <typename Real> class XXZPhantom : public ::testing::Test {};
TYPED_TEST_SUITE(XXZPhantom, test_support::RealTypes, test_support::PrecisionNames);

template <typename Scalar> std::vector<Scalar> polynomial(std::span<Scalar const> roots)
{
  std::vector<Scalar> c{Scalar{1}};
  for (auto root : roots)
  {
    c.push_back(Scalar{0});
    for (std::size_t j = c.size(); j-- > 0;)
      c[j] = (j ? c[j - 1] : Scalar{0}) - root * c[j];
  }
  c.pop_back();
  return c;
}

TYPED_TEST(XXZPhantom, RotatedRootwiseEquations)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), d = -Real{3} / Real{5};
  Real const angle = Real{2} / Real{7};
  C const rotation{std::cos(angle), std::sin(angle)}, imag{0, 1};
  std::vector<C> const roots{
      {-Real{1} / Real{4}}, {Real{1} / Real{8}, Real{1} / Real{5}}, {Real{1} / Real{8}, -Real{1} / Real{5}}};
  auto const cc = polynomial<C>(roots);
  std::vector<Real> c;
  for (auto x : cc)
    c.push_back(x.real());
  for (std::size_t n : {8, 9})
  {
    engine::PolynomialBetheSystem<Real> const system(n, roots.size());
    auto const f = system.evaluate_rotated(c, d, rotation);
    int const sign = (n - roots.size() - 1) % 2 ? -1 : 1;
    for (std::size_t i = 0; i < roots.size(); ++i)
    {
      C const z = roots[i];
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
        driving_plus *= Real{1} + imag * z;
        driving_minus *= Real{1} - imag * z;
      }
      C const expected = rotation * driving_plus * minus - Real(sign) * std::conj(rotation) * driving_minus * plus;
      C actual{};
      for (std::size_t j = c.size(); j-- > 0;)
        actual = actual * z + f.residual[j];
      actual *= sign == 1 ? C{0, 2} : C{2};
      EXPECT_LT(std::abs(actual - expected), Real{8192} * eps);
    }
    EXPECT_EQ(system.evaluate(c, d).residual, system.evaluate_rotated(c, d, C{1}).residual);
  }
}

TYPED_TEST(XXZPhantom, RotatedAnalyticJacobian)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), h = std::cbrt(eps), d = -Real{7} / Real{10};
  std::complex<Real> const rotation{std::cos(Real{1}), std::sin(Real{1})};
  for (std::size_t n : {8, 9})
  {
    engine::PolynomialBetheSystem<Real> const system(n, 3, -Real{1} / Real{7}, Real{3} / Real{5});
    std::vector<Real> const c{-Real{1} / Real{50}, Real{1} / Real{20}, Real{1} / Real{5}};
    uni20::DenseMatrix<Real> jac(3, 3);
    auto const f = system.evaluate_rotated(c, d, rotation, &jac);
    EXPECT_EQ(f.residual, system.evaluate_rotated(c, d, rotation).residual);
    for (std::size_t j = 0; j < 3; ++j)
    {
      auto plus = c, minus = c;
      plus[j] += h;
      minus[j] -= h;
      auto const fp = system.evaluate_rotated(plus, d, rotation), fm = system.evaluate_rotated(minus, d, rotation);
      for (std::size_t i = 0; i < 3; ++i)
        EXPECT_REAL_NEAR((jac[i, j]), (fp.residual[i] - fm.residual[i]) / (Real{2} * h),
                         Real{8192} * h * h * (Real{1} + std::abs(jac[i, j])));
    }
  }
}

TYPED_TEST(XXZPhantom, TwistedOneMagnon)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  Real const phi = Real{2} / Real{5}, d = -Real{3} / Real{5};
  for (std::size_t n : {7, 8})
  {
    Real const k = (Real{6} * pi + phi) / Real(n), z = Real{1} / std::tan(k / Real{2});
    std::vector<Real> const c{-z};
    std::complex<Real> const rotation{std::cos(phi / Real{2}), std::sin(phi / Real{2})};
    engine::PolynomialBetheSystem<Real> const system(n, 1);
    auto const check = engine::check_regular_polynomial<Real>(system, c, d, Real{256} * eps, rotation);
    EXPECT_EQ(check.status, engine::RegularityStatus::regular_on_shell);
    EXPECT_REAL_NEAR(system.energy(c, d), Real(n) * d / Real{4} + std::cos(k) - d, Real{64} * eps);
    EXPECT_GT(system.evaluate(c, d).norm, Real{1} / Real{10});
    EXPECT_GT(system.momentum_defect(c), Real{1} / Real{10}); // Not the twisted momentum condition.
  }
}

TYPED_TEST(XXZPhantom, ContinuedMixedClusters)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {5, 7, 9, 13, 17, 21})
  {
    auto const m = n / 2;
    for (std::size_t p = 1; p < m; ++p)
    {
      // N-2*(M-p)=2*p+1 for these minimal-|Sz| sectors.
      Real const d = p == 1 ? -Real{1} / Real{2} : -std::cos(pi / Real(2 * p + 1));
      SCOPED_TRACE(::testing::Message() << "N=" << n << " p=" << p);
      auto const branch =
          engine::continue_odd_polynomial(n, d, uni20::from_twice(std::int64_t{1}), {.max_iterations = 500});
      ASSERT_TRUE(branch.equations_converged);
      engine::PolynomialBetheSystem<Real> const system(n, m, branch.center, branch.coordinate_scale);
      auto const reduced = engine::reduce_phantom_polynomial<Real>(system, branch.coefficients, d, p, -1);
      ASSERT_TRUE(reduced.finite_regularity);
      if (p == 1) EXPECT_EQ(reduced.status, engine::PhantomReductionStatus::regular_reduced_equations);
      // Deflation and the changed residual normalization can amplify the
      // full polynomial's roundoff. The default must still be enforced.
      if (reduced.finite_regularity->residual_norm > Real{32} * eps)
      {
        EXPECT_EQ(reduced.status, engine::PhantomReductionStatus::unresolved_reduced_state);
        EXPECT_EQ(reduced.finite_regularity->status, engine::RegularityStatus::off_shell);
      }
      else
        EXPECT_EQ(reduced.status, engine::PhantomReductionStatus::regular_reduced_equations);
      // This separate, explicitly requested coefficient-residual tolerance
      // audits the larger clusters; it is NOT an automatic fallback.
      auto const audit = engine::reduce_phantom_polynomial<Real>(system, branch.coefficients, d, p, -1, Real{128} * eps,
                                                                 Real{512} * eps);
      EXPECT_EQ(audit.status, engine::PhantomReductionStatus::regular_reduced_equations)
          << "factor=" << uni20::format_scalar(audit.reconstruction_error)
          << "phase=" << uni20::format_scalar(audit.commensurability_error);
      EXPECT_EQ(reduced.finite_coefficients, audit.finite_coefficients);
      EXPECT_LT(reduced.energy_difference, Real{1024} * Real(n) * eps);
    }
  }
}

TYPED_TEST(XXZPhantom, OneFiniteMagnonWithPhantomsIsNonzeroEigenvector)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (unsigned n : {5, 7, 9})
    for (unsigned p = 1; p < n / 2; ++p)
      for (int chirality : {-1, 1})
        for (unsigned label = 0; label < n; ++label)
        {
          // The commensurability condition uses r=1, irrespective of p.
          Real const gamma = pi - pi / Real(n - 2), k0 = Real(chirality) * gamma;
          Real const d = -std::cos(pi / Real(n - 2));
          Real const k = (Real{2} * pi * Real(label) - Real{2} * Real(p) * k0) / Real(n);
          // Exclude k=0 (the polynomial coordinate is infinite) and k=+-gamma
          // (this would not be a REGULAR finite root). All other labels here
          // are separated from these exceptional values by more than 0.001.
          if (std::abs(std::sin(k / Real{2})) < Real{1} / Real{1000} ||
              std::abs(std::cos(k) - d) < Real{1} / Real{1000})
            continue;
          SCOPED_TRACE(::testing::Message()
                       << "N=" << n << " p=" << p << " chirality=" << chirality << " label=" << label);
          Real const z0 = Real(chirality) * std::sqrt((Real{1} + d) / (Real{1} - d));
          Real const z = Real{1} / std::tan(k / Real{2});
          std::vector<Real> roots(p, z0);
          roots.push_back(z);
          auto const c = polynomial<Real>(roots);
          engine::PolynomialBetheSystem<Real> const system(n, roots.size());
          auto const reduction =
              engine::reduce_phantom_polynomial<Real>(system, c, d, p, chirality, Real{1024} * eps, Real{1024} * eps);
          ASSERT_EQ(reduction.status, engine::PhantomReductionStatus::regular_reduced_equations);
          Real const energy = Real(n) * d / Real{4} + std::cos(k) - d;
          EXPECT_REAL_NEAR(system.energy(c, d), energy, Real{2048} * eps);
          auto phase = [](Real angle) { return C{std::cos(angle), std::sin(angle)}; };
          // Derived coordinate ansatz: crossing one phantom multiplies the
          // finite-magnon amplitude by exp(2*i*k0). No Bethe residual is used
          // to check the Hamiltonian action below.
          auto amplitude = [&](unsigned bits) {
            unsigned rank = 0, sum = 0;
            C terms{};
            for (unsigned j = 0; j < n; ++j)
              if (bits & (1U << j))
              {
                terms += phase(Real{2} * k0 * Real(rank) + (k - k0) * Real(j));
                sum += j;
                ++rank;
              }
            return phase(k0 * Real(sum)) * terms;
          };
          Real largest = Real{0}, residual = Real{0};
          for (unsigned bits = 0; bits < (1U << n); ++bits)
            if (std::popcount(bits) == int(p + 1))
            {
              C const wave = amplitude(bits);
              largest = std::max(largest, std::abs(wave));
              C action{};
              for (unsigned j = 0; j < n; ++j)
              {
                auto const next = (j + 1) % n;
                bool const anti = ((bits >> j) & 1U) != ((bits >> next) & 1U);
                action += d * (anti ? -Real{1} : Real{1}) * wave / Real{4};
                if (anti) action += amplitude(bits ^ (1U << j) ^ (1U << next)) / Real{2};
              }
              residual = std::max(residual, std::abs(action - energy * wave));
            }
          ASSERT_GT(largest, Real{1} / Real{100});
          EXPECT_LT(residual / largest, Real{2048} * eps);
        }
}

TYPED_TEST(XXZPhantom, FactorsPhaseAndRegularityAreSeparate)
{
  using Real = TypeParam;
  Real const d = -Real{1} / Real{2}, a = std::sqrt(Real{1} / Real{3});
  // Exact endpoint factor, but the finite root is deliberately off shell.
  std::vector<Real> roots{-a, Real{1} / Real{5}};
  auto c = polynomial<Real>(roots);
  engine::PolynomialBetheSystem<Real> const system(5, 2);
  auto const off = engine::reduce_phantom_polynomial<Real>(system, c, d, 1, -1);
  EXPECT_EQ(off.status, engine::PhantomReductionStatus::unresolved_reduced_state);
  ASSERT_TRUE(off.finite_regularity);
  EXPECT_EQ(off.finite_regularity->status, engine::RegularityStatus::off_shell);
  engine::PolynomialBetheSystem<Real> const wrong_length(7, 2);
  EXPECT_EQ(engine::reduce_phantom_polynomial<Real>(wrong_length, c, d, 1, -1).status,
            engine::PhantomReductionStatus::phase_mismatch);
  EXPECT_EQ(engine::reduce_phantom_polynomial<Real>(system, c, d, 1, 1).status,
            engine::PhantomReductionStatus::no_endpoint_factor);
  c[0] += Real{1} / Real{1000};
  EXPECT_EQ(engine::reduce_phantom_polynomial<Real>(system, c, d, 1, -1).status,
            engine::PhantomReductionStatus::no_endpoint_factor);
  roots.push_back(Real{2} / Real{5});
  auto const cubic = polynomial<Real>(roots);
  engine::PolynomialBetheSystem<Real> const larger(7, 3);
  EXPECT_EQ(engine::reduce_phantom_polynomial<Real>(larger, cubic, d, 2, -1).status,
            engine::PhantomReductionStatus::no_endpoint_factor);
}

TYPED_TEST(XXZPhantom, UnrepresentableDeflationIsUnresolved)
{
  using Real = TypeParam;
  engine::PolynomialBetheSystem<Real> const system(5, 2, -uni20::numeric_limits<Real>::max() / Real{2});
  std::vector<Real> const c(2, Real{0});
  EXPECT_EQ(engine::reduce_phantom_polynomial<Real>(system, c, -Real{1} / Real{2}, 1, -1).status,
            engine::PhantomReductionStatus::nonfinite);
}

TYPED_TEST(XXZPhantom, InvalidInputs)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  engine::PolynomialBetheSystem<Real> const system(5, 2);
  std::vector<Real> c{Real{1}, Real{0}};
  Real const d = -Real{1} / Real{2};
  for (std::size_t p : {0, 2, 3})
    EXPECT_THROW(engine::reduce_phantom_polynomial<Real>(system, c, d, p, -1), std::invalid_argument);
  EXPECT_THROW(engine::reduce_phantom_polynomial<Real>(system, c, d, 1, 0), std::invalid_argument);
  for (Real bad :
       {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    EXPECT_THROW(engine::reduce_phantom_polynomial<Real>(system, c, d, 1, -1, bad), std::invalid_argument);
    EXPECT_THROW(engine::reduce_phantom_polynomial<Real>(system, c, d, 1, -1, Real{1}, bad), std::invalid_argument);
  }
  for (C rotation : {C{0}, C{2}, C{uni20::numeric_limits<Real>::quiet_NaN(), 0}})
  {
    EXPECT_THROW(system.evaluate_rotated(c, d, rotation), std::invalid_argument);
    EXPECT_THROW(engine::check_regular_polynomial<Real>(system, c, d, Real{1}, rotation), std::invalid_argument);
  }
  c.pop_back();
  EXPECT_THROW(engine::reduce_phantom_polynomial<Real>(system, c, d, 1, -1), std::invalid_argument);
}
} // namespace
