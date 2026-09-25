// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxz.hpp>
#include <bethe/xyz.hpp>
#include <uni20/linalg/ops/self_adjoint_eigh.hpp>

namespace
{
template <typename Real> class XYZ : public ::testing::Test {};
TYPED_TEST_SUITE(XYZ, test_support::RealTypes, test_support::PrecisionNames);

TEST(XYZExact, SmallRings)
{
  for (unsigned n : {2u, 4u, 6u, 8u})
    for (double eta : {0.1, 0.4, 0.5, 0.7, 0.9})
      for (double t : {0.2, 0.5, 2.0})
      {
        SCOPED_TRACE(::testing::Message() << n << " eta=" << eta << " t=" << t);
        auto const state = bethe::xyz::ground_state(n, eta, t);
        ASSERT_TRUE(state.converged);
        auto const j = bethe::xyz::couplings(eta, t);
        unsigned const dim = 1u << n;
        uni20::DenseMatrix<double> h(dim, dim);
        for (unsigned a = 0; a < dim; ++a)
          for (unsigned b = 0; b < dim; ++b)
            h[a, b] = 0;
        for (unsigned a = 0; a < dim; ++a)
          for (unsigned b = 0; b < n; ++b)
          {
            unsigned const c = (b + 1) % n;
            bool const equal = ((a >> b) & 1) == ((a >> c) & 1);
            h[a, a] += (equal ? j.z : -j.z) / 4;
            h[a ^ (1u << b) ^ (1u << c), a] += (j.x + (equal ? -j.y : j.y)) / 4;
          }
        auto const exact = uni20::linalg::eigh(std::move(h)).eigenvalues[0];
        EXPECT_NEAR(*state.energy, exact, 2e-11 * (1 + std::abs(exact)));
      }
}

TYPED_TEST(XYZ, OriginalEquationsAndJacobian)
{
  using Real = TypeParam;
  using Complex = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real eta : {Real{0.2}, Real{0.5}, Real{0.8}})
    for (Real t : {Real{0.3}, Real{2}})
      for (std::size_t n : {4u, 6u, 12u})
      {
        auto const state = bethe::xyz::ground_state(n, eta, t);
        ASSERT_TRUE(state.converged);
        auto ratio = [&](Complex a, Complex b) {
          return bethe::detail::theta_ratio(bethe::detail::elliptic_theta(1, a, t), 0,
                                            bethe::detail::elliptic_theta(1, b, t), 0);
        };
        for (std::size_t j = 0; j < state.roots.size(); ++j)
        {
          auto const root = state.roots[j];
          auto product = std::pow(ratio(root + eta / Real{2}, root - eta / Real{2}), int(n));
          for (std::size_t k = 0; k < state.roots.size(); ++k)
            if (j != k) product *= ratio(root - state.roots[k] - eta, root - state.roots[k] + eta);
          EXPECT_REAL_NEAR(std::abs(product - Complex(1, 0)), Real{0}, Real{16384} * eps);
        }
        bethe::xyz::detail::GroundSystem<Real> system{n, eta, t};
        std::vector<Real> x;
        for (auto root : state.roots)
          if (root.imag() > Real{0}) x.push_back(root.imag());
        uni20::DenseMatrix<Real> jac(x.size(), x.size());
        system.evaluate(x, &jac);
        Real const delta = std::cbrt(eps) * t / Real{10};
        for (std::size_t col = 0; col < x.size(); ++col)
        {
          auto plus = x, minus = x;
          plus[col] += delta;
          minus[col] -= delta;
          auto const fp = system.evaluate(plus), fm = system.evaluate(minus);
          for (std::size_t row = 0; row < x.size(); ++row)
          {
            Real const numerical = (fp.residual[row] - fm.residual[row]) / (Real{2} * delta);
            EXPECT_REAL_NEAR((jac[row, col]), numerical,
                             Real{100000} * delta * delta * (Real{1} + std::abs(numerical)));
          }
        }
      }
}
TYPED_TEST(XYZ, LimitsAndFailureContracts)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  for (Real eta : {Real{0.1}, Real{0.4}, Real{0.5}, Real{0.7}, Real{0.9}})
    for (std::size_t n : {2u, 4u, 6u, 16u, 32u})
    {
      auto const state = bethe::xyz::ground_state(n, eta, Real{100});
      ASSERT_TRUE(state.converged);
      auto const xxz = bethe::xxz::ground_state(n, std::cos(pi * eta));
      ASSERT_TRUE(xxz.converged);
      EXPECT_REAL_NEAR(*state.energy, xxz.energy, Real{32768} * eps * Real(n));
    }
  for (Real eta : {Real{0.1}, Real{0.9}})
  {
    auto const two = bethe::xyz::ground_state(2, eta, Real{0.5}, {.max_iterations = 0});
    ASSERT_TRUE(two.converged);
    EXPECT_REAL_NEAR(*two.energy, -(two.exchange.x + two.exchange.y + two.exchange.z) / Real{2}, Real{4096} * eps);
    auto const failed = bethe::xyz::ground_state(8, eta, Real{1}, {.max_iterations = 0});
    EXPECT_FALSE(failed.converged);
    EXPECT_FALSE(failed.energy);
    EXPECT_EQ(failed.status, bethe::xyz::Status::iteration_limit);
  }
  EXPECT_THROW(bethe::xyz::ground_state(3, Real{0.3}, Real{1}), std::invalid_argument);
  EXPECT_THROW(bethe::xyz::ground_state(4, Real{1}, Real{1}), std::invalid_argument);
  EXPECT_THROW(bethe::xyz::ground_state(4, Real{0.3}, Real{0}), std::invalid_argument);
  EXPECT_THROW(bethe::xyz::ground_state(4, Real{0.3}, Real{1}, {.residual_tolerance = Real{0}}), std::invalid_argument);
}
} // namespace
