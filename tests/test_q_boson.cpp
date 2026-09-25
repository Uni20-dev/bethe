// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/lieb_liniger.hpp>
#include <bethe/q_boson.hpp>
#include <complex>
#include <map>
#include <uni20/linalg/ops/self_adjoint_eigh.hpp>

namespace
{
namespace model = bethe::q_boson;
template <typename Real> class QBoson : public ::testing::Test {};
TYPED_TEST_SUITE(QBoson, test_support::RealTypes, test_support::PrecisionNames);

// Occupation-space Hamiltonian, independent of any Bethe roots or phases.
double exact_ground(unsigned sites, unsigned particles, double eta)
{
  std::vector<std::vector<unsigned>> basis;
  std::vector<unsigned> occupation(sites);
  auto enumerate = [&](auto&& self, unsigned site, unsigned left) -> void {
    if (site + 1 == sites)
    {
      occupation[site] = left;
      basis.push_back(occupation);
      return;
    }
    for (unsigned n = 0; n <= left; ++n)
    {
      occupation[site] = n;
      self(self, site + 1, left - n);
    }
  };
  enumerate(enumerate, 0, particles);
  std::map<std::vector<unsigned>, std::size_t> index;
  for (std::size_t i = 0; i < basis.size(); ++i)
    index[basis[i]] = i;
  auto number = [&](unsigned n) {
    if (eta == 0) return double(n);
    // Finite geometric sum is independent of the production scattering code.
    double value = 0, power = 1;
    for (unsigned j = 0; j < n; ++j)
    {
      value += power;
      power *= std::exp(-2 * eta);
    }
    return value;
  };
  uni20::DenseMatrix<double> h(basis.size(), basis.size());
  for (std::size_t row = 0; row < basis.size(); ++row)
    for (std::size_t col = 0; col < basis.size(); ++col)
      h[row, col] = row == col ? 2 * particles : 0;
  for (std::size_t col = 0; col < basis.size(); ++col)
    for (unsigned j = 0; j < sites; ++j)
      for (unsigned direction = 0; direction < 2; ++direction)
      {
        auto const source = direction ? j : (j + 1) % sites;
        auto const target = direction ? (j + 1) % sites : j;
        auto state = basis[col];
        if (state[source] == 0) continue;
        double const hop = std::sqrt(number(state[source]) * number(state[target] + 1));
        --state[source];
        ++state[target];
        h[index.at(state), col] -= hop;
      }
  return uni20::linalg::eigh(std::move(h)).eigenvalues[0];
}

TYPED_TEST(QBoson, OccupationBasisGroundEnergies)
{
  using Real = TypeParam;
  for (unsigned sites = 2; sites <= 5; ++sites)
    for (unsigned n = 0; n <= 4; ++n)
      for (Real eta : {Real{0}, Real{0.1}, Real{1}, Real{5}, uni20::numeric_limits<Real>::infinity()})
      {
        auto const s = model::ground_state(sites, n, eta);
        ASSERT_TRUE(s.converged) << sites << ',' << n << ',' << double(eta) << " status=" << int(s.status);
        ASSERT_TRUE(s.energy);
        EXPECT_NEAR(double(*s.energy), exact_ground(sites, n, double(eta)), 1e-11);
      }
}

TYPED_TEST(QBoson, TwoSiteNativeFormulaAndWeakDeformation)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real eta : {eps * eps, Real{0.001}, Real{1}, Real{10}})
  {
    auto const s = model::ground_state(2, 2, eta);
    ASSERT_TRUE(s.converged) << int(s.status);
    Real const loss = -std::expm1(-Real{2} * eta);
    Real const expected = Real{2} * loss / (Real{1} + std::sqrt(Real{1} - loss / Real{2}));
    EXPECT_REAL_NEAR(*s.energy, expected, Real{512} * eps * expected);
    EXPECT_GT(*s.energy, Real{0});
    if constexpr (uni20::numeric_limits<Real>::digits > 53)
      if (eta == Real{1}) EXPECT_GT(std::abs(*s.energy - Real(double(*s.energy))), Real{4} * eps);
  }
  auto const many = model::ground_state(20, 12, eps * eps);
  ASSERT_TRUE(many.converged) << int(many.status);
  Real const leading = Real{2} * eps * eps * Real{12} * Real{11} / Real{20};
  EXPECT_REAL_NEAR(*many.energy, leading, Real{4096} * eps * leading);
}

TYPED_TEST(QBoson, OriginalMultiplicativeEquationsAndPhaseLimit)
{
  using Real = TypeParam;
  using Complex = std::complex<Real>;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (Real eta : {Real{0.01}, Real{0.5}, Real{2}})
    for (std::size_t sites : {2u, 7u})
    {
      auto const s = model::ground_state(sites, 9, eta);
      ASSERT_TRUE(s.converged) << int(s.status);
      Real momentum{};
      for (std::size_t j = 0; j < s.particles; ++j)
      {
        Real const p = s.momenta[j];
        Complex product = std::exp(Complex(0, Real(s.sites) * p));
        for (std::size_t k = 0; k < s.particles; ++k)
          if (k != j)
          {
            Real const half = (p - s.momenta[k]) / Real{2};
            product *= std::sin(Complex(half, -eta)) / std::sin(Complex(half, eta));
          }
        EXPECT_REAL_NEAR(std::abs(product - Complex(1, 0)), Real{0}, Real{4096} * eps);
        momentum += p;
        if (j) EXPECT_GT(p, s.momenta[j - 1]);
      }
      EXPECT_REAL_NEAR(momentum, Real{0}, Real{4096} * eps);
    }
  auto const phase = model::ground_state(7, 9, uni20::numeric_limits<Real>::infinity(), {.max_iterations = 0});
  ASSERT_TRUE(phase.converged);
  EXPECT_EQ(phase.iterations, 0u);
  Real const expected = Real{18} - Real{2} * std::sin(Real{9} * pi / Real{16}) / std::sin(pi / Real{16});
  EXPECT_REAL_NEAR(*phase.energy, expected, Real{512} * eps * expected);
}

TYPED_TEST(QBoson, ContinuumScalingAndLargerRings)
{
  using Real = TypeParam;
  auto const continuum = bethe::lieb_liniger::ground_state(3, Real{3}, Real{2});
  ASSERT_TRUE(continuum.converged);
  Real previous = Real{1};
  for (std::size_t sites : {24u, 48u, 96u})
  {
    Real const a = Real{3} / Real(sites);
    // eta=c*a/2; shifted lattice energy divided by a^2.
    auto const s = model::ground_state(sites, 3, a);
    ASSERT_TRUE(s.converged);
    Real const error = std::abs(*s.energy / (a * a) - continuum.energy);
    EXPECT_LT(error, previous);
    previous = error;
  }
  for (Real eta : {Real{0.001}, Real{0.3}, Real{3}})
    for (std::size_t sites : {2u, 32u, 128u})
    {
      auto const s = model::ground_state(sites, 32, eta);
      ASSERT_TRUE(s.converged) << int(s.status);
      EXPECT_GT(*s.energy, Real{0});
      EXPECT_LT(*s.energy, Real{64});
    }
}

TYPED_TEST(QBoson, ValidationAndBudgets)
{
  using Real = TypeParam;
  EXPECT_THROW(model::ground_state(1, 2, Real{1}), std::invalid_argument);
  EXPECT_THROW(model::ground_state(4, 2, Real{-1}), std::invalid_argument);
  EXPECT_THROW(model::ground_state(4, 2, uni20::numeric_limits<Real>::quiet_NaN()), std::invalid_argument);
  EXPECT_THROW(model::ground_state(4, 2, Real{1}, {.residual_tolerance = Real{0}}), std::invalid_argument);
  auto const failed = model::ground_state(4, 3, Real{1}, {.max_iterations = 0});
  EXPECT_FALSE(failed.converged);
  EXPECT_FALSE(failed.energy);
  EXPECT_EQ(failed.status, model::Status::iteration_limit);
  auto const tiny =
      model::ground_state(4, 2, uni20::numeric_limits<Real>::min() * uni20::numeric_limits<Real>::epsilon());
  EXPECT_FALSE(tiny.converged);
  EXPECT_FALSE(tiny.energy);
  EXPECT_EQ(tiny.status, model::Status::precision_limit);
}
} // namespace
