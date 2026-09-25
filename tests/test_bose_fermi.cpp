// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/bose_fermi.hpp>
#include <bethe/gaudin_yang.hpp>
#include <bit>
#include <map>
#include <uni20/linalg/ops/self_adjoint_eigh.hpp>

namespace
{
namespace model = bethe::bose_fermi;
template <typename Real> class BoseFermi : public ::testing::Test {};
TYPED_TEST_SUITE(BoseFermi, test_support::RealTypes, test_support::PrecisionNames);

// Independent Galerkin Hamiltonian: two bosons, three spinless fermions, P=0.
// Cutoff energies are variational upper bounds, not exact finite-basis BA values.
double plane_wave_ground(int cutoff, double c)
{
  using Key = std::pair<std::vector<int>, unsigned>;
  int const modes = 2 * cutoff + 1;
  std::vector<Key> basis;
  std::map<Key, std::size_t> index;
  for (int a = 0; a < modes; ++a)
    for (int b = a; b < modes; ++b)
      for (unsigned mask = 0; mask < (1u << modes); ++mask)
        if (std::popcount(mask) == 3)
        {
          int momentum = a + b - 2 * cutoff;
          for (int j = 0; j < modes; ++j)
            if (mask & (1u << j)) momentum += j - cutoff;
          if (momentum) continue;
          std::vector<int> bosons(modes);
          ++bosons[a];
          ++bosons[b];
          Key key{bosons, mask};
          index[key] = basis.size();
          basis.push_back(key);
        }
  uni20::DenseMatrix<double> h(basis.size(), basis.size());
  for (std::size_t i = 0; i < basis.size(); ++i)
    for (std::size_t j = 0; j < basis.size(); ++j)
      h[i, j] = 0;
  double const pi = 4 * std::atan(1.0);
  for (std::size_t column = 0; column < basis.size(); ++column)
  {
    auto const& [occupation, mask] = basis[column];
    for (int j = 0; j < modes; ++j)
      h[column, column] += 4 * pi * pi * (j - cutoff) * (j - cutoff) * (occupation[j] + bool(mask & (1u << j)));
    for (int r = 0; r < modes; ++r)
      for (int s = 0; s < modes; ++s)
        for (int p = 0; p < modes; ++p)
        {
          int const q = r + s - p;
          if (q < 0 || q >= modes) continue;
          // c b_p^dagger b_q^dagger b_s b_r (L=1).
          auto out = occupation;
          if (out[r])
          {
            double factor = std::sqrt(double(out[r]--));
            if (out[s])
            {
              factor *= std::sqrt(double(out[s]--));
              factor *= std::sqrt(double(++out[q]));
              factor *= std::sqrt(double(++out[p]));
              h[index.at({out, mask}), column] += c * factor;
            }
          }
          // 2c b_p^dagger f_q^dagger f_s b_r; preserve fermion signs.
          out = occupation;
          if (!out[r] || !(mask & (1u << s))) continue;
          double factor = std::sqrt(double(out[r]--));
          factor *= std::sqrt(double(++out[p]));
          unsigned const removed = mask ^ (1u << s);
          if (removed & (1u << q)) continue;
          int const parity = std::popcount(mask & ((1u << s) - 1)) + std::popcount(removed & ((1u << q) - 1));
          h[index.at({out, removed | (1u << q)}), column] += 2 * c * factor * (parity % 2 ? -1 : 1);
        }
  }
  for (std::size_t i = 0; i < basis.size(); ++i)
    for (std::size_t j = 0; j < basis.size(); ++j)
      EXPECT_NEAR((h[i, j]), (h[j, i]), 1e-12);
  return uni20::linalg::eigh(std::move(h)).eigenvalues[0];
}

TEST(BoseFermiED, IndependentMixedHamiltonian)
{
  for (double c : {0.5, 2.0})
  {
    auto const state = model::ground_state(2, 3, 1.0, c);
    ASSERT_TRUE(state.converged);
    double previous = 100;
    for (int cutoff : {2, 3, 4})
    {
      double const error = plane_wave_ground(cutoff, c) - *state.energy;
      EXPECT_GT(error, 0);
      EXPECT_LT(error, previous * 0.9);
      previous = error;
    }
    EXPECT_LT(previous, 0.5 * c * c);
  }
}

TYPED_TEST(BoseFermi, ReductionsAndFreeLimits)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  for (std::size_t b = 1; b <= 6; ++b)
    for (Real c : {Real{0.001}, Real{1}, Real{1000}})
    {
      auto const state = model::ground_state(b, 1, Real{1}, c);
      ASSERT_TRUE(state.converged) << b << ',' << double(c) << ',' << int(state.status) << ','
                                   << double(state.residual_norm);
      auto const ll = bethe::lieb_liniger::ground_state(b + 1, Real{1}, c);
      ASSERT_TRUE(ll.converged);
      EXPECT_REAL_NEAR(*state.energy, ll.energy, Real{2048} * eps * ll.energy);
      auto const pure = model::ground_state(b + 1, 0, Real{1}, c);
      ASSERT_TRUE(pure.converged);
      EXPECT_EQ(*pure.energy, ll.energy);
    }
  for (std::size_t f : {1u, 3u, 5u})
  {
    auto const state = model::ground_state(1, f, Real{2}, Real{3});
    auto const gy = bethe::gaudin_yang::ground_state(f, 1, Real{2}, Real{3});
    ASSERT_TRUE(state.converged && gy.converged);
    EXPECT_REAL_NEAR(*state.energy, gy.energy, Real{2048} * eps * gy.energy);
  }
  for (std::size_t b = 0; b < 4; ++b)
    for (std::size_t f = 0; f < 6; ++f)
    {
      auto const state = model::ground_state(b, f, Real{2}, Real{0}, {.max_iterations = 0});
      ASSERT_TRUE(state.converged);
      Real const exact = pi * pi * Real(f) * (Real(f) * Real(f) + (f % 2 ? -Real{1} : Real{2})) / Real{12};
      EXPECT_REAL_NEAR(*state.energy, exact, Real{64} * eps * (Real{1} + exact));
      EXPECT_EQ(state.momentum_index, f % 2 ? 0 : f / 2);
      auto const pure = model::ground_state(0, f, Real{2}, Real{5});
      EXPECT_EQ(pure.energy, state.energy);
    }
}

TYPED_TEST(BoseFermi, EquationsScalingAndSweep)
{
  using Real = TypeParam;
  using Complex = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t b : {1u, 2u, 3u, 6u, 12u})
    for (std::size_t f : {1u, 3u, 7u})
      for (Real c : {Real{0.001}, Real{1}, Real{1000}})
      {
        SCOPED_TRACE(::testing::Message() << b << ',' << f << ',' << double(c));
        auto const state = model::ground_state(b, f, Real{1}, c);
        ASSERT_TRUE(state.converged) << int(state.status) << ',' << double(state.residual_norm);
        EXPECT_EQ(state.momentum, Real{0});
        for (Real k : state.momenta)
        {
          Complex product{1, 0};
          for (Real l : state.auxiliary)
            product *= Complex(k - l, c / Real{2}) / Complex(k - l, -c / Real{2});
          EXPECT_REAL_NEAR(std::abs(product - std::exp(Complex(0, k))), Real{0}, Real{8192} * eps * Real(b + f));
        }
        for (Real l : state.auxiliary)
        {
          Complex product{1, 0};
          for (Real k : state.momenta)
            product *= Complex(k - l, c / Real{2}) / Complex(k - l, -c / Real{2});
          EXPECT_REAL_NEAR(std::abs(product - Complex(1, 0)), Real{0}, Real{8192} * eps * Real(b + f));
        }
        auto const scaled = model::ground_state(b, f, Real{2}, c / Real{2});
        ASSERT_TRUE(scaled.converged);
        EXPECT_REAL_NEAR(*scaled.energy * Real{4}, *state.energy, Real{1024} * eps * *state.energy);
      }
}

TYPED_TEST(BoseFermi, LimitsAndBudgets)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  auto const analytic = model::ground_state(1, 1, Real{1}, pi);
  ASSERT_TRUE(analytic.converged);
  EXPECT_REAL_NEAR(*analytic.energy, pi * pi / Real{2}, Real{512} * eps);
  if constexpr (uni20::numeric_limits<Real>::digits > 53)
    EXPECT_GT(std::abs(Real(double(pi * pi / Real{2})) - pi * pi / Real{2}), Real{512} * eps);
  for (std::size_t b : {1u, 2u, 5u})
  {
    auto const tiny = model::ground_state(b, 1, Real{1}, eps * eps);
    ASSERT_TRUE(tiny.converged) << b << ',' << int(tiny.status);
    Real const leading = eps * eps * Real(b) * Real(b + 1);
    EXPECT_REAL_NEAR(*tiny.energy, leading, Real{8192} * eps * leading);
    auto const free = model::ground_state(b, 3, Real{1}, Real{0});
    Real previous = Real{1};
    for (Real c : {Real{0.01}, Real{0.001}, Real{0.0001}})
    {
      auto const weak = model::ground_state(b, 3, Real{1}, c);
      ASSERT_TRUE(weak.converged);
      Real const error = std::abs((*weak.energy - *free.energy) / (c * Real(b) * Real(b - 1 + 6)) - Real{1});
      EXPECT_LT(error, previous / Real{3});
      previous = error;
    }
    auto const strong = model::ground_state(b, 3, Real{1}, Real{1000000});
    ASSERT_TRUE(strong.converged);
    Real const n = Real(b + 3), limit = pi * pi * n * (n * n - Real{1}) / Real{3};
    EXPECT_LT(*strong.energy, limit);
    EXPECT_LT((limit - *strong.energy) / limit, Real{0.0001});
  }
  EXPECT_THROW(model::ground_state(2, 2, Real{1}, Real{1}), std::invalid_argument);
  EXPECT_THROW(model::ground_state(2, 3, Real{0}, Real{1}), std::invalid_argument);
  EXPECT_THROW(model::ground_state(2, 3, Real{1}, Real{-1}), std::invalid_argument);
  EXPECT_THROW(model::ground_state(2, 3, Real{1}, Real{1}, {.residual_tolerance = Real{0}}), std::invalid_argument);
  auto const failed = model::ground_state(2, 3, Real{1}, Real{1}, {.max_iterations = 0});
  EXPECT_FALSE(failed.converged);
  EXPECT_FALSE(failed.energy);
  EXPECT_EQ(failed.status, model::Status::iteration_limit);
  EXPECT_THROW(model::ground_state(2, 3, Real{1}, uni20::numeric_limits<Real>::quiet_NaN()), std::invalid_argument);
  EXPECT_THROW(model::ground_state(std::numeric_limits<std::size_t>::max(), 3, Real{1}, Real{1}), std::length_error);
  auto const overflow = model::ground_state(0, 3, uni20::numeric_limits<Real>::min(), Real{0});
  EXPECT_FALSE(overflow.converged);
  EXPECT_FALSE(overflow.energy);
  EXPECT_EQ(overflow.status, model::Status::precision_limit);
}

TYPED_TEST(BoseFermi, AnalyticJacobian)
{
  using Real = TypeParam;
  Real const base = std::cbrt(uni20::numeric_limits<Real>::epsilon());
  for (std::size_t b : {2u, 3u})
    for (Real g : {Real{0.01}, Real{1}, Real{100}})
    {
      model::detail::System<Real> system{b + 3, b, g};
      auto const x = system.seed();
      uni20::DenseMatrix<Real> jacobian(system.order(), system.order());
      system.evaluate(x, &jacobian);
      for (std::size_t col = 0; col < x.size(); ++col)
      {
        Real const h = base * (Real{1} + std::abs(x[col]));
        auto plus = x, minus = x;
        plus[col] += h;
        minus[col] -= h;
        auto const fp = system.evaluate(plus), fm = system.evaluate(minus);
        for (std::size_t row = 0; row < x.size(); ++row)
        {
          Real const numerical = (fp.residual[row] - fm.residual[row]) / (Real{2} * h);
          EXPECT_REAL_NEAR((jacobian[row, col]), numerical, Real{10000} * h * h * (Real{1} + std::abs(numerical)));
        }
      }
    }
}
} // namespace
