// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/sutherland.hpp>
#include <complex>
#include <set>

namespace
{
namespace model = bethe::sutherland;
template <typename Real> class Sutherland : public ::testing::Test {};
TYPED_TEST_SUITE(Sutherland, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(Sutherland, FreeBosonsHardCoreAndOneParticle)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), length = Real{7}, q = Real{2} * pi / length;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n = 1; n <= 9; ++n)
  {
    std::vector<std::int64_t> labels(n);
    for (std::size_t j = 0; j < n; ++j)
      labels[j] = std::int64_t(j / 2) - 2;
    for (Real lambda : {Real{0}, Real{1}})
    {
      auto const s = model::evaluate(n, length, lambda, labels);
      Real energy{0}, momentum{0};
      for (std::size_t j = 0; j < n; ++j)
      {
        // Hard-core boson dual: integer fermion momenta for odd N,
        // half-integer (antiperiodic) for even N. Labels themselves may repeat.
        Real const k = q * (Real(labels[j]) + lambda * Real(2 * std::int64_t(j) + 1 - std::int64_t(n)) / Real{2});
        energy += k * k;
        momentum += k;
        EXPECT_REAL_NEAR(s.pseudomomenta[j], k, Real{32} * eps * (Real{1} + std::abs(k)));
      }
      EXPECT_REAL_NEAR(s.energy, energy, Real{64} * eps * (Real{1} + energy));
      EXPECT_REAL_NEAR(s.momentum, momentum, Real{64} * eps * (Real{1} + std::abs(momentum)));
      if (lambda == Real{0}) EXPECT_EQ(s.ground_energy, Real{0});
    }
  }
  for (Real lambda : {Real{0}, Real{0.25}, Real{1}, Real{100}})
  {
    auto const s = model::evaluate(1, length, lambda, {-3});
    EXPECT_EQ(s.ground_energy, Real{0});
    EXPECT_REAL_NEAR(s.energy, Real{9} * q * q, Real{32} * eps * q * q);
  }
}

TYPED_TEST(Sutherland, BoostReflectionScalingAndSmallGap)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  Real const length = Real{5}, lambda = Real{2.25}, q = Real{2} * pi / length;
  auto const s = model::evaluate(4, length, lambda, {-2, 0, 0, 3});
  auto const reflection = model::evaluate(4, length, lambda, {-3, 0, 0, 2});
  EXPECT_REAL_NEAR(s.energy, reflection.energy, Real{32} * eps * s.energy);
  EXPECT_EQ(s.momentum_index, -reflection.momentum_index);
  auto const boost = model::evaluate(4, length, lambda, {1, 3, 3, 6});
  Real const expected = s.energy + Real{6} * q * s.momentum + Real{36} * q * q;
  EXPECT_REAL_NEAR(boost.energy, expected, Real{64} * eps * expected);
  EXPECT_EQ(boost.momentum_index, s.momentum_index + 12);
  auto const scaled = model::evaluate(4, Real{2} * length, lambda, s.labels);
  EXPECT_REAL_NEAR(scaled.energy, s.energy / Real{4}, Real{32} * eps * s.energy);
  EXPECT_REAL_NEAR(scaled.momentum, s.momentum / Real{2}, Real{32} * eps * std::abs(s.momentum));
  // The boost gap survives even if addition to a huge E0 rounds it away.
  auto const huge = model::evaluate(4, length, uni20::parse_real<Real>("1e20"), {1, 1, 1, 1});
  EXPECT_REAL_NEAR(huge.gap, Real{4} * q * q, Real{32} * eps * q * q);
  EXPECT_EQ(huge.energy, huge.ground_energy);
}

TYPED_TEST(Sutherland, CompleteWindowEnumerationAndBudgets)
{
  using Real = TypeParam;
  auto const all = model::spectrum(3, Real{4}, Real{2}, 2, {}, 35);
  ASSERT_TRUE(all.complete);
  ASSERT_EQ(all.total_states, 35);
  ASSERT_EQ(all.states.size(), 35);
  EXPECT_EQ(all.states[0].labels, (std::vector<std::int64_t>{0, 0, 0}));
  std::set<std::vector<std::int64_t>> labels;
  for (std::size_t j = 0; j < all.states.size(); ++j)
  {
    auto const& s = all.states[j];
    EXPECT_TRUE(labels.insert(s.labels).second);
    EXPECT_TRUE(std::is_sorted(s.labels.begin(), s.labels.end()));
    EXPECT_GE(s.labels.front(), -2);
    EXPECT_LE(s.labels.back(), 2);
    if (j) EXPECT_LE(all.states[j - 1].gap, s.gap);
    std::vector<std::int64_t> reflected;
    for (auto it = s.labels.rbegin(); it != s.labels.rend(); ++it)
      reflected.push_back(-*it);
    auto const partner =
        std::find_if(all.states.begin(), all.states.end(), [&](auto const& r) { return r.labels == reflected; });
    ASSERT_NE(partner, all.states.end());
    EXPECT_EQ(partner->gap, s.gap);
    EXPECT_EQ(partner->momentum_index, -s.momentum_index);
  }
  auto const few = model::spectrum(3, Real{4}, Real{2}, 2, 4);
  ASSERT_EQ(few.states.size(), 4);
  EXPECT_EQ(few.states[3].labels, all.states[3].labels);
  auto const no = model::spectrum(3, Real{4}, Real{2}, 2, {}, 34);
  EXPECT_FALSE(no.complete);
  EXPECT_TRUE(no.states.empty());
  EXPECT_FALSE(model::spectrum(3, Real{4}, Real{2}, 0, {}, 0).complete);
  EXPECT_EQ(model::spectrum(3, Real{4}, Real{2}, 0).states.size(), 1);
  EXPECT_FALSE(model::spectrum(model::max_particles, Real{4}, Real{2}, 1000000).complete);
  // Independent stars-and-bars counts over a grid of small windows.
  for (std::size_t n = 1; n <= 7; ++n)
    for (std::size_t w = 0; w <= 3; ++w)
    {
      std::size_t count = 1;
      for (std::size_t j = 1; j <= n; ++j)
        count = count * (2 * w + j) / j;
      EXPECT_EQ(model::spectrum(n, Real{4}, Real{0}, w).states.size(), count);
    }
}

TYPED_TEST(Sutherland, InvalidInputsAndScalarRange)
{
  using Real = TypeParam;
  Real const inf = uni20::numeric_limits<Real>::infinity(), nan = uni20::numeric_limits<Real>::quiet_NaN();
  for (auto n : {std::size_t{0}, model::max_particles + 1})
    EXPECT_THROW(model::ground_state(n, Real{4}, Real{2}), std::invalid_argument);
  for (Real length : {Real{0}, Real{-1}, inf, nan})
    EXPECT_THROW(model::ground_state(3, length, Real{2}), std::invalid_argument);
  for (Real lambda : {Real{-1}, inf, nan})
    EXPECT_THROW(model::ground_state(3, Real{4}, lambda), std::invalid_argument);
  for (auto labels : {std::vector<std::int64_t>{}, {0, 1}, {0, 2, 1}})
    EXPECT_THROW(model::evaluate(3, Real{4}, Real{2}, labels), std::invalid_argument);
  EXPECT_THROW(model::spectrum(3, Real{4}, Real{2}, 1000001), std::invalid_argument);
  EXPECT_THROW(model::spectrum(3, Real{4}, Real{2}, 2, 0), std::invalid_argument);
  auto const lo = std::numeric_limits<std::int64_t>::min(), hi = std::numeric_limits<std::int64_t>::max();
  EXPECT_THROW(model::evaluate(2, Real{4}, Real{2}, {hi, hi}), std::overflow_error);
  EXPECT_THROW(model::evaluate(2, Real{4}, Real{2}, {lo, lo}), std::overflow_error);
  // Prefix overflow is not total overflow; the adjacent difference exceeds INT64_MAX.
  auto const cancellation = model::evaluate(4, Real{4}, Real{2}, {lo, lo, hi, hi});
  EXPECT_EQ(cancellation.momentum_index, -2);
  EXPECT_TRUE(uni20::isfinite(cancellation.energy));
  Real const max = uni20::numeric_limits<Real>::max();
  EXPECT_THROW(model::ground_state(3, Real{4}, max), std::overflow_error);
  EXPECT_THROW(model::ground_state(3, max, Real{1}), std::underflow_error);
  EXPECT_THROW(model::evaluate(1, max, Real{0}, {1}), std::underflow_error);
}

// Independent coordinate-space checks: apply the differential Hamiltonian to
// Jastrow and degree-two Jack wavefunctions, without using pseudomomentum rules.
TYPED_TEST(Sutherland, JastrowAndJackHamiltonianOracles)
{
  using Real = TypeParam;
  using Complex = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  Real const length = Real{7}, a = pi / length, q = Real{2} * a;
  for (std::size_t n = 2; n <= 6; ++n)
    for (Real lambda : {Real{0}, Real{0.25}, Real{0.5}, Real{1}, Real{2}, Real{3.25}})
      for (Real perturb : {Real{0.125}, Real{0.25}})
      {
        SCOPED_TRACE(::testing::Message() << n << ',' << uni20::format_real(lambda));
        std::vector<Real> x(n), gradient(n, Real{0});
        std::vector<Complex> z(n);
        for (std::size_t j = 0; j < n; ++j)
        {
          x[j] = length * (Real(j + 1) + perturb * Real(j * j) / Real(n + 1)) / Real(n + 2);
          z[j] = {std::cos(q * x[j]), std::sin(q * x[j])};
        }
        Real local_energy{0};
        for (std::size_t i = 0; i < n; ++i)
        {
          Real log_second{0};
          for (std::size_t j = 0; j < n; ++j)
            if (i != j)
            {
              Real const y = a * (x[i] - x[j]), sine = std::sin(y);
              gradient[i] += lambda * a * std::cos(y) / sine;
              log_second -= lambda * a * a / (sine * sine);
              if (i < j) local_energy += Real{2} * lambda * (lambda - Real{1}) * a * a / (sine * sine);
            }
          local_energy -= gradient[i] * gradient[i] + log_second;
        }
        auto const ground = model::ground_state(n, length, lambda);
        EXPECT_REAL_NEAR(local_energy, ground.energy, Real{512} * eps * (Real{1} + ground.energy) * Real(n));
        Real const c = Real{2} * lambda / (Real{1} + lambda);
        Complex polynomial{0}, kinetic{0}, drift{0};
        for (std::size_t i = 0; i < n; ++i)
        {
          Complex cross{0};
          for (std::size_t j = 0; j < n; ++j)
            if (i != j) cross += z[i] * z[j];
          polynomial += z[i] * z[i] + c * cross / Real{2};
          kinetic += q * q * (Real{4} * z[i] * z[i] + c * cross);
          drift -= Real{2} * gradient[i] * Complex{0, q} * (Real{2} * z[i] * z[i] + c * cross);
        }
        std::vector<std::int64_t> labels(n, 0);
        labels.back() = 2;
        auto const state = model::evaluate(n, length, lambda, labels);
        Complex const actual = local_energy * polynomial + kinetic + drift, expected = state.energy * polynomial;
        Real const tolerance = Real{2048} * eps * (Real{1} + state.energy) * Real(n * n);
        EXPECT_REAL_NEAR(actual.real(), expected.real(), tolerance);
        EXPECT_REAL_NEAR(actual.imag(), expected.imag(), tolerance);
      }
}

template <typename Real> Real gegenbauer(int n, Real lambda, Real x)
{
  if (n < 0) return Real{0};
  if (n == 0) return Real{1};
  Real previous{1}, current = Real{2} * lambda * x;
  for (int j = 2; j <= n; ++j)
  {
    Real next =
        (Real{2} * (Real(j - 1) + lambda) * x * current - (Real(j - 2) + Real{2} * lambda) * previous) / Real(j);
    previous = current;
    current = next;
  }
  return current;
}

TYPED_TEST(Sutherland, TwoParticleExcitedWavefunctions)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  Real const length = Real{5}, factor = Real{2} * (pi / length) * (pi / length);
  for (Real lambda : {Real{0.25}, Real{0.5}, Real{1}, Real{2}, Real{3.25}})
    for (int d = 0; d <= 8; ++d)
      for (Real y : {pi / Real{5}, pi * Real{3} / Real{7}})
      {
        Real const sine = std::sin(y), cosine = std::cos(y), cot = cosine / sine;
        Real const g = gegenbauer(d, lambda, cosine);
        Real const first = Real{2} * lambda * gegenbauer(d - 1, lambda + Real{1}, cosine);
        Real const second = Real{4} * lambda * (lambda + Real{1}) * gegenbauer(d - 2, lambda + Real{2}, cosine);
        Real const f_second = (lambda * lambda * cot * cot - lambda / (sine * sine)) * g -
                              (Real{2} * lambda + Real{1}) * cosine * first + sine * sine * second;
        // psi=exp(i*pi*M*(x1+x2)/L)*sin(y)^lambda*C_d^lambda(cos(y)).
        // Strip the nonzero Jastrow/COM factor, but do not divide by g at nodes.
        int const m = d - 6;
        Real const actual = factor * (Real(m * m) * g - f_second + lambda * (lambda - Real{1}) * g / (sine * sine));
        auto const state = model::evaluate(2, length, lambda, {-3, d - 3});
        Real const expected = state.energy * g;
        EXPECT_REAL_NEAR(actual, expected, Real{1024} * eps * (Real{1} + state.energy) * (Real{1} + std::abs(g)));
      }
}
} // namespace
