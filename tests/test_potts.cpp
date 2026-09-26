// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/potts.hpp>
#include <complex>
#include <uni20/linalg/ops/self_adjoint_eigh.hpp>

namespace
{
template <typename Real> class Potts : public ::testing::Test {};
TYPED_TEST_SUITE(Potts, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(Potts, NativePrecisionClockHamiltonianReferences)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  // Independent 75-digit diagonalization of literal clock charge blocks,
  // not Bethe equations. scripts/reference_potts.py --sites 3 --digits 75.
  auto const vacuum = bethe::potts::ground_state<Real>(3, {eps, 10000});
  auto const magnetic = bethe::potts::charged_one_hole<Real>(3, 1, 0, {eps, 10000});
  auto const descendant = bethe::potts::charged_one_hole<Real>(3, -1, 1, {eps, 10000});
  for (auto const* state : {&vacuum, &magnetic, &descendant})
  {
    ASSERT_TRUE(state->converged);
    EXPECT_LE(state->residual_norm, eps);
    EXPECT_EQ(state->auxiliary_roots.size(), 3);
  }
  EXPECT_REAL_NEAR(vacuum.energy,
                   uni20::parse_real<Real>("-7.6846584384264908247321147839611155377207988380604306515979503596"),
                   Real{128} * eps);
  EXPECT_REAL_NEAR(magnetic.energy,
                   uni20::parse_real<Real>("-6.9243439920202489313622513170042223878990963349288496338373686912"),
                   Real{128} * eps);
  EXPECT_REAL_NEAR(descendant.energy,
                   uni20::parse_real<Real>("-2.3775516419230208751905199448242467572796214788739660597038311918"),
                   Real{128} * eps);
  auto const tiny = bethe::potts::ground_state<Real>(2);
  ASSERT_TRUE(tiny.converged);
  EXPECT_REAL_NEAR(tiny.energy, -Real{2} - Real{2} * std::sqrt(Real{3}), Real{1024} * eps);
  auto const charged = bethe::potts::charged_one_hole<Real>(2, 1, 0);
  EXPECT_REAL_NEAR(charged.energy, -(Real{1} + std::sqrt(Real{57})) / Real{2}, Real{1024} * eps);
}

TYPED_TEST(Potts, ChargeReflectionMomentaAndDistinctHoleLabels)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  for (std::size_t l : {3, 4, 7})
    for (std::size_t h = 0; h < l; ++h)
    {
      auto const a = bethe::potts::charged_one_hole<Real>(l, 1, h);
      auto const b = bethe::potts::charged_one_hole<Real>(l, -1, h);
      auto const reflected = bethe::potts::charged_one_hole<Real>(l, 1, (l - h) % l);
      ASSERT_TRUE(a.converged);
      ASSERT_TRUE(b.converged);
      ASSERT_TRUE(reflected.converged);
      EXPECT_EQ(a.energy, b.energy);
      EXPECT_EQ(a.charge, 1);
      EXPECT_EQ(b.charge, -1);
      EXPECT_EQ(a.momentum_index, h);
      EXPECT_REAL_NEAR(a.momentum, Real{2} * pi * Real(h) / Real(l), Real{8} * eps);
      EXPECT_REAL_NEAR(a.energy, reflected.energy, Real{4096} * Real(l) * eps);
      std::int64_t sum = 0;
      for (auto i : a.auxiliary_numbers)
        sum += i.twice();
      EXPECT_EQ(sum, -2 * static_cast<std::int64_t>(h));
      for (std::size_t i = 1; i < l; ++i)
      {
        EXPECT_LT(a.auxiliary_numbers[i - 1], a.auxiliary_numbers[i]);
        EXPECT_LT(a.auxiliary_roots[i - 1], a.auxiliary_roots[i]);
      }
    }
}

TYPED_TEST(Potts, ExplicitBudgetAndInputFailures)
{
  using Real = TypeParam;
  auto const initial = bethe::potts::charged_one_hole<Real>(4, 1, 1, {Real{1e-6}, 0});
  EXPECT_FALSE(initial.converged);
  EXPECT_EQ(initial.iterations, 0);
  EXPECT_GT(initial.residual_norm, Real{1e-6});
  EXPECT_THROW((void)bethe::potts::ground_state<Real>(1), std::invalid_argument);
  EXPECT_THROW((void)bethe::potts::ground_state<Real>(std::numeric_limits<std::size_t>::max()), std::invalid_argument);
  EXPECT_THROW((void)bethe::potts::charged_one_hole<Real>(4, 0, 0), std::invalid_argument);
  EXPECT_THROW((void)bethe::potts::charged_one_hole<Real>(4, 2, 0), std::invalid_argument);
  EXPECT_THROW((void)bethe::potts::charged_one_hole<Real>(4, 1, 4), std::invalid_argument);
  EXPECT_THROW((void)bethe::potts::ground_state<Real>(4, {Real{0}, 10}), std::invalid_argument);
  EXPECT_THROW((void)bethe::potts::ground_state<Real>(4, {uni20::numeric_limits<Real>::infinity(), 10}),
               std::invalid_argument);
}

TYPED_TEST(Potts, ExponentiatedAuxiliaryEquations)
{
  using Real = TypeParam;
  using Complex = std::complex<Real>;
  Real const delta = std::sqrt(Real{3}) / Real{2}, eps = uni20::numeric_limits<Real>::epsilon();
  // Multiplicative complex equations independently check logarithm branches,
  // twist sign, root count and half-integer parity at both odd/even clock L.
  for (std::size_t l : {3, 4, 7})
    for (bool vacuum : {false, true})
    {
      auto const s = vacuum ? bethe::potts::ground_state<Real>(l) : bethe::potts::charged_one_hole<Real>(l, 1, l / 2);
      ASSERT_TRUE(s.converged);
      auto const& z = s.auxiliary_roots;
      for (std::size_t i = 0; i < l; ++i)
      {
        Complex equation(std::cos(s.auxiliary_twist), std::sin(s.auxiliary_twist));
        Complex const drive = Complex(z[i], Real{1}) / Complex(z[i], -Real{1});
        for (std::size_t n = 0; n < 2 * l; ++n)
          equation *= drive;
        for (std::size_t j = 0; j < l; ++j)
          if (i != j)
          {
            Real const a = delta * (z[i] - z[j]), b = Real{1} + delta - (Real{1} - delta) * z[i] * z[j];
            equation *= Complex(a, -b) / Complex(a, b);
          }
        EXPECT_LT(std::abs(equation - Complex(Real{1}, Real{0})), Real{512} * Real(l) * eps);
      }
    }
}

TYPED_TEST(Potts, ConformalScalingAndFiniteMomentumLine)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), velocity = bethe::potts::velocity<Real>();
  Real old_c_error = Real{1}, old_x_error = Real{1};
  for (std::size_t l : {16, 32, 64})
  {
    auto const vacuum = bethe::potts::ground_state<Real>(l);
    auto const magnetic = bethe::potts::charged_one_hole<Real>(l, 1, 0);
    auto const desc = bethe::potts::charged_one_hole<Real>(l, 1, 1);
    auto const edge = bethe::potts::charged_one_hole<Real>(l, 1, l / 2);
    ASSERT_TRUE(vacuum.converged);
    ASSERT_TRUE(magnetic.converged);
    ASSERT_TRUE(desc.converged);
    ASSERT_TRUE(edge.converged);
    Real const c =
        -Real{6} * Real(l) * (vacuum.energy - Real(l) * bethe::potts::bulk_energy_density<Real>()) / (pi * velocity);
    Real const x = Real(l) * (magnetic.energy - vacuum.energy) / (Real{2} * pi * velocity);
    Real const xd = Real(l) * (desc.energy - vacuum.energy) / (Real{2} * pi * velocity);
    auto const c_error = std::abs(c - Real{4} / Real{5}), x_error = std::abs(x - Real{2} / Real{15});
    EXPECT_LT(c_error, old_c_error);
    EXPECT_LT(x_error, old_x_error);
    EXPECT_LT(c_error, Real{0.005});
    EXPECT_LT(x_error, Real{0.001});
    EXPECT_LT(std::abs(xd - Real{17} / Real{15}), Real{0.004});
    EXPECT_LT(std::abs(edge.energy - vacuum.energy - Real{2} * velocity), Real{4} / Real(l));
    old_c_error = c_error;
    old_x_error = x_error;
  }
}

// Independent literal clock Hamiltonian, resolved into charge and translation
// orbits BEFORE diagonalization. No Bethe equations or CFT dimensions involved.
std::vector<double> clock_spectrum(std::size_t l, int charge, std::size_t momentum)
{
  using Complex = std::complex<double>;
  std::vector<std::size_t> powers(l + 1, 1);
  for (std::size_t j = 0; j < l; ++j)
    powers[j + 1] = 3 * powers[j];
  auto const dimension = powers[l];
  auto translate = [&](std::size_t s) { return (3 * s) % dimension + s / powers[l - 1]; };
  std::vector<int> orbit(dimension, -1), position(dimension, -1);
  std::vector<std::size_t> representatives, periods;
  std::vector<bool> seen(dimension, false);
  for (std::size_t s = 0; s < dimension; ++s)
  {
    int sum = 0;
    for (std::size_t j = 0; j < l; ++j)
      sum += int(s / powers[j] % 3);
    if (sum % 3 != (charge + 3) % 3 || seen[s]) continue;
    std::vector<std::size_t> cycle;
    auto t = s;
    do
    {
      cycle.push_back(t);
      seen[t] = true;
      t = translate(t);
    }
    while (t != s);
    if ((momentum * cycle.size()) % l) continue;
    for (std::size_t r = 0; r < cycle.size(); ++r)
    {
      orbit[cycle[r]] = int(representatives.size());
      position[cycle[r]] = int(r);
    }
    representatives.push_back(s);
    periods.push_back(cycle.size());
  }
  auto const n = representatives.size();
  uni20::DenseMatrix<Complex> h(n, n);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      h[i, j] = Complex{};
  double const p = 8 * std::atan(1.0) * double(momentum) / double(l);
  for (std::size_t b = 0; b < n; ++b)
  {
    auto const s = representatives[b];
    std::vector<int> digits(l);
    for (std::size_t j = 0; j < l; ++j)
    {
      digits[j] = int(s / powers[j] % 3);
      h[b, b] += digits[j] == 0 ? -2. : 1.;
    }
    for (std::size_t j = 0; j < l; ++j)
      for (int step : {1, 2})
      {
        auto const k = (j + 1) % l;
        auto const target = std::int64_t(s) + ((digits[j] + step) % 3 - digits[j]) * std::int64_t(powers[j]) +
                            ((digits[k] + 3 - step) % 3 - digits[k]) * std::int64_t(powers[k]);
        auto const a = orbit[target];
        if (a >= 0)
          h[a, b] -= std::sqrt(double(periods[b]) / double(periods[a])) * std::polar(1., p * double(position[target]));
      }
  }
  auto const result = uni20::linalg::eigh(h);
  std::vector<double> values(n);
  for (std::size_t i = 0; i < n; ++i)
    values[i] = result.eigenvalues[i];
  return values;
}

TEST(PottsExact, LiteralHamiltonianBothChargesEveryMomentum)
{
  for (std::size_t l = 2; l <= 6; ++l)
  {
    SCOPED_TRACE(l);
    auto const vacuum = bethe::potts::ground_state(l);
    ASSERT_TRUE(vacuum.converged);
    EXPECT_NEAR(vacuum.energy, clock_spectrum(l, 0, 0).front(), 2e-11);
    for (int charge : {-1, 1})
      for (std::size_t k = 0; k < l; ++k)
      {
        SCOPED_TRACE(k);
        auto const state = bethe::potts::charged_one_hole(l, charge, k);
        ASSERT_TRUE(state.converged);
        auto const exact = clock_spectrum(l, charge, k);
        ASSERT_FALSE(exact.empty());
        EXPECT_NEAR(state.energy, exact.front(), 2e-11);
      }
  }
}

TEST(PottsExact, ExtraXXZStatesMustNotBeReportedAsPottsLevels)
{
  double const root3 = std::sqrt(3.);
  auto const auxiliary = bethe::xxz::sector_ground_state(8, root3 / 2, uni20::half_int(1));
  ASSERT_TRUE(auxiliary.converged);
  // This W_{1,1} energy is also in the vacuum auxiliary module W_{0,q^2}:
  // add lambda=-infinity, shift the finite I by -1/2, and set twist=pi/3.
  // Its zero bare energy preserves E, but the Potts quotient excludes it.
  double const delta = root3 / 2, pi = 4 * std::atan(1.);
  auto roots = auxiliary.rapidities;
  roots.insert(roots.begin(), -(2 + root3)); // tanh(lambda)=-1, tan(pi/12)=2-sqrt(3)
  double const labels[] = {-2.5, -1.5, -.5, .5};
  for (std::size_t i = 0; i < 4; ++i)
  {
    double residual = 16 * std::atan(roots[i]) - 2 * pi * labels[i] - pi / 3;
    for (std::size_t j = 0; j < 4; ++j)
      if (i != j)
        residual -= 2 * std::atan2(delta * (roots[i] - roots[j]), 1 + delta - (1 - delta) * roots[i] * roots[j]);
    EXPECT_NEAR(residual, 0., 2e-12);
  }
  double const unphysical = 2 * root3 * auxiliary.energy + 2;
  double nearest = 100.;
  for (int charge : {-1, 0, 1})
    for (std::size_t k = 0; k < 4; ++k)
      for (double energy : clock_spectrum(4, charge, k))
        nearest = std::min(nearest, std::abs(energy - unphysical));
  EXPECT_GT(nearest, 0.1);
}

TEST(PottsExact, PublishedFourSiteRootPatternsAndClockDimensions)
{
  // Fukai et al., SciPost Phys. 16, 003 (2024), Table 3, rows 1,2,4.
  double const lambdas[][4] = {
      {-.1783, -.0222, .1171, .4927}, {-.5795, -.0068, .1306, .5258}, {-.5470, -.1381, .1381, .5470}};
  for (std::size_t k = 0; k < 3; ++k)
  {
    auto const state = bethe::potts::charged_one_hole(4, 1, k);
    ASSERT_TRUE(state.converged);
    for (std::size_t j = 0; j < 4; ++j)
      EXPECT_NEAR(std::atanh((2 - std::sqrt(3.)) * state.auxiliary_roots[j]), lambdas[k][j], 5.1e-5);
  }
  std::size_t charge_dimension = 3;
  for (std::size_t l = 2; l <= 5; ++l, charge_dimension *= 3)
    for (int q : {-1, 0, 1})
    {
      std::size_t count = 0;
      for (std::size_t k = 0; k < l; ++k)
        count += clock_spectrum(l, q, k).size();
      EXPECT_EQ(count, charge_dimension);
      EXPECT_LT(l, count); // one selected charged level per momentum is not all states
    }
}
} // namespace
