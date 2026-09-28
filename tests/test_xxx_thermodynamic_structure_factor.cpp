// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxx_thermodynamic_structure_factor.hpp>

namespace model = bethe::heisenberg;
template <typename T> class XXXThermodynamicDSF : public ::testing::Test {};
TYPED_TEST_SUITE(XXXThermodynamicDSF, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXXThermodynamicDSF, IndependentHighPrecisionKernelReferences)
{
  using Real = TypeParam;
  // Independently evaluated with mpmath, 65 digits, split tanh-sinh over
  // [0,160]; different analytic tail subtraction (see the reference script).
  char const* points[][2] = {{"0.0001", "15.52194620516074847070346727156680040926022021782178"},
                             {"0.1", "1.717064850945922729080388492434794963237911883094265"},
                             {"0.5", "-1.466564930537033766662372479033758100468791516201768"},
                             {"1", "-3.364265989986018990374997730938263686743180501328627"},
                             {"2", "-6.822049349886392169085170632085091020323612667557795"},
                             {"5", "-16.69780025126922817405884705789034373251178086517828"}};
  model::detail::XXXTransitionRate<Real> kernel;
  for (auto const& point : points)
  {
    Real const rho = uni20::parse_real<Real>(point[0]), expected = -uni20::parse_real<Real>(point[1]);
    auto const k = kernel(rho), negative = kernel(-rho);
    ASSERT_TRUE(k.log_rate) << point[0];
    EXPECT_EQ(k.log_rate, negative.log_rate);
    EXPECT_REAL_NEAR(*k.log_rate, expected, Real{4096} * uni20::numeric_limits<Real>::epsilon());
    EXPECT_LE(std::abs(*k.log_rate - expected), *k.error);
  }
  model::ThermodynamicTwoSpinonStructureFactor<Real> sf;
  auto const point = sf(uni20::parse_real<Real>("1.5"), uni20::parse_real<Real>("1.8"));
  ASSERT_TRUE(point.value);
  Real const reference = uni20::parse_real<Real>("0.8731304551163434158090722608166793299789868199172442958");
  EXPECT_REAL_NEAR(*point.value, reference, Real{4096} * uni20::numeric_limits<Real>::epsilon());
  EXPECT_LE(std::abs(*point.value - reference), *point.error);
}

TYPED_TEST(XXXThermodynamicDSF, ThresholdAsymptotics)
{
  using Real = TypeParam;
  Real const pi = bethe::detail::pi<Real>();
  model::ThermodynamicTwoSpinonStructureFactor<Real> sf;
  auto const [lo, hi] = sf.boundaries(pi / Real{2});
  Real const small = uni20::parse_real<Real>("0.000001");
  auto const a = sf(pi / Real{2}, hi - small), b = sf(pi / Real{2}, hi - small / Real{4});
  ASSERT_TRUE(a.value);
  ASSERT_TRUE(b.value);
  // Vanishing square-root cusp at the upper edge, not the divergent DOS.
  EXPECT_REAL_NEAR(*a.value / *b.value, Real{2}, Real{0.0001});
  auto const c = sf(pi / Real{2}, lo + small), d = sf(pi / Real{2}, lo + small / Real{4});
  ASSERT_TRUE(c.value);
  ASSERT_TRUE(d.value);
  EXPECT_GT(*d.value, Real{2} * *c.value); // logarithmic enhancement
}

TYPED_TEST(XXXThermodynamicDSF, SupportSymmetryScalingAndFailures)
{
  using Real = TypeParam;
  Real const pi = bethe::detail::pi<Real>(), eps = uni20::numeric_limits<Real>::epsilon();
  model::ThermodynamicTwoSpinonStructureFactor<Real> sf, scaled(Real{3});
  EXPECT_EQ(*sf(Real{0}, Real{0}).value, Real{0});
  EXPECT_EQ(*sf(Real{2} * pi, Real{0}).value, Real{0});
  auto const [lo, hi] = sf.boundaries(pi / Real{2});
  EXPECT_EQ(sf(pi / Real{2}, lo).status, model::SpectralDensityStatus::lower_threshold);
  EXPECT_FALSE(sf(pi / Real{2}, lo).value);
  EXPECT_EQ(*sf(pi / Real{2}, hi).value, Real{0});
  EXPECT_EQ(*sf(pi / Real{2}, lo / Real{2}).value, Real{0});
  EXPECT_EQ(*sf(pi, Real{-1}).value, Real{0});
  EXPECT_EQ(sf(pi, Real{0}).status, model::SpectralDensityStatus::lower_threshold);
  Real const tiny_q = uni20::parse_real<Real>("1e-30");
  auto const collapsed = sf.boundaries(tiny_q);
  EXPECT_EQ(collapsed.first, collapsed.second);
  EXPECT_FALSE(sf(tiny_q, collapsed.first).converged());
  EXPECT_EQ(*sf(tiny_q, Real{-1}).value, Real{0});
  for (Real fraction : {Real{0.01}, Real{0.5}, Real{0.99}})
  {
    Real const w = lo + fraction * (hi - lo);
    auto const a = sf(pi / Real{2}, w), b = sf(Real{3} * pi / Real{2}, w), c = scaled(pi / Real{2}, Real{3} * w);
    ASSERT_TRUE(a.value);
    ASSERT_TRUE(b.value);
    ASSERT_TRUE(c.value);
    EXPECT_GT(*a.value, Real{0});
    EXPECT_REAL_NEAR(*a.value, *b.value, Real{100000} * eps * *a.value);
    EXPECT_REAL_NEAR(*a.value, Real{3} * *c.value, Real{100000} * eps * *a.value);
  }
  model::ThermodynamicTwoSpinonStructureFactor<Real> failed(Real{1}, {.max_evaluations = 0});
  EXPECT_FALSE(failed(pi, Real{1}).converged());
  EXPECT_FALSE(failed(pi, Real{1}).value);
  EXPECT_THROW(sf(Real{-1}, Real{1}), std::invalid_argument);
  EXPECT_THROW(sf(pi, uni20::numeric_limits<Real>::quiet_NaN()), std::invalid_argument);
  EXPECT_THROW((model::ThermodynamicTwoSpinonStructureFactor<Real>(Real{0})), std::invalid_argument);
  EXPECT_THROW((model::ThermodynamicTwoSpinonStructureFactor<Real>(Real{1}, {.tolerance = eps})),
               std::invalid_argument);
}

TEST(XXXThermodynamicSumRules, PublishedIntensityAndFirstMoment)
{
  // Analytically integrate q before doing the remaining rho integral.
  // K1(q)=sqrt(U^2-L^2)/4 * integral exp(-I)*sech(pi*rho)^2 d rho.
  // Four times the integrated zz intensity is integral
  // 2*rho*exp(-I)*sech(pi*rho)^2/tanh(pi*rho) d rho.
  model::detail::XXXTransitionRate<double> kernel;
  auto const rule = bethe::detail::gauss_legendre<double>(32);
  double first = 0, total = 0, actual_first = 0, pi = bethe::detail::pi<double>();
  model::ThermodynamicTwoSpinonStructureFactor<double> sf;
  for (int p = 0; p < 32; ++p)
    for (std::size_t j = 0; j < rule.x.size(); ++j)
    {
      double rho = .5 * (p + (1 + rule.x[j]) / 2);
      auto const k = kernel(rho);
      ASSERT_TRUE(k.log_rate) << rho;
      double f = std::exp(*k.log_rate) / std::pow(std::cosh(pi * rho), 2);
      first += .25 * rule.w[j] * f;
      total += .25 * rule.w[j] * 2 * rho * f / std::tanh(pi * rho);
      // Independently exercise the public density and its energy Jacobian,
      // rather than only the integrated-kernel normalization.
      double w = pi / std::cosh(pi * rho);
      auto const density = sf(pi, w);
      ASSERT_TRUE(density.value);
      actual_first += .25 * rule.w[j] * w * w * std::tanh(pi * rho) * *density.value / 2;
    }
  EXPECT_NEAR(total, .7289, .00005);
  EXPECT_NEAR(first * 3 * pi / (16 * (std::log(2.) - .25)), .7130, .00005);
  EXPECT_NEAR(actual_first, pi * first / 4, 1e-12);
}
