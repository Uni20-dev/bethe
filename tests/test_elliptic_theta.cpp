// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/detail/elliptic_theta.hpp>
#include <bethe/xyz.hpp>

namespace
{
template <typename Real> class EllipticTheta : public ::testing::Test {};
TYPED_TEST_SUITE(EllipticTheta, test_support::RealTypes, test_support::PrecisionNames);
template <typename Real> auto actual(bethe::detail::ThetaJet<Real> const& jet, unsigned derivative = 0)
{
  return std::exp(jet.log_scale) * jet.derivative[derivative];
}

TYPED_TEST(EllipticTheta, IndependentHighPrecisionValues)
{
  using Real = TypeParam;
  using Complex = std::complex<Real>;
  // Independently generated with mpmath 90-digit jtheta/diff, decimal inputs.
  std::array<std::array<char const*, 2>, 4> const reference{
      {{"0.69219558972113418875657598555125487595689802745775963095538321277550698026515883",
        "0.30295842568440231704359734352559357616656419512939285695290802635041584770007338"},
       {"0.97370036030675749908185309462479722684202005621228971607055473603044726260643021",
        "-0.23134011872194558048036441888787309500055703385516093235047558690923750031527775"},
       {"1.0820563314307416382141902966895917672243922526193374482791519769533824502881634",
        "-0.14172351980229852150061527154356932129353427560991712023159400835264615038631135"},
       {"0.91701403643718451700049740426185963471839253409158838269591152285563434881747893",
        "0.14114932439016119087397800722902770968704101496998396921409793366468937888358458"}}};
  Real const eps = uni20::numeric_limits<Real>::epsilon(), t = uni20::parse_real<Real>("0.7");
  Complex const u(uni20::parse_real<Real>("0.2"), uni20::parse_real<Real>("0.1"));
  for (unsigned k = 1; k <= 4; ++k)
  {
    auto const jet = bethe::detail::elliptic_theta(k, u, t);
    Complex const expected(uni20::parse_real<Real>(reference[k - 1][0]), uni20::parse_real<Real>(reference[k - 1][1]));
    EXPECT_REAL_NEAR(std::abs(actual(jet) - expected), Real{0}, Real{512} * eps);
  }
  auto const jet = bethe::detail::elliptic_theta(1, u, t);
  Complex const first(
      uni20::parse_real<Real>("3.140373175101863767109126990443280790919216481608393060449146373266355474"),
      uni20::parse_real<Real>("-0.542140463819793229105053265208634636705604224096074518836382433938263608"));
  Complex const second(
      uni20::parse_real<Real>("-5.256238827195634225526959228714510493158642457037190918060930185253611353"),
      uni20::parse_real<Real>("-3.365849286153387228448025426576455363601810484091901718883225391186697882"));
  EXPECT_REAL_NEAR(std::abs(actual(jet, 1) - first), Real{0}, Real{1024} * eps);
  EXPECT_REAL_NEAR(std::abs(actual(jet, 2) - second), Real{0}, Real{4096} * eps);
}

TYPED_TEST(EllipticTheta, DualSeriesPeriodsAndDerivatives)
{
  using Real = TypeParam;
  using Complex = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  for (Real t : {Real{0.5}, Real{1}, Real{2}})
  {
    Complex const u(Real{0.125}, t / Real{8});
    auto const fourier = bethe::detail::theta_detail::theta1_reduced(u, t, false);
    auto const gaussian = bethe::detail::theta_detail::theta1_reduced(u, t, true);
    for (unsigned d = 0; d < 3; ++d)
      EXPECT_REAL_NEAR(std::abs(actual(fourier, d) - actual(gaussian, d)), Real{0},
                       Real{4096} * eps * (Real{1} + std::abs(actual(fourier, d))));
    for (unsigned k = 1; k <= 4; ++k)
    {
      auto const base = bethe::detail::elliptic_theta(k, u, t);
      for (int n : {-2, -1, 0, 1, 2})
        for (int m : {-2, -1, 0, 1, 2})
        {
          auto const shifted = bethe::detail::elliptic_theta(k, u + Complex(Real(n), Real(m) * t), t);
          int const sign = ((k == 1 ? n + m : k == 2 ? n : k == 4 ? m : 0) % 2) ? -1 : 1;
          Complex const slope(0, -Real{2} * pi * Real(m)),
              phase = Real(sign) * std::exp(Complex(0, -Real{2} * pi * Real(m) * u.real()));
          Real const exponent = pi * Real(m) * (Real(m) * t + Real{2} * u.imag());
          for (unsigned d = 0; d < 3; ++d)
          {
            auto expected = actual(base, d);
            if (d == 1) expected += slope * actual(base);
            if (d == 2) expected += Real{2} * slope * actual(base, 1) + slope * slope * actual(base);
            auto const normalized = shifted.derivative[d] * std::exp(shifted.log_scale - exponent);
            EXPECT_REAL_NEAR(std::abs(normalized - phase * expected), Real{0},
                             Real{8192} * eps * (Real{1} + std::abs(expected)));
          }
        }
      Real const h = std::cbrt(eps);
      auto const plus = bethe::detail::elliptic_theta(k, u + h, t), minus = bethe::detail::elliptic_theta(k, u - h, t);
      EXPECT_REAL_NEAR(std::abs((actual(plus) - actual(minus)) / (Real{2} * h) - actual(base, 1)), Real{0},
                       Real{10000} * h * h);
      EXPECT_REAL_NEAR(std::abs((actual(plus, 1) - actual(minus, 1)) / (Real{2} * h) - actual(base, 2)), Real{0},
                       Real{100000} * h * h);
    }
  }
}

TYPED_TEST(EllipticTheta, ZerosTinyArgumentsAndScaledLimits)
{
  using Real = TypeParam;
  using Complex = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  for (Real t : {Real{0.001}, Real{0.1}, Real{1}, Real{10}, Real{10000}})
  {
    auto const zero = bethe::detail::elliptic_theta(1, Complex{}, t);
    EXPECT_EQ(std::abs(zero.derivative[0]), Real{0});
    EXPECT_EQ(std::abs(zero.derivative[2]), Real{0});
    auto const tiny = bethe::detail::elliptic_theta(1, Complex(eps * eps, 0), t);
    EXPECT_REAL_NEAR(std::abs(tiny.derivative[0] / (eps * eps * tiny.derivative[1]) - Real{1}), Real{0},
                     Real{2048} * eps);
    for (unsigned k = 2; k <= 4; ++k)
    {
      Complex const u = k == 2   ? Complex(Real{0.5}, 0)
                        : k == 3 ? Complex(Real{0.5}, t / Real{2})
                                 : Complex(0, t / Real{2});
      auto const root = bethe::detail::elliptic_theta(k, u, t);
      EXPECT_EQ(std::abs(root.derivative[0]), Real{0});
    }
  }
  auto const trigonometric = bethe::detail::elliptic_theta(1, Complex(Real{0.25}, 0), Real{10000});
  for (Real t : {Real{0.01}, Real{0.5}, Real{1}, Real{10}})
    for (unsigned kind = 2; kind <= 4; ++kind)
    {
      auto const zero = bethe::detail::elliptic_theta(kind, Complex{}, t);
      auto const tiny = bethe::detail::elliptic_theta(kind, Complex(eps * eps, 0), t);
      EXPECT_EQ(std::abs(zero.derivative[1]), Real{0});
      auto const normalized =
          tiny.derivative[1] * std::exp(tiny.log_scale - zero.log_scale) / (eps * eps * zero.derivative[2]);
      EXPECT_REAL_NEAR(std::abs(normalized - Real{1}), Real{0}, Real{2048} * eps);
    }
  EXPECT_REAL_NEAR(trigonometric.log_scale, -pi * Real{2500}, Real{64} * eps * pi * Real{2500});
  EXPECT_REAL_NEAR(std::abs(trigonometric.derivative[1] / trigonometric.derivative[0] - pi), Real{0},
                   Real{256} * eps * pi);
  auto const narrow = bethe::detail::elliptic_theta(1, Complex(Real{0.25}, 0), Real{1} / Real{1000000});
  EXPECT_REAL_NEAR(std::abs(narrow.derivative[1] / narrow.derivative[0] - pi * Real{500000}), Real{0},
                   Real{256} * eps * pi * Real{500000});
  EXPECT_THROW(bethe::detail::elliptic_theta(0, Complex{}, Real{1}), std::invalid_argument);
  EXPECT_THROW(bethe::detail::elliptic_theta(1, Complex{}, Real{0}), std::invalid_argument);
  EXPECT_THROW(bethe::detail::elliptic_theta(1, Complex(uni20::numeric_limits<Real>::infinity(), 0), Real{1}),
               std::invalid_argument);
  EXPECT_THROW(bethe::detail::elliptic_theta(1, Complex(Real{1} / eps, 0), Real{1}), std::overflow_error);
  auto const zero = bethe::detail::elliptic_theta(1, Complex{}, Real{1});
  EXPECT_THROW(bethe::detail::theta_ratio(zero, 1, zero, 0), std::domain_error);
  EXPECT_THROW(bethe::detail::theta_ratio(zero, 3, zero, 1), std::invalid_argument);
}

TYPED_TEST(EllipticTheta, XYZCouplingNormalization)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  auto const values = bethe::xyz::couplings(Real{0.25}, uni20::parse_real<Real>("0.7"));
  EXPECT_REAL_NEAR(
      values.x, uni20::parse_real<Real>("1.28413294604596728052548431355937740558049584455828014366415979541534619403"),
      Real{2048} * eps);
  EXPECT_REAL_NEAR(
      values.y,
      uni20::parse_real<Real>("0.818012679007936827203729999038278734809654235550546858902930502983096986726"),
      Real{2048} * eps);
  EXPECT_REAL_NEAR(
      values.z, uni20::parse_real<Real>("0.68992199448261857320856457867294294680755351238831796860570067495746513562"),
      Real{2048} * eps);
  for (Real t : {Real{0.01}, Real{0.5}, Real{1}, Real{10}, Real{10000}})
  {
    auto const xxx = bethe::xyz::couplings(Real{0}, t);
    EXPECT_EQ(xxx.x, Real{1});
    EXPECT_EQ(xxx.y, Real{1});
    EXPECT_EQ(xxx.z, Real{1});
    auto const xy = bethe::xyz::couplings(Real{0.5}, t);
    EXPECT_EQ(xy.z, Real{0});
    EXPECT_REAL_NEAR(xy.x * xy.y, Real{1}, Real{4096} * eps);
  }
  auto const xxz = bethe::xyz::couplings(Real{0.25}, Real{10000});
  EXPECT_REAL_NEAR(xxz.x, Real{1}, Real{4096} * eps);
  EXPECT_REAL_NEAR(xxz.y, Real{1}, Real{4096} * eps);
  EXPECT_REAL_NEAR(xxz.z, std::cos(pi / Real{4}), Real{4096} * eps);
  EXPECT_THROW(bethe::xyz::couplings(Real{0.25}, Real{0}), std::invalid_argument);
}
TYPED_TEST(EllipticTheta, XYZTwoSiteEnergyNormalization)
{
  using Real = TypeParam;
  using Complex = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real t : {Real{0.1}, Real{0.38}, Real{1}, Real{5}})
    for (Real eta : {Real{0.2}, Real{0.4}, Real{0.7}})
    {
      auto const j = bethe::xyz::couplings(eta, t);
      // Bell eigenstates of 2*(Jx*Sx*Sx+Jy*Sy*Sy+Jz*Sz*Sz).
      std::array<Real, 4> expected{-(j.x + j.y + j.z) / Real{2}, (-j.x + j.y + j.z) / Real{2},
                                   (j.x - j.y + j.z) / Real{2}, (j.x + j.y - j.z) / Real{2}};
      std::array<Complex, 4> root{Complex{}, Complex(0, t / Real{2}), Complex(Real{0.5}, t / Real{2}),
                                  Complex(Real{0.5}, 0)};
      for (unsigned k = 0; k < 4; ++k)
      {
        auto const e = bethe::xyz::candidate_energy<Real>(2, std::span<Complex const>(&root[k], 1), eta, t);
        EXPECT_REAL_NEAR(e.real(), expected[k], Real{8192} * eps * (Real{1} + std::abs(expected[k])));
        EXPECT_REAL_NEAR(e.imag(), Real{0}, Real{8192} * eps * (Real{1} + std::abs(expected[k])));
      }
    }
  std::array<Complex, 1> pole{Complex(Real{0.25}, 0)};
  EXPECT_THROW(bethe::xyz::candidate_energy<Real>(2, pole, Real{0.5}, Real{1}), std::domain_error);
  EXPECT_THROW(bethe::xyz::candidate_energy<Real>(3, pole, Real{0.5}, Real{1}), std::invalid_argument);
}
} // namespace
