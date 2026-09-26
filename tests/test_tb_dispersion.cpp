// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/takhtajan_babujian.hpp>
#include <bethe/tb_dispersion.hpp>

namespace
{
namespace model = bethe::takhtajan_babujian;
template <typename Real> class TBDispersion : public ::testing::Test {};
TYPED_TEST_SUITE(TBDispersion, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(TBDispersion, NormalizationNativePrecisionAndDressedRapidity)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  model::SpinonDispersion<Real> const band;
  EXPECT_EQ(band.spin(), uni20::from_twice(1));
  EXPECT_EQ(band.gap(), Real{0});
  EXPECT_EQ(band.energy(Real{0}), Real{0});
  EXPECT_EQ(band.energy(pi), Real{0});
  EXPECT_EQ(band.maximum_energy(), Real{2} * pi);
  EXPECT_REAL_NEAR(
      band.energy(uni20::parse_real<Real>("0.3")),
      uni20::parse_real<Real>("1.856808220469203776013916923017469578630411818675737307724914301655265856822454147584"),
      Real{16} * eps);
  // Independent dressed two-string hole energy, obtained by dividing
  // -4*pi*(a1+a3) by 1+2*a2+a4 in Fourier space. Thus e=4*pi*rho,
  // rho=1/(2*cosh(pi*lambda)); p integrates 2*pi*rho from lambda to infinity.
  for (Real lambda : {-Real{2}, -Real{0.5}, Real{0}, Real{0.5}, Real{2}})
  {
    Real const p = Real{2} * std::atan(std::exp(-pi * lambda));
    EXPECT_REAL_NEAR(band.energy(p), Real{2} * pi / std::cosh(pi * lambda), Real{32} * eps);
  }
  // Shared XXX band kinematics, but TB's local J=1 is four times XXX J=1.
  for (unsigned j = 0; j <= 32; ++j)
  {
    Real const p = pi * Real(j) / Real{32};
    EXPECT_REAL_NEAR(band.energy(p), Real{4} * bethe::heisenberg::spinon_energy(p), Real{32} * eps);
  }
  Real const tiny = eps * eps;
  EXPECT_REAL_NEAR(band.energy(tiny) / tiny, band.velocity(), Real{8} * eps);
  Real const subnormal = uni20::numeric_limits<Real>::denorm_min();
  EXPECT_EQ(band.energy(subnormal), band.velocity() * subnormal);
  auto const edge = band.continuum(subnormal);
  EXPECT_EQ(edge.lower, band.energy(subnormal));
  EXPECT_EQ(edge.upper, band.energy(subnormal));
  EXPECT_EQ(band.continuum(subnormal, bethe::SpinonMomentum::folded).lower, edge.lower);
  EXPECT_EQ(band.four_spinon_continuum(subnormal, bethe::SpinonMomentum::folded).lower, edge.lower);
}

TYPED_TEST(TBDispersion, TwoSpinonDirectScan)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  model::SpinonDispersion<Real> const band;
  for (unsigned j = 0; j <= 32; ++j)
  {
    Real const q = Real{2} * pi * Real(j) / Real{32};
    Real const lo = std::max(Real{0}, q - pi), hi = std::min(q, pi);
    Real low = Real{100}, high = Real{0};
    for (unsigned i = 0; i <= 1024; ++i)
    {
      Real const p = lo + (hi - lo) * Real(i) / Real{1024};
      Real const value = band.energy(std::clamp(p, Real{0}, pi)) + band.energy(std::clamp(q - p, Real{0}, pi));
      low = std::min(low, value);
      high = std::max(high, value);
    }
    auto const edge = band.continuum(q);
    EXPECT_REAL_NEAR(edge.lower, low, Real{128} * eps);
    EXPECT_REAL_NEAR(edge.upper, high, Real{1} / Real{10000});
  }
}

TYPED_TEST(TBDispersion, FourSpinonConvolutionAndFolding)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  model::SpinonDispersion<Real> const band;
  for (unsigned j = 0; j <= 32; ++j)
  {
    Real const q = Real{2} * pi * Real(j) / Real{32};
    Real low = Real{100}, high = Real{0};
    for (unsigned i = 0; i <= 1024; ++i)
    {
      Real const p = Real{2} * pi * Real(i) / Real{1024};
      Real r = q - p;
      if (r < Real{0}) r += Real{2} * pi;
      auto const a = band.continuum(p), b = band.continuum(r);
      low = std::min(low, a.lower + b.lower);
      high = std::max(high, a.upper + b.upper);
    }
    auto const edge = band.four_spinon_continuum(q);
    EXPECT_REAL_NEAR(edge.lower, low, Real{256} * eps);
    EXPECT_REAL_NEAR(edge.upper, high, Real{1} / Real{1000});
    Real const shifted = q < pi ? q + pi : q - pi;
    for (bool four : {false, true})
    {
      auto const a = four ? band.four_spinon_continuum(q) : band.continuum(q);
      auto const b = four ? band.four_spinon_continuum(shifted) : band.continuum(shifted);
      auto const folded = four ? band.four_spinon_continuum(q, bethe::SpinonMomentum::folded)
                               : band.continuum(q, bethe::SpinonMomentum::folded);
      EXPECT_REAL_NEAR(folded.lower, std::min(a.lower, b.lower), Real{128} * eps);
      EXPECT_REAL_NEAR(folded.upper, std::max(a.upper, b.upper), Real{128} * eps);
    }
  }
  EXPECT_EQ(band.continuum(Real{0}).upper, Real{0});
  EXPECT_EQ(band.four_spinon_continuum(Real{0}).upper, Real{4} * band.maximum_energy());
  EXPECT_REAL_NEAR(band.four_spinon_continuum(pi).upper, Real{4} * pi * std::sqrt(Real{2}), Real{64} * eps);
}

TYPED_TEST(TBDispersion, ScalingAndNumericalFailures)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  model::SpinonDispersion<Real> const band, scaled(Real{3});
  EXPECT_REAL_NEAR(scaled.energy(Real{1}), Real{3} * band.energy(Real{1}), Real{64} * eps);
  EXPECT_REAL_NEAR(scaled.four_spinon_continuum(pi).upper, Real{3} * band.four_spinon_continuum(pi).upper,
                   Real{256} * eps);
  EXPECT_THROW(model::SpinonDispersion<Real>(Real{0}), std::invalid_argument);
  EXPECT_THROW(model::SpinonDispersion<Real>(-Real{1}), std::invalid_argument);
  EXPECT_THROW(model::SpinonDispersion<Real>{uni20::numeric_limits<Real>::infinity()}, std::invalid_argument);
  EXPECT_THROW(model::SpinonDispersion<Real>{uni20::numeric_limits<Real>::max()}, std::overflow_error);
  EXPECT_THROW(band.energy(-eps), std::invalid_argument);
  EXPECT_THROW(band.energy(Real{4}), std::invalid_argument);
  EXPECT_THROW(band.continuum(uni20::numeric_limits<Real>::quiet_NaN()), std::invalid_argument);
  EXPECT_THROW(band.four_spinon_continuum(Real{7}), std::invalid_argument);
  EXPECT_THROW(band.four_spinon_continuum(Real{0}, static_cast<bethe::SpinonMomentum>(7)), std::invalid_argument);
  model::SpinonDispersion<Real> const huge(uni20::numeric_limits<Real>::max() / Real{16});
  EXPECT_THROW(huge.four_spinon_continuum(Real{0}), std::overflow_error);
  model::SpinonDispersion<Real> const tiny(uni20::numeric_limits<Real>::min());
  EXPECT_THROW(tiny.energy(uni20::numeric_limits<Real>::min()), std::underflow_error);
}

TEST(TBDispersionOracles, FiniteRingSpinAndMomentumChecks)
{
  // scripts/reference_tb_excitations.py: independent spin-1 matrices and
  // exact translation/S^2 projection. The finite-size corrections are not
  // used as an error tolerance for the native thermodynamic formulas.
  struct Fixture
  {
      unsigned length;
      double ground, small_q_gap, staggered_gap, spin_two_gap;
  };
  Fixture const data[] = {{4, -17.403124237432877, 8.167056259933087, 3.4031242374328805, 6.5746971126866836},
                          {6, -24.877409871293818, 6.207722752576089, 2.2257393676516166, 4.686482195819234},
                          {8, -32.64273973285925, 4.852117819005649, 1.6619603613968792, 3.6541907480856963},
                          {10, -40.50829245662726, 3.94871104155294, 1.328463417081032, 2.9980204068189664}};
  double const pi = bethe::detail::pi<double>();
  model::SpinonDispersion<double> const band;
  double last_error = 2, last_gap = 10;
  for (auto const& f : data)
  {
    auto const state = model::ground_state<double>(f.length);
    ASSERT_TRUE(state.converged);
    EXPECT_NEAR(state.energy, f.ground, 2e-10); // actual roots retain string deviations
    double const q = 2 * pi / f.length, v = f.small_q_gap / q;
    EXPECT_LT(std::abs(v - band.velocity()), last_error);
    last_error = std::abs(v - band.velocity());
    EXPECT_NEAR(f.length * f.staggered_gap / (2 * pi * band.velocity()), 3.0 / 8.0, 0.045);
    EXPECT_GT(f.staggered_gap, band.continuum(pi).lower);
    EXPECT_GT(f.spin_two_gap, band.four_spinon_continuum(0.0).lower);
    EXPECT_LT(f.staggered_gap, last_gap);
    last_gap = f.staggered_gap;
    EXPECT_NEAR(f.small_q_gap, band.continuum(q).upper, 0.18 * band.continuum(q).upper);
  }
}
} // namespace
