// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/su3_dispersion.hpp>

namespace
{
using bethe::su3::Continuum;
using bethe::su3::Momentum;
using bethe::su3::Particle;
template <typename Real> class SU3Dispersion : public ::testing::Test {};
TYPED_TEST_SUITE(SU3Dispersion, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(SU3Dispersion, NativePrecisionBandsAndSoftEndpoints)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  bethe::su3::ExcitationDispersion<Real> const band;
  char const* reference[] = {"0.672943171908697176654766974738581598626465808498253705740688258859362234071373801",
                             "0.5649289750707720073545109739397314537938087372855711660759212755774816704769289641"};
  unsigned index = 0;
  for (auto particle : {Particle::fundamental, Particle::antifundamental})
  {
    Real const end = band.momentum_max(particle);
    EXPECT_REAL_NEAR(band.energy(particle, uni20::parse_real<Real>("0.3")), uni20::parse_real<Real>(reference[index++]),
                     Real{16} * eps);
    EXPECT_EQ(band.energy(particle, Real{0}), Real{0});
    EXPECT_EQ(band.energy(particle, end), Real{0});
    EXPECT_REAL_NEAR(band.energy(particle, end / Real{2}), band.maximum_energy(particle), Real{16} * eps);
    // A cosine difference would round this positive energy to zero.
    Real const tiny = eps * eps;
    EXPECT_REAL_NEAR(band.energy(particle, tiny) / tiny, band.velocity(), Real{16} * eps);
    EXPECT_REAL_NEAR(band.energy(particle, end * eps) / (end * eps), band.velocity(), Real{64} * eps);
    for (unsigned j = 0; j <= 32; ++j)
    {
      Real const p = end * Real(j) / Real{32};
      EXPECT_REAL_NEAR(band.energy(particle, p), band.energy(particle, end - p), Real{16} * eps);
    }
  }
  EXPECT_REAL_NEAR(band.velocity(), Real{2} * pi / Real{3}, Real{8} * eps);
  EXPECT_EQ(band.momentum_max(Particle::fundamental), Real{4} * pi / Real{3});
  EXPECT_EQ(band.momentum_max(Particle::antifundamental), Real{2} * pi / Real{3});
  Real const subnormal = uni20::numeric_limits<Real>::denorm_min();
  EXPECT_EQ(band.energy(Particle::fundamental, subnormal), band.velocity() * subnormal);
  EXPECT_GT(band.energy(Particle::fundamental, subnormal), Real{0});
}

TYPED_TEST(SU3Dispersion, TwoParticleBoundsAgainstDirectMomentumScan)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  bethe::su3::ExcitationDispersion<Real> const band;
  Real const narrow = band.momentum_max(Particle::antifundamental), wide = band.momentum_max(Particle::fundamental);
  for (unsigned j = 0; j <= 96; ++j)
  {
    Real const q = Real{2} * pi * Real(j) / Real{96};
    Real const lo = std::max(Real{0}, q - wide), hi = std::min(narrow, q);
    Real low = Real{100}, high = Real{0};
    for (unsigned i = 0; i <= 1024; ++i)
    {
      Real const p = lo + (hi - lo) * Real(i) / Real{1024};
      Real const value = band.energy(Particle::antifundamental, std::clamp(p, Real{0}, narrow)) +
                         band.energy(Particle::fundamental, std::clamp(q - p, Real{0}, wide));
      low = std::min(low, value);
      high = std::max(high, value);
    }
    auto const edges = band.continuum(q);
    EXPECT_LE(edges.lower, low + Real{64} * eps);
    EXPECT_GE(edges.upper, high - Real{64} * eps);
    EXPECT_REAL_NEAR(edges.lower, low, Real{1} / Real{1024});
    EXPECT_REAL_NEAR(edges.upper, high, Real{1} / Real{1024});
  }
}

TYPED_TEST(SU3Dispersion, IndependentNestedSeaDensities)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  bethe::su3::ExcitationDispersion<Real> const band;
  // Inverse Fourier transform of sinh((3-j)w/2)/sinh(3w/2),
  // Doikou--Nepomechie Eq. (2.46), and integrated hole momentum (2.50).
  // These expressions do not use the trigonometric dispersion in momentum.
  for (unsigned sea : {1u, 2u})
    for (Real lambda : {-Real{2}, -Real{0.25}, Real{0}, Real{0.25}, Real{2}})
    {
      Real const angle = Real(sea) * pi / Real{3};
      Real const density = std::sin(angle) / (Real{3} * (std::cosh(Real{2} * pi * lambda / Real{3}) - std::cos(angle)));
      Real const momentum =
          pi - angle - Real{2} * std::atan(std::tanh(pi * lambda / Real{3}) / std::tan(angle / Real{2}));
      auto const particle = sea == 1 ? Particle::fundamental : Particle::antifundamental;
      EXPECT_REAL_NEAR(band.energy(particle, momentum), Real{2} * pi * density, Real{64} * eps);
    }
}

TYPED_TEST(SU3Dispersion, FourParticleBoundsAndThreeSiteFolding)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  bethe::su3::ExcitationDispersion<Real> const band;
  for (unsigned j = 0; j <= 48; ++j)
  {
    Real const q = Real{2} * pi * Real(j) / Real{48};
    // An independent convolution of the already scan-validated two-particle
    // intervals. Include both total-momentum aliases, not just p1+p2=Q.
    Real low = Real{100}, high = Real{0};
    for (unsigned i = 0; i <= 1536; ++i)
    {
      Real const p = Real{2} * pi * Real(i) / Real{1536};
      Real r = q - p;
      if (r < Real{0}) r += Real{2} * pi;
      auto const a = band.continuum(p), b = band.continuum(r);
      low = std::min(low, a.lower + b.lower);
      high = std::max(high, a.upper + b.upper);
    }
    auto const edges = band.continuum(q, Continuum::four_soliton);
    EXPECT_REAL_NEAR(edges.lower, low, Real{1} / Real{256});
    EXPECT_REAL_NEAR(edges.upper, high, Real{1} / Real{256});
    EXPECT_LE(edges.lower, low + Real{128} * eps);
    EXPECT_GE(edges.upper, high - Real{128} * eps);
    for (auto content : {Continuum::two_soliton, Continuum::four_soliton})
    {
      Real lower = Real{100}, upper = Real{0};
      for (unsigned image = 0; image < 3; ++image)
      {
        Real p = q + Real(2 * image) * pi / Real{3};
        if (p > Real{2} * pi) p -= Real{2} * pi;
        auto const e = band.continuum(p, content);
        lower = std::min(lower, e.lower);
        upper = std::max(upper, e.upper);
      }
      auto const folded = band.continuum(q, content, Momentum::three_site_folded);
      EXPECT_REAL_NEAR(folded.lower, lower, Real{64} * eps);
      EXPECT_REAL_NEAR(folded.upper, upper, Real{64} * eps);
    }
  }
  // Crucial iMPS case: the two-particle lower edge is NOT the lowest
  // adjoint-channel kinematic threshold between the two soft momenta.
  auto const two = band.continuum(pi), four = band.continuum(pi, Continuum::four_soliton);
  EXPECT_REAL_NEAR(four.lower * Real{2}, two.lower, Real{16} * eps);
  EXPECT_GT(two.lower, four.lower);
  EXPECT_REAL_NEAR(four.upper, Real{8} * pi * std::sqrt(Real{2}) / (Real{3} * std::sqrt(Real{3})), Real{64} * eps);
}

TYPED_TEST(SU3Dispersion, ExchangeScalingAndInvalidInputs)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  bethe::su3::ExcitationDispersion<Real> const band, scaled(Real{7});
  for (auto particle : {Particle::fundamental, Particle::antifundamental})
    EXPECT_REAL_NEAR(scaled.energy(particle, Real{1}), Real{7} * band.energy(particle, Real{1}), Real{128} * eps);
  for (auto content : {Continuum::two_soliton, Continuum::four_soliton})
  {
    auto const a = band.continuum(pi, content), b = scaled.continuum(pi, content);
    EXPECT_REAL_NEAR(b.lower, Real{7} * a.lower, Real{256} * eps);
    EXPECT_REAL_NEAR(b.upper, Real{7} * a.upper, Real{256} * eps);
  }
  EXPECT_THROW(bethe::su3::ExcitationDispersion<Real>(Real{0}), std::invalid_argument);
  EXPECT_THROW(bethe::su3::ExcitationDispersion<Real>(-Real{1}), std::invalid_argument);
  EXPECT_THROW(bethe::su3::ExcitationDispersion<Real>{uni20::numeric_limits<Real>::infinity()}, std::invalid_argument);
  EXPECT_THROW(bethe::su3::ExcitationDispersion<Real>{uni20::numeric_limits<Real>::max()}, std::overflow_error);
  EXPECT_THROW(band.energy(Particle::antifundamental, pi), std::invalid_argument);
  EXPECT_THROW(band.energy(Particle::fundamental, -eps), std::invalid_argument);
  EXPECT_THROW(band.energy(static_cast<Particle>(7), Real{0}), std::invalid_argument);
  EXPECT_THROW(band.continuum(uni20::numeric_limits<Real>::quiet_NaN()), std::invalid_argument);
  EXPECT_THROW(band.continuum(Real{7}), std::invalid_argument);
  EXPECT_THROW(band.continuum(Real{0}, static_cast<Continuum>(7)), std::invalid_argument);
  EXPECT_THROW(band.continuum(Real{0}, Continuum::two_soliton, static_cast<Momentum>(7)), std::invalid_argument);
  bethe::su3::ExcitationDispersion<Real> const huge(uni20::numeric_limits<Real>::max() / Real{8});
  EXPECT_THROW(huge.continuum(Real{0}, Continuum::four_soliton), std::overflow_error);
  bethe::su3::ExcitationDispersion<Real> const tiny(uni20::numeric_limits<Real>::min());
  EXPECT_THROW(tiny.energy(Particle::fundamental, uni20::numeric_limits<Real>::min()), std::underflow_error);
}

TEST(SU3DispersionOracles, FiniteRingAdjointScalingAndFourParticleCounterexample)
{
  // Independent color-word, momentum-projected ED, including C2=3 selection.
  // Reproduce with scripts/reference_su3_excitations.py. These are fp64 ED
  // references, not exact finite-size BA results or high-precision fixtures.
  double const pi = bethe::detail::pi<double>();
  bethe::su3::ExcitationDispersion<double> const band;
  struct Fixture
  {
      unsigned length;
      double small_q_gap, soft_gap;
  };
  Fixture const data[] = {{6, 2.0367016659324833, 1.4011497251345526},
                          {9, 1.430375396369497, 0.9322358471981751},
                          {12, 1.0907287552890352, 0.7003049619050863}};
  double previous_error = 1;
  for (auto const& f : data)
  {
    double const velocity = f.small_q_gap * f.length / (2 * pi);
    double const error = std::abs(velocity - band.velocity());
    EXPECT_LT(error, previous_error);
    EXPECT_LT(error, 0.08 * band.velocity());
    previous_error = error;
    // SU(3)_1 adjoint tower x=2/3, with logarithmic finite-size corrections.
    double const dimension = f.length * f.soft_gap / (2 * pi * band.velocity());
    EXPECT_NEAR(dimension, 2.0 / 3.0, 0.04);
    EXPECT_GT(f.soft_gap, 0); // unlike the thermodynamic soft threshold
  }
  EXPECT_EQ(band.continuum(2 * pi / 3).lower, 0);
  // Already at L=12 the local-response (adjoint) sector has a level below
  // the two-soliton envelope at pi. It approaches the lower multiparticle
  // threshold, not a falsely identified two-particle spectral onset.
  double const adjoint_gap_at_pi = 2.0188788699927294;
  EXPECT_LT(adjoint_gap_at_pi, band.continuum(pi).lower);
  EXPECT_GT(adjoint_gap_at_pi, band.continuum(pi, Continuum::four_soliton).lower);
}
} // namespace
