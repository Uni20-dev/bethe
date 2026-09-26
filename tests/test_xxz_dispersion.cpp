// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxz_bulk.hpp>
#include <bethe/xxz_dispersion.hpp>

namespace
{
template <typename Real> class XXZDispersion : public ::testing::Test {};
TYPED_TEST_SUITE(XXZDispersion, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXZDispersion, IndependentThetaReferences)
{
  using Real = TypeParam;
  // mpmath 90-digit jtheta, exact binary-representable anisotropies. This
  // includes a mass below epsilon: computing sqrt(1-k*k) would give zero.
  struct Fixture
  {
      char const *delta, *maximum, *gap, *energy;
  };
  Fixture const fixtures[] = {{"1.0078125", "1.57488481565395445394845819362865977565068654844206579053926431182299",
                               "4.394403365495107781556670501508773572180778665838746595180932219401351e-17",
                               "0.46541028618986229976987350755986308590070717986873884313497766512903"},
                              {"1.125", "1.63571640876064637419378889905167016073183235090174640576579357096031",
                               "0.000305893591956058457594898867913667451120769202880443511451506459191114",
                               "0.48338733949036485155879806515678519020849380658165864154997203191258"},
                              {"2", "2.07049615480713443465916002155733565857530307500104608496194981910538",
                               "0.19490113571001470727701163515916987988002921827177622392600310333179",
                               "0.63957652647962607833811009496693871208889567824798565871139993315029"},
                              {"4", "3.04785480879093419508699977578530438582676121563206802945619617491017",
                               "1.07909680831248795103812190041789171073708239068431083876640207198874",
                               "1.36894896911584580564823854825782568984447313217071688608105066829942"},
                              {"100", "51.0024750625000155089746771518275259570445846487037805310423228804276",
                               "49.0025250624999842566307904262413198986957841196193139452782692587505",
                               "49.1804263673154190878316270259864683636413280312440806444403880489983"}};
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  for (auto const& f : fixtures)
  {
    SCOPED_TRACE(f.delta);
    bethe::xxz::SpinonDispersion<Real> const band(uni20::parse_real<Real>(f.delta));
    auto const expected_gap = uni20::parse_real<Real>(f.gap);
    EXPECT_REAL_NEAR(band.gap() / expected_gap, Real{1}, Real{256} * eps);
    EXPECT_REAL_NEAR(band.maximum_energy() / uni20::parse_real<Real>(f.maximum), Real{1}, Real{64} * eps);
    EXPECT_REAL_NEAR(band.energy(uni20::parse_real<Real>("0.3")) / uni20::parse_real<Real>(f.energy), Real{1},
                     Real{64} * eps);
    EXPECT_EQ(band.energy(Real{0}), band.gap());
    EXPECT_EQ(band.energy(pi), band.gap());
    EXPECT_EQ(band.energy(pi / Real{2}), band.maximum_energy());
    // Scaling J must preserve all energies, including a tiny mass.
    bethe::xxz::SpinonDispersion<Real> const scaled(uni20::parse_real<Real>(f.delta), Real{3});
    EXPECT_REAL_NEAR(scaled.gap() / band.gap(), Real{3}, Real{1024} * eps);
  }
}

TYPED_TEST(XXZDispersion, GaplessAndIsingLimits)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  for (Real delta : {-Real{0.75}, Real{0}, Real{0.5}, Real{1}})
  {
    bethe::xxz::SpinonDispersion<Real> const band(delta);
    EXPECT_EQ(band.gap(), Real{0});
    for (unsigned i = 0; i <= 16; ++i)
    {
      Real const q = pi * Real(i) / Real{16};
      auto const edges = band.continuum(q);
      EXPECT_REAL_NEAR(band.energy(q), bethe::xxz::spinon_energy(q, delta), Real{32} * eps);
      EXPECT_REAL_NEAR(edges.lower, band.energy(q), Real{32} * eps);
      EXPECT_REAL_NEAR(edges.upper, Real{2} * band.energy(q / Real{2}), Real{32} * eps);
    }
  }
  Real const large = Real{1000000};
  bethe::xxz::SpinonDispersion<Real> const ising(large);
  EXPECT_REAL_NEAR(ising.gap(), large / Real{2} - Real{1} + Real{1} / (Real{4} * large),
                   Real{1} / (large * large) + Real{128} * eps * large);
  EXPECT_REAL_NEAR(ising.maximum_energy(), large / Real{2} + Real{1} + Real{1} / (Real{4} * large),
                   Real{1} / (large * large) + Real{128} * eps * large);
  bethe::xxz::SpinonDispersion<Real> const near_critical(Real{1} + Real{1} / Real{4096});
  EXPECT_GT(near_critical.gap(), Real{0});
  EXPECT_REAL_NEAR(near_critical.maximum_energy(), pi / Real{2}, Real{1} / Real{4096});
  Real const enormous = uni20::parse_real<Real>("1e100");
  bethe::xxz::SpinonDispersion<Real> const flat(enormous);
  EXPECT_REAL_NEAR(flat.maximum_energy() / enormous, Real{0.5}, Real{8} * eps);
  EXPECT_REAL_NEAR(flat.gap() / enormous, Real{0.5}, Real{8} * eps);
  auto const bulk_ising = bethe::xxz::bulk_energy_density(enormous);
  ASSERT_TRUE(bulk_ising.converged);
  EXPECT_REAL_NEAR(*bulk_ising.energy / enormous, -Real{0.25}, Real{8} * eps);
}

TYPED_TEST(XXZDispersion, ContinuumAgainstIndependentMomentumScan)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  // Test the shared kinematic helper directly, with no elliptic evaluation.
  // At Q>pi/2 the lower edge is m+e(Q), not 2*e((pi-Q)/2).
  for (Real mass : {Real{0}, Real{1} / Real{1000}, Real{0.2}, Real{0.8}, Real{1}})
  {
    bethe::detail::SpinonBand<Real> const band(Real{1}, mass);
    for (unsigned j = 0; j <= 20; ++j)
    {
      Real const q = pi * Real(j) / Real{20};
      Real low = Real{3}, high = Real{0};
      for (unsigned i = 0; i <= 2048; ++i)
      {
        Real const p = q * Real(i) / Real{2048};
        Real const energy = band.energy(p) + band.energy(q - p);
        low = std::min(low, energy);
        high = std::max(high, energy);
      }
      auto const edge = band.continuum(q);
      EXPECT_LE(edge.lower, low + Real{128} * eps);
      EXPECT_GE(edge.upper, high - Real{128} * eps);
      // |de/dp|<=I: the scan is within twice its mesh spacing of an extremum.
      EXPECT_REAL_NEAR(edge.lower, low, Real{2} * pi / Real{2048});
      EXPECT_REAL_NEAR(edge.upper, high, Real{2} * pi / Real{2048});
      auto const shifted = band.continuum(pi - q);
      auto const folded = band.continuum(q, bethe::SpinonMomentum::folded);
      EXPECT_EQ(folded.lower, std::min(edge.lower, shifted.lower));
      EXPECT_EQ(folded.upper, std::max(edge.upper, shifted.upper));
      auto const reflected = band.continuum(Real{2} * pi - q);
      EXPECT_REAL_NEAR(reflected.lower, edge.lower, Real{64} * eps);
      EXPECT_REAL_NEAR(reflected.upper, edge.upper, Real{64} * eps);
    }
  }
}

TYPED_TEST(XXZDispersion, TinyMassTransitionAndFailures)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  bethe::detail::SpinonBand<Real> const tiny(Real{1}, eps * eps);
  Real const transition = Real{2} * std::atan(eps);
  auto const below = tiny.continuum(transition / Real{2});
  EXPECT_EQ(below.lower, Real{2} * tiny.energy(transition / Real{4}));
  auto const above = tiny.continuum(Real{2} * transition);
  EXPECT_REAL_NEAR(above.lower / std::sin(Real{2} * transition), Real{1}, Real{64} * eps);
  EXPECT_THROW(bethe::xxz::SpinonDispersion<Real>(-Real{1}), std::invalid_argument);
  EXPECT_THROW((bethe::xxz::SpinonDispersion<Real>(Real{2}, Real{0})), std::invalid_argument);
  EXPECT_THROW(bethe::xxz::SpinonDispersion<Real>{uni20::numeric_limits<Real>::infinity()}, std::invalid_argument);
  EXPECT_THROW(bethe::xxz::SpinonDispersion<Real>(Real{1} + eps), std::underflow_error);
  EXPECT_THROW(tiny.energy(-eps), std::invalid_argument);
  EXPECT_THROW(tiny.energy(Real{2} * pi), std::invalid_argument);
  EXPECT_THROW(tiny.continuum(-eps), std::invalid_argument);
  EXPECT_THROW(tiny.continuum(uni20::numeric_limits<Real>::quiet_NaN()), std::invalid_argument);
  bethe::detail::SpinonBand<Real> const huge(uni20::numeric_limits<Real>::max(), Real{0});
  EXPECT_THROW(huge.continuum(pi), std::overflow_error);
}

TYPED_TEST(XXZDispersion, BulkEnergyIndependentReferences)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  struct Fixture
  {
      char const *delta, *energy;
  };
  Fixture const fixtures[] = {{"1.0078125", "-0.444302644254333433552569624296847819038663213571795387869335226226313"},
                              {"1.125", "-0.461984157634790441030839697958010732438226735259084766039838672090467"},
                              {"2", "-0.617222045975865643852780656734321406704332851676111761438426802683836"},
                              {"4", "-1.06152368348350153416525074935899958726966867551680265811725009239714"},
                              {"100", "-25.0024999375000000390644531738269040374542223548331250907367689333422"},
                              {"-0.75", "-0.258815227610764176177266609981691875265553453477375827131353522390755"},
                              {"-0.5", "-0.274519052838328985072792378064702137603551970178892735520927617294475"},
                              {"0.999", "-0.442999488316368675703345472070812426605156134303499126053649792444211"},
                              {"1.0001", "-0.443161952367489168036439661170345717641228856736460233892309954427280"},
                              {"0.5", "-0.375"}};
  for (auto const& f : fixtures)
  {
    SCOPED_TRACE(f.delta);
    auto const energy = bethe::xxz::bulk_energy_density(uni20::parse_real<Real>(f.delta));
    ASSERT_TRUE(energy.converged);
    ASSERT_TRUE(energy.energy);
    EXPECT_REAL_NEAR(*energy.energy, uni20::parse_real<Real>(f.energy), *energy.error);
    auto const scaled = bethe::xxz::bulk_energy_density(uni20::parse_real<Real>(f.delta), Real{3});
    ASSERT_TRUE(scaled.converged);
    EXPECT_REAL_NEAR(*scaled.energy, Real{3} * *energy.energy,
                     Real{1024} * eps * std::max(Real{1}, std::abs(*scaled.energy)));
  }
  auto const xxx = bethe::xxz::bulk_energy_density(Real{1});
  EXPECT_EQ(*xxx.energy, bethe::heisenberg::bulk_energy_density<Real>());
  auto const xx = bethe::xxz::bulk_energy_density(Real{0});
  EXPECT_EQ(*xx.energy, -Real{1} / bethe::detail::pi<Real>());
  // The bulk integral must stay affordable even where the mass underflows.
  for (Real delta : {Real{1} - eps, Real{1} + eps})
  {
    auto const near = bethe::xxz::bulk_energy_density(delta);
    ASSERT_TRUE(near.converged);
    EXPECT_LT(near.evaluations, 10000u);
    EXPECT_REAL_NEAR(*near.energy, *xxx.energy, Real{256} * eps);
  }
  auto const ferro_limit = bethe::xxz::bulk_energy_density(-Real{1} + eps);
  ASSERT_TRUE(ferro_limit.converged);
  EXPECT_REAL_NEAR(*ferro_limit.energy, -Real{1} / Real{4}, Real{256} * eps);
}

TYPED_TEST(XXZDispersion, BulkResummationAndBudgetFailures)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real delta : {Real{1.25}, Real{1.5}, Real{1.625}})
  {
    Real const eta = std::acosh(delta);
    bethe::detail::CompensatedSum<Real> terms;
    // Independent direct Fourier sum overlaps both sides of eta=1.
    for (unsigned n = 1; n < 256; ++n)
      terms.add(Real{1} / (Real{1} + std::exp(Real{2} * Real(n) * eta)));
    Real const expected = delta / Real{4} - std::sinh(eta) * (Real{0.5} + Real{2} * terms.value());
    auto const energy = bethe::xxz::bulk_energy_density(delta);
    ASSERT_TRUE(energy.converged);
    EXPECT_REAL_NEAR(*energy.energy, expected, Real{256} * eps);
    bethe::xxz::BulkOptions<Real> controls;
    controls.max_evaluations = 0;
    auto const failed = bethe::xxz::bulk_energy_density(delta, Real{1}, controls);
    EXPECT_FALSE(failed.converged);
    EXPECT_EQ(failed.status, bethe::xxz::BulkStatus::evaluation_limit);
    EXPECT_FALSE(failed.energy);
    EXPECT_FALSE(failed.error);
    EXPECT_EQ(failed.evaluations, 0u);
    controls.max_evaluations = 200000;
    controls.max_levels = 0;
    if (eta < Real{1})
      EXPECT_EQ(bethe::xxz::bulk_energy_density(delta, Real{1}, controls).status, bethe::xxz::BulkStatus::mesh_limit);
    controls.tolerance = eps;
    EXPECT_THROW(bethe::xxz::bulk_energy_density(delta, Real{1}, controls), std::invalid_argument);
  }
}
} // namespace
