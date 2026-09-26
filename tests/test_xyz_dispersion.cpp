// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxz_dispersion.hpp>
#include <bethe/xyz.hpp>
#include <bethe/xyz_dispersion.hpp>
#include <bit>
#include <uni20/linalg/ops/self_adjoint_eigh.hpp>

namespace
{
template <typename Real> class XYZDispersion : public ::testing::Test {};
TYPED_TEST_SUITE(XYZDispersion, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XYZDispersion, IndependentEllipticReferences)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  // mpmath 65-digit evaluation in JKM's cn/dn/sn parametrization, not the
  // simplified theta derivative prefactor used by the implementation.
  struct Fixture
  {
      char const *eta, *t, *maximum, *gap, *bound1, *bound2;
  };
  Fixture const fixtures[] = {{"0.75", "1", "0.50231250248767710143496779750111542398194818522710862683886191805",
                               "0.233183237167837415865617249661776631695819568738375906775664238",
                               "0.47394351390478720711395987040107571233531207985859309211578682791",
                               "0.52548335561148757239797255701204032628025013237732926803857645261"},
                              {"0.875", "0.75", "0.28343933547384184378281333767383562373510220535323708988857414324",
                               "0.23090927993737522632365395297971984612540171015905540221778031329",
                               "0.35258192520998911237697126764789897299955950024544368557695131339",
                               "0.34343293076329186174178406726122240471569562115197293043952809318"}};
  for (auto const& f : fixtures)
  {
    SCOPED_TRACE(f.eta);
    auto parse = [](char const* s) { return uni20::parse_real<Real>(s); };
    bethe::xyz::SpinonDispersion<Real> const band(parse(f.eta), parse(f.t));
    EXPECT_REAL_NEAR(band.maximum_energy(), parse(f.maximum), Real{128} * eps);
    EXPECT_REAL_NEAR(band.gap(), parse(f.gap), Real{128} * eps);
    EXPECT_REAL_NEAR(band.bound_energy(1, parse("0.7")), parse(f.bound1), Real{128} * eps);
    EXPECT_REAL_NEAR(band.bound_energy(2, parse("0.7")), parse(f.bound2), Real{128} * eps);
    EXPECT_EQ(band.energy(Real{0}), band.gap());
    EXPECT_EQ(band.energy(bethe::detail::pi<Real>()), band.gap());
  }
}

TYPED_TEST(XYZDispersion, FreeFermionXYAndGaplessXXZLimits)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  for (Real t : {Real{0.25}, Real{1}, Real{4}})
  {
    auto const j = bethe::xyz::couplings(Real{0.5}, t);
    bethe::xyz::SpinonDispersion<Real> const band(Real{0.5}, t);
    EXPECT_REAL_NEAR(band.maximum_energy(), (j.x + j.y) / Real{2}, Real{256} * eps * j.x);
    EXPECT_REAL_NEAR(band.gap(), (j.x - j.y) / Real{2}, Real{256} * eps * j.x);
    for (unsigned i = 0; i <= 8; ++i)
    {
      Real const p = pi * Real(i) / Real{8};
      Real const exact = std::hypot((j.x - j.y) * std::cos(p), (j.x + j.y) * std::sin(p)) / Real{2};
      EXPECT_REAL_NEAR(band.energy(p), exact, Real{256} * eps * j.x);
    }
    EXPECT_FALSE(band.bound_exists(1));
  }
  for (Real eta : {Real{0.25}, Real{0.5}, Real{0.75}})
  {
    bethe::xyz::SpinonDispersion<Real> const band(eta, Real{80});
    // Jx=Jy=1, Jz=cos(pi*eta), up to corrections exp(-pi*t).
    Real const maximum = std::sin(pi * eta) / (Real{2} * eta);
    EXPECT_REAL_NEAR(band.maximum_energy(), maximum, Real{128} * eps);
    EXPECT_LT(band.gap(), eps);
    EXPECT_GT(band.gap(), Real{0});
    EXPECT_REAL_NEAR(band.energy(pi / Real{4}), maximum / std::sqrt(Real{2}), Real{128} * eps);
  }
}

TYPED_TEST(XYZDispersion, BoundExistenceParitiesAndThresholds)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  bethe::xyz::SpinonDispersion<Real> const band(Real{0.75}, Real{1});
  EXPECT_FALSE(band.bound_exists(0));
  EXPECT_TRUE(band.bound_exists(1));
  EXPECT_TRUE(band.bound_exists(2));
  EXPECT_FALSE(band.bound_exists(3)); // exactly at the third branch merger
  EXPECT_EQ(band.bound_x_parity(1), -1);
  EXPECT_EQ(band.bound_x_parity(2), 1);
  EXPECT_LT(band.bound_gap(1), band.bound_gap(2));
  EXPECT_LT(band.bound_gap(2), Real{2} * band.gap());
  for (unsigned s : {1u, 2u})
  {
    EXPECT_EQ(band.bound_gap(s), band.bound_energy(s, Real{0}));
    EXPECT_EQ(band.bound_gap(s), band.bound_energy(s, Real{2} * pi));
    for (unsigned i = 0; i <= 8; ++i)
    {
      Real const q = pi * Real(i) / Real{8};
      EXPECT_GE(band.bound_energy(s, q), band.bound_gap(s));
      EXPECT_REAL_NEAR(band.bound_energy(s, q), band.bound_energy(s, Real{2} * pi - q), Real{256} * eps);
    }
  }
  // Approach the first merger from the attractive side. At a_s=1 the curve
  // tends to the equal-momentum two-spinon line 2*e(Q/2), not the entire edge.
  Real const offset = std::sqrt(eps);
  bethe::xyz::SpinonDispersion<Real> const emerging(Real{0.5} + offset, Real{1});
  EXPECT_TRUE(emerging.bound_exists(1));
  for (Real q : {Real{0}, pi / Real{3}, pi})
    EXPECT_REAL_NEAR(emerging.bound_energy(1, q), Real{2} * emerging.energy(q / Real{2}), Real{1024} * eps);

  // A longer tower checks the general index dependence, not only s=1,2.
  // Independent Jacobi-function reference script, eta=7/8, t=3/4.
  bethe::xyz::SpinonDispersion<Real> const tower(Real{0.875}, Real{0.75});
  char const* gaps[] = {"0.13058224205719948561933329089134799802139205814300486511872317612",
                        "0.24480853474942555156885209195082721978439558510565732930863934119",
                        "0.33334194230605970477759809988069632026208397843539297023900629882",
                        "0.39516410732250358509732643630221514031955687428335085645548729133",
                        "0.43416576932568663809359693848268586052180043756447034999438847151",
                        "0.45521455909996010055873684230075515555854853887237412528321418502"};
  for (unsigned s = 1; s <= 6; ++s)
  {
    ASSERT_TRUE(tower.bound_exists(s));
    EXPECT_REAL_NEAR(tower.bound_gap(s), uni20::parse_real<Real>(gaps[s - 1]), Real{128} * eps);
  }
  EXPECT_FALSE(tower.bound_exists(7));
  EXPECT_LT(tower.bound_gap(1), tower.gap());
}

TYPED_TEST(XYZDispersion, MassiveAntiferromagneticAndFerromagneticXXZLimits)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  Real const t = Real{1} / Real{100}, lambda = std::acosh(Real{2});
  // After axis rotations and rescaling, t->0 at eta=t*lambda/pi gives
  // AF Delta=2. At eta=1-t*lambda/pi it gives ferro Delta=2 and s magnons.
  for (bool ferro : {false, true})
  {
    Real const eta = ferro ? Real{1} - t * lambda / pi : t * lambda / pi;
    auto const j = bethe::xyz::couplings(eta, t);
    bethe::xyz::SpinonDispersion<Real> const band(eta, t, Real{1} / j.y);
    EXPECT_REAL_NEAR(j.x / j.y, Real{2}, Real{4096} * eps);
    EXPECT_REAL_NEAR(j.z / j.y, ferro ? -Real{1} : Real{1}, Real{4096} * eps);
    if (!ferro)
    {
      bethe::xxz::SpinonDispersion<Real> const xxz(Real{2});
      EXPECT_REAL_NEAR(band.gap(), xxz.gap(), Real{4096} * eps);
      EXPECT_REAL_NEAR(band.maximum_energy(), xxz.maximum_energy(), Real{4096} * eps);
    }
    else
      for (unsigned s : {1u, 2u, 3u, 4u})
        for (unsigned i = 0; i <= 8; ++i)
        {
          Real const q = pi * Real(i) / Real{8};
          Real const exact =
              std::sinh(lambda) / std::sinh(Real(s) * lambda) * (std::cosh(Real(s) * lambda) - std::cos(q));
          EXPECT_REAL_NEAR(band.bound_energy(s, q), exact, Real{8192} * eps);
        }
  }
}

TEST(XYZDispersionExact, FirstBoundBranchFiniteSizeAndMomentumCopies)
{
  // Independent literal spin Hamiltonian. This is a finite-size convergence
  // check in fp64, not a test of high-precision thermodynamic exact equality.
  auto const j = bethe::xyz::couplings(0.75, 1.0);
  bethe::xyz::SpinonDispersion<double> const band(0.75, 1.0);
  double previous_error = 1;
  for (unsigned n : {6u, 8u, 10u})
  {
    SCOPED_TRACE(n);
    unsigned const dim = 1u << n, mask = dim - 1;
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
    auto const exact = uni20::linalg::eigh(std::move(h));
    auto const& v = exact.eigenvectors;
    double const vacuum_parity = n % 4 ? -1 : 1;
    std::vector<double> gaps;
    std::vector<int> z_parities;
    for (unsigned level = 0; level < 8 && gaps.size() < 2; ++level)
    {
      double rx = 0, rz = 0, translation = 0;
      for (unsigned a = 0; a < dim; ++a)
      {
        rx += v[a, level] * v[a ^ mask, level];
        rz += v[a, level] * v[a, level] * (std::popcount(a) % 2 ? -1 : 1);
        unsigned const translated = ((a << 1) & mask) | (a >> (n - 1));
        translation += v[a, level] * v[translated, level];
      }
      if (rx * vacuum_parity > 0) continue;
      EXPECT_NEAR(rx, -vacuum_parity, 1e-10);
      EXPECT_NEAR(std::abs(rz), 1, 1e-10);
      // Relative T=relative Rz at reduced Q=0: the two copies have physical
      // momenta 0 and pi relative to the same reference vacuum.
      EXPECT_NEAR(translation, rz, 1e-10);
      gaps.push_back(exact.eigenvalues[level] - exact.eigenvalues[0]);
      z_parities.push_back(rz > 0 ? 1 : -1);
    }
    ASSERT_EQ(gaps.size(), 2u);
    EXPECT_EQ(z_parities[0], -z_parities[1]);
    EXPECT_LT(gaps[0], band.bound_gap(1));
    EXPECT_GT(gaps[1], band.bound_gap(1));
    double const error = std::max(band.bound_gap(1) - gaps[0], gaps[1] - band.bound_gap(1));
    EXPECT_LT(error, previous_error);
    previous_error = error;
  }
  EXPECT_LT(previous_error, 0.014);
}

TYPED_TEST(XYZDispersion, ScalingAndExplicitFailures)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  bethe::xyz::SpinonDispersion<Real> const band(Real{0.75}, Real{1}), scaled(Real{0.75}, Real{1}, Real{3});
  EXPECT_REAL_NEAR(scaled.gap(), Real{3} * band.gap(), Real{128} * eps);
  EXPECT_REAL_NEAR(scaled.bound_energy(1, Real{0.4}), Real{3} * band.bound_energy(1, Real{0.4}), Real{256} * eps);
  for (auto eta : {Real{0}, Real{1}, -Real{1}})
    EXPECT_THROW((bethe::xyz::SpinonDispersion<Real>(eta, Real{1})), std::invalid_argument);
  EXPECT_THROW((bethe::xyz::SpinonDispersion<Real>(Real{0.5}, Real{0})), std::invalid_argument);
  EXPECT_THROW((bethe::xyz::SpinonDispersion<Real>(Real{0.5}, Real{1}, Real{0})), std::invalid_argument);
  EXPECT_THROW(band.bound_energy(0, Real{0}), std::invalid_argument);
  EXPECT_THROW(band.bound_gap(3), std::invalid_argument);
  EXPECT_THROW(band.bound_x_parity(3), std::invalid_argument);
  EXPECT_THROW(band.bound_energy(1, -Real{1}), std::invalid_argument);
  EXPECT_THROW(band.energy(-Real{1}), std::invalid_argument);
  // t pushes the strictly positive mass below the normal exponent range in
  // every supported scalar. Returning a false gapless zero is not acceptable.
  Real const large_t = -std::log(uni20::numeric_limits<Real>::min());
  EXPECT_THROW((bethe::xyz::SpinonDispersion<Real>(Real{0.25}, large_t)), std::underflow_error);
}
} // namespace
