// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "exact_spectrum.hpp"
#include "test_support.hpp"
#include <bethe/xxz_odd_sectors.hpp>

namespace
{
namespace engine = bethe::xxz::detail;
template <typename Real> class XXZOddSectors : public ::testing::Test {};
TYPED_TEST_SUITE(XXZOddSectors, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXZOddSectors, EverySectorAndGlobalMinimumAgainstED)
{
  using Real = TypeParam;
  for (unsigned n : {3, 5, 7, 9})
    for (Real d : {Real{0}, -Real{1} / Real{5}, -Real{7} / Real{10}, -Real{9} / Real{10}, -Real{999} / Real{1000}})
    {
      SCOPED_TRACE(::testing::Message() << "N=" << n << " Delta=" << uni20::format_scalar(d));
      auto const scan = engine::scan_odd_polynomial_sectors(n, d);
      ASSERT_TRUE(scan.equations_complete);
      EXPECT_TRUE(scan.regular_states_complete);
      EXPECT_TRUE(scan.state_checks_complete);
      ASSERT_TRUE(scan.lowest_index);
      ASSERT_EQ(scan.sectors.size(), n / 2 + 1);
      double lowest = std::numeric_limits<double>::infinity();
      std::size_t updates = 0;
      for (unsigned m = 0; m <= n / 2; ++m)
      {
        auto const& sector = scan.sectors[m];
        EXPECT_EQ(sector.sz, uni20::from_twice(std::int64_t(n - 2 * m)));
        EXPECT_EQ(sector.branch.coefficients.size(), m);
        EXPECT_TRUE(sector.branch.equations_converged);
        EXPECT_TRUE(sector.wronskian);
        EXPECT_TRUE(sector.regularity);
        auto const spectrum = test_support::exact_spectrum(n, m, 0, true, static_cast<double>(d));
        EXPECT_REAL_NEAR(static_cast<double>(sector.branch.energy), spectrum.front(), 3e-10);
        lowest = std::min(lowest, spectrum.front());
        updates += sector.branch.iterations;
      }
      EXPECT_EQ(scan.iterations, updates);
      EXPECT_REAL_NEAR(static_cast<double>(scan.sectors[*scan.lowest_index].branch.energy), lowest, 3e-10);
      EXPECT_NE(std::find(scan.nearby_indices.begin(), scan.nearby_indices.end(), *scan.lowest_index),
                scan.nearby_indices.end());
      for (auto const& sector : scan.sectors)
        EXPECT_LE(scan.sectors[*scan.lowest_index].branch.energy, sector.branch.energy);
    }
}

TYPED_TEST(XXZOddSectors, NativeFiveSiteSectorCrossing)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  Real const C = std::cos(pi / Real{5});
  for (Real d : {-Real{1} / Real{5}, -Real{7} / Real{10}, -Real{9} / Real{10}, -Real{1} + Real{128} * eps})
  {
    auto const scan = engine::scan_odd_polynomial_sectors(5, d);
    ASSERT_TRUE(scan.lowest_index);
    Real const vacuum = Real{5} * d / Real{4};
    Real const energies[] = {vacuum, vacuum - d - C,
                             vacuum + (-Real{3} * d - C - std::sqrt((d + C) * (d + C) + Real{4} * C * C)) / Real{2}};
    for (std::size_t m = 0; m < 3; ++m)
      EXPECT_REAL_NEAR(scan.sectors[m].branch.energy, energies[m], Real{2048} * eps);
    auto const expected = d < -C ? std::size_t{0} : std::size_t{2};
    EXPECT_EQ(*scan.lowest_index, expected);
    EXPECT_EQ(scan.nearby_indices, std::vector<std::size_t>{expected});
    // The strongest regression: returning smallest |Sz| at -0.9 would give
    // approximately -0.989803 instead of the polarized ground energy -1.125.
    EXPECT_REAL_NEAR(scan.sectors[*scan.lowest_index].branch.energy, energies[expected], Real{2048} * eps);
  }
}

TYPED_TEST(XXZOddSectors, LargerAllSectorReferences)
{
  using Real = TypeParam;
  // Independent sparse spin-basis ED, including sectors absent from the
  // earlier minimal-|Sz| audit. Regenerate with reference_xxz_odd_ed.py
  // --all-sectors --sites 13 17 --deltas -0.7 -0.99. These are fp64 oracles;
  // native-precision analytic tests above retain a separate tighter contract.
  std::vector<std::vector<double>> const references{
      {-2.2749999999999999, -2.5459418174260531, -2.7862137983442641, -2.9862027696162934, -3.1403482283902431,
       -3.245000812454792, -3.2978800209497248},
      {-3.2174999999999998, -3.1984418174260529, -3.1826571736276357, -3.1700968826370568, -3.1607180848016134,
       -3.1544861327781684, -3.1513765211381943},
      {-2.9749999999999996, -3.2579730996839023, -3.5202327880400834, -3.7537578355204055, -3.9538868341184363,
       -4.1173893086771054, -4.2419248441994482, -4.3258360207149629, -4.3680541417037979},
      {-4.2074999999999996, -4.2004730996839026, -4.1943384503462049, -4.1890908277270995, -4.184725535227467,
       -4.181238486811834, -4.1786262893819766, -4.176886319417636, -4.1760167879441754}};
  for (std::size_t i = 0; i < references.size(); ++i)
  {
    std::size_t const n = i < 2 ? 13 : 17;
    Real const d = i % 2 ? -Real{99} / Real{100} : -Real{7} / Real{10};
    SCOPED_TRACE(::testing::Message() << "N=" << n << " Delta=" << uni20::format_scalar(d));
    auto const scan = engine::scan_odd_polynomial_sectors(n, d, {.max_iterations = 500});
    ASSERT_EQ(scan.sectors.size(), references[i].size());
    for (std::size_t m = 0; m < references[i].size(); ++m)
    {
      auto const& branch = scan.sectors[m].branch;
      ASSERT_TRUE(branch.equations_converged)
          << "M=" << m << " status=" << int(branch.status) << " root_delta=" << uni20::format_scalar(branch.root_delta);
      EXPECT_REAL_NEAR(static_cast<double>(branch.energy), references[i][m], 3e-10);
    }
    ASSERT_TRUE(scan.equations_complete);
    ASSERT_TRUE(scan.lowest_index);
    EXPECT_EQ(*scan.lowest_index, i % 2 ? std::size_t{0} : n / 2);
  }
}

TYPED_TEST(XXZOddSectors, CollisionReportsNearbySectorsWithoutChangingEnergies)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1});
  for (std::size_t n : {3, 5, 7, 9})
  {
    Real const d = -std::cos(pi / Real(n));
    auto const scan = engine::scan_odd_polynomial_sectors(n, d);
    ASSERT_TRUE(scan.equations_complete);
    ASSERT_TRUE(scan.lowest_index);
    EXPECT_EQ(scan.nearby_indices.size(), scan.sectors.size());
    for (auto const& sector : scan.sectors)
    {
      EXPECT_REAL_NEAR(sector.branch.energy, Real(n) * d / Real{4}, scan.comparison_band);
      auto const direct = engine::continue_odd_polynomial(n, d, sector.sz);
      EXPECT_EQ(sector.branch.energy, direct.energy);
    }
  }
}

TYPED_TEST(XXZOddSectors, FailedSectorsNeverBecomeAnIncompleteMinimum)
{
  using Real = TypeParam;
  Real const d = -Real{9} / Real{10};
  for (std::size_t budget : {0, 1, 3})
  {
    auto const scan = engine::scan_odd_polynomial_sectors(9, d, {.max_iterations = budget});
    EXPECT_FALSE(scan.equations_complete);
    EXPECT_FALSE(scan.wronskians_consistent);
    EXPECT_FALSE(scan.regular_states_complete);
    EXPECT_FALSE(scan.state_checks_complete);
    EXPECT_FALSE(scan.lowest_index);
    EXPECT_TRUE(scan.nearby_indices.empty());
    EXPECT_EQ(scan.comparison_band, Real{0});
    ASSERT_EQ(scan.sectors.size(), std::size_t{5});
    EXPECT_EQ(scan.iterations, 3 * budget);
    for (std::size_t m = 0; m < scan.sectors.size(); ++m)
    {
      auto const& sector = scan.sectors[m];
      EXPECT_EQ(sector.branch.equations_converged, m <= 1);
      EXPECT_EQ(bool(sector.wronskian), m <= 1);
      EXPECT_EQ(bool(sector.regularity), m <= 1);
      if (m > 1)
      {
        EXPECT_EQ(sector.branch.status, engine::PolynomialContinuationStatus::iteration_limit);
        EXPECT_EQ(sector.branch.iterations, budget);
        engine::PolynomialBetheSystem<Real> const system(9, m, sector.branch.center, sector.branch.coordinate_scale);
        EXPECT_EQ(sector.branch.energy, system.energy(sector.branch.coefficients, d));
      }
    }
  }
}

TYPED_TEST(XXZOddSectors, AdmissibilityDiagnosticsRemainSeparate)
{
  using Real = TypeParam;
  // At a root of unity the generic Wronskian can be inconclusive even for
  // the vacuum. Do not turn that into a fake continuation failure or claim
  // that all equation-converged states have been physically certified.
  auto const scan = engine::scan_odd_polynomial_sectors(7, -Real{1} / Real{2});
  ASSERT_TRUE(scan.equations_complete);
  ASSERT_TRUE(scan.lowest_index);
  EXPECT_FALSE(scan.wronskians_consistent);
  ASSERT_TRUE(scan.sectors.front().wronskian);
  EXPECT_EQ(scan.sectors.front().wronskian->status, engine::WronskianStatus::ill_conditioned);
  ASSERT_TRUE(scan.sectors.front().regularity);
  EXPECT_EQ(scan.sectors.front().regularity->status, engine::RegularityStatus::regular_on_shell);
  EXPECT_FALSE(scan.regular_states_complete);
  EXPECT_FALSE(scan.state_checks_complete);

  // Merely changing the independent diagnostic tolerance cannot change the
  // continuation, selected energy, Newton budget, or near-degeneracy report.
  auto const tight = engine::scan_odd_polynomial_sectors(
      5, -Real{7} / Real{10}, {}, uni20::numeric_limits<Real>::epsilon() * uni20::numeric_limits<Real>::epsilon());
  auto const loose = engine::scan_odd_polynomial_sectors(5, -Real{7} / Real{10}, {}, Real{1});
  EXPECT_EQ(tight.lowest_index, loose.lowest_index);
  EXPECT_EQ(tight.iterations, loose.iterations);
  EXPECT_EQ(tight.nearby_indices, loose.nearby_indices);
  for (std::size_t m = 0; m < tight.sectors.size(); ++m)
    EXPECT_EQ(tight.sectors[m].branch.energy, loose.sectors[m].branch.energy);
}

TYPED_TEST(XXZOddSectors, InvalidInput)
{
  using Real = TypeParam;
  Real const d = -Real{7} / Real{10};
  for (std::size_t n : {0, 1, 2, 6})
    EXPECT_THROW(engine::scan_odd_polynomial_sectors(n, d), std::invalid_argument);
  for (Real bad :
       {-Real{1}, Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW(engine::scan_odd_polynomial_sectors(5, bad), std::invalid_argument);
  for (Real bad :
       {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    EXPECT_THROW(engine::scan_odd_polynomial_sectors(5, d, {}, bad), std::invalid_argument);
    EXPECT_THROW(engine::scan_odd_polynomial_sectors(5, d, {.residual_tolerance = bad}), std::invalid_argument);
  }
}
} // namespace
