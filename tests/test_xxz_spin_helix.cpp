// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/xxz_spin_helix.hpp>

#include "test_support.hpp"
#include <bethe/xxz_helix_check.hpp>
#include <bethe/xxz_odd_sectors.hpp>
#include <bit>

namespace
{
namespace xxz = bethe::xxz;
namespace engine = xxz::detail;
template <typename Real> class XXZSpinHelix : public ::testing::Test {};
TYPED_TEST_SUITE(XXZSpinHelix, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXZSpinHelix, EveryWindingAndSectorAgainstHamiltonian)
{
  using Real = TypeParam;
  using Complex = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  for (unsigned n = 2; n <= 7; ++n)
    for (unsigned winding = 0; winding < n; ++winding)
      for (unsigned m = 0; m <= n; ++m)
      {
        auto const state =
            xxz::periodic_spin_helix<Real>(n, winding, uni20::from_twice(std::int64_t(n) - 2 * std::int64_t(m)));
        EXPECT_EQ(state.down_spins, m);
        EXPECT_EQ(state.momentum_index, (m * winding) % n);
        EXPECT_REAL_NEAR(state.delta, std::cos(Real{2} * pi * Real(winding) / Real(n)), Real{16} * eps);
        auto amplitude = [&](unsigned bits) {
          std::vector<std::size_t> occupied;
          for (unsigned j = 0; j < n; ++j)
            if (bits & (1U << j)) occupied.push_back(j);
          return state.unnormalized_amplitude(occupied);
        };
        for (unsigned bits = 0; bits < (1U << n); ++bits)
          if (std::popcount(bits) == int(m))
          {
            Complex const wave = amplitude(bits);
            EXPECT_REAL_NEAR(std::abs(wave), Real{1}, Real{8} * eps);
            Complex action{};
            for (unsigned j = 0; j < n; ++j)
            {
              auto const k = (j + 1) % n;
              bool const anti = ((bits >> j) & 1U) != ((bits >> k) & 1U);
              action += state.delta * (anti ? -Real{1} : Real{1}) * wave / Real{4};
              if (anti) action += amplitude(bits ^ (1U << j) ^ (1U << k)) / Real{2};
            }
            EXPECT_LT(std::abs(action - state.energy * wave), Real{64} * Real(n) * eps);
            auto const translated = ((bits << 1) & ((1U << n) - 1)) | (bits >> (n - 1));
            Complex const phase{std::cos(state.momentum), std::sin(state.momentum)};
            EXPECT_LT(std::abs(amplitude(translated) - phase * wave), Real{64} * eps);
          }
      }
}

TYPED_TEST(XXZSpinHelix, LargeIntegerPhasesDoNotOverflow)
{
  using Real = TypeParam;
  auto const n = std::size_t(std::numeric_limits<std::int64_t>::max() / 4);
  auto const large = xxz::periodic_spin_helix<Real>(n, n - 1, uni20::from_twice(std::int64_t{1}));
  EXPECT_EQ(large.momentum_index, n - n / 2);
  auto const two = xxz::periodic_spin_helix<Real>(n, n - 1, uni20::from_twice(std::int64_t(n - 4)));
  EXPECT_EQ(two.momentum_index, n - 2);
  std::vector<std::size_t> occupied{n - 2, n - 1};
  auto const amplitude = two.unnormalized_amplitude(occupied);
  Real const angle = Real{24} * std::atan(Real{1}) / Real(n);
  EXPECT_REAL_NEAR(amplitude.imag() / angle, Real{1}, Real{8} * uni20::numeric_limits<Real>::epsilon());
}

TYPED_TEST(XXZSpinHelix, RepeatedRootPolynomialInAffineCoordinates)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {5, 7, 9})
    for (std::size_t winding : {n / 2, n / 2 + 1})
    {
      auto const state = xxz::periodic_spin_helix<Real>(n, winding, uni20::from_twice(std::int64_t{1}));
      Real const a = std::sqrt((Real{1} + state.delta) / (Real{1} - state.delta));
      Real const z = winding <= n / 2 ? a : -a;
      for (std::size_t m = 1; m <= n / 2; ++m)
      {
        // Q(x)=(x-r)^M, with a nontrivial coordinate shift and scale.
        Real const scale = Real{3} / Real{7}, center = z - scale / Real{3}, r = Real{1} / Real{3};
        std::vector<Real> c(m);
        Real coefficient = Real{1};
        for (std::size_t j = m; j-- > 0;)
        {
          coefficient *= -r * Real(j + 1) / Real(m - j);
          c[j] = coefficient;
        }
        engine::PolynomialBetheSystem<Real> const system(n, m, center, scale);
        auto const match = engine::check_helix_polynomial<Real>(system, c, state.delta, winding);
        EXPECT_EQ(match.status, engine::HelixMatchStatus::compatible);
        EXPECT_LT(match.coefficient_error, Real{128} * eps);
        EXPECT_EQ(match.coupling_error, Real{0});
        EXPECT_EQ(match.eigenvector_defect_bound, Real{0});
        EXPECT_REAL_NEAR(system.energy(c, state.delta), state.energy, Real{128} * Real(n) * eps);
        EXPECT_EQ(engine::check_helix_polynomial<Real>(system, c, state.delta, n - winding).status,
                  engine::HelixMatchStatus::different_polynomial);
        c[0] += Real{1} / Real{1000};
        EXPECT_EQ(engine::check_helix_polynomial<Real>(system, c, state.delta, winding).status,
                  engine::HelixMatchStatus::different_polynomial);
      }
    }
}

TYPED_TEST(XXZSpinHelix, MatchingDoesNotSnapNearbyCouplings)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  auto const state = xxz::periodic_spin_helix<Real>(7, 4, uni20::from_twice(std::int64_t{1}));
  Real const z = -std::tan(Real{2} * std::atan(Real{1}) / Real{7});
  engine::PolynomialBetheSystem<Real> const system(7, 3, z);
  std::vector<Real> const c(3, Real{0});
  Real const nearby = state.delta + Real{8} * eps;
  auto const match = engine::check_helix_polynomial<Real>(system, c, nearby, 4);
  EXPECT_EQ(match.status, engine::HelixMatchStatus::compatible);
  EXPECT_EQ(match.helix_delta, state.delta);
  EXPECT_GT(match.coupling_error, Real{0});
  EXPECT_EQ(match.eigenvector_defect_bound, Real{7} * std::abs(nearby - state.delta) / Real{4});
  EXPECT_EQ(engine::check_helix_polynomial<Real>(system, c, nearby, 4, eps).status,
            engine::HelixMatchStatus::different_coupling);
  EXPECT_EQ(engine::check_helix_polynomial<Real>(system, c, state.delta + Real{1} / Real{1000}, 4).status,
            engine::HelixMatchStatus::different_coupling);
  // Another root-of-unity coupling is not this all-phantom state.
  EXPECT_EQ(engine::check_helix_polynomial<Real>(system, c, -Real{1} / Real{2}, 4).status,
            engine::HelixMatchStatus::different_coupling);
}

TYPED_TEST(XXZSpinHelix, ContinuedPhantomSectorsAreRecognized)
{
  using Real = TypeParam;
  for (std::size_t n : {3, 5, 7, 9, 13, 17, 21})
  {
    SCOPED_TRACE(::testing::Message() << "N=" << n);
    auto const helix = xxz::periodic_spin_helix<Real>(n, (n + 1) / 2, uni20::from_twice(std::int64_t{1}));
    auto const scan = engine::scan_odd_polynomial_sectors(n, helix.delta, {.max_iterations = 500});
    ASSERT_TRUE(scan.equations_complete);
    EXPECT_FALSE(scan.regular_states_complete);
    EXPECT_TRUE(scan.state_checks_complete);
    for (std::size_t m = 1; m < scan.sectors.size(); ++m)
    {
      auto const& sector = scan.sectors[m];
      ASSERT_TRUE(sector.helix);
      EXPECT_EQ(sector.helix->status, engine::HelixMatchStatus::compatible)
          << "M=" << m << " coefficient error=" << uni20::format_scalar(sector.helix->coefficient_error);
      EXPECT_TRUE(sector.phantom_lifts.empty()); // The direct helix construction suffices.
      EXPECT_EQ(sector.branch.delta, helix.delta);
      EXPECT_EQ(sector.branch.momentum_index, (m * ((n + 1) / 2)) % n);
    }
  }
  auto const other = engine::scan_odd_polynomial_sectors(7, -Real{1} / Real{2});
  ASSERT_TRUE(other.sectors.back().helix);
  EXPECT_EQ(other.sectors.back().helix->status, engine::HelixMatchStatus::different_coupling);
  ASSERT_FALSE(other.sectors.back().phantom_lifts.empty());
  EXPECT_EQ(other.sectors.back().phantom_lifts.back().status, engine::PhantomLiftStatus::nonzero_witness);
  EXPECT_TRUE(other.state_checks_complete); // Resolved by a mixed witness, not a helix match.
  auto const failed = engine::scan_odd_polynomial_sectors(7, -Real{7} / Real{10}, {.max_iterations = 0});
  EXPECT_FALSE(failed.state_checks_complete);
  for (auto const& sector : failed.sectors)
    if (!sector.branch.equations_converged) EXPECT_FALSE(sector.helix);
}

TYPED_TEST(XXZSpinHelix, UnrepresentableAffinePolynomialIsUnresolved)
{
  using Real = TypeParam;
  auto const state = xxz::periodic_spin_helix<Real>(5, 3, uni20::from_twice(std::int64_t{1}));
  engine::PolynomialBetheSystem<Real> const system(5, 2, uni20::numeric_limits<Real>::max() / Real{2});
  std::vector<Real> const c(2, Real{0});
  auto const result = engine::check_helix_polynomial<Real>(system, c, state.delta, 3);
  EXPECT_EQ(result.status, engine::HelixMatchStatus::nonfinite);
}

TYPED_TEST(XXZSpinHelix, InvalidInputs)
{
  using Real = TypeParam;
  auto const half = uni20::from_twice(std::int64_t{1});
  EXPECT_THROW(xxz::periodic_spin_helix<Real>(1, 0, half), std::invalid_argument);
  EXPECT_THROW(xxz::periodic_spin_helix<Real>(5, 5, half), std::invalid_argument);
  EXPECT_THROW(xxz::periodic_spin_helix<Real>(5, 2, uni20::half_int{0}), std::invalid_argument);
  EXPECT_THROW(xxz::periodic_spin_helix<Real>(5, 2, uni20::half_int{3}), std::invalid_argument);
  auto const state = xxz::periodic_spin_helix<Real>(5, 2, half);
  for (std::vector<std::size_t> const& occupied : {std::vector<std::size_t>{0}, {0, 0}, {1, 0}, {0, 5}})
    EXPECT_THROW(state.unnormalized_amplitude(occupied), std::invalid_argument);
  engine::PolynomialBetheSystem<Real> const system(5, 2);
  std::vector<Real> c(2, Real{0});
  for (Real bad :
       {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW(engine::check_helix_polynomial<Real>(system, c, state.delta, 2, bad), std::invalid_argument);
  for (Real bad :
       {-Real{1}, Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW(engine::check_helix_polynomial<Real>(system, c, bad, 2), std::invalid_argument);
  EXPECT_THROW(engine::check_helix_polynomial<Real>(system, c, state.delta, 0), std::invalid_argument);
  EXPECT_THROW(engine::check_helix_polynomial<Real>(system, c, state.delta, 5), std::invalid_argument);
  c.pop_back();
  EXPECT_THROW(engine::check_helix_polynomial<Real>(system, c, state.delta, 2), std::invalid_argument);
}
} // namespace
