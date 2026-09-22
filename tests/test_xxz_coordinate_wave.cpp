// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/xxz_coordinate_wave.hpp>

#include "test_support.hpp"
#include <array>
#include <bethe/polynomial_roots.hpp>
#include <bethe/xxz.hpp>
#include <bethe/xxz_odd_continuation.hpp>
#include <numeric>

namespace
{
namespace engine = bethe::xxz::detail;
template <typename Real> class XXZCoordinateWave : public ::testing::Test {};
TYPED_TEST_SUITE(XXZCoordinateWave, test_support::RealTypes, test_support::PrecisionNames);

// Deliberately factorial, test-only reference with direct integer powers.
template <typename Real>
engine::CoordinateAmplitude<Real> permutations(Real d, std::vector<std::complex<Real>> const& v,
                                               std::vector<std::size_t> const& occupied)
{
  using C = std::complex<Real>;
  std::vector<std::size_t> permutation(v.size());
  std::iota(permutation.begin(), permutation.end(), 0);
  Real normalization{1};
  for (std::size_t i = 0; i < v.size(); ++i)
    for (std::size_t j = i + 1; j < v.size(); ++j)
      normalization *= std::max({Real{1}, std::abs(Real{1} + v[i] * v[j] - Real{2} * d * v[i]),
                                 std::abs(Real{1} + v[i] * v[j] - Real{2} * d * v[j])});
  engine::CoordinateAmplitude<Real> sum;
  do
  {
    C term{1};
    for (std::size_t a = 0; a < v.size(); ++a)
    {
      auto const i = permutation[a];
      for (std::size_t x = 0; x < occupied[a]; ++x)
        term *= v[i];
      for (std::size_t b = a + 1; b < v.size(); ++b)
      {
        auto const j = permutation[b];
        term *= Real{1} + v[i] * v[j] - Real{2} * d * v[i];
        if (i > j) term = -term;
      }
    }
    term /= normalization;
    sum.value += term;
    sum.absolute_term_sum += std::abs(term);
  }
  while (std::next_permutation(permutation.begin(), permutation.end()));
  return sum;
}

TYPED_TEST(XXZCoordinateWave, SubsetRecurrenceMatchesPermutationDefinition)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real d : {-Real{1} / Real{2}, Real{0}, Real{1}, Real{2}})
    for (unsigned r = 0; r <= 6; ++r)
    {
      std::vector<C> v;
      std::vector<std::size_t> occupied;
      for (unsigned j = 0; j < r; ++j)
      {
        v.emplace_back(Real(j + 2) / Real{7}, Real(int(j % 3) - 1) / Real{3});
        occupied.push_back(2 * j + 1);
      }
      engine::CoordinateBetheWave<Real> const wave(14, d, v);
      EXPECT_EQ(wave.subsets_per_amplitude(), std::size_t{1} << r);
      auto const actual = wave.evaluate(occupied), expected = permutations(d, v, occupied);
      ASSERT_GT(expected.absolute_term_sum, Real{0});
      Real const tolerance = Real{8192} * eps * expected.absolute_term_sum;
      EXPECT_LT(std::abs(actual.value - expected.value), tolerance);
      EXPECT_REAL_NEAR(actual.absolute_term_sum, expected.absolute_term_sum, tolerance);
      if (r >= 2)
      {
        std::swap(v[0], v[1]);
        auto const swapped = engine::CoordinateBetheWave<Real>(14, d, v).evaluate(occupied);
        EXPECT_LT(std::abs(swapped.value + actual.value), tolerance);
        EXPECT_REAL_NEAR(swapped.absolute_term_sum, actual.absolute_term_sum, tolerance);
      }
    }
}

TYPED_TEST(XXZCoordinateWave, VacuumPlaneWavesAndCancellation)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  auto const vacuum = engine::CoordinateBetheWave<Real>(8, Real{0}, {}).evaluate({});
  EXPECT_TRUE(vacuum.value == C{1});
  EXPECT_EQ(vacuum.absolute_term_sum, Real{1});
  std::array<C, 2> const duplicate{C{1}, C{1}};
  std::array<std::size_t, 2> const occupied{0, 3};
  auto const cancelled = engine::CoordinateBetheWave<Real>(8, Real{0}, duplicate).evaluate(occupied);
  EXPECT_TRUE(cancelled.value == C{});
  EXPECT_EQ(cancelled.absolute_term_sum, Real{2});
  auto const zero_pairs = engine::CoordinateBetheWave<Real>(8, Real{1}, duplicate).evaluate(occupied);
  EXPECT_TRUE(zero_pairs.value == C{});
  EXPECT_EQ(zero_pairs.absolute_term_sum, Real{0}); // No division by a zero scattering factor.
  std::array<C, 1> const phase{C{0, 1}};
  auto const n = std::size_t(std::numeric_limits<std::int64_t>::max() / 4);
  std::array<std::size_t, 1> const far{n - 1};
  auto const plane = engine::CoordinateBetheWave<Real>(n, Real{2}, phase).evaluate(far);
  std::array<C, 4> const powers{C{1}, C{0, 1}, C{-1}, C{0, -1}};
  EXPECT_TRUE(plane.value == powers[(n - 1) % 4]);
  EXPECT_EQ(plane.absolute_term_sum, Real{1});
}

TYPED_TEST(XXZCoordinateWave, TwelveParticleFreeFermionFourierDeterminant)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  constexpr unsigned r = 12, n = 3 * r;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  std::vector<C> v;
  std::vector<std::size_t> occupied;
  Real factorial{1}, determinant_magnitude{1};
  for (unsigned j = 0; j < r; ++j)
  {
    Real const angle = pi * Real(2 * int(j) - int(r - 1)) / Real(n);
    v.emplace_back(std::cos(angle), std::sin(angle));
    occupied.push_back(3 * j);
    factorial *= Real(j + 1);
    if (j % 2) determinant_magnitude *= Real(r);
  }
  // Delta=0 makes the pair factor symmetric. All |F_ij|>1 here,
  // so pair normalization leaves only a unit global phase. The wave
  // is a Fourier determinant: |det|=r^(r/2), each permutation has |term|=1.
  auto const result = engine::CoordinateBetheWave<Real>(n, Real{0}, v).evaluate(occupied);
  EXPECT_REAL_NEAR(std::abs(result.value), determinant_magnitude, Real{4096} * eps * determinant_magnitude);
  EXPECT_REAL_NEAR(result.absolute_term_sum, factorial, Real{4096} * eps * factorial);
}

template <typename Real>
Real hamiltonian_residual(unsigned n, unsigned r, Real d, Real energy, std::vector<std::complex<Real>> const& v)
{
  using C = std::complex<Real>;
  engine::CoordinateBetheWave<Real> const wave(n, d, v);
  std::vector<C> psi(1U << n);
  Real amplitude{0}, absolute_terms{0}, error{0};
  for (unsigned bits = 0; bits < (1U << n); ++bits)
    if (std::popcount(bits) == int(r))
    {
      std::vector<std::size_t> occupied;
      for (unsigned j = 0; j < n; ++j)
        if (bits & (1U << j)) occupied.push_back(j);
      auto const result = wave.evaluate(occupied);
      psi[bits] = result.value;
      absolute_terms = std::max(absolute_terms, result.absolute_term_sum);
      amplitude = std::max(amplitude, std::abs(psi[bits]));
    }
  EXPECT_GT(amplitude, Real{64} * uni20::numeric_limits<Real>::epsilon() * absolute_terms);
  for (unsigned bits = 0; bits < (1U << n); ++bits)
    if (std::popcount(bits) == int(r))
    {
      C residual = -energy * psi[bits];
      for (unsigned j = 0; j < n; ++j)
      {
        auto const k = (j + 1) % n;
        bool const anti = ((bits >> j) & 1U) != ((bits >> k) & 1U);
        residual += d * (anti ? -Real{1} : Real{1}) * psi[bits] / Real{4};
        if (anti) residual += psi[bits ^ (1U << j) ^ (1U << k)] / Real{2};
      }
      error = std::max(error, std::abs(residual));
    }
  return amplitude > Real{0} ? error / amplitude : uni20::numeric_limits<Real>::infinity();
}

TYPED_TEST(XXZCoordinateWave, RealRootStatesSatisfyDirectSpinHamiltonian)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (unsigned n : {8, 10})
    for (unsigned r : {3, 4})
      for (Real d : {-Real{1} / Real{2}, Real{0}, Real{1} / Real{2}, Real{1}, Real{2}})
      {
        SCOPED_TRACE(::testing::Message() << "N=" << n << " r=" << r << " Delta=" << uni20::format_scalar(d));
        auto const state = bethe::xxz::sector_ground_state<Real>(n, d, uni20::from_twice(std::int64_t(n - 2 * r)));
        ASSERT_TRUE(state.converged);
        std::vector<C> v;
        for (Real z : state.rapidities)
          v.push_back(-C(Real{1}, -z) / C(Real{1}, z));
        ASSERT_EQ(v.size(), r);
        EXPECT_LT(hamiltonian_residual(n, r, d, state.energy, v), Real{32768} * Real(n) * eps);
        // A finite off-shell input remains evaluable but is NOT an eigenstate.
        v[0] *= C{std::cos(Real{1} / Real{10}), std::sin(Real{1} / Real{10})};
        EXPECT_GT(hamiltonian_residual(n, r, d, state.energy, v), Real{1} / Real{1000});
      }
}

TYPED_TEST(XXZCoordinateWave, RecoveredContinuedPolynomialsGiveEigenvectors)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (unsigned n : {5, 7, 9, 11})
    for (unsigned r = 2; r <= n / 2; ++r)
      for (Real d : {-Real{1} / Real{5}, -Real{7} / Real{10}, -Real{9} / Real{10}, -Real{999} / Real{1000}})
      {
        SCOPED_TRACE(::testing::Message() << "N=" << n << " r=" << r << " Delta=" << uni20::format_scalar(d));
        auto const branch = engine::continue_odd_polynomial(n, d, uni20::from_twice(std::int64_t(n - 2 * r)));
        ASSERT_TRUE(branch.equations_converged);
        auto const result = bethe::detail::recover_polynomial_roots<Real>(branch.coefficients);
        ASSERT_EQ(result.status, bethe::detail::PolynomialRootStatus::resolved)
            << "residual=" << uni20::format_scalar(result.residual_norm)
            << " reconstruction=" << uni20::format_scalar(result.reconstruction_error);
        std::vector<C> v;
        for (C x : result.roots)
        {
          C const z = branch.center + branch.coordinate_scale * x, imaginary{0, 1};
          v.push_back(-(Real{1} - imaginary * z) / (Real{1} + imaginary * z));
        }
        EXPECT_LT(hamiltonian_residual(n, r, d, branch.energy, v), Real{32768} * Real(n) * eps);
      }
}

TYPED_TEST(XXZCoordinateWave, InputVariationEnclosesPerturbedAmplitudes)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), delta = -Real{1} / Real{2};
  for (unsigned r : {2, 3, 4})
  {
    std::vector<C> v;
    std::vector<std::size_t> occupied;
    for (unsigned j = 0; j < r; ++j)
    {
      v.emplace_back(Real(j + 3) / Real{5}, Real(int(j) - 1) / Real{7});
      occupied.push_back(2 * j);
    }
    auto normalization = [&](auto const& factors) {
      Real result{1};
      for (unsigned i = 0; i < r; ++i)
        for (unsigned j = i + 1; j < r; ++j)
          result *= std::max({Real{1}, std::abs(Real{1} + factors[i] * factors[j] - Real{2} * delta * factors[i]),
                              std::abs(Real{1} + factors[i] * factors[j] - Real{2} * delta * factors[j])});
      return result;
    };
    engine::CoordinateBetheWave<Real> const wave(12, delta, v);
    auto const nominal = wave.evaluate(occupied);
    EXPECT_EQ(nominal.input_variation, Real{0});
    for (Real radius : {Real{1} / Real{1000000}, Real{1} / Real{20}})
    {
      std::vector<Real> radii(r, radius);
      auto const bounded = wave.evaluate(occupied, radii);
      EXPECT_TRUE(bounded.value == nominal.value);
      EXPECT_EQ(bounded.absolute_term_sum, nominal.absolute_term_sum);
      EXPECT_GT(bounded.input_variation, Real{0});
      for (unsigned pattern = 0; pattern < (1U << (2 * r)); ++pattern)
      {
        auto perturbed = v;
        for (unsigned j = 0; j < r; ++j)
        {
          std::array<C, 4> const directions{C{1}, C{-1}, C{0, 1}, C{0, -1}};
          perturbed[j] += radius * directions[(pattern >> (2 * j)) & 3U];
        }
        // The uncertainty envelope holds the NOMINAL global pair scaling
        // fixed. Undo the reference's perturbed normalization accordingly.
        auto const reference = permutations(delta, perturbed, occupied);
        C const actual = reference.value * (normalization(perturbed) / normalization(v));
        EXPECT_LE(std::abs(actual - nominal.value),
                  bounded.input_variation +
                      Real{1024} * eps * std::max(nominal.absolute_term_sum, reference.absolute_term_sum));
      }
    }
  }
  // A zero nominal amplitude can have nonzero sensitivity to its inputs.
  std::array<C, 2> const identical{C{1}, C{1}};
  auto const cancelled =
      engine::CoordinateBetheWave<Real>(8, Real{0}, identical)
          .evaluate(std::array<std::size_t, 2>{0, 1}, std::array<Real, 2>{Real{1} / Real{100}, Real{1} / Real{100}});
  EXPECT_TRUE(cancelled.value == C{});
  EXPECT_GT(cancelled.input_variation, Real{1} / Real{100});
  engine::CoordinateBetheWave<Real> const wave(8, delta, std::array<C, 2>{C{1}, C{0, 1}});
  std::array<std::size_t, 2> const occupied{0, 3};
  EXPECT_THROW(wave.evaluate(occupied, std::array<Real, 1>{Real{0}}), std::invalid_argument);
  EXPECT_THROW(wave.evaluate(occupied, std::array<Real, 2>{-Real{1}, Real{0}}), std::invalid_argument);
  EXPECT_THROW(wave.evaluate(occupied, std::array<Real, 2>{Real{0}, uni20::numeric_limits<Real>::infinity()}),
               std::invalid_argument);
  EXPECT_THROW(wave.evaluate(occupied, std::array<Real, 2>{uni20::numeric_limits<Real>::max(), Real{0}}),
               std::overflow_error);
}

TYPED_TEST(XXZCoordinateWave, RejectsInvalidInputsAndExcessWork)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  using Wave = engine::CoordinateBetheWave<Real>;
  std::vector<C> v(4, C{1});
  EXPECT_THROW((Wave(8, Real{0}, v, 15)), std::length_error);
  EXPECT_EQ(Wave(8, Real{0}, v, 16).subsets_per_amplitude(), 16U);
  EXPECT_THROW((Wave(8, Real{0}, {}, 0)), std::length_error);
  v.resize(std::numeric_limits<std::size_t>::digits, C{1});
  EXPECT_THROW((Wave(100, Real{0}, v, std::numeric_limits<std::size_t>::max())), std::length_error);
  EXPECT_THROW((Wave(1, Real{0}, {})), std::invalid_argument);
  EXPECT_THROW((Wave(2, Real{0}, std::array<C, 3>{C{1}, C{1}, C{1}})), std::invalid_argument);
  EXPECT_THROW((Wave(8, uni20::numeric_limits<Real>::infinity(), {})), std::invalid_argument);
  for (auto bad : {C{}, C{uni20::numeric_limits<Real>::quiet_NaN(), 0}, C{0, uni20::numeric_limits<Real>::infinity()}})
    EXPECT_THROW((Wave(8, Real{0}, std::array<C, 1>{bad})), std::invalid_argument);
  Wave const wave(8, Real{0}, std::array<C, 2>{C{1}, C{0, 1}});
  for (auto const& occupied : {std::vector<std::size_t>{0}, {0, 0}, {3, 1}, {0, 8}})
    EXPECT_THROW(wave.evaluate(occupied), std::invalid_argument);
  Real const large = uni20::numeric_limits<Real>::max() / Real{4};
  EXPECT_THROW((Wave(8, Real{0}, std::array<C, 2>{C{large}, C{large}})), std::overflow_error);
  Wave const overflow(8, Real{0}, std::array<C, 1>{C{large}});
  EXPECT_THROW(overflow.evaluate(std::array<std::size_t, 1>{2}), std::overflow_error);
}
} // namespace
