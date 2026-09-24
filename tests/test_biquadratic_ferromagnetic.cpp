// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "tl_ed.hpp"
#include <bethe/biquadratic_ferromagnetic.hpp>

namespace
{
namespace bq = bethe::biquadratic;
namespace ferro = bq::ferromagnetic;
template <typename Real> class BiquadraticFerro : public ::testing::Test {};
TYPED_TEST_SUITE(BiquadraticFerro, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(BiquadraticFerro, AnalyticBandAndNativePrecision)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon(), pi = R{4} * std::atan(R{1});
  for (std::size_t n : {2, 3, 4, 5, 16, 127, 128})
  {
    SCOPED_TRACE(n);
    auto const ground = ferro::ground_space<R>(n);
    EXPECT_EQ(ground.energy, R(n - 1));
    EXPECT_EQ(ground.through_lines, n);
    R previous{};
    for (std::size_t j = n - 1; j > 0; --j)
    {
      auto const s = ferro::one_defect_level<R>(n, j);
      EXPECT_EQ(s.mode, j);
      EXPECT_EQ(s.through_lines, n - 2);
      EXPECT_EQ(s.multiplicity, bethe::temperley_lieb::spin_chain_multiplicity(3, n - 2));
      EXPECT_REAL_NEAR(s.wave_number, pi * R(j) / R(n), R{8} * eps);
      EXPECT_REAL_NEAR(s.gap, R{3} + R{2} * std::cos(pi * R(j) / R(n)), R{32} * eps);
      EXPECT_REAL_NEAR(s.energy, ground.energy + s.gap, R{4} * R(n) * eps);
      EXPECT_GT(s.gap, previous);
      previous = s.gap;
      if (j == n - 1) EXPECT_REAL_NEAR(s.gap, ferro::spectral_gap<R>(n), R{8} * eps);
    }
  }
  auto const four = ferro::one_defect_level<R>(4, 3);
  EXPECT_REAL_NEAR(four.gap, R{3} - std::sqrt(R{2}), R{16} * eps);
  EXPECT_EQ(four.multiplicity, 8);
  EXPECT_EQ(ferro::ground_space<R>(4).multiplicity, 55);
  if constexpr (uni20::numeric_limits<R>::digits > 53) EXPECT_GT(std::abs(R(double(four.gap)) - four.gap), R{16} * eps);
  // Extensive total energy may round away a unit-scale gap; the gap is independent.
  auto const huge =
      ferro::one_defect_level<R>(std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max() - 1);
  EXPECT_GE(huge.gap, R{1});
  EXPECT_LE(huge.gap, R{2});
  EXPECT_FALSE(huge.multiplicity);
}

TYPED_TEST(BiquadraticFerro, RealScansReverseBeforeTruncationAndKeepAuxiliaryConvention)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon();
  for (std::size_t n : {4, 5, 8, 9, 16})
  {
    auto const all = ferro::real_excitations<R>(n, n - 2, {.count = 100});
    auto const one = ferro::real_excitations<R>(n, n - 2, {.count = 1});
    ASSERT_TRUE(all.converged());
    ASSERT_TRUE(one.converged());
    ASSERT_EQ(all.levels.size(), n - 1);
    ASSERT_EQ(one.levels.size(), 1);
    EXPECT_EQ(one.levels[0].state.energy, all.levels[0].state.energy);
    EXPECT_EQ(one.levels[0].state.reference.quantum_numbers.front(), uni20::half_int(std::int64_t(n - 1)));
    EXPECT_EQ(all.ground_state.energy, R(n - 1));
    EXPECT_EQ(all.ground_state.through_lines, n);
    for (std::size_t i = 0; i < n - 1; ++i)
    {
      auto const& level = all.levels[i];
      auto const& state = level.state;
      auto const analytic = ferro::one_defect_level<R>(n, n - 1 - i);
      ASSERT_TRUE(level.gap);
      EXPECT_REAL_NEAR(*level.gap, analytic.gap, R{128} * R(n) * eps);
      EXPECT_REAL_NEAR(state.energy, analytic.energy, R{128} * R(n) * eps);
      EXPECT_EQ(state.exchange, bq::Exchange::ferromagnetic);
      auto const af = bq::solve_real<R>(n, state.reference.quantum_numbers);
      EXPECT_EQ(af.exchange, bq::Exchange::antiferromagnetic);
      EXPECT_EQ(state.energy, -af.energy);
      EXPECT_EQ(state.tl_energy, -af.tl_energy);
      EXPECT_EQ(state.reference.energy, af.reference.energy);
      EXPECT_EQ(state.reference.rapidities, af.reference.rapidities);
    }
  }
}

TYPED_TEST(BiquadraticFerro, ComplexStatesAreEssentialAtFerroModuleEdge)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon();
  auto const real = ferro::real_excitations<R>(4, 0, {.count = 100});
  auto const full = ferro::qsystem::spectrum<R>(4, 0);
  auto const original = bq::qsystem::spectrum<R>(4, 0);
  ASSERT_TRUE(real.converged());
  ASSERT_TRUE(full.complete());
  ASSERT_EQ(real.levels.size(), 1);
  ASSERT_EQ(full.states.size(), 2);
  ASSERT_TRUE(original.complete());
  EXPECT_EQ(full.attempts, original.attempts);
  EXPECT_EQ(full.failed_attempts, original.failed_attempts);
  for (std::size_t i = 0; i < 2; ++i)
  {
    EXPECT_EQ(full.states[i].reference.energy, original.states[1 - i].reference.energy);
    EXPECT_EQ(full.states[i].reference.coefficients, original.states[1 - i].reference.coefficients);
    auto const& a = full.states[i].reference.roots.roots;
    auto const& b = original.states[1 - i].reference.roots.roots;
    ASSERT_EQ(a.size(), b.size());
    for (std::size_t j = 0; j < a.size(); ++j)
    {
      test_support::expect_exact(a[j].real(), b[j].real(), "unchanged root real part");
      test_support::expect_exact(a[j].imag(), b[j].imag(), "unchanged root imaginary part");
    }
    EXPECT_EQ(original.states[i].exchange, bq::Exchange::antiferromagnetic);
  }
  R const low = (R{15} - std::sqrt(R{17})) / R{2}, high = (R{15} + std::sqrt(R{17})) / R{2};
  EXPECT_REAL_NEAR(*full.states[0].energy, low, R{2048} * eps);
  EXPECT_REAL_NEAR(*full.states[1].energy, high, R{2048} * eps);
  EXPECT_REAL_NEAR(real.levels[0].state.energy, high, R{2048} * eps);
  EXPECT_GT(real.levels[0].state.energy, *full.states[0].energy + R{4});
  EXPECT_REAL_NEAR(*full.states[0].tl_energy, low - R{3}, R{2048} * eps);
  EXPECT_GT(std::abs(full.states[0].reference.roots.roots.front().imag()), R{1});
  auto const selected = ferro::qsystem::solve<R>(4, full.states[0].reference.coefficients);
  ASSERT_TRUE(selected.reference.converged);
  EXPECT_EQ(selected.exchange, bq::Exchange::ferromagnetic);
  EXPECT_REAL_NEAR(*selected.energy, low, R{2048} * eps);
}

TYPED_TEST(BiquadraticFerro, InvalidInputsAndPartialResults)
{
  using R = TypeParam;
  for (std::size_t n : {0, 1})
  {
    EXPECT_THROW((void)ferro::ground_space<R>(n), std::invalid_argument);
    EXPECT_THROW((void)ferro::spectral_gap<R>(n), std::invalid_argument);
    EXPECT_THROW((void)ferro::one_defect_level<R>(n, 1), std::invalid_argument);
  }
  for (std::size_t j : {0, 4, 5})
    EXPECT_THROW((void)ferro::one_defect_level<R>(4, j), std::invalid_argument);
  EXPECT_FALSE(ferro::ground_space<R>(46).multiplicity);
  EXPECT_THROW((void)ferro::real_excitations<R>(5, 2), std::invalid_argument);
  EXPECT_THROW((void)ferro::real_excitations<R>(8, 2, {.count = 1, .max_candidates = 9}), std::length_error);
  EXPECT_THROW((void)ferro::real_excitations<R>(8, 2, {.count = 0}), std::invalid_argument);
  auto const failed = ferro::real_excitations<R>(8, 2, {}, {.max_iterations = 0});
  EXPECT_FALSE(failed.converged());
  EXPECT_TRUE(failed.ground_state.reference.converged);
  EXPECT_EQ(failed.ground_state.energy, R{7});
  EXPECT_TRUE(failed.first_unconverged);
  auto const vacuum = ferro::real_excitations<R>(8, 8, {}, {.max_iterations = 0});
  ASSERT_TRUE(vacuum.converged());
  ASSERT_EQ(vacuum.levels.size(), 1);
  EXPECT_EQ(vacuum.levels[0].gap, R{0});
  auto const incomplete = ferro::qsystem::spectrum<R>(6, 2, {.max_attempts = 0});
  EXPECT_FALSE(incomplete.complete());
  EXPECT_TRUE(incomplete.states.empty());
  EXPECT_EQ(incomplete.expected_count, 9);
  auto const pole = ferro::qsystem::solve<R>(4, std::vector<R>{-R{1.5}}, {.max_iterations = 0});
  EXPECT_FALSE(pole.reference.converged);
  EXPECT_FALSE(pole.energy);
  EXPECT_FALSE(pole.tl_energy);
}

TEST(BiquadraticFerroED, GroundMultiplicityAndGlobalPositiveGapOddAndEven)
{
  // Physical spin matrices, independent of the TL spectral correspondence.
  for (unsigned n : {2, 3, 4, 5, 6})
  {
    SCOPED_TRACE(n);
    auto values = bethe::test::biquadratic_ed(n);
    for (auto& e : values)
      e = -e;
    std::sort(values.begin(), values.end());
    auto const ground = ferro::ground_space(n);
    auto first = std::find_if(values.begin(), values.end(), [&](double e) { return e > ground.energy + 1e-8; });
    ASSERT_NE(first, values.end());
    EXPECT_EQ(std::size_t(first - values.begin()), *ground.multiplicity);
    EXPECT_NEAR(*first - ground.energy, ferro::spectral_gap(n), 2e-11);
    for (unsigned j = 1; j < n; ++j)
    {
      auto const s = ferro::one_defect_level(n, j);
      auto count = std::count_if(values.begin(), values.end(), [&](double e) { return std::abs(e - s.energy) < 1e-9; });
      EXPECT_GE(std::uint64_t(count), *s.multiplicity);
    }
  }
}

TEST(BiquadraticFerroED, CompleteSmallModuleSpectraAndPhysicalWeights)
{
  for (unsigned n : {2, 4, 6})
  {
    SCOPED_TRACE(n);
    std::vector<double> reconstructed;
    for (unsigned ell = 0; ell <= n; ell += 2)
    {
      auto const spectrum = ferro::qsystem::spectrum<double>(n, ell);
      ASSERT_TRUE(spectrum.complete()) << ell;
      auto expected = bethe::test::quantum_group_module_ed(n, ell, 1.5);
      for (auto& e : expected)
        e = 1.75 * (n - 1) - 2 * e;
      std::sort(expected.begin(), expected.end());
      ASSERT_EQ(spectrum.states.size(), expected.size());
      for (std::size_t i = 0; i < expected.size(); ++i)
      {
        auto const& s = spectrum.states[i];
        ASSERT_TRUE(s.reference.converged);
        ASSERT_TRUE(s.energy);
        EXPECT_NEAR(*s.energy, expected[i], 3e-10);
        reconstructed.insert(reconstructed.end(), *s.multiplicity, *s.energy);
      }
    }
    auto physical = bethe::test::biquadratic_ed(n);
    for (auto& e : physical)
      e = -e;
    std::sort(physical.begin(), physical.end());
    std::sort(reconstructed.begin(), reconstructed.end());
    ASSERT_EQ(physical.size(), reconstructed.size());
    for (std::size_t i = 0; i < physical.size(); ++i)
      EXPECT_NEAR(physical[i], reconstructed[i], 3e-10);
  }
}
} // namespace
