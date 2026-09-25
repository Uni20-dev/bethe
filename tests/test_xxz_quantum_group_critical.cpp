// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxz_quantum_group_critical.hpp>

namespace
{
namespace model = bethe::xxz::quantum_group::critical;
template <typename Real> class CriticalQuantumGroup : public ::testing::Test {};
TYPED_TEST_SUITE(CriticalQuantumGroup, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(CriticalQuantumGroup, NativeSmallChains)
{
  using Real = TypeParam;
  auto const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real d : {Real{1} / Real{100}, Real{1} / Real{4}, Real{3} / Real{5}, Real{9} / Real{10}, Real{99} / Real{100}})
  {
    auto const two = model::sea_state<Real>(2, d);
    ASSERT_TRUE(two.converged);
    EXPECT_REAL_NEAR(*two.energy, -Real{3} * d / Real{4}, Real{128} * eps);
    auto const four = model::sea_state<Real>(4, d);
    ASSERT_TRUE(four.converged);
    Real const exact = -Real{3} * d / Real{4} - std::sqrt(d * d + Real{2}) / Real{2};
    EXPECT_REAL_NEAR(*four.energy, exact, Real{512} * eps);
    auto const three = model::sea_state<Real>(3, d, 1);
    ASSERT_TRUE(three.converged);
    EXPECT_REAL_NEAR(*three.energy, -(Real{1} + d) / Real{2}, Real{128} * eps);
    test_support::expect_exact(uni20::parse_real<Real>(uni20::format_real(*four.energy)), *four.energy,
                               "critical QG I/O");
  }
}

TYPED_TEST(CriticalQuantumGroup, OriginalComplexBetheEquations)
{
  using Real = TypeParam;
  using C = uni20::complex<Real>;
  auto const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {4, 7, 16, 32})
    for (Real d : {Real{1} / Real{4}, Real{3} / Real{5}, Real{9} / Real{10}})
      for (bool excited : {false, true})
      {
        auto labels = bethe::xxz::quantum_group::detail::consecutive(n / 2 - (excited ? 1 : 0));
        if (excited) labels.back() += uni20::half_int{1};
        auto const s = model::solve_real<Real>(n, d, labels);
        ASSERT_TRUE(s.converged) << n << ' ' << uni20::format_real(d);
        Real const gamma = std::acos(d);
        auto ratio = [&](Real x, Real width) { return std::sinh(C{x, width}) / std::sinh(C{x, -width}); };
        Real energy = Real(n - 1) * d / Real{4};
        for (std::size_t i = 0; i < s.rapidities.size(); ++i)
        {
          auto const lambda = s.rapidities[i];
          ASSERT_TRUE(uni20::isfinite(lambda));
          EXPECT_GT(lambda, Real{0});
          C lhs{1, 0}, rhs{1, 0};
          auto const drive = ratio(lambda, gamma / Real{2});
          for (std::size_t k = 0; k < 2 * n; ++k)
            lhs *= drive;
          for (std::size_t j = 0; j < s.rapidities.size(); ++j)
            if (i != j) rhs *= ratio(lambda - s.rapidities[j], gamma) * ratio(lambda + s.rapidities[j], gamma);
          EXPECT_REAL_NEAR(lhs.real(), rhs.real(), Real{1024} * Real(n) * eps);
          EXPECT_REAL_NEAR(lhs.imag(), rhs.imag(), Real{1024} * Real(n) * eps);
          energy -= (Real{1} - d * d) / (std::cosh(Real{2} * lambda) - d);
        }
        EXPECT_REAL_NEAR(*s.energy, energy, Real{512} * Real(n) * eps);
      }
}

TYPED_TEST(CriticalQuantumGroup, InputsBudgetsAndExplicitLabels)
{
  using Real = TypeParam;
  Real const d = Real{1} / Real{4};
  auto const failed = model::sea_state<Real>(8, d, 0, {.max_iterations = 0});
  EXPECT_FALSE(failed.converged);
  EXPECT_EQ(failed.status, model::Status::iteration_limit);
  EXPECT_FALSE(failed.energy);
  EXPECT_FALSE(failed.energy_shift);
  auto const vacuum = model::sea_state<Real>(8, d, 8, {.max_iterations = 0});
  ASSERT_TRUE(vacuum.converged);
  EXPECT_EQ(*vacuum.energy, Real{7} * d / Real{4});
  std::vector<uni20::half_int> labels{uni20::half_int{2}};
  auto const one = model::solve_real<Real>(8, d, labels);
  ASSERT_TRUE(one.converged);
  Real const pi = Real{4} * std::atan(Real{1});
  EXPECT_REAL_NEAR(*one.energy, Real{7} * d / Real{4} - d - std::cos(pi / Real{4}),
                   Real{128} * uni20::numeric_limits<Real>::epsilon());
  labels[0] = uni20::half_int{7};
  EXPECT_THROW(model::solve_real<Real>(8, d, labels), std::invalid_argument);
  for (Real invalid : {Real{0}, Real{1}, -d, uni20::numeric_limits<Real>::infinity()})
    EXPECT_THROW(model::sea_state<Real>(4, invalid), std::invalid_argument);
  EXPECT_THROW(model::sea_state<Real>(5, d), std::invalid_argument);
  EXPECT_THROW(model::sea_state<Real>(4, d, 1), std::invalid_argument);
  EXPECT_THROW(model::sea_state<Real>(4, d, 0, {.residual_tolerance = Real{0}}), std::invalid_argument);
}

TYPED_TEST(CriticalQuantumGroup, IndependentNonHermitianSpinBasis)
{
  using Real = TypeParam;
  // Direct complex spin-basis eigensolver, not a Hermitian diagonalization
  // or the logarithmic equations. See scripts/reference_xxz_nonhermitian.py.
  struct Reference
  {
      std::size_t n;
      int numerator;
      double energy;
  };
  for (auto const& r : {Reference{5, 25, -1.3016421401931328},
                        {5, 60, -1.582686814421488},
                        {5, 90, -1.8397334893284407},
                        {6, 25, -1.605519747797113},
                        {6, 60, -1.9906369970263755},
                        {6, 90, -2.362192244915388},
                        {8, 25, -2.3007582273337226},
                        {8, 60, -2.7645096966847658},
                        {8, 90, -3.2148081413167406}})
  {
    auto const s = model::sea_state<Real>(r.n, Real(r.numerator) / Real{100}, r.n % 2);
    ASSERT_TRUE(s.converged);
    EXPECT_REAL_NEAR(*s.energy, Real(r.energy), Real{1e-11});
  }
}
TYPED_TEST(CriticalQuantumGroup, ScanWindowAndAnalyticOneParticle)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon(), pi = R{4} * std::atan(R{1});
  // A strict threshold at I=5: the root at infinity is not admitted.
  EXPECT_EQ(model::real_quantum_number_window(8, R{1} / R{2}, 4).slots, 4u);
  EXPECT_EQ(model::real_quantum_number_window(8, R{49} / R{100}, 4).slots, 4u);
  EXPECT_EQ(model::real_quantum_number_window(8, R{51} / R{100}, 4).slots, 5u);
  EXPECT_EQ(model::real_excitation_count(8, R{1} / R{2}, 4), 6u);
  EXPECT_EQ(model::real_excitation_count(8, R{51} / R{100}, 4), 10u);
  for (R d : {R{1} / R{4}, R{3} / R{5}, R{9} / R{10}})
  {
    auto const w = model::real_quantum_number_window(8, d, 6);
    auto const all = model::real_excitations(8, d, 6, {.count = 100}, {.max_iterations = 0});
    ASSERT_TRUE(all.converged());
    ASSERT_EQ(all.levels.size(), w.slots);
    for (std::size_t i = 0; i < w.slots; ++i)
    {
      auto const& level = all.levels[i];
      EXPECT_EQ(level.state.quantum_numbers.front(), uni20::half_int(static_cast<int>(i + 1)));
      EXPECT_REAL_NEAR(*level.state.energy_shift, -d - std::cos(pi * R(i + 1) / R{8}), R{128} * eps);
      EXPECT_REAL_NEAR(*level.gap, std::cos(pi / R{8}) - std::cos(pi * R(i + 1) / R{8}), R{128} * eps);
    }
    auto const lowest = model::real_excitations(8, d, 6, {.count = 2});
    ASSERT_EQ(lowest.levels.size(), 2u);
    EXPECT_EQ(lowest.candidate_count, all.candidate_count);
    EXPECT_EQ(lowest.converged_count, all.converged_count);
    for (std::size_t i = 0; i < 2; ++i)
      EXPECT_EQ(lowest.levels[i].state.quantum_numbers, all.levels[i].state.quantum_numbers);
  }
  auto const vacuum = model::real_excitations(7, R{1} / R{4}, 7, {}, {.max_iterations = 0});
  ASSERT_TRUE(vacuum.converged());
  ASSERT_EQ(vacuum.levels.size(), 1u);
  EXPECT_EQ(vacuum.candidate_count, 1u);
  EXPECT_EQ(*vacuum.levels.front().gap, R{0});
}

TYPED_TEST(CriticalQuantumGroup, ScanAllCandidatesAndFailureAccounting)
{
  using R = TypeParam;
  using C = uni20::complex<R>;
  R const eps = uni20::numeric_limits<R>::epsilon();
  for (std::size_t n : {7, 8, 12})
    for (R d : {R{1} / R{4}, R{3} / R{5}, R{9} / R{10}})
    {
      auto const ell = n - 4; // Two roots; independently enumerate all label pairs.
      auto const w = model::real_quantum_number_window(n, d, ell);
      auto const result = model::real_excitations(n, d, ell, {.count = 1000});
      std::vector<model::RealState<R>> explicit_states;
      std::size_t failures = 0;
      for (std::size_t i = 1; i < w.slots; ++i)
        for (std::size_t j = i + 1; j <= w.slots; ++j)
        {
          std::vector<uni20::half_int> labels{uni20::half_int(static_cast<int>(i)),
                                              uni20::half_int(static_cast<int>(j))};
          auto s = model::solve_real<R>(n, d, labels);
          if (s.converged)
            explicit_states.push_back(std::move(s));
          else
            ++failures;
        }
      EXPECT_EQ(result.candidate_count, w.slots * (w.slots - 1) / 2);
      EXPECT_EQ(result.converged_count, explicit_states.size());
      EXPECT_EQ(result.family_converged(), failures == 0);
      EXPECT_EQ(result.first_unconverged.has_value(), failures != 0);
      std::sort(explicit_states.begin(), explicit_states.end(), [](auto const& a, auto const& b) {
        return a.energy_shift != b.energy_shift ? a.energy_shift < b.energy_shift
                                                : a.quantum_numbers < b.quantum_numbers;
      });
      ASSERT_EQ(result.levels.size(), explicit_states.size());
      auto const gamma = std::acos(d);
      auto ratio = [&](R x, R width) { return std::sinh(C{x, width}) / std::sinh(C{x, -width}); };
      for (std::size_t i = 0; i < result.levels.size(); ++i)
      {
        auto const& level = result.levels[i];
        EXPECT_EQ(level.state.quantum_numbers, explicit_states[i].quantum_numbers);
        ASSERT_TRUE(level.state.energy);
        ASSERT_TRUE(level.gap);
        EXPECT_EQ(*level.gap, *level.state.energy_shift - *result.ground_state.energy_shift);
        // Original multiplicative equations, not the scanner/log residual.
        auto const& roots = level.state.rapidities;
        for (std::size_t r = 0; r < 2; ++r)
        {
          C lhs{1, 0};
          auto const drive = ratio(roots[r], gamma / R{2});
          for (std::size_t power = 0; power < 2 * n; ++power)
            lhs *= drive;
          auto const rhs = ratio(roots[r] - roots[1 - r], gamma) * ratio(roots[r] + roots[1 - r], gamma);
          EXPECT_REAL_NEAR(lhs.real(), rhs.real(), R{2048} * R(n) * eps);
          EXPECT_REAL_NEAR(lhs.imag(), rhs.imag(), R{2048} * R(n) * eps);
        }
      }
    }
  auto const failed = model::real_excitations(8, R{3} / R{5}, 4, {}, {.max_iterations = 0});
  EXPECT_FALSE(failed.converged());
  EXPECT_FALSE(failed.ground_state.energy);
  EXPECT_TRUE(failed.first_unconverged);
  EXPECT_EQ(failed.converged_count, 0u);
  EXPECT_TRUE(failed.levels.empty());
  // A partially converged scan must retain valid levels without disguising
  // the failed candidates. Choose a native iteration budget from a full run
  // rather than hard-coding a double-specific convergence iteration count.
  auto const full = model::real_excitations(12, R{3} / R{5}, 8, {.count = 1000});
  ASSERT_GT(full.levels.size(), 1u);
  auto minimum = full.levels.front().state.iterations;
  auto maximum = minimum;
  for (auto const& level : full.levels)
  {
    minimum = std::min(minimum, level.state.iterations);
    maximum = std::max(maximum, level.state.iterations);
  }
  ASSERT_LT(minimum, maximum);
  auto const partial = model::real_excitations(12, R{3} / R{5}, 8, {.count = 1000}, {.max_iterations = minimum});
  EXPECT_GT(partial.converged_count, 0u);
  EXPECT_LT(partial.converged_count, partial.candidate_count);
  EXPECT_FALSE(partial.family_converged());
  ASSERT_TRUE(partial.first_unconverged);
  EXPECT_FALSE(partial.first_unconverged->energy);
  for (auto const& level : partial.levels)
  {
    EXPECT_TRUE(level.state.converged);
    EXPECT_TRUE(level.state.energy);
    EXPECT_EQ(level.gap.has_value(), partial.ground_state.converged);
  }
}

TYPED_TEST(CriticalQuantumGroup, ScanInputAndPreflight)
{
  using R = TypeParam;
  R const d = R{1} / R{2};
  EXPECT_THROW(model::real_excitations(8, d, 4, {.count = 0}), std::invalid_argument);
  EXPECT_THROW(model::real_excitations(8, d, 4, {.max_candidates = 0}), std::invalid_argument);
  EXPECT_THROW(model::real_excitations(8, d, 4, {.max_candidates = 5}), std::length_error);
  EXPECT_THROW(model::real_excitation_count(8, d, 4, 5), std::length_error);
  EXPECT_THROW(model::real_quantum_number_window(8, R{0}, 4), std::invalid_argument);
  EXPECT_THROW(model::real_quantum_number_window(8, R{1}, 4), std::invalid_argument);
  EXPECT_THROW(model::real_quantum_number_window(7, d, 4), std::invalid_argument);
  EXPECT_THROW(model::real_quantum_number_window(8, d, 10), std::invalid_argument);
  EXPECT_THROW(model::real_excitations(8, d, 4, {}, {.residual_tolerance = R{0}}), std::invalid_argument);
  EXPECT_THROW(model::real_excitation_count(1000000000, d, 500000000, 100), std::length_error);
  // Near Delta=0 the finite precision margin can remove the last sea label.
  auto const tiny = uni20::numeric_limits<R>::epsilon();
  EXPECT_EQ(model::real_excitation_count(4, tiny, 0), 0u);
  EXPECT_THROW(model::real_excitations(4, tiny, 0), std::invalid_argument);
}
} // namespace
