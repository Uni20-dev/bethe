// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "tl_ed.hpp"
#include <bethe/biquadratic.hpp>

namespace
{
namespace tl = bethe::temperley_lieb;
namespace qg = bethe::xxz::quantum_group;
template <typename Real> class TemperleyLieb : public ::testing::Test {};
TYPED_TEST_SUITE(TemperleyLieb, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(TemperleyLieb, AnalyticTwoAndFourSites)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real loop : {Real{2} + Real{256} * eps, Real{3}, Real{4}, Real{10}, Real{10000}})
  {
    SCOPED_TRACE(uni20::format_real(loop));
    auto const two = tl::open_ground_state(2, loop);
    ASSERT_TRUE(two.reference.converged);
    EXPECT_REAL_NEAR(two.energy, -loop, Real{256} * loop * eps);
    auto const four = tl::open_ground_state(4, loop);
    ASSERT_TRUE(four.reference.converged);
    Real const exact = -(Real{3} * loop + std::sqrt(loop * loop + Real{8})) / Real{2};
    EXPECT_REAL_NEAR(four.energy, exact, Real{1024} * loop * eps);
  }
  auto const two = bethe::biquadratic::ground_state<Real>(2);
  EXPECT_REAL_NEAR(two.energy, -Real{4}, Real{64} * eps);
  auto const four = bethe::biquadratic::ground_state<Real>(4);
  ASSERT_TRUE(four.reference.converged);
  Real const exact = -(Real{15} + std::sqrt(Real{17})) / Real{2};
  EXPECT_REAL_NEAR(four.energy, exact, Real{256} * eps);
  if constexpr (uni20::numeric_limits<Real>::digits > 53)
    EXPECT_GT(std::abs(Real(double(exact)) - exact), Real{256} * eps);
  test_support::expect_exact(uni20::parse_real<Real>(uni20::format_real(four.energy)), four.energy, "biquadratic I/O");
}

TYPED_TEST(TemperleyLieb, NativePrecisionOriginalBetheEquations)
{
  using Real = TypeParam;
  using C = uni20::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  for (std::size_t n : {4, 8, 32, 128})
    for (Real delta : {Real{101} / Real{100}, Real{3} / Real{2}, Real{4}})
      for (std::size_t ell : {0, 2, 4})
      {
        SCOPED_TRACE(::testing::Message() << "N=" << n << " Delta=" << uni20::format_real(delta) << " ell=" << ell);
        auto const m = (n - ell) / 2;
        auto labels = qg::detail::consecutive(m);
        // Ground state at ell=0; nontrivial high-label excited configurations otherwise.
        for (auto& label : labels)
          label += uni20::half_int(std::int64_t(ell));
        auto const s = qg::solve_real<Real>(n, delta, labels);
        ASSERT_TRUE(s.converged) << int(s.status) << " " << uni20::format_real(s.residual_norm);
        EXPECT_LE(s.residual_norm, Real{32} * eps);
        ASSERT_EQ(s.rapidities.size(), m);
        Real const eta = std::acosh(delta);
        auto ratio = [&](Real alpha, Real width) {
          return std::sinh(C{width, alpha / Real{2}}) / std::sinh(C{width, -alpha / Real{2}});
        };
        Real energy = Real(n - 1) * delta / Real{4};
        for (std::size_t i = 0; i < m; ++i)
        {
          Real const a = s.rapidities[i];
          EXPECT_GT(a, Real{0});
          EXPECT_LT(a, pi);
          if (i) EXPECT_GT(a, s.rapidities[i - 1]);
          EXPECT_EQ(s.quantum_numbers[i], labels[i]);
          C lhs{Real{1}, Real{0}}, rhs{Real{1}, Real{0}};
          auto const drive = ratio(a, eta / Real{2});
          for (std::size_t k = 0; k < 2 * n; ++k)
            lhs *= drive;
          for (std::size_t j = 0; j < m; ++j)
            if (i != j) rhs *= ratio(a - s.rapidities[j], eta) * ratio(a + s.rapidities[j], eta);
          EXPECT_REAL_NEAR(lhs.real(), rhs.real(), Real{512} * Real(n) * eps);
          EXPECT_REAL_NEAR(lhs.imag(), rhs.imag(), Real{512} * Real(n) * eps);
          energy -= (delta * delta - Real{1}) / (delta - std::cos(a));
        }
        EXPECT_REAL_NEAR(energy, s.energy, Real{4096} * Real(n) * eps);
      }
}

TYPED_TEST(TemperleyLieb, AnalyticJacobian)
{
  using Real = TypeParam;
  Real const h = std::cbrt(uni20::numeric_limits<Real>::epsilon());
  for (Real delta : {Real{101} / Real{100}, Real{3} / Real{2}, Real{10}})
  {
    qg::detail::GroundSystem<Real> const system(8, delta);
    for (auto const& x : {system.seed(), qg::ground_state(8, delta).angles})
    {
      std::vector<Real> jac;
      (void)system.evaluate(x, &jac);
      for (std::size_t j = 0; j < 4; ++j)
      {
        auto p = x, m = x;
        p[j] += h;
        m[j] -= h;
        auto const plus = system.evaluate(p), minus = system.evaluate(m);
        for (std::size_t i = 0; i < 4; ++i)
          EXPECT_REAL_NEAR((plus.residual[i] - minus.residual[i]) / (Real{2} * h), jac[i * 4 + j], Real{4096} * h * h);
      }
    }
  }
}

TYPED_TEST(TemperleyLieb, InvalidInputsAndBudgetDiagnostics)
{
  using Real = TypeParam;
  for (auto n : {0u, 1u, 3u, 5u})
    EXPECT_THROW((void)bethe::biquadratic::ground_state<Real>(n), std::invalid_argument);
  for (Real value :
       {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    EXPECT_THROW((void)bethe::biquadratic::ground_state<Real>(4, {.residual_tolerance = value}), std::invalid_argument);
    EXPECT_THROW((void)tl::open_ground_state(4, value), std::invalid_argument);
  }
  EXPECT_THROW((void)tl::open_ground_state(4, Real{2}), std::invalid_argument);
  EXPECT_THROW((void)qg::ground_state(4, Real{1}), std::invalid_argument);
  EXPECT_THROW((void)bethe::biquadratic::ground_state<Real>(std::numeric_limits<std::size_t>::max()),
               std::invalid_argument);
  EXPECT_THROW((void)bethe::biquadratic::ground_state<Real>(1000000000000000000ULL), std::length_error);
  for (std::size_t budget : {0, 1, 2})
  {
    auto const s = bethe::biquadratic::ground_state<Real>(8, {.max_iterations = budget});
    auto const& r = s.reference;
    EXPECT_FALSE(r.converged);
    EXPECT_EQ(r.status, qg::SolveStatus::iteration_limit);
    EXPECT_EQ(r.iterations, budget);
    qg::detail::GroundSystem<Real> const system(8, Real{3} / Real{2});
    EXPECT_EQ(r.residual_norm, system.evaluate(r.angles).norm);
    test_support::expect_exact(s.energy, Real{2} * r.energy_shift - Real{7},
                               "direct TL contribution at failed iterate");
    EXPECT_REAL_NEAR(s.energy, Real{2} * r.energy - Real{49} / Real{4},
                     Real{128} * uni20::numeric_limits<Real>::epsilon());
  }
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  auto const strict =
      bethe::biquadratic::ground_state<Real>(8, {.residual_tolerance = eps * eps, .max_iterations = 30});
  EXPECT_FALSE(strict.reference.converged);
  EXPECT_NE(strict.reference.status, qg::SolveStatus::converged);
}

TEST(TemperleyLiebRepresentation, ExactMultiplicitiesAndOverflow)
{
  std::vector<std::uint64_t> const expected{1, 3, 8, 21, 55, 144, 377};
  for (std::size_t ell = 0; ell < expected.size(); ++ell)
  {
    EXPECT_EQ(tl::spin_chain_multiplicity(3, ell), expected[ell]);
    EXPECT_EQ(tl::spin_chain_multiplicity(2, ell), ell + 1);
  }
  EXPECT_EQ(tl::spin_chain_multiplicity(3, 45), 7540113804746346429ULL);
  EXPECT_FALSE(tl::spin_chain_multiplicity(3, 46));
  EXPECT_EQ(tl::spin_chain_multiplicity(4294967296ULL, 2), std::numeric_limits<std::uint64_t>::max());
  EXPECT_FALSE(tl::spin_chain_multiplicity(4294967297ULL, 2));
  EXPECT_THROW((void)tl::spin_chain_multiplicity(1, 0), std::invalid_argument);
}

TEST(TemperleyLiebED, FullSpinOneSpectrumAndRepresentationWeights)
{
  for (unsigned n : {2, 4, 6})
  {
    SCOPED_TRACE(n);
    auto const physical = bethe::test::biquadratic_ed(n);
    std::vector<double> reconstructed, previous;
    for (unsigned m = 0; m <= n / 2; ++m)
    {
      auto const spectrum = bethe::test::quantum_group_xxz_ed(n, m, 1.5);
      auto module = spectrum;
      // The Sz=N/2-m space contains modules ell=N,N-2,...,N-2m once.
      // Subtract the previous magnetization space to isolate the new module.
      for (double e : previous)
      {
        auto found = std::find_if(module.begin(), module.end(), [&](double f) { return std::abs(e - f) < 1e-10; });
        ASSERT_NE(found, module.end());
        module.erase(found);
      }
      auto const multiplicity = tl::spin_chain_multiplicity(3, n - 2 * m);
      ASSERT_TRUE(multiplicity);
      for (double e : module)
        for (std::uint64_t k = 0; k < *multiplicity; ++k)
          reconstructed.push_back(2 * e - 1.75 * (n - 1));
      previous = spectrum;
    }
    std::sort(reconstructed.begin(), reconstructed.end());
    ASSERT_EQ(reconstructed.size(), physical.size());
    for (std::size_t j = 0; j < physical.size(); ++j)
      EXPECT_NEAR(reconstructed[j], physical[j], 3e-11);
    auto const state = bethe::biquadratic::ground_state(n);
    ASSERT_TRUE(state.reference.converged);
    EXPECT_NEAR(state.energy, physical.front(), 3e-11);
    EXPECT_GT(physical[1] - physical[0], 1e-3); // Unique finite even-chain ground state.
    auto const wrong = bethe::test::quantum_group_xxz_ed(n, n / 2, 1.5, false);
    EXPECT_GT(std::abs(2 * wrong.front() - 1.75 * (n - 1) - physical.front()), 0.01);
  }
}

TEST(TemperleyLiebRepresentation, PhysicalSpinMultipletsAndOverflow)
{
  std::vector<std::vector<std::uint64_t>> const expected{{1}, {0, 1}, {0, 1, 1}, {1, 1, 2, 1}, {1, 3, 3, 3, 1}};
  for (std::size_t ell = 0; ell < expected.size(); ++ell)
    EXPECT_EQ(tl::spin_one_multiplets(ell), expected[ell]);
  for (std::size_t ell = 0; ell <= 45; ++ell)
  {
    auto const counts = tl::spin_one_multiplets(ell);
    ASSERT_TRUE(counts);
    ASSERT_EQ(counts->size(), ell + 1);
    EXPECT_EQ(counts->back(), 1); // Unique fully polarized multiplet.
    std::uint64_t dimension = 0;
    for (std::size_t spin = 0; spin <= ell; ++spin)
    {
      ASSERT_LE((*counts)[spin], std::numeric_limits<std::uint64_t>::max() / (2 * spin + 1));
      auto const states = (2 * spin + 1) * (*counts)[spin];
      ASSERT_LE(states, std::numeric_limits<std::uint64_t>::max() - dimension);
      dimension += states;
    }
    EXPECT_EQ(dimension, tl::spin_chain_multiplicity(3, ell));
  }
  EXPECT_FALSE(tl::spin_one_multiplets(46));
  EXPECT_FALSE(tl::spin_one_multiplets(100000));
  EXPECT_FALSE(tl::spin_one_multiplets(std::numeric_limits<std::size_t>::max()));
}

TEST(TemperleyLiebED, EveryPhysicalSpinSpectrum)
{
  // Resolve irreps independently by subtracting the physical Sz=S+1 spectrum
  // from Sz=S. Compare full energy multisets, not just total dimensions.
  // This catches assigning auxiliary spin ell/2 or physical spin ell alone.
  for (unsigned n : {2, 3, 4, 5, 6})
  {
    std::vector<std::vector<double>> modules(n / 2 + 1);
    for (unsigned m = 0; m <= n / 2; ++m)
      modules[m] = bethe::test::quantum_group_module_ed(n, n - 2 * m, 1.5);
    std::vector<double> higher;
    for (int spin = int(n); spin >= 0; --spin)
    {
      SCOPED_TRACE(::testing::Message() << "N=" << n << " S=" << spin);
      auto const sz = bethe::test::biquadratic_ed(n, spin);
      auto physical = sz;
      for (double e : higher)
      {
        auto found = std::find_if(physical.begin(), physical.end(), [e](double f) { return std::abs(e - f) < 1e-10; });
        ASSERT_NE(found, physical.end());
        physical.erase(found);
      }
      higher = sz;
      std::vector<double> reconstructed;
      for (unsigned m = 0; m <= n / 2; ++m)
      {
        auto const counts = tl::spin_one_multiplets(n - 2 * m);
        ASSERT_TRUE(counts);
        if (std::size_t(spin) >= counts->size()) continue;
        for (double e : modules[m])
          for (std::uint64_t k = 0; k < (*counts)[spin]; ++k)
            reconstructed.push_back(2 * e - 1.75 * (n - 1));
      }
      std::sort(reconstructed.begin(), reconstructed.end());
      ASSERT_EQ(reconstructed.size(), physical.size());
      for (std::size_t i = 0; i < physical.size(); ++i)
        EXPECT_NEAR(reconstructed[i], physical[i], 3e-11);
      // Reversing the Hamiltonian negates these same spin-resolved levels;
      // representation content cannot depend on the energy sign.
    }
  }
}

TEST(TemperleyLiebED, GenericLoopWeightGroundStates)
{
  for (unsigned n : {2, 4, 6, 8})
    for (double delta : {1.01, 1.5, 3.0, 10.0})
    {
      auto const s = qg::ground_state(n, delta);
      ASSERT_TRUE(s.converged);
      EXPECT_NEAR(s.energy, bethe::test::quantum_group_xxz_ed(n, n / 2, delta).front(), 3e-11);
    }
}

TYPED_TEST(TemperleyLieb, ExcitationNativePrecisionAndWeights)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  auto const scan = bethe::biquadratic::real_excitations<Real>(4, 2, {.count = 100});
  ASSERT_TRUE(scan.converged());
  ASSERT_EQ(scan.candidate_count, 3);
  ASSERT_EQ(scan.levels.size(), 3);
  EXPECT_EQ(scan.ground_state.through_lines, 0);
  EXPECT_EQ(scan.ground_state.multiplicity, 1);
  Real const root2 = std::sqrt(Real{2});
  std::array<Real, 3> const energies{-Real{6} - root2, -Real{6}, -Real{6} + root2};
  Real const e0 = -(Real{15} + std::sqrt(Real{17})) / Real{2};
  for (std::size_t j = 0; j < 3; ++j)
  {
    auto const& level = scan.levels[j];
    EXPECT_EQ(level.state.through_lines, 2);
    EXPECT_EQ(level.state.multiplicity, 8);
    EXPECT_REAL_NEAR(level.state.energy, energies[j], Real{128} * eps);
    ASSERT_TRUE(level.gap);
    EXPECT_REAL_NEAR(*level.gap, energies[j] - e0, Real{256} * eps);
    EXPECT_EQ(level.state.reference.quantum_numbers, qg::QuantumNumbers({uni20::half_int(std::int64_t(j + 1))}));
  }
  if constexpr (uni20::numeric_limits<Real>::digits > 53)
    EXPECT_GT(std::abs(Real(double(energies[0])) - energies[0]), Real{128} * eps);
  auto const vacuum = bethe::biquadratic::sector_ground_state<Real>(4, 4);
  EXPECT_TRUE(vacuum.reference.converged);
  EXPECT_EQ(vacuum.reference.iterations, 0);
  EXPECT_EQ(vacuum.energy, -Real{3});
  EXPECT_EQ(vacuum.multiplicity, 55);
  auto const overflow = bethe::biquadratic::sector_ground_state<Real>(46, 46);
  EXPECT_TRUE(overflow.reference.converged);
  EXPECT_FALSE(overflow.multiplicity);
  EXPECT_EQ(overflow.energy, -Real{45});
}

TYPED_TEST(TemperleyLieb, ExcitationScanContracts)
{
  using Real = TypeParam;
  namespace bq = bethe::biquadratic;
  EXPECT_EQ(bq::real_excitation_count(8, 2), 10);
  EXPECT_EQ(bq::real_excitation_count(8, 0), 1);
  EXPECT_EQ(bq::real_excitation_count(8, 8), 1);
  EXPECT_THROW((void)bq::real_excitations<Real>(8, 2, {.count = 0}), std::invalid_argument);
  EXPECT_THROW((void)bq::real_excitations<Real>(8, 2, {.max_candidates = 9}), std::length_error);
  EXPECT_THROW((void)bq::real_excitations<Real>(1000000000000000000ULL, 2), std::length_error);
  EXPECT_THROW((void)bq::real_excitation_count(128, 64), std::length_error);
  for (std::size_t ell : {1, 3, 10})
  {
    EXPECT_THROW((void)bq::real_excitations<Real>(8, ell), std::invalid_argument);
    EXPECT_THROW((void)bq::sector_ground_state<Real>(8, ell), std::invalid_argument);
  }
  for (auto const& raw : std::vector<std::vector<int>>{{0}, {1, 1}, {2, 1}, {1, 4}, {1, 2, 3}, {4}})
  {
    qg::QuantumNumbers labels;
    for (int i : raw)
      labels.emplace_back(i);
    EXPECT_THROW((void)bq::solve_real<Real>(4, labels), std::invalid_argument);
  }
  EXPECT_THROW((void)bq::solve_real<Real>(4, qg::QuantumNumbers{uni20::from_twice(std::int64_t{1})}),
               std::invalid_argument);
  auto const all = bq::real_excitations<Real>(8, 2, {.count = 100});
  auto const two = bq::real_excitations<Real>(8, 2, {.count = 2});
  ASSERT_TRUE(all.converged());
  ASSERT_TRUE(two.converged());
  ASSERT_EQ(all.levels.size(), 10);
  ASSERT_EQ(two.levels.size(), 2);
  EXPECT_EQ(two.candidate_count, 10);
  for (std::size_t j = 0; j < 2; ++j)
  {
    EXPECT_EQ(two.levels[j].state.energy, all.levels[j].state.energy);
    EXPECT_EQ(two.levels[j].state.reference.quantum_numbers, all.levels[j].state.reference.quantum_numbers);
  }
  auto const failed = bq::real_excitations<Real>(8, 2, {}, {.max_iterations = 0});
  EXPECT_FALSE(failed.converged());
  EXPECT_FALSE(failed.family_converged());
  EXPECT_EQ(failed.converged_count, 0);
  EXPECT_TRUE(failed.levels.empty());
  ASSERT_TRUE(failed.first_unconverged);
  EXPECT_EQ(failed.first_unconverged->multiplicity, 8);
  auto const no_ground = bq::real_excitations<Real>(8, 8, {}, {.max_iterations = 0});
  EXPECT_FALSE(no_ground.converged());
  EXPECT_TRUE(no_ground.family_converged());
  ASSERT_EQ(no_ground.levels.size(), 1);
  EXPECT_FALSE(no_ground.levels[0].gap);
}

TEST(TemperleyLiebED, AllRealFamiliesBelongToTheirModules)
{
  namespace bq = bethe::biquadratic;
  for (unsigned n : {2, 4, 6, 8, 10})
    for (unsigned ell = 0; ell <= n; ell += 2)
    {
      SCOPED_TRACE(::testing::Message() << "N=" << n << " ell=" << ell);
      auto const module = bethe::test::quantum_group_module_ed(n, ell, 1.5);
      auto available = module;
      auto const scan = bq::real_excitations(n, ell, {.count = 10000});
      ASSERT_TRUE(scan.converged());
      EXPECT_EQ(scan.levels.size(), bq::real_excitation_count(n, ell));
      auto const minimum = bq::sector_ground_state(n, ell);
      ASSERT_TRUE(minimum.reference.converged);
      EXPECT_NEAR(minimum.energy, 2 * module.front() - 1.75 * (n - 1), 4e-11);
      EXPECT_EQ(minimum.energy, scan.levels.front().state.energy);
      for (auto const& level : scan.levels)
      {
        auto const found = std::find_if(available.begin(), available.end(), [&](double e) {
          return std::abs(2 * e - 1.75 * (n - 1) - level.state.energy) < 4e-11;
        });
        ASSERT_NE(found, available.end()); // Also excludes duplicate/descendant solutions.
        available.erase(found);
        EXPECT_EQ(level.state.through_lines, ell);
        EXPECT_EQ(level.state.multiplicity, tl::spin_chain_multiplicity(3, ell));
      }
      if (ell == n || ell == n - 2) EXPECT_TRUE(available.empty()); // Complete zero-/one-root modules.
      if (n == 4 && ell == 0) EXPECT_EQ(available.size(), 1);       // Explicitly missing complex-root singlet.
    }
  auto const physical = bethe::test::biquadratic_ed(4);
  auto const scan = bq::real_excitations(4, 2);
  for (auto const& level : scan.levels)
    EXPECT_EQ(std::count_if(physical.begin(), physical.end(),
                            [&](double e) { return std::abs(e - level.state.energy) < 1e-10; }),
              8);
  EXPECT_NEAR(scan.levels.front().state.energy, physical[1], 1e-11); // First excited energy, eightfold.
}
} // namespace
