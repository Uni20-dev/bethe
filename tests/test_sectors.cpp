// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include "exact_spectrum.hpp"
#include "test_support.hpp"
#include <bethe/heisenberg.hpp>

#include <bit>
#include <string>
#include <string_view>

namespace
{

using test_support::exact_spectrum;

template <uni20::Real Real> void check_state(std::size_t n, bethe::heisenberg::RealState<Real> const& state)
{
  using std::abs;
  using std::atan;
  using std::cos;
  using std::sin;
  Real const pi = Real{4} * atan(Real{1});
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real residual = Real{0};
  Real energy = static_cast<Real>(n) / Real{4};
  Real momentum = Real{0};
  ASSERT_TRUE(state.quantum_numbers.size() == state.rapidities.size()) << "root/quantum-number count";
  for (std::size_t i = 0; i < state.rapidities.size(); ++i)
  {
    SCOPED_TRACE(::testing::Message() << " i=" << i);
    auto const z = state.rapidities[i];
    Real f = Real{2} * static_cast<Real>(n) * atan(z) - pi * static_cast<Real>(state.quantum_numbers[i].twice());
    for (std::size_t j = 0; j < state.rapidities.size(); ++j)
      if (i != j) f -= Real{2} * atan((z - state.rapidities[j]) / Real{2});
    residual = std::max(residual, abs(f) / static_cast<Real>(n));
    energy -= Real{2} / (Real{1} + z * z);
    momentum += pi - Real{2} * atan(z);
    if (i > 0 && state.converged)
    {
      ASSERT_TRUE(state.rapidities[i - 1] < z) << "ordered converged roots";
    }
  }
  Real const allowance = Real{64} * static_cast<Real>(n) * eps;
  EXPECT_REAL_NEAR(energy, state.energy, allowance) << "returned energy matches roots";
  EXPECT_REAL_NEAR(residual, state.residual_norm, allowance) << "returned residual matches roots";
  ASSERT_TRUE(state.momentum_index < n) << "momentum index in range";
  if (state.converged)
  {
    // Summing individual momentum errors costs O(N) times the phase residual.
    Real const momentum_allowance =
        Real{4} * static_cast<Real>(n * n) * eps + static_cast<Real>(n) * state.residual_norm;
    EXPECT_REAL_NEAR(cos(momentum), cos(state.momentum), momentum_allowance) << "cos(P) from rapidities";
    EXPECT_REAL_NEAR(sin(momentum), sin(state.momentum), momentum_allowance) << "sin(P) from rapidities";
  }
}

template <typename Real> class Sectors : public ::testing::Test {};
TYPED_TEST_SUITE(Sectors, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(Sectors, EnergiesAgainstED)
{
  using Real = TypeParam;

  using namespace bethe::heisenberg;
  auto half = [](std::int64_t twice) { return uni20::from_twice(twice); };
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const tolerance = Real{4096} * eps;

  for (unsigned n = 2; n <= 9; ++n)
  {
    SCOPED_TRACE(::testing::Message() << " n=" << n);
    auto const sectors = sector_ground_states<Real>(n);
    ASSERT_TRUE(sectors.size() == n + 1) << "all Sz sectors present";
    for (unsigned m = 0; m <= n / 2; ++m)
    {
      SCOPED_TRACE(::testing::Message() << " m=" << m);
      auto const ed = exact_spectrum(n, m);
      auto const& state = sectors[n - m];

      ASSERT_TRUE(state.converged) << "sector minimum disagrees with exact diagonalization";
      ASSERT_REAL_NEAR(state.energy, static_cast<Real>(ed.front()), Real{1e-11})
          << "sector minimum disagrees with exact diagonalization";

      ASSERT_TRUE(state.sz == half(n - 2 * m) && state.rapidities.size() == m) << "sector labels";
      auto const& reversed = sectors[m];
      ASSERT_TRUE(reversed.energy == state.energy && reversed.sz == -state.sz) << "spin-reversed sectors";
      ASSERT_TRUE(reversed.momentum_index == state.momentum_index) << "spin reversal preserves momentum";
      ASSERT_TRUE(reversed.spin_reversed == (reversed.sz.twice() < 0)) << "spin-reversed reference label";
      ASSERT_NO_FATAL_FAILURE(check_state(n, state));
      auto const direct = sector_ground_state<Real>(n, -state.sz);
      ASSERT_TRUE(direct.energy == state.energy) << "direct negative Sz solve";
      if (n % 2 != 0 && m != 0)
      {
        auto reflected_numbers = state.quantum_numbers;
        std::reverse(reflected_numbers.begin(), reflected_numbers.end());
        for (auto& number : reflected_numbers)
          number = -number;
        auto const reflected = solve_real<Real>(n, reflected_numbers);

        ASSERT_TRUE(reflected.converged) << "odd-N reflection-related sector minima";
        ASSERT_REAL_NEAR(reflected.energy, state.energy, tolerance) << "odd-N reflection-related sector minima";

        ASSERT_TRUE((reflected.momentum_index + state.momentum_index) % n == 0) << "reflected momentum";
      }
    }
    auto const gs = ground_state<Real>(n);
    ASSERT_TRUE(gs.converged && gs.energy == sectors[(n + 1) / 2].energy) << "ground-state wrapper";
  }
}

TYPED_TEST(Sectors, OneMagnonPrecision)
{
  using Real = TypeParam;

  using namespace bethe::heisenberg;
  using std::abs;
  using std::sqrt;
  auto half = [](std::int64_t twice) { return uni20::from_twice(twice); };
  Real const eps = uni20::numeric_limits<Real>::epsilon();

  // Precision-sensitive one-magnon sector oracle: E=N/4-1-cos(pi/N), odd N.
  auto const magnon = sector_ground_state<Real>(5, half(3));
  Real const exact_magnon = -sqrt(Real{5}) / Real{4};
  EXPECT_REAL_NEAR(magnon.energy, exact_magnon, Real{64} * eps) << "one-magnon irrational energy precision";
  if constexpr (uni20::numeric_limits<Real>::digits > uni20::numeric_limits<double>::digits)
  {
    ASSERT_TRUE(abs(static_cast<Real>(static_cast<double>(exact_magnon)) - exact_magnon) > Real{64} * eps)
        << "one-magnon oracle discriminates double narrowing";
  }
}

TYPED_TEST(Sectors, SpecifiedStatesAgainstED)
{
  using Real = TypeParam;

  using namespace bethe::heisenberg;
  using std::abs;
  using std::cos;
  auto half = [](std::int64_t twice) { return uni20::from_twice(twice); };

  // Exercise general specified states, not only the sector-ground or spinon
  // generators: every supported two-root configuration for these small chains.
  for (unsigned n : {6, 7, 8})
  {
    SCOPED_TRACE(::testing::Message() << " n=" << n);
    auto const ed = exact_spectrum(n, 2, 0.371);
    auto const bound = static_cast<std::int64_t>(n) - 3;
    for (std::int64_t first = -bound; first <= bound; first += 2)
      for (std::int64_t second = first + 2; second <= bound; second += 2)
      {
        SCOPED_TRACE(::testing::Message() << " first=" << first);
        SCOPED_TRACE(::testing::Message() << " second=" << second);
        QuantumNumbers const numbers{half(first), half(second)};
        auto const state = solve_real<Real>(n, numbers);
        ASSERT_TRUE(state.converged && state.quantum_numbers == numbers) << "specified two-root state";
        ASSERT_NO_FATAL_FAILURE(check_state(n, state));
        double const target = static_cast<double>(state.energy) + 0.371 * std::cos(static_cast<double>(state.momentum));
        ASSERT_TRUE(std::any_of(ed.begin(), ed.end(), [&](double e) { return std::abs(e - target) < 1e-11; }))
            << "specified real-root energy/momentum disagrees with ED";
      }
  }
}

TYPED_TEST(Sectors, SpinonBranchAgainstED)
{
  using Real = TypeParam;

  using namespace bethe::heisenberg;
  using std::abs;
  using std::atan;
  using std::cos;
  using std::sin;
  auto half = [](std::int64_t twice) { return uni20::from_twice(twice); };
  Real const pi = Real{4} * atan(Real{1});
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const tolerance = Real{4096} * eps;

  for (unsigned n : {3, 5, 7, 9})
  {
    SCOPED_TRACE(::testing::Message() << " n=" << n);
    auto const branch = one_spinon_branch<Real>(n);
    auto const ed = exact_spectrum(n, (n - 1) / 2, 0.371);
    ASSERT_TRUE(branch.size() == (n + 1) / 2) << "one-spinon state count";
    std::vector<std::size_t> momenta;
    for (std::size_t i = 0; i < branch.size(); ++i)
    {
      SCOPED_TRACE(::testing::Message() << " i=" << i);
      auto const& point = branch[i];
      auto const& state = point.state;
      ASSERT_TRUE(state.converged && state.sz == half(1)) << "one-spinon convergence and Sz";
      ASSERT_TRUE(point.hole.twice() == static_cast<std::int64_t>((n - 1) / 2) - 2 * static_cast<std::int64_t>(i))
          << "hole enumeration";
      ASSERT_TRUE(std::find(state.quantum_numbers.begin(), state.quantum_numbers.end(), point.hole) ==
                  state.quantum_numbers.end())
          << "specified hole is absent";
      ASSERT_TRUE(point.spinon_momentum > Real{0} && point.spinon_momentum < pi) << "spinon k range";
      if (i > 0)
      {
        ASSERT_TRUE(branch[i - 1].spinon_momentum < point.spinon_momentum) << "ascending spinon k";
      }
      auto const& mirror = branch[branch.size() - 1 - i];
      EXPECT_REAL_NEAR(point.spinon_momentum + mirror.spinon_momentum, pi, tolerance) << "spinon k reflection";
      EXPECT_REAL_NEAR(state.energy, mirror.state.energy, tolerance) << "spinon energy reflection";
      EXPECT_REAL_NEAR(point.bulk_subtracted_energy, (state.energy - Real(n) * bulk_energy_density<Real>()), tolerance)
          << "bulk energy subtraction";
      Real const p = pi * static_cast<Real>((n - 1) / 2) + pi / Real{2} - point.spinon_momentum;

      ASSERT_REAL_NEAR(cos(p), cos(state.momentum), tolerance) << "spinon/lattice momentum offset";
      ASSERT_REAL_NEAR(sin(p), sin(state.momentum), tolerance) << "spinon/lattice momentum offset";

      double const target = static_cast<double>(state.energy) + 0.371 * std::cos(static_cast<double>(state.momentum));
      ASSERT_TRUE(std::any_of(ed.begin(), ed.end(), [&](double e) { return std::abs(e - target) < 1e-11; }))
          << "one-spinon energy/momentum disagrees with ED translation oracle";
      momenta.push_back(state.momentum_index);
      ASSERT_NO_FATAL_FAILURE(check_state(n, state));
    }
    std::sort(momenta.begin(), momenta.end());
    ASSERT_TRUE(std::adjacent_find(momenta.begin(), momenta.end()) == momenta.end()) << "unique lattice momenta";
    EXPECT_REAL_NEAR(branch.front().state.energy, ground_state<Real>(n).energy, tolerance) << "branch endpoint is GS";
  }
}

TYPED_TEST(Sectors, SpinonFiniteSize)
{
  using Real = TypeParam;

  using namespace bethe::heisenberg;
  using std::abs;
  using std::atan;
  auto half = [](std::int64_t twice) { return uni20::from_twice(twice); };
  Real const pi = Real{4} * atan(Real{1});
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const tolerance = Real{4096} * eps;

  // N=5 central-hole state has roots +/-1, E=-3/4.
  auto const center = one_spinon_state<Real>(5, half(0));
  ASSERT_EQ(center.state.rapidities.size(), 2u);
  EXPECT_REAL_NEAR(center.state.energy + Real{3} / Real{4}, Real{0}, tolerance) << "N=5 central spinon exact energy";
  EXPECT_REAL_NEAR(center.state.rapidities[1], Real{1}, tolerance) << "N=5 central spinon exact root";
  ASSERT_TRUE(center.state.momentum_index == 0) << "N=5 central spinon momentum";
  // At fixed k=pi/2, bulk-subtracted finite-size energies approach pi/2.
  auto const c17 = one_spinon_state<Real>(17, half(0));
  auto const c65 = one_spinon_state<Real>(65, half(0));
  ASSERT_TRUE(c17.state.converged && c65.state.converged) << "larger one-spinon convergence";
  EXPECT_REAL_NEAR(c65.bulk_subtracted_energy, pi / Real{2}, abs(c17.bulk_subtracted_energy - pi / Real{2}))
      << "one-spinon thermodynamic convergence";
  EXPECT_REAL_NEAR(c65.bulk_subtracted_energy, pi / Real{2}, Real{0.02}) << "one-spinon thermodynamic normalization";
}

TYPED_TEST(Sectors, BudgetAndRestart)
{
  using Real = TypeParam;

  using namespace bethe::heisenberg;
  using Options = SolverOptions<Real>;
  auto half = [](std::int64_t twice) { return uni20::from_twice(twice); };

  auto const center = one_spinon_state<Real>(5, half(0));
  auto const no_updates = one_spinon_state<Real>(7, half(1), Options{.max_iterations = 0});
  ASSERT_TRUE(!no_updates.state.converged && no_updates.state.iterations == 0) << "spinon zero budget";
  ASSERT_NO_FATAL_FAILURE(check_state(7, no_updates.state));
  auto const one_update = sector_ground_state<Real>(7, half(1), Options{.max_iterations = 1});
  ASSERT_TRUE(!one_update.converged && one_update.iterations == 1) << "sector one-update budget";
  ASSERT_NO_FATAL_FAILURE(check_state(7, one_update));
  auto const vacuum = solve_real<Real>(8, QuantumNumbers{}, Options{.max_iterations = 0});
  ASSERT_TRUE(vacuum.converged && vacuum.energy == Real{2} && vacuum.momentum_index == 0) << "empty-root vacuum";
  auto const restarted = solve_real<Real>(5, center.state.quantum_numbers, {}, center.state.rapidities);
  ASSERT_TRUE(restarted.converged && restarted.iterations == 0 && restarted.energy == center.state.energy)
      << "initial roots remain separate from state identity";
}

TYPED_TEST(Sectors, AnalyticDispersion)
{
  using Real = TypeParam;

  using namespace bethe::heisenberg;
  using std::atan;
  using std::sin;
  using std::sqrt;
  Real const pi = Real{4} * atan(Real{1});
  Real const eps = uni20::numeric_limits<Real>::epsilon();

  ASSERT_TRUE(spinon_energy(Real{0}) == Real{0} && spinon_energy(pi) == Real{0}) << "exact dispersion endpoints";
  EXPECT_REAL_NEAR(spinon_energy(pi / Real{2}), pi / Real{2}, Real{8} * eps) << "XXX dispersion maximum";
  EXPECT_REAL_NEAR(spinon_energy(pi / Real{6}), pi / Real{4}, Real{8} * eps) << "XXX dispersion interior";
  ASSERT_TRUE(spinon_energy(pi / Real{3}, Real{2}) == Real{2} * spinon_energy(pi / Real{3})) << "exchange scaling";
  for (Real const k : {Real{0}, pi / Real{7}, pi / Real{2}, pi})
  {
    SCOPED_TRACE(::testing::Message() << " k=" << uni20::format_scalar(k));
    ASSERT_TRUE(bethe::xxz::spinon_energy(k, Real{1}) == spinon_energy(k)) << "XXZ isotropic limit";
    Real const xx = k == pi ? Real{0} : sin(k);
    EXPECT_REAL_NEAR(bethe::xxz::spinon_energy(k, Real{0}), xx, Real{8} * eps) << "XX free-fermion dispersion";
    EXPECT_REAL_NEAR(bethe::xxz::spinon_energy(k, Real{1} - Real{8} * eps), spinon_energy(k), Real{32} * eps)
        << "XXZ stable near-isotropic limit";
    EXPECT_REAL_NEAR(bethe::xxz::spinon_energy(k, Real{1} / Real{2}), Real{3} * sqrt(Real{3}) / Real{4} * xx,
                     Real{16} * eps)
        << "XXZ Delta=1/2 analytic coefficient";
    EXPECT_REAL_NEAR(bethe::xxz::spinon_energy(k, -Real{1} / Real{2}), Real{3} * sqrt(Real{3}) / Real{8} * xx,
                     Real{16} * eps)
        << "XXZ Delta=-1/2 analytic coefficient";
  }
}

TYPED_TEST(Sectors, InvalidInputs)
{
  using Real = TypeParam;

  using namespace bethe::heisenberg;
  using std::atan;
  auto half = [](std::int64_t twice) { return uni20::from_twice(twice); };
  Real const pi = Real{4} * atan(Real{1});

  auto const center = one_spinon_state<Real>(5, half(0));
  for (auto const sz : {half(1), half(10), half(-10), half(std::numeric_limits<std::int64_t>::min())})
    EXPECT_THROW(([&] { (void)sector_ground_state<Real>(4, sz); })(), std::invalid_argument);
  EXPECT_THROW(([&] { (void)sector_ground_state<Real>(5, half(0)); })(), std::invalid_argument);
  for (QuantumNumbers const& numbers :
       {QuantumNumbers{half(-1), half(-1)}, QuantumNumbers{half(1), half(-1)}, QuantumNumbers{half(-2), half(2)},
        QuantumNumbers{half(-3), half(3)}, QuantumNumbers{half(-2), half(0), half(2)}})
    EXPECT_THROW(([&] { (void)solve_real<Real>(4, numbers); })(), std::invalid_argument);
  EXPECT_THROW(([&] { (void)one_spinon_branch<Real>(4); })(), std::invalid_argument);
  EXPECT_THROW(([&] { (void)one_spinon_state<Real>(5, half(1)); })(), std::invalid_argument);
  EXPECT_THROW(([&] { (void)one_spinon_state<Real>(5, half(4)); })(), std::invalid_argument);
  EXPECT_THROW(([&] { (void)solve_real<Real>(5, center.state.quantum_numbers, {}, std::vector<Real>{Real{0}}); })(),
               std::invalid_argument);
  EXPECT_THROW(([&] {
                 (void)solve_real<Real>(5, center.state.quantum_numbers, {},
                                        std::vector<Real>{Real{0}, uni20::numeric_limits<Real>::infinity()});
               })(),
               std::invalid_argument);
  for (Real const bad : {-Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    SCOPED_TRACE(::testing::Message() << " bad=" << uni20::format_scalar(bad));
    EXPECT_THROW(([&] { (void)spinon_energy(bad); })(), std::invalid_argument);
    EXPECT_THROW(([&] { (void)spinon_energy(Real{1}, bad); })(), std::invalid_argument);
    EXPECT_THROW(([&] { (void)bethe::xxz::spinon_energy(Real{1}, bad); })(), std::invalid_argument);
  }
  EXPECT_THROW(([&] { (void)spinon_energy(pi + Real{1}); })(), std::invalid_argument);
  EXPECT_THROW(([&] { (void)spinon_energy(Real{1}, Real{0}); })(), std::invalid_argument);
  EXPECT_THROW(([&] { (void)bethe::xxz::spinon_energy(Real{1}, Real{2}); })(), std::invalid_argument);
}

} // namespace
