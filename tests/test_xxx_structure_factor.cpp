// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "exact_spectrum.hpp"
#include "test_support.hpp"
#include <bethe/xxx_structure_factor.hpp>
#include <bethe/xxz_coordinate_wave.hpp>

namespace model = bethe::heisenberg;
template <typename T> class XXXStructureFactor : public ::testing::Test {};
TYPED_TEST_SUITE(XXXStructureFactor, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXXStructureFactor, AnalyticDimerAndFourSites)
{
  using Real = TypeParam;
  Real const tol = Real{65536} * uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {2, 4})
  {
    auto const result = model::two_spinon_structure_factor<Real>(n);
    ASSERT_TRUE(result.converged());
    EXPECT_REAL_NEAR(result.weight_fraction, Real{1}, tol);
    EXPECT_REAL_NEAR(*result.first_moment_fraction, Real{1}, tol);
    for (auto const& line : result.lines)
    {
      Real const expected = n == 2 ? Real{1} : line.momentum_index == 2 ? Real{4} / Real{3} : Real{1} / Real{3};
      EXPECT_REAL_NEAR(line.weight, expected, tol);
      EXPECT_REAL_NEAR(result.moments[line.momentum_index].first_moment,
                       model::raising_first_moment(n, line.momentum_index, result.scan.ground_state.energy), tol);
    }
    EXPECT_EQ(result.moments[0].weight, Real{0});
  }
}

TYPED_TEST(XXXStructureFactor, IndependentExactDiagonalization)
{
  using Real = TypeParam;
  using C = std::complex<double>;
  for (unsigned n : {4, 6, 8, 10})
  {
    SCOPED_TRACE(n);
    auto const result = model::two_spinon_structure_factor<Real>(n);
    ASSERT_TRUE(result.converged());
    auto const gs = test_support::exact_eigensystem(n, n / 2);
    auto const ex = test_support::exact_eigensystem(n, n / 2 - 1);
    auto const g = std::min_element(gs.energies.begin(), gs.energies.end()) - gs.energies.begin();
    auto const dim = ex.basis.size();
    std::vector<std::size_t> index(1U << n);
    for (std::size_t i = 0; i < dim; ++i)
      index[ex.basis[i]] = i;
    double total = 0;
    for (unsigned q = 0; q < n; ++q)
    {
      std::vector<C> raised(dim);
      for (std::size_t i = 0; i < gs.basis.size(); ++i)
        for (unsigned j = 0; j < n; ++j)
          if ((gs.basis[i] >> j) & 1U)
            raised[index[gs.basis[i] ^ (1U << j)]] +=
                gs.vectors[i * gs.basis.size() + g] *
                std::polar(1 / std::sqrt(double(n)), -2 * std::acos(-1.0) * q * j / n);
      std::vector<double> weights(dim);
      double first = 0;
      for (std::size_t k = 0; k < dim; ++k)
      {
        C overlap{};
        for (std::size_t i = 0; i < dim; ++i)
          overlap += ex.vectors[i * dim + k] * raised[i];
        weights[k] = std::norm(overlap);
        total += weights[k];
        first += (ex.energies[k] - gs.energies[g]) * weights[k];
      }
      EXPECT_NEAR(first, static_cast<double>(model::raising_first_moment(n, q, result.scan.ground_state.energy)),
                  2e-11);
      // Compare the weight into the full degenerate energy subspace. Arbitrary
      // Jacobi eigenvectors mix +/- momenta; individual ED vector weights do not
      // correspond to individual Bethe states.
      for (auto const& line : result.lines)
        if (line.momentum_index == q)
        {
          double oracle = 0, bethe = 0;
          for (std::size_t k = 0; k < dim; ++k)
            if (std::abs(ex.energies[k] - gs.energies[g] - static_cast<double>(line.gap)) < 1e-8) oracle += weights[k];
          for (auto const& other : result.lines)
            if (other.momentum_index == q && std::abs(static_cast<double>(other.gap - line.gap)) < 1e-8)
              bethe += static_cast<double>(other.weight);
          EXPECT_NEAR(bethe, oracle, 3e-11);
        }
      EXPECT_LE(static_cast<double>(result.moments[q].first_moment), first + 2e-11);
    }
    EXPECT_NEAR(total / n, .5, 3e-12);
    if (n >= 6) EXPECT_LT(result.weight_fraction, Real{1});
  }
}

TYPED_TEST(XXXStructureFactor, MomentumConventionAndNativeReflection)
{
  using Real = TypeParam;
  auto const result = model::two_spinon_structure_factor<Real>(6);
  ASSERT_TRUE(result.converged());
  EXPECT_EQ(result.scan.ground_state.momentum_index, 3U);
  using C = std::complex<Real>;
  Real const tol = Real{65536} * uni20::numeric_limits<Real>::epsilon();
  for (auto const& line : result.lines)
  {
    auto const& state = result.scan.levels[line.state_id - 1].state;
    EXPECT_EQ(line.momentum_index, (3U + 6 - state.momentum_index) % 6);
    std::vector<C> factors;
    for (Real z : state.rapidities)
      factors.push_back(C(z, Real{1}) / C(z, Real{-1}));
    bethe::xxz::detail::CoordinateBetheWave<Real> wave(6, Real{1}, factors);
    std::vector<std::size_t> sites{0, 1}, translated{1, 2};
    C ratio = wave.evaluate(translated).value / wave.evaluate(sites).value;
    C expected = std::polar(Real{1}, state.momentum);
    EXPECT_REAL_NEAR(Real(std::abs(ratio - expected)), Real{0}, tol);
    EXPECT_REAL_NEAR(result.moments[line.momentum_index].weight, result.moments[(6 - line.momentum_index) % 6].weight,
                     tol);
  }
}

TYPED_TEST(XXXStructureFactor, GuardrailsAndFailures)
{
  using Real = TypeParam;
  EXPECT_THROW((void)model::two_spinon_structure_factor<Real>(5), std::invalid_argument);
  EXPECT_THROW((void)model::two_spinon_structure_factor<Real>(64, {}, 527), std::length_error);
  auto failed = model::two_spinon_structure_factor<Real>(8, {.max_iterations = 0});
  EXPECT_FALSE(failed.converged());
  EXPECT_TRUE(failed.lines.empty());
  EXPECT_FALSE(failed.first_moment_fraction);
  auto result = model::two_spinon_structure_factor<Real>(8);
  auto ground = result.scan.ground_state;
  auto excited = result.scan.levels.front().state;
  ground.converged = false;
  EXPECT_EQ(model::raising_form_factor(8, ground, excited).status, model::FormFactorStatus::roots_unconverged);
  ground = result.scan.ground_state;
  ground.rapidities[0] += Real{0.001};
  auto const stale = model::raising_form_factor(8, ground, excited);
  EXPECT_EQ(stale.status, model::FormFactorStatus::precision_limit);
  EXPECT_FALSE(stale.weight);
}

TYPED_TEST(XXXStructureFactor, NativeCoordinateWaveWeights)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  constexpr unsigned n = 8;
  auto const result = model::two_spinon_structure_factor<Real>(n);
  ASSERT_TRUE(result.converged());
  auto wave = [&](auto const& state) {
    std::vector<C> factors;
    for (Real z : state.rapidities)
      factors.push_back(C(z, Real{1}) / C(z, Real{-1}));
    bethe::xxz::detail::CoordinateBetheWave<Real> evaluator(n, Real{1}, factors);
    std::vector<C> vector(1U << n);
    for (unsigned bits = 0; bits < (1U << n); ++bits)
      if (std::popcount(bits) == int(factors.size()))
      {
        std::vector<std::size_t> sites;
        for (unsigned j = 0; j < n; ++j)
          if ((bits >> j) & 1U) sites.push_back(j);
        vector[bits] = evaluator.evaluate(sites).value;
      }
    return vector;
  };
  auto const ground = wave(result.scan.ground_state);
  Real ground_norm{0};
  for (auto x : ground)
    ground_norm += std::norm(x);
  for (auto const& line : result.lines)
  {
    auto const excited = wave(result.scan.levels[line.state_id - 1].state);
    Real excited_norm{0};
    C overlap{};
    for (unsigned bits = 0; bits < (1U << n); ++bits)
    {
      excited_norm += std::norm(excited[bits]);
      if (bits & 1U) overlap += std::conj(excited[bits ^ 1U]) * ground[bits];
    }
    EXPECT_REAL_NEAR(line.weight, Real(n) * std::norm(overlap) / (ground_norm * excited_norm),
                     Real{1048576} * uni20::numeric_limits<Real>::epsilon());
  }
}

TYPED_TEST(XXXStructureFactor, LogDeterminantScalingAndSingularRejection)
{
  using Real = TypeParam;
  using std::log;
  std::vector<Real> matrix{Real{1e150}, Real{0}, Real{0}, Real{1e-150}};
  auto const d = bethe::detail::log_absolute_determinant<Real>(matrix, 2);
  ASSERT_TRUE(d.log_absolute);
  EXPECT_REAL_NEAR(*d.log_absolute, Real(log(matrix[0]) + log(matrix[3])),
                   Real{1024} * uni20::numeric_limits<Real>::epsilon());
  EXPECT_FALSE(bethe::detail::log_absolute_determinant<Real>(std::vector<Real>{1, 2, 2, 4}, 2).log_absolute);
  EXPECT_EQ(*bethe::detail::log_absolute_determinant<Real>(std::vector<Real>{}, 0).log_absolute, Real{0});
}
