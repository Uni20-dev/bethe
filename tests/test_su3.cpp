// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/heisenberg.hpp>
#include <bethe/su3.hpp>
#include <uni20/linalg/ops/self_adjoint_eigh.hpp>

namespace
{
namespace model = bethe::su3;
template <typename Real> class SU3 : public ::testing::Test {};
TYPED_TEST_SUITE(SU3, test_support::RealTypes, test_support::PrecisionNames);

template <uni20::Real Real> void check_state(model::State<Real> const& state)
{
  using std::abs;
  using std::atan;
  Real const pi = Real{4} * atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  std::array<Real, 2> residuals{};
  Real energy = Real(state.sites), phase = Real{0};
  EXPECT_EQ(state.converged, state.status == model::SolveStatus::converged);
  EXPECT_EQ(state.populations[0] + state.populations[1] + state.populations[2], state.sites);
  for (std::size_t a = 0; a < 2; ++a)
  {
    auto const& roots = state.rapidities[a];
    ASSERT_EQ(roots.size(), state.quantum_numbers[a].size());
    ASSERT_EQ(roots.size(), (2 - a) * (state.sites / 3));
    for (std::size_t j = 0; j < roots.size(); ++j)
    {
      SCOPED_TRACE(::testing::Message() << a << "," << j);
      Real const lambda = roots[j];
      if (j) EXPECT_GT(lambda, roots[j - 1]);
      EXPECT_EQ(lambda, -roots[roots.size() - 1 - j]);
      Real f = -pi * Real(state.quantum_numbers[a][j].twice());
      if (a == 0)
      {
        f += Real(state.sites) * Real{2} * atan(Real{2} * lambda);
        energy -= Real{1} / (lambda * lambda + Real{1} / Real{4});
        phase += pi - Real{2} * atan(Real{2} * lambda);
      }
      for (std::size_t k = 0; k < roots.size(); ++k)
        if (k != j) f -= Real{2} * atan(lambda - roots[k]);
      for (Real other : state.rapidities[1 - a])
        f += Real{2} * atan(Real{2} * (lambda - other));
      residuals[a] = std::max(residuals[a], abs(f) / Real(state.sites));
    }
    EXPECT_REAL_NEAR(residuals[a], state.level_residuals[a], Real{64} * eps);
    if (state.converged) EXPECT_LT(residuals[a], Real{128} * eps);
  }
  EXPECT_EQ(state.residual_norm, std::max(state.level_residuals[0], state.level_residuals[1]));
  EXPECT_REAL_NEAR(energy, state.energy, Real{32} * Real(state.sites) * eps);
  EXPECT_EQ(state.momentum_index, 0u);
  EXPECT_EQ(state.momentum, Real{0});
  EXPECT_REAL_NEAR(phase, Real(state.rapidities[0].size()) * pi, Real{32} * Real(state.sites) * eps);
}

TYPED_TEST(SU3, GroundCountsAndDiagnostics)
{
  using Real = TypeParam;
  for (std::size_t n : {3, 6, 9, 12, 24, 48, 96, 192})
  {
    SCOPED_TRACE(n);
    auto const state = model::ground_state<Real>(n);
    ASSERT_TRUE(state.converged) << state.iterations;
    EXPECT_EQ(state.quantum_numbers, model::ground_quantum_numbers(n));
    for (auto count : state.populations)
      EXPECT_EQ(count, n / 3);
    EXPECT_LE(state.residual_norm, Real{32} * uni20::numeric_limits<Real>::epsilon());
    ASSERT_NO_FATAL_FAILURE(check_state(state));
  }
}

TYPED_TEST(SU3, AnalyticThreeSites)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  auto const state = model::ground_state<Real>(3);
  ASSERT_TRUE(state.converged);
  Real const root = Real{1} / std::sqrt(Real{12});
  EXPECT_REAL_NEAR(state.energy, -Real{3}, Real{256} * eps);
  EXPECT_REAL_NEAR(state.rapidities[0][1], root, Real{64} * eps);
  EXPECT_EQ(state.rapidities[1][0], Real{0});
  if constexpr (uni20::numeric_limits<Real>::digits > 53)
    EXPECT_GT(std::abs(Real(double(root)) - root), Real{64} * eps);
  test_support::expect_exact(uni20::parse_real<Real>(uni20::format_real(state.energy)), state.energy, "energy I/O");
}

TYPED_TEST(SU3, OriginalMultiplicativeEquations)
{
  using Real = TypeParam;
  using Complex = uni20::complex<Real>;
  auto const e = [](Real x, Real width) { return Complex{x, width / Real{2}} / Complex{x, -width / Real{2}}; };
  for (std::size_t n : {3, 6, 9, 12})
  {
    auto const state = model::ground_state<Real>(n);
    ASSERT_TRUE(state.converged);
    Real const tolerance = Real{256} * Real(n) * uni20::numeric_limits<Real>::epsilon();
    for (std::size_t a = 0; a < 2; ++a)
      for (Real root : state.rapidities[a])
      {
        // Original Eq. (2.17), including the SAME-level self factor e_2(0)=-1.
        // This also checks logarithmic parity/sign conventions, without atan.
        Complex product{-Real{1}, Real{0}};
        if (a == 0)
          for (std::size_t site = 0; site < n; ++site)
            product *= e(root, -Real{1});
        for (Real other : state.rapidities[a])
          product *= e(root - other, Real{2});
        for (Real other : state.rapidities[1 - a])
          product *= e(root - other, -Real{1});
        EXPECT_REAL_NEAR(product.real(), Real{1}, tolerance);
        EXPECT_REAL_NEAR(product.imag(), Real{0}, tolerance);
      }
  }
}

TYPED_TEST(SU3, SU2Reduction)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {2, 4, 6, 12, 32, 64})
  {
    SCOPED_TRACE(n);
    // Delete the auxiliary sea in the SAME equation/Jacobian implementation.
    model::detail::GroundSystem<Real> const system(n, {n / 2, 0});
    auto const state = model::detail::solve(system, model::SolverOptions<Real>{});
    auto const xxx = bethe::heisenberg::ground_state<Real>(n);
    ASSERT_TRUE(state.converged && xxx.converged);
    EXPECT_REAL_NEAR(state.energy, Real{2} * xxx.energy + Real(n) / Real{2}, Real{256} * Real(n) * eps);
    for (std::size_t j = 0; j < n / 2; ++j)
      EXPECT_REAL_NEAR(Real{2} * state.rapidities[0][j], xxx.rapidities[j], Real{1024} * eps);
  }
}

TYPED_TEST(SU3, AnalyticSixSites)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  // In the five-dimensional singlet space, formed from pairs of three-color
  // Levi-Civita tensors, the characteristic polynomial is
  // E*(E+2)^2*(E^2+2E-12). See the guide for an explicit reduced matrix.
  Real const exact = -Real{1} - std::sqrt(Real{13}), tolerance = Real{256} * eps;
  auto const state = model::ground_state<Real>(6);
  ASSERT_TRUE(state.converged);
  EXPECT_REAL_NEAR(state.energy, exact, tolerance);
  if constexpr (uni20::numeric_limits<Real>::digits > 53) EXPECT_GT(std::abs(Real(double(exact)) - exact), tolerance);
}

struct ExactGround
{
    double energy, translation, gap;
};
ExactGround exact_ground(unsigned n, bool balanced = true)
{
  // Independent color-word Hamiltonian: swap adjacent trits, including the
  // closing bond. This oracle intentionally uses fp64, not the solver's roots.
  std::vector<unsigned> powers(n + 1, 1);
  for (unsigned j = 0; j < n; ++j)
    powers[j + 1] = 3 * powers[j];
  std::vector<unsigned> basis;
  std::vector<std::size_t> index(powers[n]);
  for (unsigned word = 0; word < powers[n]; ++word)
  {
    std::array<unsigned, 3> counts{};
    for (unsigned j = 0; j < n; ++j)
      ++counts[(word / powers[j]) % 3];
    if (!balanced || (counts[0] == n / 3 && counts[1] == n / 3 && counts[2] == n / 3))
    {
      index[word] = basis.size();
      basis.push_back(word);
    }
  }
  uni20::DenseMatrix<double> h(basis.size(), basis.size());
  for (std::size_t i = 0; i < basis.size(); ++i)
    for (std::size_t j = 0; j < basis.size(); ++j)
      h[i, j] = 0;
  for (std::size_t col = 0; col < basis.size(); ++col)
    for (unsigned j = 0; j < n; ++j)
    {
      unsigned const k = (j + 1) % n, a = (basis[col] / powers[j]) % 3, b = (basis[col] / powers[k]) % 3;
      unsigned const swapped = basis[col] - a * powers[j] - b * powers[k] + b * powers[j] + a * powers[k];
      h[index[swapped], col] += 1;
    }
  auto const eig = uni20::linalg::eigh(std::move(h));
  double translation = 0;
  for (std::size_t j = 0; j < basis.size(); ++j)
  {
    unsigned const translated = (3 * basis[j]) % powers[n] + basis[j] / powers[n - 1];
    translation += eig.eigenvectors[index[translated], 0] * eig.eigenvectors[j, 0];
  }
  return {eig.eigenvalues[0], translation, eig.eigenvalues[1] - eig.eigenvalues[0]};
}

TYPED_TEST(SU3, EnergyAndMomentumAgainstED)
{
  using Real = TypeParam;
  for (unsigned n : {3, 6, 9})
  {
    SCOPED_TRACE(n);
    auto const exact = exact_ground(n);
    auto const state = model::ground_state<Real>(n);
    ASSERT_TRUE(state.converged);
    EXPECT_REAL_NEAR(state.energy, Real(exact.energy), Real{1} / Real{10000000000});
    ASSERT_GT(exact.gap, 1e-6); // unique state, so its translation expectation is an eigenvalue
    EXPECT_NEAR(exact.translation, 1.0, 1e-12);
    EXPECT_EQ(state.momentum_index, 0u);
  }
}

TEST(SU3Oracles, FullHilbertSpaceAndULSIdentity)
{
  EXPECT_NEAR(exact_ground(3, false).energy, -3.0, 1e-13);
  EXPECT_NEAR(exact_ground(6, false).energy, -1.0 - std::sqrt(13.0), 1e-12);
  // S1.S2 in the spin-1 Sz basis. Each allowed flip-flop matrix element is
  // (sqrt(2)*sqrt(2))/2=1; its square is a true matrix product, not elementwise.
  double dot[9][9]{};
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
    {
      int const col = 3 * a + b;
      dot[col][col] = (a - 1) * (b - 1);
      if (a < 2 && b > 0) dot[3 * (a + 1) + b - 1][col] = 1;
      if (a > 0 && b < 2) dot[3 * (a - 1) + b + 1][col] = 1;
    }
  for (int row = 0; row < 9; ++row)
    for (int col = 0; col < 9; ++col)
    {
      double value = dot[row][col] - double(row == col);
      for (int k = 0; k < 9; ++k)
        value += dot[row][k] * dot[k][col];
      EXPECT_EQ(value, double(row == 3 * (col % 3) + col / 3));
    }
}

TYPED_TEST(SU3, Jacobian)
{
  using Real = TypeParam;
  Real const h = std::cbrt(uni20::numeric_limits<Real>::epsilon());
  for (std::size_t n : {3, 6, 9, 12})
  {
    model::detail::GroundSystem<Real> const system(n, {2 * (n / 3), n / 3});
    auto const x = system.seed();
    uni20::DenseMatrix<Real> jacobian(system.order, system.order);
    (void)system.evaluate(x, &jacobian);
    for (std::size_t col = 0; col < system.order; ++col)
    {
      auto plus = x, minus = x;
      plus[col] += h;
      minus[col] -= h;
      auto const fp = system.evaluate(plus), fm = system.evaluate(minus);
      for (std::size_t row = 0; row < system.order; ++row)
      {
        SCOPED_TRACE(::testing::Message() << n << "," << row << "," << col);
        Real const numerical = (fp.residual[row] - fm.residual[row]) / (Real{2} * h);
        EXPECT_REAL_NEAR((jacobian[row, col]), numerical, Real{1000} * h * h * (Real{1} + std::abs(numerical)));
      }
    }
  }
}

TYPED_TEST(SU3, BudgetsAndInvalidInputs)
{
  using Real = TypeParam;
  for (std::size_t budget : {0, 1})
  {
    auto const state = model::ground_state<Real>(12, {.max_iterations = budget});
    EXPECT_FALSE(state.converged);
    EXPECT_EQ(state.iterations, budget);
    EXPECT_EQ(state.status, model::SolveStatus::iteration_limit);
    ASSERT_NO_FATAL_FAILURE(check_state(state));
  }
  for (std::size_t n : {0, 1, 2, 4, 5, 7, 8})
  {
    EXPECT_THROW((model::ground_state<Real>(n)), std::invalid_argument);
    EXPECT_THROW((model::ground_quantum_numbers(n)), std::invalid_argument);
  }
  EXPECT_THROW((model::ground_state<Real>(std::numeric_limits<std::size_t>::max())), std::invalid_argument);
  EXPECT_THROW((model::ground_state<Real>(6000000000ULL)), std::length_error);
  for (Real bad :
       {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW((model::ground_state<Real>(3, {.residual_tolerance = bad})), std::invalid_argument);
}

TYPED_TEST(SU3, BulkLimit)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1});
  // From the digamma expression in Doikou/Nepomechie (2.49), scaled to sum P.
  Real const bulk = Real{1} - std::log(Real{3}) - pi / (Real{3} * std::sqrt(Real{3}));
  Real previous = Real{1};
  for (std::size_t n : {12, 24, 48, 96, 192})
  {
    auto const state = model::ground_state<Real>(n);
    ASSERT_TRUE(state.converged);
    Real const error = std::abs(state.energy / Real(n) - bulk);
    EXPECT_LT(error, previous / Real{3});
    previous = error;
  }
  EXPECT_LT(previous, Real{1} / Real{10000});
}
} // namespace
