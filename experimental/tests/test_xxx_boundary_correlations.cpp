// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "exact_spectrum.hpp"
#include "test_support.hpp"
#include "xxx_open_wave_oracle.hpp"
#include <bethe/detail/newton.hpp>
#include <bethe/experimental/xxx_boundary_correlations.hpp>
#include <bethe/heisenberg_excitations.hpp>

namespace model = bethe::heisenberg::open;
namespace boundary = bethe::experimental::xxx::open;
template <typename T> class ExperimentalBoundaryCorrelations : public ::testing::Test {};
TYPED_TEST_SUITE(ExperimentalBoundaryCorrelations, test_support::RealTypes, test_support::PrecisionNames);

namespace
{
unsigned reverse_bits(unsigned bits, unsigned n)
{
  unsigned result = 0;
  for (unsigned i = 0; i < n; ++i)
    result = (result << 1) | ((bits >> i) & 1U);
  return result;
}

template <uni20::Real Real>
std::array<uni20::complex<Real>, 16> trace_first_two(unsigned n, std::span<uni20::complex<Real> const> wave)
{
  auto const zero = bethe::detail::precision_zero(wave.front().real());
  std::array<uni20::complex<Real>, 16> rho;
  rho.fill(uni20::complex<Real>(zero, zero));
  for (unsigned rest = 0; rest < (1U << (n - 2)); ++rest)
    for (unsigned i = 0; i < 4; ++i)
      for (unsigned j = 0; j < 4; ++j)
      {
        auto const ket_i = (rest << 2) | (i >> 1) | ((i & 1U) << 1);
        auto const ket_j = (rest << 2) | (j >> 1) | ((j & 1U) << 1);
        rho[i * 4 + j] += wave[ket_i] * uni20::conj(wave[ket_j]);
      }
  return rho;
}

// Test-only perturbation of the driving term. Newton's Jacobian affects the
// iteration, not the target equations. No production correlation formula is
// used to evaluate the independently reconstructed wave/Hamiltonian below.
template <uni20::Real Real> std::vector<Real> impurity_roots(unsigned n, Real theta)
{
  using uni20::math::abs;
  using uni20::math::atan;
  Real const zero = bethe::detail::precision_zero(theta);
  auto options = [&] {
#if UNI20_HAS_MPREAL
    if constexpr (uni20::has_runtime_precision_v<Real>)
      return bethe::SolverOptions<Real>(zero.precision());
    else
#endif
      return bethe::SolverOptions<Real>{};
  }();
  options.residual_tolerance = 4 * uni20::numeric_limits<Real>::epsilon(zero);
  auto state = model::ground_state<Real>(n, options);
  std::vector<Real> u;
  for (Real z : state.rapidities)
    u.push_back(z / Real{2});
  auto const m = u.size();
  for (unsigned iteration = 0; iteration < 30; ++iteration)
  {
    std::vector<Real> residual(m, zero);
    Real largest = zero;
    auto matrix = boundary::boundary_detail::gaudin_matrix<Real>(n, u, zero);
    for (std::size_t i = 0; i < m; ++i)
    {
      Real const v = u[i];
      bethe::detail::CompensatedSum<Real> value(zero);
      value.add(Real{4} * Real(n - 1) * atan(Real{2} * v));
      value.add(Real{2} * atan(Real{2} * (v - theta)) + Real{2} * atan(Real{2} * (v + theta)));
      value.add(-Real{2} * bethe::detail::pi_like(zero) * Real(i + 1));
      for (std::size_t j = 0; j < m; ++j)
        if (i != j) value.add(-Real{2} * (atan(v - u[j]) + atan(v + u[j])));
      residual[i] = -value.value();
      largest = std::max(largest, abs(residual[i]));
      matrix[i * m + i] += Real{4} / (Real{1} + Real{4} * (v - theta) * (v - theta)) +
                           Real{4} / (Real{1} + Real{4} * (v + theta) * (v + theta)) -
                           Real{8} / (Real{1} + Real{4} * v * v);
    }
    if (largest <= Real{64} * Real(n) * uni20::numeric_limits<Real>::epsilon(zero)) return u;
    if (!bethe::detail::newton_step(std::move(matrix), residual, options.precision()))
      throw std::runtime_error("test impurity solve failed");
    for (std::size_t i = 0; i < m; ++i)
      u[i] += residual[i];
  }
  throw std::runtime_error("test impurity roots did not converge");
}

template <uni20::Real Real> Real impurity_energy(unsigned n, Real theta, std::span<Real const> u)
{
  auto const zero = bethe::detail::precision_zero(theta);
  bethe::detail::CompensatedSum<Real> energy(zero);
  energy.add((Real(n - 2) + Real{1} / (Real{1} + theta * theta)) / Real{4});
  for (Real v : u)
    energy.add(-Real{0.5} / (v * v + Real{0.25}));
  return energy.value();
}
} // namespace

TYPED_TEST(ExperimentalBoundaryCorrelations, ReflectionParityAgainstNativeWaves)
{
  using Real = TypeParam;
  using C = uni20::complex<Real>;
  Real const tolerance = Real{524288} * uni20::numeric_limits<Real>::epsilon();
  for (unsigned n : {2, 3, 4, 5, 6, 8})
    for (unsigned m = 0; m <= n / 2; ++m)
    {
      auto const scan =
          model::real_excitations<Real>(n, uni20::from_twice(static_cast<std::int64_t>(n - 2 * m)), {.count = 1000});
      ASSERT_TRUE(scan.converged());
      for (auto const& level : scan.levels)
      {
        int const parity = boundary::real_root_reflection_parity(n, level.state.quantum_numbers);
        auto const wave = test_support::open_aba_wave<Real>(n, level.state.rapidities);
        Real largest{};
        for (unsigned bits = 0; bits < (1U << n); ++bits)
          largest = std::max(largest, Real(uni20::math::abs(wave[reverse_bits(bits, n)] - Real(parity) * wave[bits])));
        EXPECT_REAL_NEAR(largest, Real{0}, tolerance);
      }
    }
}

TYPED_TEST(ExperimentalBoundaryCorrelations, AnalyticBoundaryDensity)
{
  using Real = TypeParam;
  using uni20::math::sqrt;
  Real const tolerance = Real{65536} * uni20::numeric_limits<Real>::epsilon();
  for (unsigned n : {2, 4})
  {
    auto const result = boundary::boundary_two_site_density(n, model::ground_state<Real>(n));
    ASSERT_TRUE(result.converged());
    ASSERT_TRUE(result.solve_info && result.solve_info->succeeded());
    ASSERT_TRUE(result.matrix());
    Real const expected = n == 2 ? Real{-0.75} : -(Real{1} + sqrt(Real{3})) / Real{4};
    EXPECT_REAL_NEAR(*result.bond_correlation, expected, tolerance);
    EXPECT_REAL_NEAR(*result.zz(), expected / Real{3}, tolerance);
    EXPECT_REAL_NEAR(*result.singlet_probability() + Real{3} * *result.triplet_component_probability(), Real{1},
                     tolerance);
    EXPECT_LT(*result.solve_backward_error, tolerance);
  }
}

TYPED_TEST(ExperimentalBoundaryCorrelations, CompleteDensityAgainstABAAndED)
{
  using Real = TypeParam;
  using C = uni20::complex<Real>;
  Real const tolerance = Real{524288} * uni20::numeric_limits<Real>::epsilon();
  for (unsigned n : {2, 4, 6, 8, 10})
  {
    SCOPED_TRACE(n);
    auto const ground = model::ground_state<Real>(n);
    auto const result = boundary::boundary_two_site_density(n, ground);
    ASSERT_TRUE(result.converged());
    auto const rho = *result.matrix();
    auto const native = test_support::open_aba_wave<Real>(n, ground.rapidities);
    auto const native_rho = trace_first_two<Real>(n, native);
    for (std::size_t i = 0; i < 16; ++i)
      EXPECT_REAL_NEAR(Real(uni20::math::abs(native_rho[i] - rho[i])), Real{0}, tolerance);

    auto const ed = test_support::exact_eigensystem(n, n / 2, 0, false);
    auto const g = std::min_element(ed.energies.begin(), ed.energies.end()) - ed.energies.begin();
    std::vector<uni20::complex<double>> wave(1U << n);
    for (std::size_t i = 0; i < ed.basis.size(); ++i)
      wave[ed.basis[i]] = ed.vectors[i * ed.basis.size() + g];
    auto const ed_rho = trace_first_two<double>(n, wave);
    for (std::size_t i = 0; i < 16; ++i)
      EXPECT_NEAR(static_cast<double>(rho[i]), ed_rho[i].real(), 3e-11);
    EXPECT_GT(*result.singlet_probability(), Real{0});
    EXPECT_GE(*result.triplet_component_probability(), -tolerance);
    EXPECT_REAL_NEAR(rho[0] + rho[5] + rho[10] + rho[15], Real{1}, tolerance);
    EXPECT_REAL_NEAR(rho[0] + rho[5], Real{0.5}, tolerance);
    EXPECT_REAL_NEAR(rho[10] + rho[15], Real{0.5}, tolerance);
  }
}

TYPED_TEST(ExperimentalBoundaryCorrelations, ImpurityFamilyAndHellmannFeynman)
{
  using Real = TypeParam;
  using C = uni20::complex<Real>;
  using uni20::math::abs;
  using uni20::math::pow;
  Real const epsilon = uni20::numeric_limits<Real>::epsilon();
  for (unsigned n : {2, 4, 6, 8})
  {
    SCOPED_TRACE(n);
    // Check a finite deformation against the independent bit-basis Hamiltonian.
    Real const theta = Real{0.2};
    auto const u = impurity_roots<Real>(n, theta);
    std::vector<Real> z;
    for (Real v : u)
      z.push_back(Real{2} * v);
    auto const wave = test_support::open_aba_wave<Real>(n, z, theta);
    Real const energy = impurity_energy<Real>(n, theta, u);
    std::vector<C> applied(1U << n);
    for (unsigned bits = 0; bits < (1U << n); ++bits)
      for (unsigned j = 0; j + 1 < n; ++j)
      {
        bool const opposite = ((bits >> j) & 1U) != ((bits >> (j + 1)) & 1U);
        Real const coupling = j == 0 ? Real{1} / (Real{1} + theta * theta) : Real{1};
        applied[bits] += coupling * (opposite ? Real{-0.25} : Real{0.25}) * wave[bits];
        if (opposite) applied[bits ^ (3U << j)] += coupling * Real{0.5} * wave[bits];
      }
    Real largest{};
    for (unsigned bits = 0; bits < (1U << n); ++bits)
      largest = std::max(largest, Real(uni20::math::abs(applied[bits] - energy * wave[bits])));
    EXPECT_REAL_NEAR(largest, Real{0}, Real{524288} * epsilon);

    // Precision-scaled Richardson test of the energy derivative at theta^2=0.
    Real const step = pow(epsilon, Real{1} / Real{6});
    auto const base = impurity_roots<Real>(n, Real{0});
    auto const first = impurity_roots<Real>(n, step);
    auto const second = impurity_roots<Real>(n, Real{2} * step);
    Real const e0 = impurity_energy<Real>(n, Real{0}, base);
    Real const d1 = (impurity_energy<Real>(n, step, first) - e0) / (step * step);
    Real const d2 = (impurity_energy<Real>(n, Real{2} * step, second) - e0) / (Real{4} * step * step);
    Real const derivative = (Real{4} * d1 - d2) / Real{3};
    auto const density = boundary::boundary_two_site_density(n, model::ground_state<Real>(n));
    ASSERT_TRUE(density.converged());
    EXPECT_REAL_NEAR(-derivative, *density.bond_correlation, Real{4096} * Real(n) * pow(epsilon, Real{2} / Real{3}));
  }
}

TEST(ExperimentalBoundaryCorrelationsGuardrails, FailedRootsAndInvalidLabels)
{
  auto ground = model::ground_state<double>(4);
  EXPECT_THROW((void)boundary::boundary_two_site_density(3, ground), std::invalid_argument);
  auto excited = model::sector_ground_state<double>(4, uni20::half_int{1});
  EXPECT_THROW((void)boundary::boundary_two_site_density(4, excited), std::invalid_argument);
  ground.converged = false;
  auto failed = boundary::boundary_two_site_density(4, ground);
  EXPECT_EQ(failed.status, boundary::BoundaryCorrelationStatus::roots_unconverged);
  EXPECT_FALSE(failed.solve_info);
  EXPECT_FALSE(failed.matrix());
  ground.converged = true;
  ground.rapidities[0] += 0.01;
  auto stale = boundary::boundary_two_site_density(4, ground);
  EXPECT_FALSE(stale.converged());
  EXPECT_TRUE(stale.ground_residual);
  EXPECT_FALSE(stale.solve_info);
  EXPECT_FALSE(stale.bond_correlation);
  EXPECT_THROW((void)boundary::real_root_reflection_parity(4, bethe::heisenberg::QuantumNumbers{uni20::half_int{0}}),
               std::invalid_argument);
}

#if UNI20_HAS_MPREAL
TEST(ExperimentalBoundaryCorrelationsMpfr, ExplicitPrecisionDensityAndResponse)
{
  using R = uni20::mpreal;
  for (int bits : {128, 256})
  {
    auto p = uni20::Precision::bits(bits);
    R const zero(0, p), one(1, p), eps = uni20::numeric_limits<R>::epsilon(zero);
    bethe::SolverOptions<R> options(p);
    for (unsigned n : {2, 4, 6})
    {
      auto ground = model::ground_state<R>(n, options);
      auto result = boundary::boundary_two_site_density(n, ground);
      ASSERT_TRUE(result.converged());
      EXPECT_EQ(result.bond_correlation->precision(), p);
      EXPECT_EQ(result.solve_backward_error->precision(), p);
      auto wave = test_support::open_aba_wave<R>(n, ground.rapidities, zero);
      auto expected = trace_first_two<R>(n, wave);
      auto actual = *result.matrix();
      for (std::size_t i = 0; i < 16; ++i)
      {
        EXPECT_EQ(actual[i].precision(), p);
        EXPECT_REAL_NEAR(uni20::math::abs(expected[i] - actual[i]), zero, 524288 * eps);
      }
      auto shifted = impurity_roots<R>(n, one / 5);
      ASSERT_EQ(shifted.size(), n / 2);
      for (auto const& u : shifted)
        EXPECT_EQ(u.precision(), p);
      // N=2 has one root: promoting it alone is a valid uniform-precision set.
      // Use an exact root there; for larger sets retain a mixed-precision case.
      ground.rapidities[0] = n == 2 ? R{1} : ground.rapidities[0].at(uni20::Precision::bits(bits + 16));
      EXPECT_THROW((void)boundary::boundary_two_site_density(n, ground), std::invalid_argument);
    }
  }
}
#endif
