// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "exact_spectrum.hpp"
#include "test_support.hpp"
#include "xxx_open_wave_oracle.hpp"
#include <bethe/experimental/xxx_boundary_form_factors.hpp>
#include <bethe/heisenberg_excitations.hpp>

namespace f = bethe::experimental::xxx::open;
namespace model = bethe::heisenberg::open;
namespace
{
template <uni20::Real R> void check_fixed_root(R const& zero)
{
  using C = uni20::complex<R>;
  R const one = uni20::scalar_like(zero, 1);
  // An operator identity on every bit-basis state, independent of root solving.
  for (unsigned n : {2, 3, 4, 5})
    for (unsigned bits = 0; bits < (1U << n); ++bits)
    {
      std::vector<C> basis(1U << n, C(zero, zero));
      basis[bits] = C(one, zero);
      auto const result = test_support::open_aba_create<R>(n, basis, C(zero, -one / 2), zero);
      for (unsigned i = 0; i < (1U << n); ++i)
      {
        R const expected =
            !(bits & (1U << (n - 1))) && i == (bits | (1U << (n - 1))) ? (n % 2 ? 2 * one : -2 * one) : zero;
        EXPECT_EQ(result[i].real(), expected);
        EXPECT_EQ(result[i].imag(), zero);
      }
    }
}

template <uni20::Real R> void check_determinants(R const& zero)
{
  using C = uni20::complex<R>;
  R const one = uni20::scalar_like(zero, 1), eps = uni20::numeric_limits<R>::epsilon(zero), tol = 1024 * eps;
  std::vector<R> matrix{zero, 2 * one, 3 * one, 4 * one};
  auto d = bethe::detail::log_determinant<R>(matrix, 2, zero);
  ASSERT_TRUE(d.log_absolute);
  ASSERT_TRUE(d.phase);
  EXPECT_EQ(*d.phase, -one);
  EXPECT_REAL_NEAR(*d.log_absolute, uni20::math::log(6 * one), tol);
  auto magnitude = bethe::detail::log_absolute_determinant<R>(matrix, 2, zero);
  EXPECT_EQ(magnitude.log_absolute, d.log_absolute);
  EXPECT_EQ(magnitude.pivot_margin, d.pivot_margin);
  std::vector<C> complex{C(zero, zero), C(zero, 2 * one), C(3 * one, zero), C(4 * one, zero)};
  auto cd = bethe::detail::log_determinant<R>(complex, 2, zero);
  ASSERT_TRUE(cd.phase);
  ASSERT_TRUE(cd.log_absolute);
  EXPECT_REAL_NEAR(cd.phase->real(), zero, tol);
  EXPECT_REAL_NEAR(cd.phase->imag(), -one, tol);
  EXPECT_REAL_NEAR(*cd.log_absolute, *d.log_absolute, tol);
  auto empty = bethe::detail::log_determinant<R>(std::vector<C>{}, 0, zero);
  ASSERT_TRUE(empty.phase);
  EXPECT_EQ(empty.phase->real(), one);
  EXPECT_EQ(empty.phase->imag(), zero);
  EXPECT_EQ(*empty.log_absolute, zero);
  auto singular = bethe::detail::log_determinant<R>(std::vector<R>{one, one, one, one}, 2, zero);
  EXPECT_FALSE(singular.log_absolute);
  EXPECT_FALSE(singular.phase);
  auto rejected = bethe::detail::log_determinant<R>(std::vector<R>{one, one, one, one + eps}, 2, zero);
  EXPECT_FALSE(rejected.log_absolute);
  EXPECT_FALSE(rejected.phase);
  EXPECT_THROW(bethe::detail::log_determinant<R>(matrix, 3, zero), std::invalid_argument);
  if constexpr (uni20::has_runtime_precision_v<R>)
  {
    EXPECT_EQ(d.phase->precision(), zero.precision());
    EXPECT_EQ(empty.phase->real().precision(), zero.precision());
    EXPECT_EQ(empty.log_absolute->precision(), zero.precision());
  }
}

template <uni20::Real R> void check_wave_weights(unsigned n, bethe::SolverOptions<R> const& options)
{
  using C = uni20::complex<R>;
  R const zero = bethe::detail::solver_zero(options), one = uni20::scalar_like(zero, 1);
  auto scan = model::real_excitations<R>(n, uni20::half_int{1}, {.count = 1000}, options);
  ASSERT_TRUE(scan.converged());
  f::BoundaryFormFactors<R> prepared(n, scan.ground_state);
  ASSERT_TRUE(prepared.ready());
  ASSERT_TRUE(prepared.ground_energy());
  EXPECT_REAL_NEAR(*prepared.ground_energy(), scan.ground_state.energy, 4096 * uni20::numeric_limits<R>::epsilon(zero));
  auto const ground_wave = test_support::open_aba_wave<R>(n, scan.ground_state.rapidities, zero);
  R total = zero;
  for (auto const& level : scan.levels)
  {
    auto factor = prepared(level.state);
    ASSERT_TRUE(factor.converged());
    ASSERT_TRUE(factor.weight);
    ASSERT_TRUE(factor.log_weight);
    EXPECT_NE(factor.determinant_sign, 0);
    auto wave = test_support::open_aba_wave<R>(n, level.state.rapidities, zero);
    C overlap(zero, zero);
    for (unsigned bits = 0; bits < (1U << n); ++bits)
      if (bits & 1U) overlap += uni20::conj(wave[bits ^ 1U]) * ground_wave[bits];
    R const expected = overlap.real() * overlap.real() + overlap.imag() * overlap.imag();
    EXPECT_REAL_NEAR(*factor.weight, expected, 524288 * uni20::numeric_limits<R>::epsilon(zero));
    total += *factor.weight;
    if constexpr (uni20::has_runtime_precision_v<R>)
    {
      EXPECT_EQ(factor.weight->precision(), zero.precision());
      EXPECT_EQ(factor.log_weight->precision(), zero.precision());
      EXPECT_EQ(factor.pivot_margin->precision(), zero.precision());
    }
  }
  EXPECT_LE(total, one / 2 + 524288 * uni20::numeric_limits<R>::epsilon(zero));
  if (n <= 4) EXPECT_REAL_NEAR(total, one / 2, 524288 * uni20::numeric_limits<R>::epsilon(zero));
}

template <uni20::Real R> void check_lifetime(bethe::SolverOptions<R> const& options)
{
  R const zero = bethe::detail::solver_zero(options), one = uni20::scalar_like(zero, 1);
  auto ground = model::ground_state<R>(4, options);
  auto excited = model::sector_ground_state<R>(4, uni20::half_int{1}, options);
  f::BoundaryFormFactors<R> source(4, ground);
  EXPECT_THROW((f::BoundaryFormFactors<R>(3, ground)), std::invalid_argument);
  EXPECT_THROW((void)source(ground), std::invalid_argument);
  auto malformed = ground;
  malformed.spin_reversed = true;
  EXPECT_THROW((f::BoundaryFormFactors<R>(4, malformed)), std::invalid_argument);
  malformed = ground;
  malformed.quantum_numbers[0] = uni20::half_int{0};
  EXPECT_THROW((f::BoundaryFormFactors<R>(4, malformed)), std::invalid_argument);
  malformed = ground;
  malformed.rapidities[0] = zero;
  EXPECT_FALSE((f::BoundaryFormFactors<R>(4, malformed)).ready());
  malformed = ground;
  std::reverse(malformed.rapidities.begin(), malformed.rapidities.end());
  EXPECT_FALSE((f::BoundaryFormFactors<R>(4, malformed)).ready());
  malformed = ground;
  malformed.converged = false;
  f::BoundaryFormFactors<R> unconverged(4, malformed);
  EXPECT_EQ(unconverged.status(), f::BoundaryFormFactorStatus::roots_unconverged);
  EXPECT_FALSE(unconverged.ground_residual());
  EXPECT_FALSE(unconverged(excited).weight);
  auto expected = source(excited);
  ASSERT_TRUE(expected.converged());
  auto residual = source.ground_residual();
  ASSERT_TRUE(residual);
  auto moved = std::move(source);
  EXPECT_FALSE(source.ready());
  EXPECT_FALSE(source.ground_residual());
  EXPECT_FALSE(source(excited).weight);
  EXPECT_EQ(moved.ground_residual(), residual);
  auto copy = moved;
  source = std::move(moved);
  EXPECT_FALSE(moved.ready());
  EXPECT_FALSE(moved.ground_residual());
  auto& alias = source;
  source = std::move(alias);
  EXPECT_TRUE(source.ready());
  ground.rapidities.front() += one / 8;
  f::BoundaryFormFactors<R> bad(4, ground);
  EXPECT_FALSE(bad.ready());
  EXPECT_TRUE(bad.ground_residual());
  EXPECT_FALSE(bad.ground_energy());
  EXPECT_FALSE(bad(excited).weight);
  ground = {};
  EXPECT_EQ(source(excited).weight, expected.weight);
  EXPECT_EQ(copy(excited).weight, expected.weight);
  excited.rapidities.front() += one / 8;
  excited.residual_norm = zero;
  EXPECT_FALSE(source(excited).weight);
  excited.converged = false;
  EXPECT_EQ(source(excited).status, f::BoundaryFormFactorStatus::roots_unconverged);
}
} // namespace

template <typename R> class ExperimentalBoundary : public ::testing::Test {};
TYPED_TEST_SUITE(ExperimentalBoundary, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(ExperimentalBoundary, DeterminantPhasesAndFailure) { check_determinants(TypeParam{0}); }
TYPED_TEST(ExperimentalBoundary, FixedImaginaryRootIsLocalOperator) { check_fixed_root(TypeParam{0}); }
TYPED_TEST(ExperimentalBoundary, IndependentABAVectors)
{
  for (unsigned n : {2, 4, 6, 8})
  {
    SCOPED_TRACE(n);
    check_wave_weights<TypeParam>(n, {});
  }
}
TYPED_TEST(ExperimentalBoundary, OwnershipMovesAndStaleRoots) { check_lifetime<TypeParam>({}); }

TEST(ExperimentalBoundaryED, IndividualLocalWeights)
{
  for (unsigned n : {2, 4, 6, 8})
  {
    SCOPED_TRACE(n);
    auto scan = model::real_excitations<double>(n, uni20::half_int{1}, {.count = 1000});
    ASSERT_TRUE(scan.converged());
    f::BoundaryFormFactors<double> prepared(n, scan.ground_state);
    auto gs = test_support::exact_eigensystem(n, n / 2, 0, false);
    auto ex = test_support::exact_eigensystem(n, n / 2 - 1, 0, false);
    auto g = std::min_element(gs.energies.begin(), gs.energies.end()) - gs.energies.begin();
    std::vector<double> raised(1U << n, 0);
    for (std::size_t i = 0; i < gs.basis.size(); ++i)
      if (gs.basis[i] & 1U) raised[gs.basis[i] ^ 1U] = gs.vectors[i * gs.basis.size() + g];
    for (auto const& level : scan.levels)
    {
      auto factor = prepared(level.state);
      ASSERT_TRUE(factor.weight);
      double weight = 0;
      for (std::size_t k = 0; k < ex.energies.size(); ++k)
        if (std::abs(ex.energies[k] - level.state.energy) < 1e-10)
        {
          double overlap = 0;
          for (std::size_t i = 0; i < ex.basis.size(); ++i)
            overlap += ex.vectors[i * ex.basis.size() + k] * raised[ex.basis[i]];
          weight += overlap * overlap;
        }
      EXPECT_NEAR(*factor.weight, weight, 2e-10);
    }
  }
}

#if UNI20_HAS_MPREAL
TEST(ExperimentalBoundaryMpfr, ExplicitPrecisionAndInputRejection)
{
  using R = uni20::mpreal;
  for (int bits : {128, 256})
  {
    auto p = uni20::Precision::bits(bits);
    R const zero(0, p);
    bethe::SolverOptions<R> options(p);
    check_determinants(zero);
    check_fixed_root(zero);
    check_lifetime(options);
    for (unsigned n : {2, 4})
      check_wave_weights(n, options);
    auto ground = model::ground_state<R>(4, options);
    auto excited = model::sector_ground_state<R>(4, uni20::half_int{1}, options);
    f::BoundaryFormFactors<R> prepared(4, ground);
    excited.rapidities.front() = excited.rapidities.front().at(uni20::Precision::bits(bits + 32));
    EXPECT_THROW((void)prepared(excited), std::invalid_argument);
    for (auto bad : {R(1), R(uni20::uninitialized), R(1, uni20::Precision::bits(bits + 32))})
    {
      ground.rapidities.front() = bad;
      EXPECT_THROW((f::BoundaryFormFactors<R>(4, ground)), std::invalid_argument);
    }
  }
}
#endif
