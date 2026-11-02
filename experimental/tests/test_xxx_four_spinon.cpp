// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "xxx_complex_fixture.hpp"
#include "xxx_ground_fixture.hpp"
#include <bethe/experimental/xxx_four_spinon.hpp>

namespace f = bethe::experimental::xxx;
namespace
{
template <uni20::Real R> R parse(std::string_view value, R const& zero)
{
#if UNI20_HAS_MPREAL
  if constexpr (uni20::has_runtime_precision_v<R>)
    return uni20::parse_real<R>(value, zero.precision());
  else
#endif
    return uni20::parse_real<R>(value);
}

template <uni20::Real R> void solved_reference(R const& zero)
{
  R const one = uni20::scalar_like(zero, 1), eps = uni20::numeric_limits<R>::epsilon(zero);
  auto result = f::solve_four_spinon<R>({12, {-3, -1, 0, 2}, 1}, zero);
  ASSERT_TRUE(result.converged()) << result.diagnostic;
  auto const& state = *result.state;
  ASSERT_EQ(state.real_rapidities().size(), 3);
  EXPECT_EQ(state.momentum_index(), 9);
  for (std::size_t i = 0; i < 3; ++i)
    EXPECT_REAL_NEAR(state.real_rapidities()[i], parse(complex_fixture::roots[i][0], zero), 8192 * eps);
  EXPECT_REAL_NEAR(state.centre(), parse(complex_fixture::roots[3][0], zero), 8192 * eps);
  EXPECT_REAL_NEAR(state.deviation() / (parse(complex_fixture::roots[3][1], zero) - one / 2), one, 32768 * eps);
  EXPECT_REAL_NEAR(state.energy(), parse(complex_fixture::energy, zero), 8192 * eps);
  EXPECT_LT(state.rational_residual(), 512 * 12 * eps);
  EXPECT_LT(state.quantum_residual(), 32 * eps);
  auto const canonical = f::solve_four_spinon<R>({12, {2, 0, -3, -1}, 1}, zero);
  ASSERT_TRUE(canonical.converged()) << canonical.diagnostic;
  EXPECT_EQ(canonical.request, result.request);
  EXPECT_EQ(canonical.state->deviation(), state.deviation());
  auto const mirror = f::solve_four_spinon<R>({12, {-2, 0, 1, 3}, -1}, zero);
  ASSERT_TRUE(mirror.converged()) << mirror.diagnostic;
  EXPECT_EQ(mirror.state->momentum_index(), 3);
  EXPECT_REAL_NEAR(mirror.state->centre(), -state.centre(), 8192 * eps);
  EXPECT_REAL_NEAR(mirror.state->energy(), state.energy(), 8192 * eps);
  auto copied = state;
  auto moved = std::move(copied);
  EXPECT_FALSE(copied.ready());
  EXPECT_THROW((void)copied.roots(), std::logic_error);
  result.state.reset();
  EXPECT_TRUE(moved.ready());
  EXPECT_EQ(moved.request().sites, 12);
}

template <uni20::Real R> void tiny_deviations(R const& zero)
{
  R const one = uni20::scalar_like(zero, 1), eps = uni20::numeric_limits<R>::epsilon(zero);
  auto const negative = f::solve_four_spinon<R>({16, {-1, 0, 2, 3}, 0}, zero);
  ASSERT_TRUE(negative.converged()) << negative.diagnostic;
  EXPECT_REAL_NEAR(negative.state->deviation() / parse("-7.56318282038341941997737176833694899404308444e-12", zero),
                   one, std::max(100000 * eps, parse("1e-40", zero)));
  auto const tiny = f::solve_four_spinon<R>({64, {-14, -12, -11, -10}, 1}, zero);
  ASSERT_TRUE(tiny.converged()) << tiny.diagnostic;
  EXPECT_REAL_NEAR(tiny.state->deviation() / parse("1.41905199057999036182880353262087001032346067e-94", zero), one,
                   std::max(1000000 * eps, parse("1e-40", zero)));
  EXPECT_NE(tiny.state->deviation(), zero);
  EXPECT_EQ(tiny.state->roots().back().imag(), -one / 2);
  EXPECT_EQ(tiny.state->momentum_index(), 16);
  EXPECT_LT(tiny.state->rational_residual(), 512 * 64 * eps);
  EXPECT_REAL_NEAR(tiny.state->log_abs_deviation(), uni20::math::log(tiny.state->deviation()), 1024 * eps);
}

template <uni20::Real R> void physical_weights(R const& zero)
{
  R const one = uni20::scalar_like(zero, 1), eps = uni20::numeric_limits<R>::epsilon(zero);
  for (auto const request :
       {f::FourSpinonRequest{12, {-3, -1, 0, 2}, 1}, {16, {-1, 0, 2, 3}, 0}, {64, {-14, -12, -11, -10}, 1}})
  {
    bethe::heisenberg::RealState<R> ground_state;
    ground_state.converged = true;
    ground_state.sz = uni20::half_int{0};
    ground_state.quantum_numbers = bethe::heisenberg::sector_ground_quantum_numbers(request.sites, ground_state.sz);
    std::span<std::string_view const> const background =
        request.sites == 12   ? std::span<std::string_view const>(complex_fixture::ground)
        : request.sites == 16 ? std::span<std::string_view const>(ground_fixture::n16)
                              : std::span<std::string_view const>(ground_fixture::n64);
    for (auto value : background)
      ground_state.rapidities.push_back(2 * parse(value, zero));
    f::PreparedSinglet<R> const ground(request.sites, ground_state);
    ASSERT_TRUE(ground.ready());
    auto const result = f::solve_four_spinon<R>(request, zero);
    ASSERT_TRUE(result.converged()) << result.diagnostic;
    auto const weight = f::four_spinon_raising_weight(ground, *result.state);
    ASSERT_TRUE(weight.converged());
    ASSERT_TRUE(weight.weight);
    auto const expected = request.sites == 12   ? complex_fixture::weight
                          : request.sites == 16 ? "8.28523087152927669547401134062756665438561737e-8"
                                                : "1.33985918781507196786052358985966929800159504e-11";
    R const tolerance = request.sites == 12 ? 1000000 * eps : std::max(10000000 * eps, parse("1e-40", zero));
    EXPECT_REAL_NEAR(*weight.weight / parse(expected, zero), one, tolerance);
    EXPECT_REAL_NEAR(*weight.log_weight, uni20::math::log(*weight.weight), 1024 * eps);
    if (request.sites == 12)
    {
      auto const direct = f::complex_roots::raising_weight<R>(ground, result.state->roots());
      ASSERT_TRUE(direct.weight);
      EXPECT_REAL_NEAR(*weight.weight / *direct.weight, one, 32768 * eps / uni20::math::abs(result.state->deviation()));
      auto const mirror = f::solve_four_spinon<R>({12, {-2, 0, 1, 3}, -1}, zero);
      ASSERT_TRUE(mirror.converged());
      auto const mirrored = f::four_spinon_raising_weight(ground, *mirror.state);
      ASSERT_TRUE(mirrored.weight);
      EXPECT_REAL_NEAR(*mirrored.weight / *weight.weight, one, 1000000 * eps);
      auto moved = *result.state;
      auto live = std::move(moved);
      EXPECT_THROW((void)f::four_spinon_raising_weight(ground, moved), std::invalid_argument);
      EXPECT_TRUE(live.ready());
#if UNI20_HAS_MPREAL
      if constexpr (uni20::has_runtime_precision_v<R>)
      {
        R const other_zero(0, uni20::Precision::bits(uni20::numeric_limits<R>::digits(zero) + 32));
        auto const other = f::solve_four_spinon<R>(request, other_zero);
        ASSERT_TRUE(other.converged());
        EXPECT_THROW((void)f::four_spinon_raising_weight(ground, *other.state), std::invalid_argument);
      }
#endif
    }
    if (request.sites == 64)
    {
      auto const direct = f::complex_roots::raising_weight<R>(ground, result.state->roots());
      EXPECT_EQ(direct.status, bethe::heisenberg::FormFactorStatus::precision_limit);
      EXPECT_FALSE(direct.weight);
      auto const matrix = f::four_detail::gaudin_matrix(*result.state);
      auto const s = result.state->real_rapidities().size() + 2;
      // This transformed row has no unit offset: even sub-epsilon delta
      // must survive. An ideal-string substitution would zero it.
      EXPECT_NE(matrix[(s - 1) * s], zero);
    }
  }
}

template <uni20::Real R> void symmetry_zero(R const& zero)
{
  bethe::heisenberg::RealState<R> ground_state;
  ground_state.converged = true;
  ground_state.sz = uni20::half_int{0};
  ground_state.quantum_numbers = bethe::heisenberg::sector_ground_quantum_numbers(12, ground_state.sz);
  for (auto value : complex_fixture::ground)
    ground_state.rapidities.push_back(2 * parse(value, zero));
  f::PreparedSinglet<R> const ground(12, ground_state);
  ASSERT_TRUE(ground.ready());
  for (auto const request : {f::FourSpinonRequest{12, {-2, -1, 0, 2}, -1}, {12, {-3, 0, 1, 2}, 0}})
  {
    auto const solved = f::solve_four_spinon<R>(request, zero);
    ASSERT_TRUE(solved.converged()) << solved.diagnostic;
    ASSERT_EQ(solved.state->momentum_index(), 0);
    auto const factor = f::four_spinon_raising_weight(ground, *solved.state);
    ASSERT_TRUE(factor.converged());
    ASSERT_TRUE(factor.weight);
    EXPECT_EQ(*factor.weight, zero);
    EXPECT_FALSE(factor.log_weight); // A symmetry zero has no finite logarithm.
    EXPECT_TRUE(uni20::isfinite(factor.pivot_margin));
    EXPECT_GT(factor.pivot_margin, zero); // Ground and excited norms still checked.
    EXPECT_EQ(uni20::precision_of(*factor.weight), uni20::precision_of(zero));
  }
}

template <uni20::Real R> void jacobian_and_failures(R const& zero)
{
  R const one = uni20::scalar_like(zero, 1), eps = uni20::numeric_limits<R>::epsilon(zero);
  std::array<int, 3> const numbers{-3, 1, 2};
  std::vector<R> const x{-6 * one / 10, one / 10, 8 * one / 10, 2 * one / 10, uni20::math::log(4 * one / 1000)};
  R const h = uni20::math::cbrt(eps);
  for (int sign : {-1, 1})
  {
    auto const jac = f::four_detail::jacobian<R>(12, sign, x, zero);
    for (std::size_t j = 0; j < x.size(); ++j)
    {
      auto plus = x, minus = x;
      plus[j] += h;
      minus[j] -= h;
      auto const fp = f::four_detail::equations<R>(12, numbers, sign, plus, zero);
      auto const fm = f::four_detail::equations<R>(12, numbers, sign, minus, zero);
      for (std::size_t i = 0; i < x.size(); ++i)
        EXPECT_REAL_NEAR((fp[i] - fm[i]) / (2 * h), jac[i * x.size() + j], 1000 * h * h);
    }
  }
  auto const singular = f::solve_four_spinon<R>({12, {-3, -1, 1, 3}, 0}, zero);
  EXPECT_EQ(singular.status, f::FourSpinonStatus::unsupported_geometry);
  EXPECT_FALSE(singular.state);
  EXPECT_FALSE(singular.diagnostic.empty());
  EXPECT_EQ(singular.request.sites, 12);
  for (auto request : {f::FourSpinonRequest{6, {-1, 0, 1, 2}, 0},
                       {12, {-3, -1, -1, 2}, 0},
                       {12, {-4, -1, 0, 2}, 0},
                       {12, {-3, -1, 0, 2}, 2}})
    EXPECT_THROW((void)f::solve_four_spinon<R>(request, zero), std::invalid_argument);
}
} // namespace
template <typename T> class ExperimentalFourSpinon : public ::testing::Test {};
TYPED_TEST_SUITE(ExperimentalFourSpinon, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(ExperimentalFourSpinon, SelectedStateAndLifetime) { solved_reference(TypeParam{0}); }
TYPED_TEST(ExperimentalFourSpinon, ExplicitTinyDeviations) { tiny_deviations(TypeParam{0}); }
TYPED_TEST(ExperimentalFourSpinon, JacobianAndFailures) { jacobian_and_failures(TypeParam{0}); }
TYPED_TEST(ExperimentalFourSpinon, NormalizedRaisingWeights) { physical_weights(TypeParam{0}); }
TYPED_TEST(ExperimentalFourSpinon, ZeroMomentumSelectionRule) { symmetry_zero(TypeParam{0}); }
#if UNI20_HAS_MPREAL
TEST(ExperimentalFourSpinonMpfr, PrecisionGrowthAndValidation)
{
  using R = uni20::mpreal;
  for (int bits : {128, 256})
  {
    R const zero(0, uni20::Precision::bits(bits));
    solved_reference(zero);
    tiny_deviations(zero);
    jacobian_and_failures(zero);
    physical_weights(zero);
    symmetry_zero(zero);
    auto const result = f::solve_four_spinon<R>({12, {-3, -1, 0, 2}, 1}, zero);
    ASSERT_TRUE(result.converged());
    EXPECT_EQ(result.state->deviation().precision(), zero.precision());
    EXPECT_EQ(result.state->energy().precision(), zero.precision());
  }
  EXPECT_THROW((void)f::solve_four_spinon<R>({12, {-3, -1, 0, 2}, 1}), std::invalid_argument);
  EXPECT_THROW((void)f::solve_four_spinon<R>({12, {-3, -1, 0, 2}, 1}, R(uni20::uninitialized)), std::invalid_argument);
}
#endif
