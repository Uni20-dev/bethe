// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "xxx_complex_fixture.hpp"
#include "xxx_zero_momentum_fixture.hpp"
#include <bethe/experimental/xxx_complex_roots.hpp>

namespace roots = bethe::experimental::xxx::complex_roots;
namespace xxx = bethe::heisenberg;
namespace exp_xxx = bethe::experimental::xxx;

namespace
{
template <uni20::Real R> R parse(std::string_view text, R const& zero)
{
#if UNI20_HAS_MPREAL
  if constexpr (uni20::has_runtime_precision_v<R>)
    return uni20::parse_real<R>(text, zero.precision());
  else
#endif
    return uni20::parse_real<R>(text);
}

template <uni20::Real R> void equation_conditioning(R const& zero)
{
  using C = uni20::complex<R>;
  R const one = uni20::scalar_like(zero, 1), eps = uni20::numeric_limits<R>::epsilon(zero);
  R const spacing = 128 * uni20::math::sqrt(eps);
  std::vector<C> lambda{{one, zero}, {one + spacing, zero}};
  // A small real-real difference does not make log(d-i)-log(d+i)
  // ill-conditioned. It must not loosen the original-equation check.
  EXPECT_EQ(roots::residual_tolerance<R>(20, lambda, zero), 256 * eps * 20);
  lambda[1] = C(one + 2 * spacing, zero);
  EXPECT_EQ(roots::residual_tolerance<R>(20, lambda, zero), 256 * eps * 20);
  lambda[1] = lambda[0];
  EXPECT_THROW((void)roots::residual_tolerance<R>(20, lambda, zero), std::runtime_error);
  // An actual near-string pole still has its original allowance and gate.
  lambda = {{one, zero}, {one, one - spacing}};
  R const pole = spacing * (2 * one - spacing);
  EXPECT_REAL_NEAR(roots::residual_tolerance<R>(20, lambda, zero) * pole / (256 * eps * 20), one, 8 * eps / spacing);
}

template <uni20::Real R> void real_reference(R const& zero)
{
  using C = uni20::complex<R>;
  R const one = uni20::scalar_like(zero, 1), eps = uni20::numeric_limits<R>::epsilon(zero);
  R const z = one / uni20::math::sqrt(3 * one);
  xxx::RealState<R> state;
  state.converged = true;
  state.sz = uni20::half_int{0};
  state.quantum_numbers = xxx::sector_ground_quantum_numbers(4, state.sz);
  state.rapidities = {-z, z};
  exp_xxx::PreparedSinglet<R> const ground(4, state);
  ASSERT_TRUE(ground.ready());
  std::vector<C> const lambda{C(zero, zero)};
  auto const result = roots::raising_weight<R>(ground, lambda);
  ASSERT_TRUE(result.converged());
  ASSERT_TRUE(result.weight);
  EXPECT_REAL_NEAR(*result.weight, 4 * one / 3, 256 * eps);
  EXPECT_REAL_NEAR(roots::energy<R>(4, lambda, zero), -one, 4 * eps);
}

template <uni20::Real R> void complex_reference(R const& zero)
{
  using C = uni20::complex<R>;
  R const one = uni20::scalar_like(zero, 1), eps = uni20::numeric_limits<R>::epsilon(zero);
  xxx::RealState<R> state;
  state.converged = true;
  state.sz = uni20::half_int{0};
  state.quantum_numbers = xxx::sector_ground_quantum_numbers(12, state.sz);
  for (auto value : complex_fixture::ground)
    state.rapidities.push_back(2 * parse(value, zero));
  exp_xxx::PreparedSinglet<R> const ground(12, state);
  ASSERT_TRUE(ground.ready());
  std::vector<C> lambda;
  for (auto const& value : complex_fixture::roots)
    lambda.emplace_back(parse(value[0], zero), parse(value[1], zero));
  auto const weight = roots::raising_weight<R>(ground, lambda);
  ASSERT_TRUE(weight.converged());
  ASSERT_TRUE(weight.weight);
  R const relative_tolerance = 32768 * eps / roots::separation<R>(lambda, zero);
  EXPECT_REAL_NEAR(*weight.weight / parse(complex_fixture::weight, zero), one, relative_tolerance);
  auto const labelled = roots::raising_weight<R>(ground, lambda, 9);
  ASSERT_TRUE(labelled.weight);
  EXPECT_EQ(*labelled.weight, *weight.weight);
  auto const wrong_label = roots::raising_weight<R>(ground, lambda, 0);
  EXPECT_EQ(wrong_label.status, xxx::FormFactorStatus::roots_unconverged);
  EXPECT_FALSE(wrong_label.weight);
  EXPECT_REAL_NEAR(roots::energy<R>(12, lambda, zero), parse(complex_fixture::energy, zero), 4096 * eps);
  EXPECT_LT(roots::rational_residual<R>(12, lambda, zero), roots::residual_tolerance<R>(12, lambda, zero));
  for (auto& z : lambda)
    z = -uni20::math::conj(z);
  auto const mirrored = roots::raising_weight<R>(ground, lambda);
  ASSERT_TRUE(mirrored.weight);
  EXPECT_REAL_NEAR(*mirrored.weight / *weight.weight, one, relative_tolerance);
  lambda.front() += one / 1000;
  auto const offshell = roots::raising_weight<R>(ground, lambda);
  EXPECT_EQ(offshell.status, xxx::FormFactorStatus::roots_unconverged);
  EXPECT_FALSE(offshell.weight);
  lambda[3] = C(one, one / 2);
  lambda[4] = C(one, -one / 2);
  EXPECT_THROW((void)roots::residual_tolerance<R>(12, lambda, zero), std::runtime_error);
  auto const unresolved = roots::raising_weight<R>(ground, lambda);
  EXPECT_EQ(unresolved.status, xxx::FormFactorStatus::precision_limit);
  EXPECT_FALSE(unresolved.weight);
}

template <uni20::Real R> void selection_zero(R const& zero)
{
  using C = uni20::complex<R>;
  R const one = uni20::scalar_like(zero, 1);
  xxx::RealState<R> state;
  state.converged = true;
  state.sz = uni20::half_int{0};
  state.quantum_numbers = xxx::sector_ground_quantum_numbers(12, state.sz);
  for (auto value : complex_fixture::ground)
    state.rapidities.push_back(2 * parse(value, zero));
  exp_xxx::PreparedSinglet<R> const ground(12, state);
  ASSERT_TRUE(ground.ready());
  std::vector<C> lambda;
  for (auto const& value : zero_momentum_fixture::roots)
    lambda.emplace_back(parse(value[0], zero), parse(value[1], zero));
  auto const result = roots::raising_weight<R>(ground, lambda, 0);
  ASSERT_TRUE(result.converged());
  ASSERT_TRUE(result.weight);
  EXPECT_EQ(*result.weight, zero);
  EXPECT_EQ(uni20::precision_of(*result.weight), uni20::precision_of(zero));
  EXPECT_FALSE(result.log_weight);
  EXPECT_TRUE(uni20::isfinite(result.pivot_margin));
  EXPECT_GT(result.pivot_margin, zero);
  auto const wrong_label = roots::raising_weight<R>(ground, lambda, 1);
  EXPECT_EQ(wrong_label.status, xxx::FormFactorStatus::roots_unconverged);
  EXPECT_FALSE(wrong_label.weight);
  EXPECT_THROW((void)roots::raising_weight<R>(ground, lambda, -1), std::invalid_argument);
  EXPECT_THROW((void)roots::raising_weight<R>(ground, lambda, 12), std::invalid_argument);
  auto offshell = lambda;
  offshell.front() += one / 1000;
  auto const failed = roots::raising_weight<R>(ground, offshell, 0);
  EXPECT_EQ(failed.status, xxx::FormFactorStatus::roots_unconverged);
  EXPECT_FALSE(failed.weight);
  lambda[3] = C(one, one / 2);
  lambda[4] = C(one, -one / 2);
  auto const unresolved = roots::raising_weight<R>(ground, lambda, 0);
  EXPECT_EQ(unresolved.status, xxx::FormFactorStatus::precision_limit);
  EXPECT_FALSE(unresolved.weight);
}

template <uni20::Real R> void six_site_momentum(R const& zero)
{
  using C = uni20::complex<R>;
  auto options = [&] {
    if constexpr (uni20::has_runtime_precision_v<R>)
      return bethe::SolverOptions<R>(uni20::precision_of(zero));
    else
      return bethe::SolverOptions<R>{};
  }();
  auto const state = xxx::sector_ground_state<R>(6, uni20::half_int{0}, options);
  ASSERT_TRUE(state.converged);
  exp_xxx::PreparedSinglet<R> const ground(6, state);
  ASSERT_TRUE(ground.ready());
  // Analytic +/- lambda triplet: total absolute momentum zero, but the
  // N=6 singlet has momentum pi. This is an allowed nonzero-transfer weight.
  R const root = uni20::math::tan(bethe::detail::pi_like(zero) / 10) / 2;
  std::vector<C> const lambda{C(-root, zero), C(root, zero)};
  auto const raw = roots::raising_weight<R>(ground, lambda);
  auto const labelled = roots::raising_weight<R>(ground, lambda, 0);
  ASSERT_TRUE(raw.weight);
  ASSERT_TRUE(labelled.weight);
  EXPECT_EQ(*raw.weight, *labelled.weight);
  EXPECT_GT(*labelled.weight, zero);
  auto const wrong_label = roots::raising_weight<R>(ground, lambda, 3);
  EXPECT_EQ(wrong_label.status, xxx::FormFactorStatus::roots_unconverged);
  EXPECT_FALSE(wrong_label.weight);
}

template <uni20::Real R> void derivative_and_cut(R const& zero)
{
  using C = uni20::complex<R>;
  R const one = uni20::scalar_like(zero, 1), eps = uni20::numeric_limits<R>::epsilon(zero);
  for (R const& signed_zero : {zero, -zero})
  {
    C const upper = roots::bethe_atan<R>(C(signed_zero, 2 * one));
    C const lower = roots::bethe_atan<R>(C(signed_zero, -2 * one));
    EXPECT_EQ(upper.real(), -lower.real());
    EXPECT_EQ(upper.imag(), -lower.imag());
    EXPECT_REAL_NEAR(upper.real(), bethe::detail::pi_like(zero) / 2, 4 * eps);
    EXPECT_REAL_NEAR(upper.imag(), uni20::math::log(3 * one) / 2, 4 * eps);
  }
  std::vector<C> const lambda{{one / 10, zero}, {one / 5, 3 * one / 10}, {one / 5, -3 * one / 10}};
  auto const matrix = roots::gaudin<R>(8, lambda, zero);
  R const h = uni20::math::cbrt(eps);
  for (std::size_t j = 0; j < lambda.size(); ++j)
  {
    auto plus = lambda, minus = lambda;
    plus[j] += h;
    minus[j] -= h;
    for (std::size_t i = 0; i < lambda.size(); ++i)
    {
      C const difference =
          (roots::log_equation<R>(8, plus, i, zero) - roots::log_equation<R>(8, minus, i, zero)) / (2 * h);
      C const expected = C(zero, -one) * matrix[i * lambda.size() + j];
      EXPECT_REAL_NEAR(difference.real(), expected.real(), 100000 * h * h);
      EXPECT_REAL_NEAR(difference.imag(), expected.imag(), 100000 * h * h);
    }
  }
  EXPECT_EQ(roots::rational_residual<R>(8, {}, zero), zero);
  EXPECT_EQ(roots::energy<R>(8, {}, zero), 2 * one);
  EXPECT_TRUE(roots::gaudin<R>(8, {}, zero).empty());
}
} // namespace

template <typename T> class ExperimentalComplexRoots : public ::testing::Test {};
TYPED_TEST_SUITE(ExperimentalComplexRoots, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(ExperimentalComplexRoots, FourSiteAnalyticRaisingWeight) { real_reference(TypeParam{0}); }
TYPED_TEST(ExperimentalComplexRoots, IndependentComplexPairReference) { complex_reference(TypeParam{0}); }
TYPED_TEST(ExperimentalComplexRoots, AnalyticDerivativeAndAtanCut) { derivative_and_cut(TypeParam{0}); }
TYPED_TEST(ExperimentalComplexRoots, ZeroTransferAfterStateAndNormValidation) { selection_zero(TypeParam{0}); }
TYPED_TEST(ExperimentalComplexRoots, SixSiteGroundMomentumIsPi) { six_site_momentum(TypeParam{0}); }
TYPED_TEST(ExperimentalComplexRoots, CoincidenceDoesNotInflateEquationAllowance)
{
  equation_conditioning(TypeParam{0});
}

#if UNI20_HAS_FLOAT128
TEST(ExperimentalComplexRootsBinary128, RejectCollapsedCandidateBeforeFormFactor)
{
  using R = uni20::float128;
  using C = uni20::complex<R>;
  R const zero{0};
  // The former N=20 three-string solver accepted this off-shell iterate.
  // Keep the roots directly so bypassing solve() exercises the common
  // form-factor validation too. The close pair is resolved arithmetically,
  // but its original-equation residual is far above the pole-only allowance.
  std::vector<C> lambda;
  for (auto value : {"-1.27023803666736062908322553098408307", "-0.169234877615986008057963141028373053",
                     "-0.0704526355542760254770468857633289889", "0.226135650644512053040391776309456972",
                     "0.361959344646961266509367488125691848", "0.616462100468790789020457681214401906",
                     "0.616462100468791790395183024690336766"})
    lambda.emplace_back(parse(value, zero), zero);
  R const centre = parse("0.665651412276037647799913525849466774", zero);
  R const height = parse("0.958396298694217976753921955944205754", zero);
  lambda.emplace_back(centre, height);
  lambda.emplace_back(centre, -height);
  ASSERT_GT(roots::separation<R>(lambda, zero), 64 * uni20::math::sqrt(uni20::numeric_limits<R>::epsilon(zero)));
  auto const state = xxx::sector_ground_state<R>(20, uni20::half_int{0});
  exp_xxx::PreparedSinglet<R> const ground(20, state);
  ASSERT_TRUE(ground.ready());
  auto const factor = roots::raising_weight<R>(ground, lambda, 6);
  EXPECT_EQ(factor.status, xxx::FormFactorStatus::roots_unconverged);
  EXPECT_FALSE(factor.weight);
}
#endif

#if UNI20_HAS_MPREAL
TEST(ExperimentalComplexRootsMpfr, ExplicitPrecisionAndValidation)
{
  using R = uni20::mpreal;
  using C = uni20::complex<R>;
  for (int bits : {128, 256})
  {
    auto const p = uni20::Precision::bits(bits);
    R const zero(0, p), one(1, p);
    real_reference(zero);
    complex_reference(zero);
    derivative_and_cut(zero);
    selection_zero(zero);
    six_site_momentum(zero);
    equation_conditioning(zero);
    EXPECT_EQ(roots::energy<R>(8, {}, zero).precision(), p);
    EXPECT_THROW((void)roots::energy<R>(8, {}), std::invalid_argument);
    for (C const& bad : {C{}, C(uni20::uninitialized), C(R{1}, R{0}),
                         C(R(1, uni20::Precision::bits(bits + 32)), R(0, uni20::Precision::bits(bits + 32)))})
    {
      std::vector<C> const values{bad};
      EXPECT_THROW((void)roots::gaudin<R>(8, values, zero), std::invalid_argument);
      EXPECT_THROW((void)roots::energy<R>(8, values, zero), std::invalid_argument);
      EXPECT_THROW((void)roots::rational_residual<R>(8, values, zero), std::invalid_argument);
      std::vector<C> invalid(5, bad);
      xxx::RealState<R> gs;
      gs.converged = true;
      gs.sz = uni20::half_int{0};
      gs.quantum_numbers = xxx::sector_ground_quantum_numbers(12, gs.sz);
      for (auto value : complex_fixture::ground)
        gs.rapidities.push_back(2 * parse(value, zero));
      exp_xxx::PreparedSinglet<R> const prepared(12, gs);
      EXPECT_THROW((void)roots::raising_weight<R>(prepared, invalid, 0), std::invalid_argument);
    }
    EXPECT_THROW((void)roots::bethe_atan<R>(C{}), std::invalid_argument);
    EXPECT_THROW((void)roots::bethe_atan<R>(C(uni20::uninitialized)), std::invalid_argument);
    EXPECT_THROW((void)roots::bethe_atan<R>(C(R{1}, R{0})), std::invalid_argument);
  }
}
#endif
