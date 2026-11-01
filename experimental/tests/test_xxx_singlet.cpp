// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/experimental/xxx_singlet.hpp>

namespace model = bethe::heisenberg;
namespace experiment = bethe::experimental::xxx;

namespace
{
template <uni20::Real R> model::RealState<R> analytic_singlet(R const& zero)
{
  R const root = uni20::scalar_like(zero, 1) / uni20::math::sqrt(uni20::scalar_like(zero, 3));
  model::RealState<R> state;
  state.rapidities = {-root, root};
  state.quantum_numbers = model::sector_ground_quantum_numbers(4, uni20::half_int{0});
  state.sz = uni20::half_int{0};
  state.converged = true;
  // Deliberately false metadata: roots, not the caller's energy, define this snapshot.
  state.energy = uni20::scalar_like(zero, 123);
  return state;
}

template <uni20::Real R> model::RealState<R> analytic_triplet(R const& zero)
{
  model::RealState<R> state;
  state.rapidities = {zero};
  state.quantum_numbers = model::sector_ground_quantum_numbers(4, uni20::half_int{1});
  state.sz = uni20::half_int{1};
  state.converged = true;
  return state;
}

template <uni20::Real R> void check_analytic_singlet(R const& zero)
{
  auto ground = analytic_singlet(zero);
  experiment::PreparedSinglet<R> const prepared(4, ground);
  ASSERT_TRUE(prepared.ready());
  ASSERT_TRUE(prepared.energy());
  ASSERT_TRUE(prepared.working_zero());
  EXPECT_EQ(*prepared.working_zero(), zero);
  ASSERT_TRUE(prepared.ground_log_prefactor());
  ASSERT_TRUE(prepared.ground_pivot_margin());
  R const tolerance = uni20::scalar_like(zero, 4096) * uni20::numeric_limits<R>::epsilon(zero);
  EXPECT_REAL_NEAR(*prepared.energy(), uni20::scalar_like(zero, -2), tolerance);
  // mu=+/-1/(2 sqrt(3)): det G=108; the normalized ground prefactor is 1/324.
  EXPECT_REAL_NEAR(*prepared.ground_log_prefactor(), -uni20::math::log(uni20::scalar_like(zero, 324)), tolerance);
  EXPECT_GT(*prepared.ground_pivot_margin(), zero);
  EXPECT_LE(*prepared.ground_pivot_margin(), uni20::scalar_like(zero, 1));
  ASSERT_EQ(prepared.aba_rapidities().size(), 2);
  EXPECT_EQ(prepared.aba_rapidities()[0], ground.rapidities[0] / 2);
  // Independent N4 Fourier raising weight from the singlet to the pi triplet.
  auto const factor = prepared.real_raising_form_factor(analytic_triplet(zero));
  ASSERT_TRUE(factor.converged());
  ASSERT_TRUE(factor.weight);
  EXPECT_REAL_NEAR(*factor.weight, uni20::scalar_like(zero, 4) / 3, tolerance);
}
} // namespace

template <typename T> class ExperimentalSinglet : public ::testing::Test {};
TYPED_TEST_SUITE(ExperimentalSinglet, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(ExperimentalSinglet, AnalyticNormalizationAndEnergy) { check_analytic_singlet(TypeParam{0}); }

TYPED_TEST(ExperimentalSinglet, RechecksRootsAndOwnsSnapshot)
{
  using R = TypeParam;
  auto ground = analytic_singlet(R{0});
  experiment::PreparedSinglet<R> const saved(4, ground);
  ASSERT_TRUE(saved.ready());
  ground.rapidities.front() += R{0.125};
  ground.residual_norm = R{0};
  experiment::PreparedSinglet<R> bad(4, ground);
  EXPECT_FALSE(bad.ready());
  EXPECT_FALSE(bad.energy());
  EXPECT_FALSE(bad.working_zero());
  EXPECT_FALSE(bad.ground_log_prefactor());
  EXPECT_FALSE(bad.ground_pivot_margin());
  EXPECT_TRUE(bad.aba_rapidities().empty());
  EXPECT_FALSE(bad.real_raising_form_factor(analytic_triplet(R{0})).weight);
  model::PreparedRaisingFormFactors<R> rejected(4, ground);
  EXPECT_FALSE(rejected.ground_log_prefactor());
  EXPECT_FALSE(rejected.ground_pivot_margin());
  ground = {};
  ASSERT_TRUE(saved.real_raising_form_factor(analytic_triplet(R{0})).converged());
  EXPECT_TRUE(saved.energy());
}

TYPED_TEST(ExperimentalSinglet, MoveAndFailedPreparationDoNotExposeUsableFactors)
{
  using R = TypeParam;
  auto ground = analytic_singlet(R{0});
  experiment::PreparedSinglet<R> source(4, ground);
  auto const expected = source.energy();
  auto target = std::move(source);
  EXPECT_TRUE(target.ready());
  EXPECT_EQ(target.energy(), expected);
  EXPECT_FALSE(source.ready());
  EXPECT_FALSE(source.energy());
  EXPECT_FALSE(source.working_zero());
  EXPECT_FALSE(source.ground_log_prefactor());
  EXPECT_TRUE(source.aba_rapidities().empty());
  ground.converged = false;
  experiment::PreparedSinglet<R> failed(4, ground);
  EXPECT_EQ(failed.status(), model::FormFactorStatus::roots_unconverged);
  EXPECT_FALSE(failed.energy());
  auto copy = failed;
  copy = target;
  EXPECT_TRUE(copy.ready());
  EXPECT_EQ(copy.energy(), expected);
  failed = std::move(target);
  EXPECT_TRUE(failed.ready());
  EXPECT_FALSE(target.ready());
  EXPECT_FALSE(target.energy());
  EXPECT_TRUE(copy.ready());
  auto& same = failed;
  failed = std::move(same);
  EXPECT_TRUE(failed.ready());
}

#if UNI20_HAS_MPREAL
TEST(ExperimentalSingletMpfr, ExplicitPrecisionAndMixedPrecisionRejection)
{
  using R = uni20::mpreal;
  for (auto bits : {128, 256})
  {
    R const zero(0, uni20::Precision::bits(bits));
    check_analytic_singlet(zero);
    experiment::PreparedSinglet<R> const prepared(4, analytic_singlet(zero));
    ASSERT_TRUE(prepared.energy());
    EXPECT_EQ(prepared.energy()->precision(), zero.precision());
    EXPECT_EQ(prepared.working_zero()->precision(), zero.precision());
    auto ground = analytic_singlet(zero);
    ground.rapidities.back() = ground.rapidities.back().at(uni20::Precision::bits(bits + 32));
    EXPECT_THROW((experiment::PreparedSinglet<R>(4, ground)), std::invalid_argument);
    ground = analytic_singlet(zero);
    ground.rapidities.back() = R(1); // Exact is not the prepared numerical precision.
    EXPECT_THROW((experiment::PreparedSinglet<R>(4, ground)), std::invalid_argument);
    ground.rapidities.back() = R(uni20::uninitialized);
    EXPECT_THROW((experiment::PreparedSinglet<R>(4, ground)), std::invalid_argument);
  }
}
#endif
