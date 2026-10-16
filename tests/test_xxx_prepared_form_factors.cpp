// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "xxx_unprepared_form_factor.hpp"
#include <bethe/heisenberg_excitations.hpp>

namespace model = bethe::heisenberg;
template <typename T> class XXXPreparedFormFactors : public ::testing::Test {};
TYPED_TEST_SUITE(XXXPreparedFormFactors, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXXPreparedFormFactors, MatchesUnpreparedAcrossFamilies)
{
  using Real = TypeParam;
  using std::abs;
  Real const tol = Real{8192} * uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {2, 4, 6, 8, 12, 16, 32})
  {
    SCOPED_TRACE(n);
    auto const count = model::real_excitation_count(n, uni20::half_int{1});
    auto const scan = model::real_excitations<Real>(n, uni20::half_int{1}, {count});
    ASSERT_TRUE(scan.converged());
    model::PreparedRaisingFormFactors<Real> const prepared(n, scan.ground_state);
    ASSERT_TRUE(prepared.ready());
    EXPECT_EQ(prepared.sites(), n);
    ASSERT_TRUE(prepared.ground_residual());
    EXPECT_LE(*prepared.ground_residual(), prepared.root_tolerance());
    ASSERT_EQ(prepared.aba_rapidities().size(), n / 2);
    for (std::size_t j = 0; j < n / 2; ++j)
      test_support::expect_exact(prepared.aba_rapidities()[j], scan.ground_state.rapidities[j] / Real{2}, "ABA lambda");
    for (auto const& level : scan.levels)
    {
      auto const expected = test_support::unprepared_raising_form_factor(n, scan.ground_state, level.state);
      auto const actual = prepared(level.state);
      ASSERT_EQ(actual.status, expected.status);
      ASSERT_TRUE(actual.converged());
      ASSERT_TRUE(actual.weight);
      ASSERT_TRUE(expected.weight);
      ASSERT_TRUE(actual.log_weight);
      ASSERT_TRUE(expected.log_weight);
      EXPECT_REAL_NEAR(*actual.weight / *expected.weight, Real{1}, tol);
      EXPECT_REAL_NEAR(*actual.log_weight, *expected.log_weight, Real(tol * (Real{1} + abs(*expected.log_weight))));
      test_support::expect_exact(actual.pivot_margin, expected.pivot_margin, "unchanged determinant kernels");
    }
  }
}

TYPED_TEST(XXXPreparedFormFactors, OwnsGroundSnapshot)
{
  using Real = TypeParam;
  auto ground = model::ground_state<Real>(8);
  auto const excited = model::sector_ground_state<Real>(8, uni20::half_int{1});
  model::PreparedRaisingFormFactors<Real> const prepared(8, ground);
  auto const expected = prepared(excited);
  ASSERT_TRUE(expected.converged());
  ground.rapidities[0] += Real{0.01};
  ground.residual_norm = Real{0}; // Deliberately stale convergence metadata.
  model::PreparedRaisingFormFactors<Real> const stale(8, ground);
  EXPECT_FALSE(stale.ready());
  EXPECT_EQ(stale.status(), model::FormFactorStatus::precision_limit);
  ASSERT_TRUE(stale.ground_residual());
  EXPECT_GT(*stale.ground_residual(), stale.root_tolerance());
  EXPECT_FALSE(stale(excited).weight);
  ground = {}; // Release caller-owned roots; the prepared snapshot remains usable.
  auto const actual = prepared(excited);
  ASSERT_TRUE(actual.converged());
  test_support::expect_exact(*actual.weight, *expected.weight, "input mutation and destruction");
  auto copied = prepared;
  auto moved = std::move(copied);
  ASSERT_TRUE(moved.ready());
  test_support::expect_exact(*moved(excited).weight, *expected.weight, "value-owned context move");
  // A moved-from context is still safe to inspect/evaluate, but need not be ready.
  auto const after_move = copied(excited);
  EXPECT_EQ(after_move.converged(), copied.ready());
  if (copied.ready())
    test_support::expect_exact(*after_move.weight, *expected.weight, "valid moved-from context");
  else
    EXPECT_FALSE(after_move.weight);
}

TYPED_TEST(XXXPreparedFormFactors, RevalidatesEachExcitedState)
{
  using Real = TypeParam;
  model::PreparedRaisingFormFactors<Real> const prepared(8, model::ground_state<Real>(8));
  auto const valid = model::sector_ground_state<Real>(8, uni20::half_int{1});
  ASSERT_TRUE(prepared(valid).converged());
  auto excited = valid;
  excited.rapidities[0] += Real{0.01};
  excited.residual_norm = Real{0};
  auto check_rejected = [&](auto const& state, model::FormFactorStatus status) {
    auto const result = prepared(state);
    EXPECT_EQ(result.status, status);
    EXPECT_FALSE(result.weight);
    EXPECT_FALSE(result.log_weight);
  };
  check_rejected(excited, model::FormFactorStatus::precision_limit);
  excited = valid;
  excited.converged = false;
  check_rejected(excited, model::FormFactorStatus::roots_unconverged);
  for (Real bad : {uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    excited = valid;
    excited.rapidities[0] = bad;
    check_rejected(excited, model::FormFactorStatus::precision_limit);
  }
  EXPECT_TRUE(prepared(valid).converged()); // Rejections never poison the context.
}

TYPED_TEST(XXXPreparedFormFactors, RejectsInvalidGroundAndStateIdentity)
{
  using Real = TypeParam;
  using Prepared = model::PreparedRaisingFormFactors<Real>;
  auto const ground = model::ground_state<Real>(8);
  auto const excited = model::sector_ground_state<Real>(8, uni20::half_int{1});
  EXPECT_THROW((void)Prepared(0, ground), std::invalid_argument);
  EXPECT_THROW((void)Prepared(7, ground), std::invalid_argument);
  EXPECT_THROW((void)Prepared(10, ground), std::invalid_argument);
  auto bad = ground;
  bad.sz = uni20::half_int{1};
  EXPECT_THROW((void)Prepared(8, bad), std::invalid_argument);
  bad = ground;
  bad.spin_reversed = true;
  EXPECT_THROW((void)Prepared(8, bad), std::invalid_argument);
  bad = ground;
  bad.quantum_numbers[0] = bad.quantum_numbers[1];
  EXPECT_THROW((void)Prepared(8, bad), std::invalid_argument);
  bad = ground;
  bad.converged = false;
  Prepared const unconverged(8, bad);
  EXPECT_EQ(unconverged.status(), model::FormFactorStatus::roots_unconverged);
  EXPECT_FALSE(unconverged.ground_residual());
  EXPECT_FALSE(unconverged(excited).weight);
  bad = ground;
  bad.rapidities[0] = uni20::numeric_limits<Real>::quiet_NaN();
  Prepared const nonfinite(8, bad);
  EXPECT_EQ(nonfinite.status(), model::FormFactorStatus::precision_limit);
  EXPECT_FALSE(nonfinite.ground_residual());
  EXPECT_FALSE(nonfinite(excited).weight);
  bad = ground;
  for (auto& z : bad.rapidities)
    z /= Real{2};
  EXPECT_FALSE(Prepared(8, bad).ready()); // ABA roots cannot be passed as solver z.

  Prepared const prepared(8, ground);
  bad = excited;
  bad.sz = uni20::half_int{0};
  EXPECT_THROW((void)prepared(bad), std::invalid_argument);
  bad = excited;
  bad.spin_reversed = true;
  EXPECT_THROW((void)prepared(bad), std::invalid_argument);
  bad = excited;
  bad.quantum_numbers[0] = bad.quantum_numbers[1];
  EXPECT_THROW((void)prepared(bad), std::invalid_argument);
  bad = excited;
  bad.rapidities.pop_back();
  EXPECT_THROW((void)prepared(bad), std::invalid_argument);
  bad = excited;
  bad.quantum_numbers.pop_back();
  EXPECT_THROW((void)prepared(bad), std::invalid_argument);
}

TYPED_TEST(XXXPreparedFormFactors, GaudinCoordinateChainRule)
{
  using Real = TypeParam;
  using std::log;
  for (std::size_t n : {2, 8, 16})
  {
    SCOPED_TRACE(n);
    auto const ground = model::ground_state<Real>(n);
    model::PreparedRaisingFormFactors<Real> const prepared(n, ground);
    ASSERT_TRUE(prepared.ready());
    auto const m = n / 2;
    auto const& z = ground.rapidities;
    // Independently differentiate F_i=2N atan(z_i)-2 sum_j atan((z_i-z_j)/2)-2pi I_i.
    // No division by N: G_lambda = 2 G_z for these SAME unscaled equations.
    std::vector<Real> gz(m * m);
    for (std::size_t a = 0; a < m; ++a)
    {
      gz[a * m + a] = Real{2} * Real(n) / (Real{1} + z[a] * z[a]);
      for (std::size_t b = 0; b < m; ++b)
        if (a != b)
        {
          Real const d = z[a] - z[b];
          gz[a * m + b] = Real{4} / (Real{4} + d * d);
          gz[a * m + a] -= gz[a * m + b];
        }
    }
    auto const dz = bethe::detail::log_absolute_determinant<Real>(std::move(gz), m);
    auto const dl = model::detail::gaudin_log_determinant<Real>(n, prepared.aba_rapidities());
    ASSERT_TRUE(dz.log_absolute);
    ASSERT_TRUE(dl.log_absolute);
    EXPECT_REAL_NEAR(Real(*dl.log_absolute - *dz.log_absolute), Real(Real(m) * log(Real{2})),
                     Real(Real{128} * Real(m) * uni20::numeric_limits<Real>::epsilon()));
  }
}
