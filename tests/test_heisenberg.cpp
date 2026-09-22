// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include "test_support.hpp"
#include <bethe/heisenberg.hpp>

#include <string>
#include <string_view>

namespace
{

template <uni20::Real Real>
void check_diagnostics(std::size_t sites, bethe::heisenberg::GroundState<Real> const& result)
{
  using std::abs;
  using std::atan;
  Real const pi = Real{4} * atan(Real{1});
  Real residual = Real{0};
  Real energy = static_cast<Real>(sites) / Real{4};
  for (std::size_t i = 0; i < result.rapidities.size(); ++i)
  {
    SCOPED_TRACE(::testing::Message() << " i=" << i);
    Real const z = result.rapidities[i];
    Real const quantum_number = static_cast<Real>(i) - static_cast<Real>(sites) / Real{4} + Real{1} / Real{2};
    // Direct evaluation of N*phi - 2*pi*I - sum(phi), independently of the
    // solver's fixed-point angle calculation and compensated accumulation.
    Real equation = static_cast<Real>(sites) * Real{2} * atan(z) - Real{2} * pi * quantum_number;
    for (std::size_t j = 0; j < result.rapidities.size(); ++j)
      if (i != j) equation -= Real{2} * atan((z - result.rapidities[j]) / Real{2});
    residual = std::max(residual, abs(equation) / static_cast<Real>(sites));
    energy -= Real{2} / (Real{1} + z * z);
  }
  Real const rounding = Real{4} * static_cast<Real>(sites) * uni20::numeric_limits<Real>::epsilon();
  ASSERT_TRUE(abs(residual - result.residual_norm) <= rounding) << "residual does not describe the returned roots";
  ASSERT_TRUE(abs(energy - result.energy) <= rounding) << "energy does not describe the returned roots";
}

template <typename Real> class Heisenberg : public ::testing::Test {};
TYPED_TEST_SUITE(Heisenberg, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(Heisenberg, SmallChains)
{
  using Real = TypeParam;

  using bethe::heisenberg::ground_state;
  using std::sqrt;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const tolerance = Real{2048} * eps;

  auto const two = ground_state<Real>(2);
  ASSERT_TRUE(two.converged && two.iterations == 0) << "N=2 initial roots should already be exact";
  ASSERT_TRUE(two.rapidities.size() == 1 && two.rapidities[0] == Real{0}) << "N=2 rapidity";
  ASSERT_TRUE(two.energy == -Real{3} / Real{2} && two.residual_norm == Real{0}) << "N=2 double-bond normalization";

  auto const four = ground_state<Real>(4);
  ASSERT_TRUE(four.converged) << "N=4 convergence";
  ASSERT_EQ(four.rapidities.size(), 2u);
  EXPECT_REAL_NEAR(four.energy + Real{2}, Real{0}, tolerance) << "N=4 exact energy";
  EXPECT_REAL_NEAR(four.rapidities[0] + Real{1} / sqrt(Real{3}), Real{0}, tolerance) << "N=4 exact negative root";
  EXPECT_REAL_NEAR(four.rapidities[1], Real{1} / sqrt(Real{3}), tolerance) << "N=4 exact positive root";
}

TYPED_TEST(Heisenberg, NativePrecisionSix)
{
  using Real = TypeParam;

  using bethe::heisenberg::ground_state;
  using std::abs;
  using std::sqrt;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const tolerance = Real{2048} * eps;

  // At N=6 the roots are {-a,0,a}; a^2 solves 3*a^4+10*a^2-9=0.
  // This gives E=-(2+sqrt(13))/2, an irrational precision-sensitive oracle.
  auto const six = ground_state<Real>(6);
  ASSERT_EQ(six.rapidities.size(), 3u);
  Real const exact_six = -(Real{2} + sqrt(Real{13})) / Real{2};
  Real const a = sqrt((Real{2} * sqrt(Real{13}) - Real{5}) / Real{3});
  Real const precision_tolerance = Real{128} * eps;

  ASSERT_TRUE(six.converged) << "N=6 exact energy at selected precision";
  ASSERT_REAL_NEAR(six.energy, exact_six, precision_tolerance) << "N=6 exact energy at selected precision";

  ASSERT_REAL_NEAR(six.rapidities[0] + a, Real{0}, tolerance) << "N=6 exact outer roots";
  ASSERT_REAL_NEAR(six.rapidities[2], a, tolerance) << "N=6 exact outer roots";

  EXPECT_REAL_NEAR(six.rapidities[1], Real{0}, tolerance) << "N=6 central root";
  if constexpr (uni20::numeric_limits<Real>::digits > uni20::numeric_limits<double>::digits)
  {
    ASSERT_TRUE(abs(static_cast<Real>(static_cast<double>(exact_six)) - exact_six) > precision_tolerance)
        << "high-precision oracle must distinguish a double-only implementation";
  }
  test_support::expect_exact(uni20::parse_real<Real>(uni20::format_real(six.energy)), six.energy,
                             "precision-preserving text round trip");
}

TYPED_TEST(Heisenberg, GroundDiagnostics)
{
  using Real = TypeParam;

  using bethe::heisenberg::ground_state;
  using Options = bethe::heisenberg::SolverOptions<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const tolerance = Real{2048} * eps;

  for (std::size_t const sites : {4, 6, 16, 64})
  {
    SCOPED_TRACE(::testing::Message() << " sites=" << sites);
    auto const result = ground_state<Real>(sites);
    ASSERT_TRUE(result.converged && result.residual_norm <= Options{}.residual_tolerance) << "equation convergence";
    ASSERT_TRUE(result.iterations <= Options{}.max_iterations) << "iteration budget";
    ASSERT_TRUE(result.rapidities.size() == sites / 2) << "ground-state root count";
    ASSERT_NO_FATAL_FAILURE(check_diagnostics(sites, result));
    for (std::size_t i = 0; i < result.rapidities.size(); ++i)
    {
      SCOPED_TRACE(::testing::Message() << " i=" << i);
      if (i != 0)
      {
        ASSERT_TRUE(result.rapidities[i - 1] < result.rapidities[i]) << "roots must be ordered";
      }
      EXPECT_REAL_NEAR(result.rapidities[i] + result.rapidities[result.rapidities.size() - 1 - i], Real{0}, tolerance)
          << "ground-state reflection symmetry";
    }
  }
}

TYPED_TEST(Heisenberg, PublishedEnergy)
{
  using Real = TypeParam;

  using bethe::heisenberg::ground_state;

  // Independent published finite-size value: Table I of cond-mat/9809163,
  // relative to the ferromagnetic energy. Its decimal precision limits the test.
  auto const sixteen = ground_state<Real>(16);
  Real const reference = uni20::parse_real<Real>("-0.696393522538549");
  EXPECT_REAL_NEAR(sixteen.energy / Real{16} - Real{1} / Real{4} - reference, Real{0}, uni20::parse_real<Real>("2e-14"))
      << "N=16 published energy convention";
}

TYPED_TEST(Heisenberg, CrossPrecision)
{
  using Real = TypeParam;

  using bethe::heisenberg::ground_state;

  auto const baseline = ground_state<double>(32);
  auto const selected = ground_state<Real>(32);
  ASSERT_TRUE(baseline.converged && selected.converged) << "cross-precision convergence";
  EXPECT_REAL_NEAR(selected.energy, static_cast<Real>(baseline.energy),
                   Real{2048} * static_cast<Real>(uni20::numeric_limits<double>::epsilon()))
      << "cross-precision energy agreement";
}

TYPED_TEST(Heisenberg, IterationBudget)
{
  using Real = TypeParam;

  using bethe::heisenberg::ground_state;
  using std::sqrt;
  using Options = bethe::heisenberg::SolverOptions<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const tolerance = Real{2048} * eps;

  auto const no_updates = ground_state<Real>(4, Options{.max_iterations = 0});
  ASSERT_TRUE(!no_updates.converged && no_updates.iterations == 0) << "zero update budget";
  ASSERT_TRUE(no_updates.energy == -Real{3}) << "initial-guess energy";
  ASSERT_NO_FATAL_FAILURE(check_diagnostics(4, no_updates));
  auto const one_update = ground_state<Real>(4, Options{.max_iterations = 1});
  ASSERT_TRUE(!one_update.converged && one_update.iterations == 1) << "iteration exhaustion must not look converged";
  EXPECT_REAL_NEAR(one_update.energy + Real{1} + sqrt(Real{2}), Real{0}, tolerance)
      << "one simultaneous update, not two";
  ASSERT_NO_FATAL_FAILURE(check_diagnostics(4, one_update));
}

TYPED_TEST(Heisenberg, InvalidInputs)
{
  using Real = TypeParam;

  using bethe::heisenberg::ground_state;
  using Options = bethe::heisenberg::SolverOptions<Real>;

  for (std::size_t const sites : {std::size_t{0}, std::size_t{1}, std::numeric_limits<std::size_t>::max()})
    EXPECT_THROW(([&] { (void)ground_state<Real>(sites); })(), std::invalid_argument);
  for (Real const bad_tolerance :
       {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW(([&] { (void)ground_state<Real>(4, Options{.residual_tolerance = bad_tolerance}); })(),
                 std::invalid_argument);
}

} // namespace
