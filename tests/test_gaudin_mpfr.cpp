// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/central_spin.hpp>
#include <bethe/richardson.hpp>

namespace
{
using Real = uni20::mpreal;
namespace pairing = bethe::richardson;
namespace central = bethe::central_spin;

class GaudinMpfr : public ::testing::TestWithParam<int> {
  protected:
    uni20::Precision precision() const { return uni20::Precision::bits(GetParam()); }
    Real value(int n) const { return Real(n, precision()); }
    Real epsilon() const { return uni20::numeric_limits<Real>::epsilon(precision()); }
    std::vector<Real> levels() const { return {value(0), value(1), value(2), value(3)}; }
};
INSTANTIATE_TEST_SUITE_P(WorkingPrecisions, GaudinMpfr, ::testing::Values(128, 256, 333));

void check_precision(Real const& v, uni20::Precision p)
{
  EXPECT_TRUE(v.initialized());
  EXPECT_FALSE(v.is_exact());
  EXPECT_EQ(v.precision(), p);
  EXPECT_TRUE(uni20::isfinite(v));
}
void check_state(pairing::State<Real> const& s, uni20::Precision p)
{
  for (auto const& v : s.levels)
    check_precision(v, p);
  for (auto const& v : s.eigenvalue_variables)
    check_precision(v, p);
  for (auto const& v :
       {s.coupling, s.reached_coupling, s.energy, s.residual_norm, s.target_residual_norm, s.particle_number_error})
    check_precision(v, p);
}
void check_state(central::State<Real> const& s, uni20::Precision p)
{
  for (auto const& v : s.couplings)
    check_precision(v, p);
  for (auto const& v : s.eigenvalue_variables)
    check_precision(v, p);
  for (auto const& v : {s.field, s.continuation_parameter, s.residual_norm, s.number_error})
    check_precision(v, p);
  if (s.energy) check_precision(*s.energy, p);
  if (s.reached_field) check_precision(*s.reached_field, p);
}

TEST_P(GaudinMpfr, PairingAnalyticEnergyAndRegularPairCollision)
{
  pairing::SolverOptions<Real> options(precision());
  auto const one = pairing::ground_state<Real>(std::vector<Real>{value(0), value(1)}, 1, value(1), {}, options);
  ASSERT_TRUE(one.converged);
  check_state(one, precision());
  EXPECT_REAL_NEAR(one.energy, -uni20::math::sqrt(value(2)), value(4096) * epsilon());

  auto const collision = pairing::ground_state<Real>(levels(), 2, value(2) / value(3), {}, options);
  ASSERT_TRUE(collision.converged);
  check_state(collision, precision());
  EXPECT_REAL_NEAR(collision.energy, value(0), value(4096) * epsilon());
  std::vector<Real> const expected{value(7) / value(9), value(2) / value(3), value(1) / value(3), value(2) / value(9)};
  for (std::size_t i = 0; i < expected.size(); ++i)
    EXPECT_REAL_NEAR(collision.eigenvalue_variables[i], expected[i], value(1024) * epsilon());
}

TEST_P(GaudinMpfr, CentralSpinAnalyticEnergyIncludingExtremeFields)
{
  central::SolverOptions<Real> options(precision());
  // The last field makes t round to one, but the separately supplied p and
  // reached physical field must not disappear. No binary64 reference used.
  for (auto const& field : {value(1), value(-1), value(0), Real("1e10", precision()), epsilon() * epsilon()})
  {
    auto const s = central::sector_ground_state<Real>(std::vector<Real>{value(1)}, field, uni20::half_int{0}, options);
    ASSERT_TRUE(s.converged);
    check_state(s, precision());
    ASSERT_TRUE(s.energy);
    ASSERT_TRUE(s.reached_field);
    EXPECT_EQ(*s.reached_field, field);
    Real const expected = -value(1) / value(4) - uni20::math::hypot(value(1), field) / value(2);
    EXPECT_REAL_NEAR(*s.energy, expected, value(8192) * epsilon() * (value(1) + uni20::math::abs(expected)));
  }
}

TEST_P(GaudinMpfr, AnalyticBranchesAndIncompleteResultsKeepWorkingPrecision)
{
  pairing::SolverOptions<Real> po(precision());
  auto check_pairing = [&](std::size_t pairs, Real g, std::vector<std::size_t> const& blocked = {}) {
    auto const s = pairing::ground_state<Real>(levels(), pairs, g, blocked, po);
    check_state(s, precision());
    return s;
  };
  EXPECT_EQ(check_pairing(0, value(1)).energy, value(0));
  EXPECT_EQ(check_pairing(4, value(1)).energy, value(8));
  EXPECT_EQ(check_pairing(1, value(0), {0, 2}).energy, value(4));
  EXPECT_EQ(check_pairing(0, value(1), {0, 1, 2, 3}).energy, value(6));
  po.max_iterations = 0;
  auto const seed = check_pairing(2, value(1));
  EXPECT_FALSE(seed.converged);
  EXPECT_EQ(seed.reached_coupling, value(0));
  EXPECT_EQ(seed.energy, value(2));
  EXPECT_EQ(seed.iterations, 0u);
  po.max_iterations = 10000;
  po.max_stages = 1;
  auto const partial = check_pairing(2, value(1));
  EXPECT_FALSE(partial.converged);
  EXPECT_EQ(partial.stages, 1u);
  EXPECT_EQ(partial.status, pairing::SolveStatus::stage_limit);

  central::SolverOptions<Real> co(precision());
  auto const isolated = central::sector_ground_state<Real>({}, value(3), uni20::half_int{-0.5}, co);
  check_state(isolated, precision());
  EXPECT_EQ(*isolated.energy, -value(3) / value(2));
  std::vector<Real> const a{value(1), value(7) / value(10), value(3) / value(10)};
  auto const polarized = central::sector_ground_state<Real>(a, value(1), uni20::half_int{2}, co);
  check_state(polarized, precision());
  EXPECT_EQ(*polarized.energy, value(1));
  co.max_iterations = 0;
  auto const cs = central::sector_ground_state<Real>(a, value(1), uni20::half_int{0}, co);
  check_state(cs, precision());
  EXPECT_FALSE(cs.converged);
  EXPECT_FALSE(cs.energy);
  EXPECT_FALSE(cs.reached_field);
  EXPECT_EQ(cs.iterations, 0u);
  co.max_iterations = 10000;
  co.max_stages = 1;
  auto const cp = central::sector_ground_state<Real>(a, value(-1), uni20::half_int{0}, co);
  check_state(cp, precision());
  EXPECT_FALSE(cp.converged);
  EXPECT_TRUE(cp.energy);
  EXPECT_TRUE(cp.reached_field);
  EXPECT_EQ(cp.status, central::SolveStatus::stage_limit);
}

TEST_P(GaudinMpfr, RejectsImplicitConversionsAndInvalidControls)
{
  pairing::SolverOptions<Real> po(precision());
  central::SolverOptions<Real> co(precision());
  auto const other = uni20::Precision::bits(GetParam() + 32);
  for (auto const& invalid : {Real{1}, Real(1, other), Real(uni20::uninitialized)})
  {
    EXPECT_THROW((void)pairing::ground_state<Real>(levels(), 2, invalid, {}, po), std::invalid_argument);
    auto bad_levels = levels();
    bad_levels[2] = invalid;
    EXPECT_THROW((void)pairing::ground_state<Real>(bad_levels, 2, value(1), {}, po), std::invalid_argument);
    EXPECT_THROW((void)central::sector_ground_state<Real>(std::vector<Real>{value(1)}, invalid, uni20::half_int{0}, co),
                 std::invalid_argument);
    EXPECT_THROW((void)central::sector_ground_state<Real>(std::vector<Real>{invalid}, value(1), uni20::half_int{0}, co),
                 std::invalid_argument);
    po.residual_tolerance = co.residual_tolerance = invalid;
    EXPECT_THROW((void)pairing::ground_state<Real>(levels(), 2, value(1), {}, po), std::invalid_argument);
    EXPECT_THROW(
        (void)central::sector_ground_state<Real>(std::vector<Real>{value(1)}, value(1), uni20::half_int{0}, co),
        std::invalid_argument);
    po = pairing::SolverOptions<Real>(precision());
    co = central::SolverOptions<Real>(precision());
  }
  EXPECT_THROW((void)pairing::SolverOptions<Real>(uni20::Precision::exact()), std::invalid_argument);
  EXPECT_THROW((void)central::SolverOptions<Real>(uni20::Precision::exact()), std::invalid_argument);
  static_assert(!std::is_default_constructible_v<pairing::SolverOptions<Real>>);
  static_assert(!std::is_default_constructible_v<central::SolverOptions<Real>>);
}

TEST(GaudinMpfrShared, NontrivialSolutionsAgreeAcrossPrecisionsWithoutChangingDefaults)
{
  auto const ambient = mpfr_get_default_prec();
  auto solve = [](int bits) {
    auto const p = uni20::Precision::bits(bits);
    auto v = [p](int n) { return Real(n, p); };
    auto pairing_state = pairing::ground_state<Real>(std::vector<Real>{v(0), v(1), v(2), v(3), v(4), v(5)}, 3,
                                                     v(3) / v(7), {}, pairing::SolverOptions<Real>(p));
    auto central_state =
        central::sector_ground_state<Real>(std::vector<Real>{v(3) / v(10), v(1), -v(7) / v(10)}, -v(2) / v(3),
                                           uni20::half_int{0}, central::SolverOptions<Real>(p));
    EXPECT_TRUE(pairing_state.converged);
    EXPECT_TRUE(central_state.converged);
    check_state(pairing_state, p);
    check_state(central_state, p);
    return std::pair{pairing_state.energy, *central_state.energy};
  };
  auto const low = solve(256), high = solve(333), repeated = solve(256);
  EXPECT_EQ(low, repeated);
  auto const p = uni20::Precision::bits(333);
  auto const tolerance = Real("1e-70", p);
  EXPECT_REAL_NEAR(low.first.at(p), high.first, tolerance);
  EXPECT_REAL_NEAR(low.second.at(p), high.second, tolerance);
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(GaudinMpfrShared, BacktrackingUsesRequestedBitsAndDoesNotCommitRejections)
{
  auto const p = uni20::Precision::bits(256);
  std::vector<Real> x{Real(0, p)};
  unsigned attempts = 0;
  ASSERT_TRUE(bethe::detail::backtrack_newton(
      x, p, [&](std::size_t) { return Real(1, p); },
      [&](auto const&, Real damping) {
        check_precision(damping, p);
        return ++attempts == 201;
      }));
  EXPECT_EQ(attempts, 201u);
  EXPECT_EQ(x[0], uni20::math::ldexp(Real(1, p), -200));
  auto const saved = x;
  attempts = 0;
  EXPECT_FALSE(bethe::detail::backtrack_newton(
      x, p, [&](std::size_t) { return Real(1, p); },
      [&](auto const&, Real) {
        ++attempts;
        return false;
      }));
  EXPECT_EQ(attempts, 257u);
  EXPECT_EQ(x, saved);
}
} // namespace
