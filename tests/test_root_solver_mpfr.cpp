// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/heisenberg_excitations.hpp>
#include <bethe/lieb_liniger_open.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using Real = uni20::mpreal;
namespace xxx = bethe::heisenberg::open;
namespace ll = bethe::lieb_liniger;
namespace math = uni20::math;

class RootSolverMpfr : public ::testing::TestWithParam<int> {
  protected:
    uni20::Precision precision() const { return uni20::Precision::bits(GetParam()); }
    Real value(int n) const { return Real(n, precision()); }
    Real epsilon() const { return uni20::numeric_limits<Real>::epsilon(precision()); }
    bethe::SolverOptions<Real> options() const { return bethe::SolverOptions<Real>(precision()); }
};
INSTANTIATE_TEST_SUITE_P(WorkingPrecisions, RootSolverMpfr, ::testing::Values(128, 256, 333));

void check_precision(Real const& value, uni20::Precision p)
{
  ASSERT_TRUE(value.initialized());
  EXPECT_FALSE(value.is_exact());
  EXPECT_EQ(value.precision(), p);
  EXPECT_TRUE(uni20::isfinite(value));
}
template <typename State> void check_state(State const& state, uni20::Precision p)
{
  check_precision(state.energy, p);
  check_precision(state.residual_norm, p);
  if constexpr (requires { state.rapidities; })
    for (auto const& root : state.rapidities)
      check_precision(root, p);
  else
  {
    check_precision(state.length, p);
    check_precision(state.interaction, p);
    for (auto const& root : state.momenta)
      check_precision(root, p);
  }
  if constexpr (requires { state.momentum; }) check_precision(state.momentum, p);
}

TEST_P(RootSolverMpfr, OpenXXXAnalyticEnergySectorsAndEmptyRoots)
{
  auto const s = xxx::ground_state<Real>(4, options());
  ASSERT_TRUE(s.converged);
  check_state(s, precision());
  EXPECT_REAL_NEAR(s.energy, -value(3) / 4 - math::sqrt(value(3)) / 2, value(128) * epsilon());
  auto const three = xxx::ground_state<Real>(3, options());
  ASSERT_TRUE(three.converged);
  EXPECT_REAL_NEAR(three.rapidities.front(), 1 / math::sqrt(value(3)), value(32) * epsilon());
  for (std::size_t n : {4, 5})
  {
    auto const sectors = xxx::sector_ground_states<Real>(n, options());
    ASSERT_EQ(sectors.size(), n + 1);
    for (std::size_t i = 0; i <= n; ++i)
    {
      ASSERT_TRUE(sectors[i].converged);
      check_state(sectors[i], precision());
      EXPECT_EQ(sectors[i].energy, sectors[n - i].energy);
      EXPECT_EQ(sectors[i].sz, -sectors[n - i].sz);
    }
  }
  auto const vacuum = xxx::solve_real<Real>(4, {}, options());
  check_state(vacuum, precision());
  EXPECT_EQ(vacuum.energy, value(3) / 4);
  EXPECT_EQ(vacuum.residual_norm, value(0));
  auto limited = options();
  limited.max_iterations = 0;
  auto const seed = xxx::ground_state<Real>(4, limited);
  EXPECT_FALSE(seed.converged);
  EXPECT_EQ(seed.iterations, 0u);
  check_state(seed, precision());
}

#if UNI20_HAS_MPLAPACK_MPREAL
// Independent scalar bisection of the TWO-body wall equations for the root
// sum and difference: x+2 atan(x/g)=label*pi. No Newton matrix or model seed.
Real wall_two_body(Real const& g, int label)
{
  auto const p = g.precision();
  auto const pi = uni20::pi<Real>.at(p);
  Real lo = (label - 1) * pi, hi = label * pi;
  for (int j = 0; j < uni20::numeric_limits<Real>::digits(p) + 4; ++j)
  {
    Real const mid = (lo + hi) / 2;
    if (mid == lo || mid == hi) break;
    if (mid + 2 * math::atan(mid / g) > label * pi)
      hi = mid;
    else
      lo = mid;
  }
  return (lo + hi) / 2;
}

TEST_P(RootSolverMpfr, LiebLinigerAnalyticAndIndependentTwoBodyReferences)
{
  auto const pi = uni20::pi<Real>.at(precision());
  auto const ring = ll::ground_state(2, value(1), pi, options());
  ASSERT_TRUE(ring.converged);
  check_state(ring, precision());
  EXPECT_REAL_NEAR(ring.energy, pi * pi / 2, value(2048) * epsilon());
  EXPECT_REAL_NEAR(ring.momenta[1], pi / 2, value(256) * epsilon());
  auto const a = 3 * pi / 2 - 2 * math::atan(value(2));
  auto const three = ll::ground_state(3, value(1), a, options());
  ASSERT_TRUE(three.converged);
  check_state(three, precision());
  EXPECT_REAL_NEAR(three.energy, 2 * a * a, value(2048) * epsilon());
  for (auto const& g : {value(1) / 100, value(1), value(100)})
  {
    auto const walls = ll::open::ground_state(2, value(1), g, options());
    ASSERT_TRUE(walls.converged);
    check_state(walls, precision());
    auto const sum = wall_two_body(g, 3), difference = wall_two_body(g, 1);
    EXPECT_REAL_NEAR(walls.momenta[0], (sum - difference) / 2, value(4096) * epsilon());
    EXPECT_REAL_NEAR(walls.momenta[1], (sum + difference) / 2, value(4096) * epsilon());
    EXPECT_REAL_NEAR(walls.energy, (sum * sum + difference * difference) / 2, value(65536) * epsilon());
  }
  for (std::size_t n : {0, 1, 4})
  {
    auto limited = options();
    limited.max_iterations = 0;
    auto const pbc = ll::ground_state(n, value(3), value(2), limited);
    auto const obc = ll::open::ground_state(n, value(3), value(2), limited);
    check_state(pbc, precision());
    check_state(obc, precision());
    EXPECT_EQ(pbc.converged, n <= 1);
    EXPECT_EQ(obc.converged, n <= 1);
    EXPECT_EQ(pbc.iterations, 0u);
    EXPECT_EQ(obc.iterations, 0u);
    if (n == 0) EXPECT_EQ(obc.energy, value(0));
    if (n == 1) EXPECT_REAL_NEAR(obc.energy, pi * pi / 9, value(16) * epsilon());
  }
}

#endif

TEST_P(RootSolverMpfr, RejectsImplicitPrecisionChanges)
{
  auto const other = uni20::Precision::bits(GetParam() + 32);
  for (auto const& invalid : {Real{1}, Real(1, other), Real(uni20::uninitialized)})
  {
    auto o = options();
#if UNI20_HAS_MPLAPACK_MPREAL
    EXPECT_THROW((void)ll::ground_state(2, invalid, value(1), o), std::invalid_argument);
    EXPECT_THROW((void)ll::open::ground_state(2, value(1), invalid, o), std::invalid_argument);
    auto ring = ll::ground_state(2, value(1), value(1), o);
    ring.momenta[0] = invalid;
    EXPECT_THROW((void)ll::solve_real<Real>(value(1), value(1), ring.quantum_numbers, o, ring.momenta),
                 std::invalid_argument);
    auto wall = ll::open::ground_state(2, value(1), value(1), o);
    wall.momenta[0] = invalid;
    EXPECT_THROW((void)ll::open::solve_real<Real>(value(1), value(1), wall.quantum_numbers, o, wall.momenta),
                 std::invalid_argument);
#endif
    auto chain = xxx::ground_state<Real>(4, o);
    chain.rapidities[0] = invalid;
    EXPECT_THROW((void)xxx::solve_real<Real>(4, chain.quantum_numbers, o, chain.rapidities), std::invalid_argument);
    o.residual_tolerance = invalid;
    EXPECT_THROW((void)xxx::ground_state<Real>(4, o), std::invalid_argument);
#if UNI20_HAS_MPLAPACK_MPREAL
    EXPECT_THROW((void)ll::ground_state(2, value(1), value(1), o), std::invalid_argument);
    EXPECT_THROW((void)ll::open::ground_state(2, value(1), value(1), o), std::invalid_argument);
#endif
  }
}

template <typename Scan> void compare_scan(Scan const& serial, Scan const& parallel, uni20::Precision p)
{
  ASSERT_TRUE(serial.converged());
  ASSERT_TRUE(parallel.converged());
  EXPECT_EQ(serial.candidate_count, parallel.candidate_count);
  ASSERT_EQ(serial.levels.size(), parallel.levels.size());
  for (std::size_t i = 0; i < serial.levels.size(); ++i)
  {
    auto const& a = serial.levels[i];
    auto const& b = parallel.levels[i];
    check_state(b.state, p);
    EXPECT_EQ(a.state.quantum_numbers, b.state.quantum_numbers);
    EXPECT_EQ(a.state.energy, b.state.energy);
    EXPECT_EQ(a.gap, b.gap);
    ASSERT_TRUE(b.gap);
    check_precision(*b.gap, p);
  }
}

TEST(RootSolverMpfrShared, AlternatingPrecisionsAndThreadedScans)
{
  auto const ambient = mpfr_get_default_prec();
  uni20::async::TbbScheduler scheduler(2);
  bethe::RealExcitationOptions serial{.count = 1000, .max_candidates = 1000, .execution = {.batch_size = 3}};
  auto parallel = serial;
  parallel.execution.scheduler = &scheduler;
  for (int bits : {128, 256, 128})
  {
    auto const p = uni20::Precision::bits(bits);
    bethe::SolverOptions<Real> options(p);
    compare_scan(xxx::real_excitations<Real>(6, uni20::half_int{1}, serial, options),
                 xxx::real_excitations<Real>(6, uni20::half_int{1}, parallel, options), p);
#if UNI20_HAS_MPLAPACK_MPREAL
    compare_scan(ll::real_excitations(2, Real(2, p), Real(1, p), 1, serial, options),
                 ll::real_excitations(2, Real(2, p), Real(1, p), 1, parallel, options), p);
    compare_scan(ll::open::real_excitations(2, Real(2, p), Real(1, p), 2, serial, options),
                 ll::open::real_excitations(2, Real(2, p), Real(1, p), 2, parallel, options), p);
#endif
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

#if UNI20_HAS_MPLAPACK_MPREAL
TEST(RootSolverMpfrShared, WeakWallCouplingNeedsMoreThanBinary128)
{
  // sqrt(c) is below binary64 root spacing; the O(c) energy shift is below
  // binary128 spacing around the free energy.
  // Check that 333 bits actually resolves it, not just the O(1) free energy.
  auto const p = uni20::Precision::bits(333);
  Real const c("1e-50", p), one(1, p);
  auto const pi = uni20::pi<Real>.at(p);
  auto const s = ll::open::ground_state(2, one, c, bethe::SolverOptions<Real>(p));
  ASSERT_TRUE(s.converged);
  check_state(s, p);
  EXPECT_GT(s.momenta[1], s.momenta[0]);
  EXPECT_REAL_NEAR((s.energy - 2 * pi * pi) / c, Real(3, p), Real("1e-45", p));
  auto const high = uni20::Precision::bits(512);
  auto const reference = ll::open::ground_state(2, one.at(high), c.at(high), bethe::SolverOptions<Real>(high));
  ASSERT_TRUE(reference.converged);
  EXPECT_REAL_NEAR(s.energy, reference.energy.at(p), 4096 * uni20::numeric_limits<Real>::epsilon(p));
  auto const ring = ll::ground_state(2, one, c, bethe::SolverOptions<Real>(p));
  ASSERT_TRUE(ring.converged);
  EXPECT_REAL_NEAR(ring.energy / c, Real(2, p), Real("1e-48", p));
  EXPECT_THROW((void)ll::open::ground_state(2, 1.0, 1e-50), std::invalid_argument);
}

TEST(RootSolverMpfrShared, NontrivialStatesImproveWithPrecisionAndRestart)
{
  auto const high = uni20::Precision::bits(512);
  for (auto const& text : {"1e-12", "1", "1e12"})
  {
    SCOPED_TRACE(text);
    auto const ring = ll::ground_state(8, Real(8, high), Real(text, high), bethe::SolverOptions<Real>(high));
    auto const walls = ll::open::ground_state(8, Real(8, high), Real(text, high), bethe::SolverOptions<Real>(high));
    ASSERT_TRUE(ring.converged);
    ASSERT_TRUE(walls.converged);
    for (int bits : {128, 256})
    {
      auto const p = uni20::Precision::bits(bits);
      bethe::SolverOptions<Real> options(p);
      auto const a = ll::ground_state(8, Real(8, p), Real(text, p), options);
      auto const b = ll::open::ground_state(8, Real(8, p), Real(text, p), options);
      ASSERT_TRUE(a.converged);
      ASSERT_TRUE(b.converged);
      auto const bound = 65536 * uni20::numeric_limits<Real>::epsilon(p);
      EXPECT_REAL_NEAR(a.energy / ring.energy.at(p), Real(1, p), bound);
      EXPECT_REAL_NEAR(b.energy / walls.energy.at(p), Real(1, p), bound);
      options.max_iterations = 0;
      auto const restarted_a = ll::solve_real<Real>(a.length, a.interaction, a.quantum_numbers, options, a.momenta);
      auto const restarted_b =
          ll::open::solve_real<Real>(b.length, b.interaction, b.quantum_numbers, options, b.momenta);
      EXPECT_TRUE(restarted_a.converged);
      EXPECT_TRUE(restarted_b.converged);
      EXPECT_EQ(a.energy, restarted_a.energy);
      EXPECT_EQ(b.energy, restarted_b.energy);
    }
  }
}
#endif
} // namespace
