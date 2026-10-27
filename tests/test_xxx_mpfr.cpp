// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/heisenberg.hpp>
#include <bethe/heisenberg_excitations.hpp>
#include <gtest/gtest.h>
#include <uni20/async/tbb_scheduler.hpp>
#if UNI20_HAS_MPREAL && UNI20_HAS_MPLAPACK_MPREAL
#include <bethe/xxx_structure_factor.hpp>
#endif

namespace
{
using Real = uni20::mpreal;
using uni20::Precision;
namespace model = bethe::heisenberg;
void near(Real const& actual, Real const& expected, int units = 4096)
{
  EXPECT_EQ(actual.precision(), expected.precision());
  EXPECT_TRUE(uni20::math::abs(actual - expected) <= units * uni20::numeric_limits<Real>::epsilon(expected));
}
} // namespace

TEST(XxxMpfr, ExplicitPrecisionAndAnalyticRootsAndEnergies)
{
  static_assert(!std::is_default_constructible_v<bethe::SolverOptions<Real>>);
  EXPECT_THROW((void)bethe::SolverOptions<Real>{Precision::exact()}, std::invalid_argument);
  for (int bits : {128, 256})
  {
    auto p = Precision::bits(bits);
    bethe::SolverOptions<Real> options(p);
    auto state = model::ground_state<Real>(4, options);
    ASSERT_TRUE(state.converged);
    near(state.energy, Real(-2, p));
    ASSERT_EQ(state.rapidities.size(), 2);
    auto root = 1 / uni20::math::sqrt(Real(3, p));
    near(state.rapidities[0], -root);
    near(state.rapidities[1], root);
    EXPECT_EQ(state.residual_norm.precision(), p);
    EXPECT_EQ(state.momentum.precision(), p);
    auto polarized = model::solve_real<Real>(4, {}, options);
    near(polarized.energy, Real(1, p));
    near(polarized.residual_norm, Real(0, p));
    near(polarized.momentum, Real(0, p));
    // Irrational energy independently detects any hidden narrowing to double.
    auto six = model::ground_state<Real>(6, options);
    ASSERT_TRUE(six.converged);
    near(six.energy, -(2 + uni20::math::sqrt(Real(13, p))) / 2, 128);
  }
}

TEST(XxxMpfr, RejectsImplicitPrecisionChanges)
{
  auto p = Precision::bits(128);
  bethe::SolverOptions<Real> options(p);
  options.residual_tolerance = Real(1);
  EXPECT_THROW((void)model::ground_state<Real>(4, options), std::invalid_argument);
  options = bethe::SolverOptions<Real>(p);
  auto state = model::ground_state<Real>(4, options);
  state.rapidities[0] = state.rapidities[0].at(Precision::bits(256));
  EXPECT_THROW((void)model::solve_real<Real>(4, state.quantum_numbers, options, state.rapidities),
               std::invalid_argument);
  options.residual_tolerance = options.residual_tolerance.at(Precision::bits(256));
  EXPECT_THROW((void)model::ground_state<Real>(4, options), std::invalid_argument);
  options = bethe::SolverOptions<Real>(p);
  auto taken = std::move(state.rapidities[0]);
  EXPECT_THROW((void)model::solve_real<Real>(4, state.quantum_numbers, options, state.rapidities),
               std::invalid_argument);
}

TEST(XxxMpfr, SectorMinimaAndSpinonBranch)
{
  auto const reference = model::one_spinon_branch<Real>(9, bethe::SolverOptions<Real>(Precision::bits(512)));
  for (int bits : {128, 256})
  {
    auto p = Precision::bits(bits);
    bethe::SolverOptions<Real> options(p);
    auto const sectors = model::sector_ground_states<Real>(5, options);
    ASSERT_EQ(sectors.size(), 6);
    for (std::size_t i = 0; i < sectors.size(); ++i)
    {
      ASSERT_TRUE(sectors[i].converged);
      EXPECT_EQ(sectors[i].energy.precision(), p);
      EXPECT_EQ(sectors[i].momentum.precision(), p);
      EXPECT_EQ(sectors[i].residual_norm.precision(), p);
      EXPECT_EQ(sectors[i].sz, -sectors[5 - i].sz);
      EXPECT_EQ(sectors[i].energy, sectors[5 - i].energy);
    }
    near(sectors[4].energy, -uni20::math::sqrt(Real(5, p)) / 4, 128);
    auto const negative = model::sector_ground_state<Real>(5, uni20::from_twice(-3), options);
    EXPECT_TRUE(negative.spin_reversed);
    EXPECT_EQ(negative.energy, sectors[1].energy);
    auto const pi = uni20::pi<Real>.at(p);
    auto const log2 = uni20::math::log(Real(2, p));
    auto const center = model::one_spinon_state<Real>(5, uni20::half_int(0), options);
    ASSERT_TRUE(center.state.converged);
    near(center.state.energy, Real(-3, p) / 4);
    near(center.spinon_momentum, pi / 2);
    near(center.bulk_subtracted_energy, 5 * log2 - 2);
    near(model::spinon_energy(center.spinon_momentum), pi / 2);
    near(model::spinon_energy(Real(0, p)), Real(0, p));
    near(model::spinon_energy(pi), Real(0, p));
    near(model::bulk_energy_density(Real(1, p)), Real(1, p) / 4 - log2);
    EXPECT_THROW((void)model::bulk_energy_density<Real>(), std::invalid_argument);
    EXPECT_THROW((void)model::spinon_energy(Real(1)), std::invalid_argument);
    EXPECT_THROW((void)model::spinon_energy(pi / 2, Real(1, Precision::bits(bits + 1))), std::invalid_argument);
    auto const branch = model::one_spinon_branch<Real>(9, options);
    ASSERT_EQ(branch.size(), reference.size());
    for (std::size_t i = 0; i < branch.size(); ++i)
    {
      ASSERT_TRUE(branch[i].state.converged);
      ASSERT_TRUE(reference[i].state.converged);
      near(branch[i].state.energy, reference[i].state.energy.at(p));
      near(branch[i].spinon_momentum, reference[i].spinon_momentum.at(p));
      near(branch[i].bulk_subtracted_energy, reference[i].bulk_subtracted_energy.at(p));
      near(model::spinon_energy(branch[i].spinon_momentum), model::spinon_energy(reference[i].spinon_momentum).at(p));
    }
  }
}

TEST(XxxMpfr, RootOnlyConcurrentPrecisionsAndThreadedExcitations)
{
  uni20::async::TbbScheduler scheduler(2);
  std::vector<model::RealExcitationScan<model::RealState<Real>>> serial(2);
  scheduler.execute_batch(2, [&](std::size_t i) {
    auto const default_before = mpfr_get_default_prec();
    // Each concurrent reference calculation needs its own serial scheduler;
    // the process-wide default DebugScheduler is single-threaded.
    uni20::async::DebugScheduler serial_scheduler;
    serial[i] = model::real_excitations<Real>(
        8, uni20::half_int(1), {.count = 100, .execution = {&serial_scheduler}},
        bethe::SolverOptions<Real>(Precision::bits(128 + 128 * i)));
    EXPECT_EQ(mpfr_get_default_prec(), default_before);
  });
  for (std::size_t i = 0; i < serial.size(); ++i)
  {
    auto p = Precision::bits(128 + 128 * i);
    auto const parallel = model::real_excitations<Real>(
        8, uni20::half_int(1), {.count = 100, .execution = {&scheduler, 3}}, bethe::SolverOptions<Real>(p));
    ASSERT_TRUE(parallel.converged());
    ASSERT_TRUE(serial[i].converged());
    ASSERT_EQ(parallel.levels.size(), 10);
    ASSERT_EQ(parallel.levels.size(), serial[i].levels.size());
    for (std::size_t j = 0; j < parallel.levels.size(); ++j)
    {
      auto const& level = parallel.levels[j];
      EXPECT_EQ(level.state.quantum_numbers, serial[i].levels[j].state.quantum_numbers);
      EXPECT_EQ(level.state.energy.precision(), p);
      EXPECT_EQ(level.state.energy, serial[i].levels[j].state.energy);
      EXPECT_EQ(level.state.rapidities, serial[i].levels[j].state.rapidities);
      ASSERT_TRUE(level.gap);
      EXPECT_EQ(level.gap->precision(), p);
      EXPECT_EQ(level.gap, serial[i].levels[j].gap);
    }
  }
}

#if UNI20_HAS_MPREAL && UNI20_HAS_MPLAPACK_MPREAL
TEST(XxxMpfr, AnalyticSpectralWeightsAndEmptyMomenta)
{
  for (int bits : {128, 256})
    for (std::size_t n : {2, 4})
    {
      auto p = Precision::bits(bits);
      auto spectrum = model::two_spinon_structure_factor<Real>(n, bethe::SolverOptions<Real>(p));
      ASSERT_TRUE(spectrum.converged());
      near(spectrum.weight_fraction, Real(1, p));
      ASSERT_TRUE(spectrum.first_moment_fraction);
      near(*spectrum.first_moment_fraction, Real(1, p));
      near(spectrum.moments[0].weight, Real(0, p));
      near(spectrum.moments[0].first_moment, Real(0, p));
      for (auto const& line : spectrum.lines)
        near(line.weight, Real(n == 2 ? 1 : line.momentum_index == 2 ? 4 : 1, p) / (n == 2 ? 1 : 3));
    }
}

TEST(XxxMpfr, AccuracyImprovesWithPrecision)
{
  auto reference = model::two_spinon_structure_factor<Real>(8, bethe::SolverOptions<Real>(Precision::bits(512)));
  ASSERT_TRUE(reference.converged());
  for (int bits : {128, 256})
  {
    auto p = Precision::bits(bits);
    auto spectrum = model::two_spinon_structure_factor<Real>(8, bethe::SolverOptions<Real>(p));
    ASSERT_TRUE(spectrum.converged());
    ASSERT_EQ(spectrum.lines.size(), reference.lines.size());
    near(spectrum.scan.ground_state.energy, reference.scan.ground_state.energy.at(p));
    for (std::size_t i = 0; i < reference.lines.size(); ++i)
      near(spectrum.lines[i].weight, reference.lines[i].weight.at(p), 32768);
  }
}

TEST(XxxMpfr, PreparedSnapshotRejectsMixedPrecisionAndSurvivesMoves)
{
  auto p = Precision::bits(128);
  bethe::SolverOptions<Real> options(p);
  auto ground = model::ground_state<Real>(8, options);
  auto excited = model::sector_ground_state<Real>(8, uni20::half_int{1}, options);
  model::PreparedRaisingFormFactors<Real> prepared(8, ground);
  ASSERT_TRUE(prepared.ready());
  auto original = prepared(excited);
  ASSERT_TRUE(original.converged());
  auto moved = std::move(prepared);
  EXPECT_FALSE(prepared.ready());
  EXPECT_FALSE(prepared(excited).weight);
  ASSERT_TRUE(moved(excited).converged());
  EXPECT_EQ(moved(excited).weight, original.weight);
  excited.rapidities[0] = excited.rapidities[0].at(Precision::bits(256));
  EXPECT_THROW((void)moved(excited), std::invalid_argument);
  ground.rapidities[0] = Real(0);
  EXPECT_THROW((model::PreparedRaisingFormFactors<Real>(8, ground)), std::invalid_argument);
  ground = model::ground_state<Real>(8, options);
  auto taken = std::move(ground.rapidities[1]);
  EXPECT_THROW((model::PreparedRaisingFormFactors<Real>(8, ground)), std::invalid_argument);
}

TEST(XxxMpfr, LogDeterminantScalingAndPrecisionDependentRejection)
{
  for (int bits : {128, 256})
  {
    auto p = Precision::bits(bits);
    Real zero(0, p), one(1, p), huge("1e1000", p), tiny("1e-1000", p);
    auto scaled = bethe::detail::log_absolute_determinant<Real>(std::vector<Real>{huge, zero, zero, tiny}, 2, zero);
    ASSERT_TRUE(scaled.log_absolute);
    // Cancellation of logs of size ~2300 loses a few bits, but never narrows.
    near(*scaled.log_absolute, zero, 65536);
    auto empty = bethe::detail::log_absolute_determinant<Real>(std::vector<Real>{}, 0, zero);
    ASSERT_TRUE(empty.log_absolute);
    near(*empty.log_absolute, zero);
    auto const epsilon = uni20::numeric_limits<Real>::epsilon(p);
    auto rejected =
        bethe::detail::log_absolute_determinant<Real>(std::vector<Real>{one, one, one, one + 32 * epsilon}, 2, zero);
    EXPECT_FALSE(rejected.log_absolute);
    auto resolved = bethe::detail::log_absolute_determinant<Real>(
        std::vector<Real>{one, one, one, one + uni20::math::sqrt(epsilon)}, 2, zero);
    EXPECT_TRUE(resolved.log_absolute);
  }
}

TEST(XxxMpfr, ConcurrentPrecisionsAndThreadedScanAreDeterministic)
{
  uni20::async::TbbScheduler scheduler(2);
  std::vector<model::TwoSpinonStructureFactor<Real>> concurrent(2);
  scheduler.execute_batch(2, [&](std::size_t i) {
    auto const default_before = mpfr_get_default_prec();
    // DebugScheduler is single-threaded; concurrent references need separate instances.
    uni20::async::DebugScheduler serial_scheduler;
    concurrent[i] = model::two_spinon_structure_factor<Real>(
        8, bethe::SolverOptions<Real>(Precision::bits(128 + 128 * i)), 100, {&serial_scheduler});
    EXPECT_EQ(mpfr_get_default_prec(), default_before);
  });
  for (std::size_t i = 0; i < concurrent.size(); ++i)
  {
    auto p = Precision::bits(128 + 128 * i);
    auto parallel = model::two_spinon_structure_factor<Real>(8, bethe::SolverOptions<Real>(p), 100, {&scheduler, 3});
    ASSERT_TRUE(concurrent[i].converged());
    ASSERT_TRUE(parallel.converged());
    ASSERT_EQ(parallel.lines.size(), concurrent[i].lines.size());
    for (std::size_t j = 0; j < parallel.lines.size(); ++j)
    {
      EXPECT_EQ(parallel.lines[j].weight.precision(), p);
      EXPECT_EQ(parallel.lines[j].weight, concurrent[i].lines[j].weight);
    }
  }
}
#endif
