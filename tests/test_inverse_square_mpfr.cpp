// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include <bethe/haldane_shastry.hpp>
#include <bethe/observable_grid.hpp>
#include <bethe/sutherland.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
namespace hs = bethe::haldane_shastry;
namespace su = bethe::sutherland;
class InverseSquareMpfr : public ::testing::TestWithParam<int> {
  protected:
    uni20::Precision p() const { return uni20::Precision::bits(GetParam()); }
    R v(int n) const { return R(n, p()); }
};
INSTANTIATE_TEST_SUITE_P(WorkingPrecisions, InverseSquareMpfr, ::testing::Values(128, 256, 333));
void finite(R const& value, uni20::Precision p)
{
  EXPECT_FALSE(value.is_exact());
  EXPECT_EQ(value.precision(), p);
  EXPECT_TRUE(uni20::isfinite(value));
}
void finite(hs::Level<R> const& s, uni20::Precision p)
{
  for (auto const& x : {s.energy, s.gap, s.momentum})
    finite(x, p);
}
void finite(su::State<R> const& s, uni20::Precision p)
{
  for (auto const& x : {s.length, s.lambda, s.energy, s.ground_energy, s.gap, s.momentum})
    finite(x, p);
  for (auto const& x : s.pseudomomenta)
    finite(x, p);
}

TEST_P(InverseSquareMpfr, SpinHamiltonianAndExactMotifEnumeration)
{
  R const pi = uni20::pi<R>.at(p()), eps = uni20::numeric_limits<R>::epsilon(p());
  // Direct spin Hamiltonian on the polarized state, using chord couplings.
  for (std::size_t n : {2u, 3u, 7u, 12u})
  {
    R energy = v(0);
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = i + 1; j < n; ++j)
      {
        R distance = uni20::math::sin(pi * R(j - i) / R(n));
        energy += (pi / R(n)) * (pi / R(n)) / (4 * distance * distance);
      }
    auto state = hs::evaluate<R>(n, {}, p());
    finite(state, p());
    EXPECT_REAL_NEAR(state.energy, energy, 64 * eps * energy);
    EXPECT_EQ(state.momentum, 0);
    auto ground = hs::ground_levels<R>(n, p());
    ASSERT_EQ(ground.size(), n % 2 ? 2u : 1u);
    for (auto const& g : ground)
    {
      finite(g, p());
      EXPECT_EQ(g.gap, 0);
      EXPECT_EQ(g.spinons, n % 2);
    }
  }
  auto all = hs::spectrum<R>(8, p());
  ASSERT_TRUE(all.complete);
  ASSERT_EQ(all.total_motifs, 34u);
  std::uint64_t dimension = 0;
  for (std::size_t i = 0; i < all.levels.size(); ++i)
  {
    auto const& state = all.levels[i];
    finite(state, p());
    ASSERT_TRUE(state.degeneracy);
    dimension += *state.degeneracy;
    if (i) EXPECT_LE(all.levels[i - 1].energy_numerator, state.energy_numerator);
  }
  EXPECT_EQ(dimension, 256u);
  auto few = hs::spectrum<R>(8, p(), 3);
  ASSERT_EQ(few.levels.size(), 3u);
  EXPECT_EQ(few.levels.back().motif, all.levels[2].motif);
  auto refused = hs::spectrum<R>(8, p(), 1, 33);
  EXPECT_FALSE(refused.complete);
  EXPECT_TRUE(refused.levels.empty());
  auto sector = hs::sector_ground_levels<R>(7, uni20::half_int::parse("3/2"), p());
  ASSERT_EQ(sector.size(), 2u);
  EXPECT_EQ(sector[0].energy, sector[1].energy);
  EXPECT_EQ(sector[0].spinons, 3u);
}

R gegenbauer(int n, R lambda, R x)
{
  if (n < 0) return uni20::scalar_like(x, 0);
  if (n == 0) return uni20::scalar_like(x, 1);
  R previous = uni20::scalar_like(x, 1), current = 2 * lambda * x;
  for (int j = 2; j <= n; ++j)
  {
    R next = (2 * (j - 1 + lambda) * x * current - (j - 2 + 2 * lambda) * previous) / j;
    previous = current;
    current = next;
  }
  return current;
}
TEST_P(InverseSquareMpfr, DifferentialHamiltonianOnTwoBodyWavefunctions)
{
  R const pi = uni20::pi<R>.at(p()), eps = uni20::numeric_limits<R>::epsilon(p());
  R const length = v(5), factor = 2 * (pi / length) * (pi / length);
  for (R lambda : {v(1) / 4, v(1), v(13) / 4})
    for (int d = 0; d <= 8; ++d)
      for (R y : {pi / 5, 3 * pi / 7})
      {
        R sine = uni20::math::sin(y), cosine = uni20::math::cos(y), cot = cosine / sine;
        R g = gegenbauer(d, lambda, cosine);
        R first = 2 * lambda * gegenbauer(d - 1, lambda + 1, cosine);
        R second = 4 * lambda * (lambda + 1) * gegenbauer(d - 2, lambda + 2, cosine);
        R f_second = (lambda * lambda * cot * cot - lambda / (sine * sine)) * g - (2 * lambda + 1) * cosine * first +
                     sine * sine * second;
        int m = d - 6;
        // Apply H to exp(i*pi*M*(x1+x2)/L)*sin(y)^lambda*C_d^lambda(cos(y)).
        R actual = factor * (m * m * g - f_second + lambda * (lambda - 1) * g / (sine * sine));
        auto state = su::evaluate(2, length, lambda, {-3, d - 3});
        finite(state, p());
        EXPECT_REAL_NEAR(actual, state.energy * g, 1024 * eps * (1 + state.energy) * (1 + uni20::math::abs(g)));
      }
}

TEST_P(InverseSquareMpfr, SmallGapsTypedZerosAndLargeExponentRange)
{
  R const pi = uni20::pi<R>.at(p()), eps = uni20::numeric_limits<R>::epsilon(p()), unit = 2 * pi / 5;
  auto boosted = su::evaluate(4, v(5), R("1e1000", p()), {1, 1, 1, 1});
  finite(boosted, p());
  EXPECT_EQ(boosted.energy, boosted.ground_energy);
  EXPECT_REAL_NEAR(boosted.gap, 4 * unit * unit, 32 * eps * unit * unit);
  for (R lambda : {v(0), v(1), v(7) / 3})
    for (std::size_t n : {1u, 3u})
    {
      auto ground = su::ground_state(n, v(5), lambda);
      finite(ground, p());
      EXPECT_EQ(ground.gap, 0);
      EXPECT_EQ(ground.momentum, 0);
      EXPECT_EQ(ground.pseudomomenta[n / 2], 0);
      if (n == 1 || lambda == 0) EXPECT_EQ(ground.energy, 0);
    }
  auto tiny = su::evaluate(1, R("1e7000", p()), v(0), {1});
  auto huge = su::evaluate(1, R("1e-7000", p()), v(0), {1});
  finite(tiny, p());
  finite(huge, p());
  EXPECT_GT(tiny.energy, 0);
  EXPECT_REAL_NEAR(tiny.energy * huge.energy, 16 * pi * pi * pi * pi, 4096 * eps);
  auto lo = std::numeric_limits<std::int64_t>::min(), hi = std::numeric_limits<std::int64_t>::max();
  auto cancelling = su::evaluate(4, v(4), v(2), {lo, lo, hi, hi});
  finite(cancelling, p());
  EXPECT_EQ(cancelling.momentum_index, -2);
  R sum = v(0);
  for (auto k : cancelling.pseudomomenta)
    sum += k * k;
  EXPECT_REAL_NEAR(cancelling.energy, sum, 32 * eps * sum);
  EXPECT_THROW(su::evaluate(2, v(4), v(2), {hi, hi}), std::overflow_error);
}

TEST_P(InverseSquareMpfr, SutherlandWindowAndPrecisionValidation)
{
  auto all = su::spectrum(3, v(4), v(2), 2);
  ASSERT_TRUE(all.complete);
  ASSERT_EQ(all.states.size(), 35u);
  for (std::size_t i = 0; i < all.states.size(); ++i)
  {
    auto const& s = all.states[i];
    finite(s, p());
    if (i) EXPECT_LE(all.states[i - 1].gap, s.gap);
    auto labels = s.labels;
    std::reverse(labels.begin(), labels.end());
    for (auto& j : labels)
      j = -j;
    auto reflected = su::evaluate(3, v(4), v(2), labels);
    EXPECT_EQ(reflected.gap, s.gap);
    EXPECT_EQ(reflected.momentum_index, -s.momentum_index);
  }
  auto refused = su::spectrum(3, v(4), v(2), 2, {}, 34);
  EXPECT_FALSE(refused.complete);
  EXPECT_TRUE(refused.states.empty());
  for (R bad : {R(uni20::uninitialized), R(4), R(4, uni20::Precision::bits(GetParam() + 1))})
  {
    EXPECT_THROW(su::ground_state(3, v(4), bad), std::invalid_argument);
    EXPECT_THROW(su::ground_state(3, bad, v(2)), std::invalid_argument);
    EXPECT_THROW(su::spectrum(3, v(4), bad, 2, {}, 0), std::invalid_argument);
  }
  for (R bad : {v(-1), R("nan", p()), R("inf", p())})
    EXPECT_THROW(su::ground_state(3, v(4), bad), std::invalid_argument);
  EXPECT_THROW(su::ground_state(3, v(0), v(2)), std::invalid_argument);
  auto exact = R(0).precision();
  EXPECT_THROW(hs::evaluate<R>(4, {1, 3}, exact), std::invalid_argument);
  EXPECT_THROW(hs::ground_levels<R>(4, exact), std::invalid_argument);
  EXPECT_THROW(hs::spectrum<R>(4, exact, {}, 0), std::invalid_argument);
}

TEST(InverseSquareMpfrShared, ExactOrderingDespiteRoundedEnergyTies)
{
  auto coarse = hs::spectrum<R>(10, uni20::Precision::bits(4));
  auto fine = hs::spectrum<R>(10, uni20::Precision::bits(256));
  ASSERT_TRUE(coarse.complete);
  ASSERT_EQ(coarse.levels.size(), fine.levels.size());
  bool collapsed = false;
  for (std::size_t i = 0; i < coarse.levels.size(); ++i)
  {
    EXPECT_EQ(coarse.levels[i].motif, fine.levels[i].motif);
    EXPECT_EQ(coarse.levels[i].energy_numerator, fine.levels[i].energy_numerator);
    if (i && coarse.levels[i].energy_numerator != coarse.levels[i - 1].energy_numerator &&
        coarse.levels[i].energy == coarse.levels[i - 1].energy)
      collapsed = true;
  }
  EXPECT_TRUE(collapsed);
}

TEST(InverseSquareMpfrShared, AccuracyGrowsBeyondBinary128)
{
  auto reference_precision = uni20::Precision::bits(512);
  R pi = uni20::pi<R>.at(reference_precision);
  // The dimer singlet is -3/4 times its single exchange coupling.
  R hs_reference = -3 * pi * pi / 16;
  // N=2, lambda=9/4, labels (-3,2): relative quantum number d=5, COM=-1.
  // The two-body relative differential operator has eigenvalue (d+lambda)^2.
  R su_reference = 2 * (pi / R("7.25", reference_precision)) * (pi / R("7.25", reference_precision)) *
                   (1 + R("7.25", reference_precision) * R("7.25", reference_precision));
  R previous_hs(1, reference_precision), previous_su(1, reference_precision);
  for (int bits : {128, 256, 333})
  {
    auto p = uni20::Precision::bits(bits);
    R e_hs = uni20::math::abs(hs::ground_levels<R>(2, p)[0].energy.at(reference_precision) - hs_reference);
    R e_su = uni20::math::abs(su::evaluate(2, R("7.25", p), R("2.25", p), {-3, 2}).energy.at(reference_precision) -
                              su_reference);
    EXPECT_LT(e_hs, previous_hs / R("1e15", reference_precision));
    EXPECT_LT(e_su, previous_su / R("1e15", reference_precision));
    EXPECT_LT(e_hs, (64 * uni20::numeric_limits<R>::epsilon(p)).at(reference_precision));
    EXPECT_LT(e_su, (1024 * uni20::numeric_limits<R>::epsilon(p)).at(reference_precision));
    previous_hs = e_hs;
    previous_su = e_su;
  }
}

TEST(InverseSquareMpfrShared, OwningMixedPrecisionSpectraAndFailures)
{
  struct Result
  {
      hs::Spectrum<R> spin;
      su::Spectrum<R> gas;
  };
  auto const ambient = mpfr_get_default_prec();
  auto submit = [](uni20::async::IAsyncScheduler& scheduler, std::vector<Result>& rows, bool fail) {
    std::vector<R> lengths;
    for (int bits : {128, 256, 333, 192})
      lengths.emplace_back(4, uni20::Precision::bits(bits));
    if (fail) lengths.back() = R(4);
    return bethe::evaluate_observable_grid_async(
        lengths.size(),
        [lengths](std::size_t i) {
          auto p = lengths[i].precision();
          return Result{hs::spectrum<R>(8, p, {}, i == 2 ? 0 : 1000),
                        su::spectrum(3, lengths[i], uni20::scalar_like(lengths[i], 2), 2, {}, i == 2 ? 0 : 1000)};
        },
        [&rows](std::size_t i, Result value) {
          EXPECT_EQ(i, rows.size());
          rows.push_back(std::move(value));
        },
        {&scheduler, 128, 2});
  };
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  std::vector<Result> reference;
  auto baseline = submit(deferred, reference, false);
  EXPECT_EQ(baseline.get_wait(deferred), 4u);
  ASSERT_EQ(reference.size(), 4u);
  EXPECT_FALSE(reference[2].spin.complete);
  EXPECT_FALSE(reference[2].gas.complete);
  for (int workers : {1, 4})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<Result> actual, failed;
    auto bad = submit(scheduler, failed, true);
    auto good = submit(scheduler, actual, false);
    EXPECT_THROW((void)bad.get_wait(scheduler), std::invalid_argument);
    EXPECT_EQ(good.get_wait(scheduler), 4u);
    ASSERT_EQ(actual.size(), reference.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      auto const& a = actual[i];
      auto const& b = reference[i];
      EXPECT_EQ(a.spin.complete, b.spin.complete);
      EXPECT_EQ(a.spin.total_motifs, b.spin.total_motifs);
      ASSERT_EQ(a.spin.levels.size(), b.spin.levels.size());
      for (std::size_t j = 0; j < a.spin.levels.size(); ++j)
      {
        auto const& s = a.spin.levels[j];
        auto const& r = b.spin.levels[j];
        finite(s, r.energy.precision());
        EXPECT_EQ(s.sites, r.sites);
        EXPECT_EQ(s.motif, r.motif);
        EXPECT_EQ(s.energy_numerator, r.energy_numerator);
        EXPECT_EQ(s.momentum_index, r.momentum_index);
        EXPECT_EQ(s.spinons, r.spinons);
        EXPECT_EQ(s.maximum_spin, r.maximum_spin);
        EXPECT_EQ(s.degeneracy, r.degeneracy);
        EXPECT_EQ(s.energy, r.energy);
        EXPECT_EQ(s.gap, r.gap);
        EXPECT_EQ(s.momentum, r.momentum);
      }
      EXPECT_EQ(a.gas.complete, b.gas.complete);
      EXPECT_EQ(a.gas.total_states, b.gas.total_states);
      ASSERT_EQ(a.gas.states.size(), b.gas.states.size());
      for (std::size_t j = 0; j < a.gas.states.size(); ++j)
      {
        auto const& s = a.gas.states[j];
        auto const& r = b.gas.states[j];
        finite(s, r.energy.precision());
        EXPECT_EQ(s.particles, r.particles);
        EXPECT_EQ(s.length, r.length);
        EXPECT_EQ(s.lambda, r.lambda);
        EXPECT_EQ(s.labels, r.labels);
        EXPECT_EQ(s.pseudomomenta, r.pseudomomenta);
        EXPECT_EQ(s.momentum_index, r.momentum_index);
        EXPECT_EQ(s.momentum, r.momentum);
        EXPECT_EQ(s.energy, r.energy);
        EXPECT_EQ(s.gap, r.gap);
        EXPECT_EQ(s.ground_energy, r.ground_energy);
      }
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
} // namespace
