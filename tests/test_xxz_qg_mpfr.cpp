// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/xxz_quantum_group_critical.hpp>
#include <bethe/xxz_quantum_group_free.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
namespace model = bethe::xxz::quantum_group::critical;
namespace free = bethe::xxz::quantum_group::free;
using R = uni20::mpreal;
class QuantumGroupMpfr : public ::testing::TestWithParam<int> {
  protected:
    uni20::Precision p() const { return uni20::Precision::bits(GetParam()); }
    R v(int value) const { return R(value, p()); }
};
INSTANTIATE_TEST_SUITE_P(WorkingPrecisions, QuantumGroupMpfr, ::testing::Values(128, 256, 333));

void precision(R const& x, uni20::Precision p)
{
  ASSERT_TRUE(x.initialized());
  EXPECT_FALSE(x.is_exact());
  EXPECT_EQ(x.precision(), p);
}
void state_precision(model::RealState<R> const& s, uni20::Precision p)
{
  precision(s.delta, p);
  precision(s.residual_norm, p);
  if (s.energy) precision(*s.energy, p);
  if (s.energy_shift) precision(*s.energy_shift, p);
  for (auto const& z : s.scaled_roots)
    precision(z, p);
  for (auto const& x : s.rapidities)
    precision(x, p);
}
void same_state(model::RealState<R> const& a, model::RealState<R> const& b)
{
  EXPECT_EQ(std::tie(a.sites, a.iterations, a.delta, a.residual_norm, a.quantum_numbers, a.scaled_roots, a.rapidities,
                     a.energy, a.energy_shift, a.converged, a.status),
            std::tie(b.sites, b.iterations, b.delta, b.residual_norm, b.quantum_numbers, b.scaled_roots, b.rapidities,
                     b.energy, b.energy_shift, b.converged, b.status));
  state_precision(a, b.delta.precision());
}

TEST_P(QuantumGroupMpfr, SmallSpinHamiltoniansAndUnreducedEquations)
{
  bethe::SolverOptions<R> options(p());
  R const eps = uni20::numeric_limits<R>::epsilon(p());
  for (R d : {v(1) / 100, v(1) / 4, v(3) / 5, v(99) / 100})
  {
    auto two = model::sea_state(2, d, 0, options);
    ASSERT_TRUE(two.converged);
    EXPECT_REAL_NEAR(*two.energy, -3 * d / 4, v(512) * eps);
    auto three = model::sea_state(3, d, 1, options);
    ASSERT_TRUE(three.converged);
    EXPECT_REAL_NEAR(*three.energy, -(1 + d) / 2, v(512) * eps);
    auto four = model::sea_state(4, d, 0, options);
    ASSERT_TRUE(four.converged);
    EXPECT_REAL_NEAR(*four.energy, -3 * d / 4 - uni20::math::sqrt(d * d + 2) / 2, v(2048) * eps);
    state_precision(four, p());
  }
  using C = uni20::complex<R>;
  for (R d : {v(1) / 4, v(3) / 5})
  {
    auto state = model::solve_real<R>(8, d, std::vector<uni20::half_int>{1, 3}, options);
    ASSERT_TRUE(state.converged);
    state_precision(state, p());
    R const gamma = uni20::math::acos(d);
    auto ratio = [&](R x, R w) {
      C numerator(uni20::math::sinh(x) * uni20::math::cos(w), uni20::math::cosh(x) * uni20::math::sin(w));
      return numerator / uni20::math::conj(numerator);
    };
    R energy = 7 * d / 4;
    for (std::size_t i = 0; i < state.rapidities.size(); ++i)
    {
      R const x = state.rapidities[i];
      C lhs(v(1), v(0)), rhs(v(1), v(0));
      auto drive = ratio(x, gamma / 2);
      for (int k = 0; k < 16; ++k)
        lhs *= drive;
      for (std::size_t j = 0; j < state.rapidities.size(); ++j)
        if (i != j) rhs *= ratio(x - state.rapidities[j], gamma) * ratio(x + state.rapidities[j], gamma);
      EXPECT_REAL_NEAR(lhs.real(), rhs.real(), v(65536) * eps);
      EXPECT_REAL_NEAR(lhs.imag(), rhs.imag(), v(65536) * eps);
      energy -= (1 - d * d) / (uni20::math::cosh(2 * x) - d);
    }
    EXPECT_REAL_NEAR(*state.energy, energy, v(4096) * eps);
  }
}

TEST_P(QuantumGroupMpfr, FreeBlocksCompleteSectorsAndBudgets)
{
  auto b = free::block<R>(4, std::vector<std::size_t>{3}, 1, p());
  EXPECT_EQ(b.block_size, 2u);
  EXPECT_REAL_NEAR(b.energy, -uni20::math::sqrt(v(2)) / 2, v(32) * uni20::numeric_limits<R>::epsilon(p()));
  for (std::size_t n = 2; n <= 8; ++n)
    for (std::size_t down = 0; down <= n; ++down)
    {
      auto blocks = free::sector<R>(n, down, p());
      std::size_t dimension = 0, pairs = 0;
      for (auto const& block : blocks)
      {
        precision(block.energy, p());
        EXPECT_EQ(block.down, down);
        dimension += block.block_size;
        pairs += block.block_size == 2;
      }
      EXPECT_EQ(dimension, bethe::detail::bounded_binomial(n, down, 10000));
      EXPECT_EQ(pairs, n % 2 || down == 0 ? 0 : bethe::detail::bounded_binomial(n - 2, down - 1, 10000));
    }
  auto zero = free::block<R>(4, std::vector<std::size_t>{1, 3}, 0, p());
  EXPECT_EQ(zero.energy, 0);
  precision(zero.energy, p());
  EXPECT_EQ(free::block<R>(4, {}, 2, p()).energy, 0);
  EXPECT_THROW((void)free::sector<R>(4, 2, p(), {.max_blocks = 3}), std::length_error);
  EXPECT_THROW((void)free::sector<R>(4, 2, p(), {.max_mode_entries = 3}), std::length_error);
  EXPECT_THROW((void)free::block<R>(4, {}, 0, uni20::Precision::exact()), std::invalid_argument);
  EXPECT_THROW((void)free::sector<R>(4, 2, uni20::Precision::exact()), std::invalid_argument);
  EXPECT_THROW((void)free::sector<double>(4, 2, p()), std::invalid_argument);
}

TEST_P(QuantumGroupMpfr, OwningFreeSectorQueries)
{
  using Sector = std::vector<free::Block<R>>;
  auto submit = [p = p()](uni20::async::IAsyncScheduler& scheduler, auto& output) {
    std::vector<std::pair<std::size_t, std::size_t>> requests{{8, 0}, {8, 4}, {9, 4}, {9, 9}};
    return bethe::evaluate_observable_grid_async(
        requests.size(),
        [p, requests](std::size_t i) {
          auto const [n, down] = requests[i];
          return free::sector<R>(n, down, p);
        },
        [&output](std::size_t i, Sector result) {
          EXPECT_EQ(i, output.size());
          output.push_back(std::move(result));
        },
        {&scheduler, 4, 1});
  };
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  std::vector<Sector> expected;
  auto final = submit(deferred, expected);
  EXPECT_EQ(final.get_wait(deferred), 4u);
  for (int workers : {1, 4})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<Sector> actual;
    auto complete = submit(scheduler, actual);
    EXPECT_EQ(complete.get_wait(scheduler), expected.size());
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      ASSERT_EQ(actual[i].size(), expected[i].size());
      for (std::size_t j = 0; j < actual[i].size(); ++j)
      {
        auto const& a = actual[i][j];
        auto const& b = expected[i][j];
        EXPECT_EQ(std::tie(a.energy, a.modes, a.zero_occupation, a.down, a.block_size),
                  std::tie(b.energy, b.modes, b.zero_occupation, b.down, b.block_size));
        precision(a.energy, p());
      }
    }
  }
}

TEST_P(QuantumGroupMpfr, WindowsFailuresAndCompleteConcurrentScans)
{
  bethe::SolverOptions<R> options(p());
  for (R bad : {R(1), R(uni20::uninitialized), R(1, uni20::Precision::bits(GetParam() + 1)), R("nan", p()),
                R("inf", p()), v(0), v(1), -v(1)})
    EXPECT_THROW((void)model::sea_state(4, bad, 0, options), std::invalid_argument);
  EXPECT_EQ(model::real_quantum_number_window(8, v(1) / 2, 4).slots, 4u);
  EXPECT_EQ(model::real_quantum_number_window(8, v(51) / 100, 4).slots, 5u);
  auto window = model::real_quantum_number_window(8, v(3) / 5, 4);
  precision(window.exclusive_threshold, p());
  auto vacuum = model::sea_state(8, v(1) / 4, 8, options);
  ASSERT_TRUE(vacuum.converged);
  EXPECT_EQ(*vacuum.energy_shift, 0);
  state_precision(vacuum, p());
  auto reference = model::real_excitations(16, v(3) / 5, 12, {.count = 100}, options);
  ASSERT_TRUE(reference.converged());
  ASSERT_EQ(reference.levels.size(), 45u);
  auto min_updates = reference.levels[0].state.iterations;
  auto max_updates = min_updates;
  for (auto const& level : reference.levels)
  {
    min_updates = std::min(min_updates, level.state.iterations);
    max_updates = std::max(max_updates, level.state.iterations);
  }
  ASSERT_LT(min_updates, max_updates);
  for (auto budget : {std::size_t{0}, min_updates, std::size_t{10000}})
  {
    options.max_iterations = budget;
    auto expected = model::real_excitations(16, v(3) / 5, 12, {.count = 100}, options);
    if (budget == min_updates)
    {
      EXPECT_GT(expected.converged_count, 0u);
      EXPECT_LT(expected.converged_count, expected.candidate_count);
    }
    for (int workers : {1, 4})
    {
      uni20::async::TbbScheduler scheduler(workers);
      auto actual = model::real_excitations(16, v(3) / 5, 12, {.count = 100, .execution = {&scheduler, 8, 3}}, options);
      EXPECT_EQ(actual.candidate_count, expected.candidate_count);
      EXPECT_EQ(actual.converged_count, expected.converged_count);
      same_state(actual.ground_state, expected.ground_state);
      ASSERT_EQ(bool(actual.first_unconverged), bool(expected.first_unconverged));
      if (actual.first_unconverged) same_state(*actual.first_unconverged, *expected.first_unconverged);
      ASSERT_EQ(actual.levels.size(), expected.levels.size());
      for (std::size_t i = 0; i < actual.levels.size(); ++i)
      {
        same_state(actual.levels[i].state, expected.levels[i].state);
        EXPECT_EQ(actual.levels[i].gap, expected.levels[i].gap);
        if (actual.levels[i].gap) precision(*actual.levels[i].gap, p());
      }
    }
  }
}

TEST(QuantumGroupMpfrShared, PrecisionGrowthAndEndpointResolution)
{
  auto const high = uni20::Precision::bits(512);
  R const d = R(1, high) / 4;
  R const exact = -3 * d / 4 - uni20::math::sqrt(d * d + 2) / 2;
  R previous(1, high);
  for (int bits : {128, 256, 333})
  {
    auto const p = uni20::Precision::bits(bits);
    auto s = model::sea_state(4, R(1, p) / 4, 0, bethe::SolverOptions<R>(p));
    ASSERT_TRUE(s.converged);
    R const error = uni20::math::abs(s.energy->at(high) - exact);
    EXPECT_LT(error, previous / R("1e15", high));
    EXPECT_LT(error, R(2048, high) * uni20::numeric_limits<R>::epsilon(p).at(high));
    previous = error;
  }
  auto const p = uni20::Precision::bits(333);
  R const tiny("1e-60", p);
  EXPECT_EQ(model::real_quantum_number_window(4, tiny, 0).slots, 2u);
  EXPECT_EQ(model::real_quantum_number_window(4, tiny.at(uni20::Precision::bits(128)), 0).slots, 1u);
  for (R delta : {tiny, R(1, p) - tiny})
  {
    auto s = model::sea_state(4, delta, 0, bethe::SolverOptions<R>(p));
    ASSERT_TRUE(s.converged);
    state_precision(s, p);
    R const reference = -3 * delta / 4 - uni20::math::sqrt(delta * delta + 2) / 2;
    EXPECT_REAL_NEAR(*s.energy, reference, R("1e-95", p));
  }
}

TEST(QuantumGroupMpfrShared, OwningMixedPrecisionQueries)
{
  struct Request
  {
      std::size_t sites;
      R delta;
      bethe::SolverOptions<R> options;
  };
  auto submit = [](uni20::async::IAsyncScheduler& scheduler, auto& output) {
    std::vector<Request> requests;
    for (int bits : {128, 256, 333})
    {
      auto p = uni20::Precision::bits(bits);
      bethe::SolverOptions<R> options(p);
      requests.push_back({4, R(1, p) / 4, options});
      options.max_iterations = 0;
      requests.push_back({6, R(3, p) / 5, options});
    }
    return bethe::evaluate_observable_grid_async(
        requests.size(),
        [requests](std::size_t i) {
          auto const& r = requests[i];
          return model::sea_state(r.sites, r.delta, 0, r.options);
        },
        [&output](std::size_t i, auto result) {
          EXPECT_EQ(i, output.size());
          output.push_back(std::move(result));
        },
        {&scheduler, 8, 2});
  };
  auto const ambient = mpfr_get_default_prec();
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  std::vector<model::RealState<R>> expected;
  auto final = submit(deferred, expected);
  EXPECT_EQ(final.get_wait(deferred), 6u);
  EXPECT_FALSE(expected.back().converged);
  for (int workers : {1, 4})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<model::RealState<R>> actual;
    auto complete = submit(scheduler, actual);
    EXPECT_EQ(complete.get_wait(scheduler), expected.size());
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
      same_state(actual[i], expected[i]);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
} // namespace
