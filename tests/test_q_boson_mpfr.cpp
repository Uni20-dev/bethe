// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/q_boson.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
namespace model = bethe::q_boson;
using R = uni20::mpreal;
class QBosonMpfr : public ::testing::TestWithParam<int> {
  protected:
    uni20::Precision p() const { return uni20::Precision::bits(GetParam()); }
    R v(int n) const { return R(n, p()); }
};
INSTANTIATE_TEST_SUITE_P(WorkingPrecisions, QBosonMpfr, ::testing::Values(128, 256, 333));

void precision(R const& x, uni20::Precision p)
{
  ASSERT_TRUE(x.initialized());
  EXPECT_FALSE(x.is_exact());
  EXPECT_EQ(x.precision(), p);
}
void state_precision(model::State<R> const& s, uni20::Precision p)
{
  for (auto const& x : {s.eta, s.residual_norm, s.momentum})
    precision(x, p);
  if (s.energy) precision(*s.energy, p);
  for (auto const& x : s.momenta)
    precision(x, p);
  for (auto const& x : s.deviations)
    precision(x, p);
}
void same_state(model::State<R> const& a, model::State<R> const& b)
{
  EXPECT_EQ(std::tie(a.sites, a.particles, a.eta, a.residual_norm, a.quantum_numbers, a.momenta, a.modes, a.deviations,
                     a.momentum_index, a.momentum, a.energy, a.iterations, a.converged, a.status),
            std::tie(b.sites, b.particles, b.eta, b.residual_norm, b.quantum_numbers, b.momenta, b.modes, b.deviations,
                     b.momentum_index, b.momentum, b.energy, b.iterations, b.converged, b.status));
  state_precision(a, b.eta.precision());
}

TEST_P(QBosonMpfr, FreePhaseAndInitialStateWithScalars)
{
  model::SolverOptions<R> options(p());
  options.max_iterations = 0;
  R const pi = uni20::pi<R>.at(p()), eps = uni20::numeric_limits<R>::epsilon(p());
  for (std::size_t n : {0, 1, 4})
  {
    auto s = model::ground_state(5, n, v(0), options);
    ASSERT_TRUE(s.converged && s.energy);
    EXPECT_EQ(*s.energy, 0);
    EXPECT_EQ(s.momentum, 0);
    state_precision(s, p());
  }
  auto s = model::ground_state(7, 9, R("inf", p()), options);
  ASSERT_TRUE(s.converged && s.energy);
  state_precision(s, p());
  EXPECT_EQ(s.iterations, 0u);
  R expected = 18 - 2 * uni20::math::sin(9 * pi / 16) / uni20::math::sin(pi / 16);
  EXPECT_REAL_NEAR(*s.energy, expected, v(1024) * eps);
  auto one = model::solve_modes<R>(5, std::vector<std::size_t>{4}, v(1), options);
  ASSERT_TRUE(one.converged && one.energy);
  EXPECT_REAL_NEAR(*one.energy, (5 - uni20::math::sqrt(v(5))) / 2, v(128) * eps);
  state_precision(one, p());
  auto failed = model::ground_state(5, 3, v(1), options);
  EXPECT_FALSE(failed.converged);
  EXPECT_FALSE(failed.energy);
  EXPECT_EQ(failed.status, model::Status::iteration_limit);
  EXPECT_EQ(failed.iterations, 0u);
  state_precision(failed, p());
}

TEST_P(QBosonMpfr, LargePhaseModeSumRetainsExactLabels)
{
  std::size_t const sites = std::numeric_limits<std::int64_t>::max() / 4;
  // Nine near-zone-boundary labels overflow int64_t when summed. Only the
  // labels are large: the solver still has nine roots and a 9x9 Jacobian.
  std::vector<std::size_t> modes(9, sites - 1);
  modes[0] = 0;
  model::SolverOptions<R> options(p());
  options.max_iterations = 0;
  auto const state = model::solve_modes<R>(sites, modes, R("inf", p()), options);
  ASSERT_TRUE(state.converged && state.energy);
  state_precision(state, p());
  EXPECT_EQ(state.momentum_index, sites - 8);
  R const pi = uni20::pi<R>.at(p());
  R const lifted = 2 * pi * (R(8) * R(sites - 1)) / R(sites);
  R const tolerance = v(1024) * uni20::numeric_limits<R>::epsilon(p());
  for (std::size_t j = 0; j < modes.size(); ++j)
  {
    R const label = 2 * R(modes[j]) + 2 * R(j) - 8;
    R const expected = (pi * label + lifted) / (R(sites) + 9);
    EXPECT_REAL_NEAR(state.momenta[j], expected, tolerance);
  }
}

TEST_P(QBosonMpfr, RejectsMixedExactUnsetAndInvalidInputs)
{
  model::SolverOptions<R> options(p());
  for (R bad : {R(uni20::uninitialized), R(), R(1), R(1, uni20::Precision::bits(GetParam() + 1)), -v(1), R("nan", p()),
                R("-inf", p())})
    EXPECT_THROW((void)model::ground_state(4, 0, bad, options), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state(1, 0, v(1), options), std::invalid_argument);
  options.residual_tolerance = R(1);
  EXPECT_THROW((void)model::ground_state(4, 0, v(1), options), std::invalid_argument);
}

TEST_P(QBosonMpfr, IndependentTwoSiteSpectrumAndMultiplicativeEquations)
{
  model::SolverOptions<R> options(p());
  R const eps = uni20::numeric_limits<R>::epsilon(p());
  for (R eta : {v(1) / 100, v(1), v(10)})
  {
    auto scan = model::real_excitations(2, 2, eta, {.count = 3}, options);
    ASSERT_TRUE(scan.converged());
    ASSERT_EQ(scan.levels.size(), 3u);
    // Eigenvalues of the three-dimensional occupation-space Hamiltonian.
    R const loss = -uni20::math::expm1(-2 * eta);
    R const low = 2 * loss / (1 + uni20::math::sqrt(1 - loss / 2));
    std::vector<R> expected{low, v(4), 8 - low};
    for (std::size_t j = 0; j < 3; ++j)
    {
      auto const& level = scan.levels[j];
      state_precision(level.state, p());
      ASSERT_TRUE(level.gap);
      precision(*level.gap, p());
      EXPECT_REAL_NEAR(*level.state.energy, expected[j], v(1024) * eps);
      EXPECT_REAL_NEAR(*level.gap, expected[j] - low, v(2048) * eps);
    }
    auto state = model::solve_modes<R>(5, std::vector<std::size_t>{0, 0, 2, 4, 4}, eta, options);
    ASSERT_TRUE(state.converged);
    state_precision(state, p());
    using C = uni20::complex<R>;
    for (std::size_t j = 0; j < state.particles; ++j)
    {
      C product = uni20::math::exp(C(v(0), 5 * state.momenta[j]));
      for (std::size_t k = 0; k < state.particles; ++k)
        if (k != j)
        {
          R const half = (state.momenta[j] - state.momenta[k]) / 2;
          product *= uni20::math::sin(C(half, -eta)) / uni20::math::sin(C(half, eta));
        }
      EXPECT_REAL_NEAR(product.real(), v(1), v(32768) * eps);
      EXPECT_REAL_NEAR(product.imag(), v(0), v(32768) * eps);
    }
  }
}

TEST(QBosonMpfrShared, PrecisionGrowthAndWeakDeformation)
{
  auto const high = uni20::Precision::bits(512);
  R const loss = -uni20::math::expm1(R(-2, high));
  R const reference = 2 * loss / (1 + uni20::math::sqrt(1 - loss / 2));
  R previous(1, high);
  for (int bits : {128, 256, 333})
  {
    auto const p = uni20::Precision::bits(bits);
    auto state = model::ground_state(2, 2, R(1, p), model::SolverOptions<R>(p));
    ASSERT_TRUE(state.converged && state.energy);
    R error = uni20::math::abs(state.energy->at(high) - reference);
    EXPECT_LT(error, previous / R("1e15", high));
    EXPECT_LT(error, R(1024, high) * uni20::numeric_limits<R>::epsilon(p).at(high));
    previous = error;
  }
  auto const p = uni20::Precision::bits(333);
  for (auto text : {"1e-60", "1e-7000"})
  {
    R eta(text, p);
    auto state = model::ground_state(2, 2, eta, model::SolverOptions<R>(p));
    ASSERT_TRUE(state.converged && state.energy);
    EXPECT_GT(*state.energy, 0);
    EXPECT_REAL_NEAR(*state.energy / (2 * eta), R(1, p), R("1e-55", p));
    if (eta == R("1e-60", p)) EXPECT_REAL_NEAR((*state.energy / (2 * eta) - 1) / eta, -R(3, p) / 4, R("1e-35", p));
    state_precision(state, p);
  }
  auto clustered =
      model::solve_modes<R>(5, std::vector<std::size_t>{1, 1}, R("1e-7000", p), model::SolverOptions<R>(p));
  ASSERT_TRUE(clustered.converged);
  EXPECT_EQ(clustered.momenta[0], clustered.momenta[1]);
  EXPECT_LT(clustered.deviations[0], clustered.deviations[1]);
}

TEST_P(QBosonMpfr, CompleteAndPartialParallelScans)
{
  model::SolverOptions<R> options(p());
  for (std::size_t budget : {0, 100})
  {
    options.max_iterations = budget;
    auto reference = model::real_excitations(4, 3, v(1), {.count = 100}, options);
    EXPECT_EQ(reference.converged(), budget != 0);
    for (int workers : {1, 4})
    {
      uni20::async::TbbScheduler scheduler(workers);
      auto actual = model::real_excitations(4, 3, v(1), {.count = 100, .execution = {&scheduler, 8, 3}}, options);
      EXPECT_EQ(actual.candidate_count, reference.candidate_count);
      EXPECT_EQ(actual.converged_count, reference.converged_count);
      ASSERT_EQ(actual.levels.size(), reference.levels.size());
      same_state(actual.ground_state, reference.ground_state);
      ASSERT_EQ(bool(actual.first_unconverged), bool(reference.first_unconverged));
      if (actual.first_unconverged) same_state(*actual.first_unconverged, *reference.first_unconverged);
      for (std::size_t i = 0; i < actual.levels.size(); ++i)
      {
        same_state(actual.levels[i].state, reference.levels[i].state);
        EXPECT_EQ(actual.levels[i].gap, reference.levels[i].gap);
      }
    }
  }
}

TEST(QBosonMpfrShared, OwningMixedPrecisionQueries)
{
  struct Request
  {
      std::size_t sites, particles;
      R eta;
      model::SolverOptions<R> options;
  };
  auto submit = [](uni20::async::IAsyncScheduler& scheduler, auto& output) {
    std::vector<Request> requests;
    for (int bits : {128, 256, 333})
    {
      auto const p = uni20::Precision::bits(bits);
      model::SolverOptions<R> options(p);
      requests.push_back({4, 3, R(0, p), options});
      requests.push_back({6, 4, R(1, p), options});
      options.max_iterations = 0;
      requests.push_back({5, 3, R(1, p), options});
    }
    return bethe::evaluate_observable_grid_async(
        requests.size(),
        [requests](std::size_t i) {
          auto const& r = requests[i];
          return model::ground_state(r.sites, r.particles, r.eta, r.options);
        },
        [&output](std::size_t i, auto state) {
          EXPECT_EQ(i, output.size());
          output.push_back(std::move(state));
        },
        {&scheduler, 8, 2});
  };
  auto const ambient = mpfr_get_default_prec();
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  std::vector<model::State<R>> expected;
  auto serial = submit(deferred, expected);
  EXPECT_EQ(serial.get_wait(deferred), 9u);
  ASSERT_FALSE(expected.back().converged);
  for (int workers : {1, 4})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<model::State<R>> actual;
    auto complete = submit(scheduler, actual);
    EXPECT_EQ(complete.get_wait(scheduler), expected.size());
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
      same_state(actual[i], expected[i]);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
} // namespace
