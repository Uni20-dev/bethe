// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/xxz_excitations.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using uni20::Precision;
namespace model = bethe::xxz;
namespace math = uni20::math;
void finite(R const& x, Precision p)
{
  ASSERT_TRUE(x.initialized());
  EXPECT_FALSE(x.is_exact());
  EXPECT_EQ(x.precision(), p);
  EXPECT_TRUE(uni20::isfinite(x));
}
template <class State> void check(State const& s, Precision p)
{
  finite(s.delta, p);
  finite(s.energy, p);
  finite(s.residual_norm, p);
  for (auto const& x : s.rapidities)
    finite(x, p);
  if constexpr (requires { s.momentum; }) finite(s.momentum, p);
  if constexpr (requires { s.log_rapidities; })
    for (auto const& x : s.log_rapidities)
      finite(x, p);
  if constexpr (requires { s.root_delta; }) finite(s.root_delta, p);
  if constexpr (requires { s.boundary_root; })
    if (s.boundary_root)
    {
      finite(s.boundary_root->log_distance, p);
      finite(s.boundary_root->inverse_square, p);
    }
}
void near(R const& got, R const& expected, int units = 2048)
{
  EXPECT_EQ(got.precision(), expected.precision());
  EXPECT_REAL_NEAR(got, expected, units * uni20::numeric_limits<R>::epsilon(expected) * (1 + math::abs(expected)));
}
// Independent N=4 open-chain Hamiltonian, even spin-reversal block:
// [d/4, 1/2, 0; 1/2, -3d/4, 1; 0, 1, -d/4]. Its smallest
// eigenvalue lies below every diagonal entry. Bisect its characteristic cubic;
// no Bethe roots, continuation, Newton matrix or library eigensolver is used.
R open_four_reference(R d)
{
  auto polynomial = [&](R e) { return (e - d / 4) * ((e + 3 * d / 4) * (e + d / 4) - 1) - (e + d / 4) / 4; };
  R low = -2 * (1 + math::abs(d)), high = std::min({d / 4, -3 * d / 4, -d / 4});
  for (int i = 0; i < uni20::numeric_limits<R>::digits(d.precision()) + 4; ++i)
  {
    R mid = (low + high) / 2;
    if (mid == low || mid == high) break;
    if (polynomial(mid) > 0)
      high = mid;
    else
      low = mid;
  }
  return (low + high) / 2;
}
class XxzMpfr : public ::testing::TestWithParam<int> {};
INSTANTIATE_TEST_SUITE_P(WorkingPrecisions, XxzMpfr, ::testing::Values(128, 256, 333));
} // namespace

TEST_P(XxzMpfr, ExplicitRealStatesAndEndpointZeros)
{
  auto p = Precision::bits(GetParam());
  model::SolverOptions<R> options(p);
  for (R d : {R(0, p), R(1, p) / 2, R(1, p)})
  {
    auto ring = model::solve_real<R>(4, d, model::sector_ground_quantum_numbers(4, uni20::half_int(0)), options);
    auto wall =
        model::open::solve_real<R>(3, d, model::open::sector_ground_quantum_numbers(3, uni20::from_twice(1)), options);
    ASSERT_TRUE(ring.converged);
    ASSERT_TRUE(wall.converged);
    check(ring, p);
    check(wall, p);
    near(ring.energy, -d / 2 - math::sqrt(d * d / 4 + 2));
    near(wall.energy, -d / 4 - math::sqrt(d * d / 16 + R(1, p) / 2));
    auto vacuum = model::solve_real<R>(4, d, {}, options);
    auto open_vacuum = model::open::solve_real<R>(4, d, {}, options);
    check(vacuum, p);
    check(open_vacuum, p);
    EXPECT_EQ(vacuum.residual_norm, 0);
    EXPECT_EQ(vacuum.momentum, 0);
    EXPECT_EQ(vacuum.energy, d);
    EXPECT_EQ(open_vacuum.energy, 3 * d / 4);
    auto failed_options = options;
    failed_options.max_iterations = 0;
    auto failed = model::solve_real<R>(4, d, ring.quantum_numbers, failed_options);
    EXPECT_FALSE(failed.converged);
    check(failed, p);
    auto window = model::real_quantum_number_window(12, d, uni20::half_int(1));
    auto native_window = model::real_quantum_number_window(12, static_cast<double>(d), uni20::half_int(1));
    EXPECT_EQ(window.first, native_window.first);
    EXPECT_EQ(window.last, native_window.last);
  }
}

TEST(XxzMpfrControls, RejectsMixedExactAndUnsetInputs)
{
  auto p = Precision::bits(128), q = Precision::bits(256);
  model::SolverOptions<R> options(p);
  auto numbers = model::sector_ground_quantum_numbers(4, uni20::half_int(0));
  for (R d : {R(1), R(1, q), R(uni20::uninitialized)})
  {
    EXPECT_THROW((void)model::solve_real<R>(4, d, numbers, options), std::invalid_argument);
#if UNI20_HAS_MPLAPACK_MPREAL
    EXPECT_THROW((void)model::ground_state(4, d, options), std::invalid_argument);
    EXPECT_THROW((void)model::open::ground_state(4, d, options), std::invalid_argument);
    EXPECT_THROW((void)model::open::sector_ground_states(4, d, options), std::invalid_argument);
#endif
  }
  options.residual_tolerance = R(1);
  EXPECT_THROW((void)model::solve_real<R>(4, R(1, p), numbers, options), std::invalid_argument);
  options = model::SolverOptions<R>(p);
  auto roots = model::solve_real<R>(4, R(1, p) / 2, numbers, options).rapidities;
  roots.front() = roots.front().at(q);
  EXPECT_THROW((void)model::solve_real<R>(4, R(1, p) / 2, numbers, options, roots), std::invalid_argument);
  EXPECT_THROW((void)model::open::solve_real<R>(
                   4, R(1, p) / 2, model::QuantumNumbers{uni20::half_int(1), uni20::half_int(2)}, options, roots),
               std::invalid_argument);
}

#if UNI20_HAS_MPLAPACK_MPREAL
TEST(XxzMpfrControls, MassiveContinuationRequiresRepresentableProgress)
{
  auto p = Precision::bits(4);
  R d(2, p);
  model::SolverOptions<R> options(p);
  // The loose default tolerance accepts the seed without a Newton update;
  // 1.125 + (1.125-1)/4 rounds back to 1.125 at four bits.
  for (std::size_t budget : {1u, 10000u})
  {
    options.max_iterations = budget;
    auto state = model::open::ground_state(3, d, options);
    EXPECT_FALSE(state.converged);
    EXPECT_EQ(state.status, model::open::GroundSolveStatus::continuation_stalled);
    EXPECT_EQ(state.iterations, 0u);
    EXPECT_EQ(state.delta, d);
    EXPECT_EQ(state.root_delta, R(9, p) / 8);
    EXPECT_LT(state.root_delta, state.delta);
    check(state, p);
    ASSERT_EQ(state.rapidities.size(), 1u);
    // Returned energy still refers to the requested Hamiltonian, not the
    // unfinished continuation stage, even when no update was possible.
    R theta = math::atan(state.rapidities.front());
    EXPECT_EQ(state.energy, -d / 2 - math::cos(2 * theta));
  }
  options.max_iterations = 0;
  auto budget = model::open::ground_state(3, d, options);
  EXPECT_EQ(budget.status, model::open::GroundSolveStatus::iteration_limit);
  // Reaching the requested stage needs no continuation step, even with zero
  // remaining updates; do not reject low precision as a blanket policy.
  auto reached = model::open::ground_state(3, R(9, p) / 8, options);
  EXPECT_TRUE(reached.converged);
  auto q = Precision::bits(3);
  auto boundary = model::open::ground_state(4, R(2, q), model::SolverOptions<R>(q));
  EXPECT_FALSE(boundary.converged);
  EXPECT_EQ(boundary.status, model::open::GroundSolveStatus::continuation_stalled);
  check(boundary, q);
  q = Precision::bits(2);
  // The near-isotropic initial point itself rounds to one. Reject before
  // constructing the even-chain boundary seed, which would divide by zero.
  EXPECT_THROW((void)model::open::ground_state(4, R(2, q), model::SolverOptions<R>(q)), std::invalid_argument);
  auto vacuum = model::open::sector_ground_state(4, R(2, q), uni20::half_int(2), model::SolverOptions<R>(q));
  EXPECT_TRUE(vacuum.converged);
  check(vacuum, q);
}

TEST_P(XxzMpfr, AllGroundRepresentationsAgainstSmallChainHamiltonians)
{
  auto p = Precision::bits(GetParam());
  model::SolverOptions<R> options(p);
  for (R d : {R(-1, p) / 2, R(0, p), R(1, p) / 2, R(1, p), R(2, p), R(8, p)})
  {
    SCOPED_TRACE(static_cast<double>(d));
    auto ring = model::ground_state(4, d, options);
    ASSERT_TRUE(ring.converged);
    check(ring, p);
    near(ring.energy, -d / 2 - math::sqrt(d * d / 4 + 2));
    for (std::size_t n : {2, 3, 4})
    {
      auto s = model::open::ground_state(n, d, options);
      ASSERT_TRUE(s.converged);
      check(s, p);
      R expected = n == 2   ? -R(1, p) / 2 - d / 4
                   : n == 3 ? -d / 4 - math::sqrt(d * d / 16 + R(1, p) / 2)
                            : open_four_reference(d.at(Precision::bits(512))).at(p);
      near(s.energy, expected);
      if (n % 2 == 0 && d > 1) ASSERT_TRUE(s.boundary_root);
    }
    auto sectors = model::open::sector_ground_states(4, d, options);
    auto rings = model::sector_ground_states(4, d, options);
    for (std::size_t i = 0; i <= 4; ++i)
    {
      ASSERT_TRUE(sectors[i].converged);
      ASSERT_TRUE(rings[i].converged);
      check(sectors[i], p);
      check(rings[i], p);
      EXPECT_EQ(sectors[i].energy, sectors[4 - i].energy);
      EXPECT_EQ(rings[i].energy, rings[4 - i].energy);
    }
    options.max_iterations = 0;
    auto failed = model::open::ground_state(4, d, options);
    EXPECT_FALSE(failed.converged);
    check(failed, p);
    options.max_iterations = 10000;
  }
  EXPECT_THROW((void)model::ground_state(3, R(-1, p) / 2, options), std::invalid_argument);
  EXPECT_THROW((void)model::solve_real<R>(4, R(-1, p) / 2, {}, options), std::invalid_argument);
}

TEST(XxzMpfrAccuracy, GrowthAndNearbyEndpointsAreNotSnapped)
{
  auto high = Precision::bits(512);
  for (int twice_delta : {-1, 1, 4})
  {
    R reference = open_four_reference(R(twice_delta, high) / 2), previous(1, high);
    for (int bits : {128, 256, 333})
    {
      auto p = Precision::bits(bits);
      auto s = model::open::ground_state(4, R(twice_delta, p) / 2, model::SolverOptions<R>(p));
      ASSERT_TRUE(s.converged);
      R error = math::abs(s.energy.at(high) - reference);
      EXPECT_LT(error, previous * R("1e-18", high));
      EXPECT_LT(error, 2048 * uni20::numeric_limits<R>::epsilon(p).at(high));
      previous = error;
    }
  }
  auto p = Precision::bits(333);
  model::SolverOptions<R> options(p);
  R tiny("1e-60", p);
  auto wall = model::open::ground_state(2, tiny, options);
  ASSERT_TRUE(wall.converged);
  near(wall.energy, -R(1, p) / 2 - tiny / 4);
  EXPECT_NE(wall.energy, -R(1, p) / 2);
  for (int sign : {-1, 1})
  {
    R d = R(1, p) + sign * tiny;
    auto s = model::open::ground_state(4, d, options);
    ASSERT_TRUE(s.converged);
    near(s.energy, open_four_reference(d.at(high)).at(p));
    EXPECT_NE(s.energy, model::open::ground_state(4, R(1, p), options).energy);
  }
}

TEST(XxzMpfrAsync, OwningMixedPrecisionSweepAndExistingContinuousScans)
{
  for (int workers : {1, 3})
  {
    uni20::async::TbbScheduler scheduler(workers);
    struct Request
    {
        R delta;
        model::SolverOptions<R> options;
    };
    std::vector<Request> requests;
    for (int bits : {128, 192, 256})
      requests.push_back({R(-1, Precision::bits(bits)) / 2, model::SolverOptions<R>(Precision::bits(bits))});
    requests.back().options.max_iterations = 0;
    auto reference = requests;
    std::vector<model::open::GroundState<R>> states;
    auto const count = requests.size();
    auto done = bethe::evaluate_observable_grid_async(
        count,
        [requests = std::move(requests)](std::size_t i) {
          return model::open::ground_state(4, requests[i].delta, requests[i].options);
        },
        [&](std::size_t i, auto state) {
          EXPECT_EQ(i, states.size());
          states.push_back(std::move(state));
        },
        {.scheduler = &scheduler, .max_in_flight = 2});
    EXPECT_EQ(done.get_wait(scheduler), 3u);
    ASSERT_EQ(states.size(), 3u);
    for (std::size_t i = 0; i < states.size(); ++i)
    {
      auto expected = model::open::ground_state(4, reference[i].delta, reference[i].options);
      check(states[i], reference[i].options.precision());
      EXPECT_EQ(states[i].energy, expected.energy);
      EXPECT_EQ(states[i].rapidities, expected.rapidities);
      EXPECT_EQ(states[i].residual_norm, expected.residual_norm);
      EXPECT_EQ(states[i].status, expected.status);
    }
    EXPECT_FALSE(states.back().converged);
    auto p = Precision::bits(160);
    model::SolverOptions<R> options(p);
    auto compare = [&](auto const& serial, auto const& parallel) {
      EXPECT_TRUE(serial.converged());
      EXPECT_TRUE(parallel.converged());
      EXPECT_EQ(serial.candidate_count, parallel.candidate_count);
      EXPECT_EQ(serial.converged_count, parallel.converged_count);
      EXPECT_EQ(serial.ground_state.energy, parallel.ground_state.energy);
      ASSERT_EQ(serial.levels.size(), parallel.levels.size());
      for (std::size_t i = 0; i < serial.levels.size(); ++i)
      {
        auto const& a = serial.levels[i];
        auto const& b = parallel.levels[i];
        check(b.state, p);
        EXPECT_EQ(a.state.energy, b.state.energy);
        EXPECT_EQ(a.state.rapidities, b.state.rapidities);
        EXPECT_EQ(a.state.quantum_numbers, b.state.quantum_numbers);
        EXPECT_EQ(a.state.residual_norm, b.state.residual_norm);
        EXPECT_EQ(a.state.iterations, b.state.iterations);
        EXPECT_EQ(a.gap, b.gap);
        ASSERT_TRUE(b.gap);
        finite(*b.gap, p);
      }
    };
    model::RealExcitationOptions serial{.count = 100};
    auto parallel = serial;
    parallel.execution = {.scheduler = &scheduler, .max_in_flight = 2};
    compare(model::real_excitations(6, R(1, p) / 2, uni20::half_int(1), serial, options),
            model::real_excitations(6, R(1, p) / 2, uni20::half_int(1), parallel, options));
    compare(model::open::real_excitations(6, R(1, p) / 2, uni20::half_int(1), serial, options),
            model::open::real_excitations(6, R(1, p) / 2, uni20::half_int(1), parallel, options));
  }
}
#endif
