// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include <bethe/hubbard.hpp>
#include <bethe/hubbard_open.hpp>
#include <bethe/observable_grid.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using uni20::Precision;
namespace model = bethe::hubbard;
namespace math = uni20::math;
void finite(R const& value, Precision p)
{
  ASSERT_TRUE(value.initialized());
  EXPECT_FALSE(value.is_exact());
  EXPECT_EQ(value.precision(), p);
  EXPECT_TRUE(uni20::isfinite(value));
}
template <class State> void check(State const& s, Precision p)
{
  for (auto const* value : {&s.interaction, &s.root_interaction, &s.energy_offset, &s.energy, &s.charge_residual,
                            &s.spin_residual, &s.residual_norm})
    finite(*value, p);
  for (auto const& value : s.charge_momenta)
    finite(value, p);
  for (auto const& value : s.spin_rapidities)
    finite(value, p);
  if constexpr (requires { s.momentum; }) finite(s.momentum, p);
  EXPECT_EQ(s.residual_norm, std::max(s.charge_residual, s.spin_residual));
  EXPECT_EQ(s.converged, s.status == model::SolveStatus::converged);
}
void near(R const& actual, R const& expected, int units = 4096)
{
  EXPECT_EQ(actual.precision(), expected.precision());
  EXPECT_REAL_NEAR(actual, expected, units * uni20::numeric_limits<R>::epsilon(expected) * (1 + math::abs(expected)));
}
R dimer(R const& u, bool periodic)
{
  R const squared_hopping = uni20::scalar_like(u, periodic ? 64 : 16);
  R const root = math::sqrt(u * u + squared_hopping);
  return u >= 0 ? -squared_hopping / (2 * (u + root)) : (u - root) / 2;
}
// Reconstruct the unreduced logarithmic equations from conventional returned
// k and Lambda, rather than using the solver's reduced system or Jacobian.
template <class State> void full_equations(State const& s)
{
  auto p = s.interaction.precision();
  constexpr bool periodic = requires { s.momentum; };
  R zero(0, p), pi = uni20::pi<R>.at(p), energy = s.energy_offset;
  for (auto const& k : s.charge_momenta)
    energy -= 2 * math::cos(k);
  near(energy, s.energy);
  if (s.free_fermion) return;
  R const u = s.root_interaction / 4, norm = uni20::scalar_like(u, periodic ? s.sites : 2 * (s.sites + 1));
  R charge = zero, spin = zero;
  std::vector<R> sine;
  for (auto const& k : s.charge_momenta)
    sine.push_back(periodic && (k == 0 || k == pi) ? zero : math::sin(k));
  for (std::size_t j = 0; j < s.charge_momenta.size(); ++j)
  {
    R f = norm * s.charge_momenta[j] - pi * R(s.quantum_numbers.charge[j].twice());
    for (auto const& l : s.spin_rapidities)
    {
      f += 2 * math::atan((sine[j] - l) / u);
      if constexpr (!periodic) f += 2 * math::atan((sine[j] + l) / u);
    }
    charge = std::max(charge, math::abs(f) / norm);
  }
  for (std::size_t a = 0; a < s.spin_rapidities.size(); ++a)
  {
    R const l = s.spin_rapidities[a];
    R f = -pi * R(s.quantum_numbers.spin[a].twice());
    for (auto const& v : sine)
    {
      f += 2 * math::atan((l - v) / u);
      if constexpr (!periodic) f += 2 * math::atan((l + v) / u);
    }
    for (std::size_t b = 0; b < s.spin_rapidities.size(); ++b)
    {
      if (a == b) continue;
      f -= 2 * math::atan((l - s.spin_rapidities[b]) / (2 * u));
      if constexpr (!periodic) f -= 2 * math::atan((l + s.spin_rapidities[b]) / (2 * u));
    }
    spin = std::max(spin, math::abs(f) / norm);
  }
  near(charge, s.charge_residual);
  near(spin, s.spin_residual);
}
class HubbardMpfr : public ::testing::TestWithParam<int> {};
INSTANTIATE_TEST_SUITE_P(WorkingPrecisions, HubbardMpfr, ::testing::Values(128, 256, 333));
} // namespace

TEST_P(HubbardMpfr, ScalarFreePathsAndUnreducedEquations)
{
  auto p = Precision::bits(GetParam());
  for (R u : {R(0, p), R(4, p), R(-4, p)})
    for (std::size_t n : {0, 1, 4})
    {
      auto ring = model::detail::free_state<R>(4, n, 0, u);
      auto wall = model::open::detail::free_state<R>(4, n, 0, u);
      check(ring, p);
      check(wall, p);
      full_equations(ring);
      full_equations(wall);
      EXPECT_EQ(ring.residual_norm, 0);
      EXPECT_EQ(wall.residual_norm, 0);
      if (n == 0 || n == 4)
      {
        EXPECT_EQ(ring.energy, 0);
        EXPECT_EQ(wall.energy, 0);
      }
    }
  // This configuration needs only MPFR scalars, not an MPFR LU provider.
  model::detail::GroundSystem<R> periodic(6, 6, 3, p);
  model::open::detail::GroundSystem<R> open(3, 3, 1, p);
  for (R u : {R(1, p) / 4, R(2, p)})
  {
    auto a = periodic.evaluate(periodic.seed(), u);
    auto b = open.evaluate(open.seed(), u);
    finite(a.charge, p);
    finite(a.spin, p);
    finite(b.charge, p);
    finite(b.spin, p);
    for (auto const& v : a.residual)
      finite(v, p);
    for (auto const& v : b.residual)
      finite(v, p);
  }
}

TEST(HubbardMpfrControls, RejectsMixedExactAndUnsetInputs)
{
  auto p = Precision::bits(128), q = Precision::bits(256);
  model::SolverOptions<R> controls(p);
  for (R u : {R(4), R(4, q), R(uni20::uninitialized)})
    EXPECT_THROW((void)model::detail::validate_sector(4, 4, uni20::half_int(0), u, controls), std::invalid_argument);
  controls.residual_tolerance = R(1);
  EXPECT_THROW((void)model::detail::validate_sector(4, 4, uni20::half_int(0), R(4, p), controls),
               std::invalid_argument);
#if UNI20_HAS_MPLAPACK_MPREAL
  EXPECT_THROW((void)model::ground_state(4, R(4, p), controls), std::invalid_argument);
  EXPECT_THROW((void)model::open::ground_state(4, R(4, p), controls), std::invalid_argument);
#endif
}

#if UNI20_HAS_MPLAPACK_MPREAL
TEST_P(HubbardMpfr, DimersMappedSectorsAndFailureDiagnostics)
{
  auto p = Precision::bits(GetParam());
  model::SolverOptions<R> options(p);
  for (R u : {R(-16, p), R(-1, p) / 4, R(0, p), R(1, p) / 4, R(4, p), R(32, p)})
  {
    auto ring = model::ground_state(2, u, options);
    auto wall = model::open::ground_state(2, u, options);
    ASSERT_TRUE(ring.converged);
    ASSERT_TRUE(wall.converged);
    check(ring, p);
    check(wall, p);
    near(ring.energy, dimer(u, true));
    near(wall.energy, dimer(u, false));
    full_equations(ring);
    full_equations(wall);
    // Mapping and exact-zero fields, including empty and completely filled bands.
    for (std::size_t n : {0, 2, 4, 6, 8})
    {
      auto a = model::sector_ground_state(4, n, uni20::half_int(0), u, options);
      auto b = model::open::sector_ground_state(4, n, uni20::half_int(0), u, options);
      ASSERT_TRUE(a.converged);
      ASSERT_TRUE(b.converged);
      check(a, p);
      check(b, p);
      full_equations(a);
      full_equations(b);
      if (n == 0)
      {
        EXPECT_EQ(a.energy, 0);
        EXPECT_EQ(b.energy, 0);
      }
      if (n == 8)
      {
        EXPECT_EQ(a.energy, 4 * u);
        EXPECT_EQ(b.energy, 4 * u);
      }
    }
  }
  for (R u : {R(-4, p), R(0, p), R(4, p)})
    for (std::size_t n : {0, 1, 2})
    {
      auto site =
          model::open::sector_ground_state(1, n, uni20::from_twice(static_cast<std::int64_t>(n % 2)), u, options);
      ASSERT_TRUE(site.converged);
      check(site, p);
      EXPECT_EQ(site.energy, n == 2 ? u : R(0, p));
    }
  for (int sign : {-1, 1})
  {
    R u = sign * R(1, p) / 4;
    auto a = model::sector_ground_state(6, 6, uni20::half_int(-2), u, options);
    auto b = model::open::sector_ground_state(5, 5, uni20::from_twice(-1), u, options);
    ASSERT_TRUE(a.converged);
    ASSERT_TRUE(b.converged);
    check(a, p);
    check(b, p);
    full_equations(a);
    full_equations(b);
    options.max_iterations = 0;
    auto fa = model::ground_state(4, u, options);
    auto fb = model::open::ground_state(3, u, options);
    EXPECT_FALSE(fa.converged);
    EXPECT_FALSE(fb.converged);
    check(fa, p);
    check(fb, p);
    full_equations(fa);
    full_equations(fb);
    EXPECT_EQ(fa.root_interaction, math::abs(u));
    EXPECT_EQ(fb.root_interaction, math::abs(u));
    options.max_iterations = 10000;
  }
  EXPECT_THROW((void)model::sector_ground_state(6, 3, uni20::from_twice(1), R(4, p), options), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state(3, R(4, p), options), std::invalid_argument);
}

TEST(HubbardMpfrAccuracy, GrowthAgainstIndependentDimerAndTinyInteraction)
{
  auto high = Precision::bits(512);
  for (int sign : {-1, 1})
  {
    R previous_ring(1, high), previous_wall(1, high), u = sign * R(1, high) / 4;
    for (int bits : {128, 256, 333})
    {
      auto p = Precision::bits(bits);
      model::SolverOptions<R> options(p);
      auto ring = model::ground_state(2, u.at(p), options);
      auto wall = model::open::ground_state(2, u.at(p), options);
      ASSERT_TRUE(ring.converged);
      ASSERT_TRUE(wall.converged);
      R er = math::abs(ring.energy.at(high) - dimer(u, true));
      R ew = math::abs(wall.energy.at(high) - dimer(u, false));
      EXPECT_LT(er, previous_ring * R("1e-18", high));
      EXPECT_LT(ew, previous_wall * R("1e-18", high));
      previous_ring = er;
      previous_wall = ew;
    }
  }
  auto p = Precision::bits(333);
  model::SolverOptions<R> options(p);
  R tiny("1e-60", p);
  for (int sign : {-1, 1})
  {
    auto ring = model::ground_state(2, sign * tiny, options);
    auto wall = model::open::ground_state(2, sign * tiny, options);
    ASSERT_TRUE(ring.converged);
    ASSERT_TRUE(wall.converged);
    near(ring.energy, dimer(sign * tiny, true));
    near(wall.energy, dimer(sign * tiny, false));
    EXPECT_NE(ring.energy, -4);
    EXPECT_NE(wall.energy, -2);
  }
}

TEST(HubbardMpfrAsync, OwningMixedPrecisionGroundSweeps)
{
  struct Request
  {
      R interaction;
      model::SolverOptions<R> options;
  };
  for (int workers : {1, 3})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<Request> inputs;
    for (int bits : {128, 192, 256})
      inputs.push_back({R(-1, Precision::bits(bits)) / 4, model::SolverOptions<R>(Precision::bits(bits))});
    inputs.back().options.max_iterations = 0;
    auto reference = inputs;
    using Result = std::pair<model::State<R>, model::open::State<R>>;
    std::vector<Result> results;
    auto const count = inputs.size();
    auto done = bethe::evaluate_observable_grid_async(
        count,
        [inputs = std::move(inputs)](std::size_t i) {
          auto const& item = inputs[i];
          return Result{model::ground_state(4, item.interaction, item.options),
                        model::open::ground_state(3, item.interaction, item.options)};
        },
        [&](std::size_t i, Result state) {
          EXPECT_EQ(i, results.size());
          results.push_back(std::move(state));
        },
        {.scheduler = &scheduler, .max_in_flight = 2});
    EXPECT_EQ(done.get_wait(scheduler), 3u);
    ASSERT_EQ(results.size(), 3u);
    auto compare = [](auto const& a, auto const& b) {
      EXPECT_EQ(a.energy, b.energy);
      EXPECT_EQ(a.charge_momenta, b.charge_momenta);
      EXPECT_EQ(a.spin_rapidities, b.spin_rapidities);
      EXPECT_EQ(a.charge_residual, b.charge_residual);
      EXPECT_EQ(a.spin_residual, b.spin_residual);
      EXPECT_EQ(a.status, b.status);
      EXPECT_EQ(a.iterations, b.iterations);
      EXPECT_EQ(a.continuation_steps, b.continuation_steps);
    };
    for (std::size_t i = 0; i < results.size(); ++i)
    {
      auto const& item = reference[i];
      check(results[i].first, item.options.precision());
      check(results[i].second, item.options.precision());
      compare(results[i].first, model::ground_state(4, item.interaction, item.options));
      compare(results[i].second, model::open::ground_state(3, item.interaction, item.options));
    }
    EXPECT_FALSE(results.back().first.converged);
    EXPECT_FALSE(results.back().second.converged);
  }
}
#endif
