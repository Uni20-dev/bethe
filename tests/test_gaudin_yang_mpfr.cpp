// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include <bethe/gaudin_yang.hpp>
#include <bethe/observable_grid.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
namespace model = bethe::gaudin_yang;
using R = uni20::mpreal;
using C = uni20::mpcomplex;
class GaudinYangMpfr : public ::testing::TestWithParam<int> {
  protected:
    uni20::Precision p() const { return uni20::Precision::bits(GetParam()); }
    R v(int n) const { return R(n, p()); }
};
INSTANTIATE_TEST_SUITE_P(WorkingPrecisions, GaudinYangMpfr, ::testing::Values(128, 256, 333));

void finite(R const& x, uni20::Precision p)
{
  EXPECT_FALSE(x.is_exact());
  EXPECT_EQ(x.precision(), p);
  EXPECT_TRUE(uni20::isfinite(x));
}
void finite_state(model::State<R> const& s, uni20::Precision p)
{
  for (auto const& x : {s.length, s.interaction, s.energy, s.momentum, s.charge_residual, s.spin_residual,
                        s.residual_norm, s.root_interaction})
    finite(x, p);
  for (auto const& x : s.momenta)
    finite(x, p);
  for (auto const& x : s.spin_rapidities)
    finite(x, p);
}
void same_state(model::State<R> const& a, model::State<R> const& b)
{
  EXPECT_EQ(std::tie(a.up, a.down, a.length, a.interaction, a.quantum_numbers.charge, a.quantum_numbers.spin, a.momenta,
                     a.spin_rapidities, a.free_modes, a.energy, a.momentum, a.momentum_index, a.charge_residual,
                     a.spin_residual, a.residual_norm, a.root_interaction, a.iterations, a.free, a.spin_reversed,
                     a.converged, a.status),
            std::tie(b.up, b.down, b.length, b.interaction, b.quantum_numbers.charge, b.quantum_numbers.spin, b.momenta,
                     b.spin_rapidities, b.free_modes, b.energy, b.momentum, b.momentum_index, b.charge_residual,
                     b.spin_residual, b.residual_norm, b.root_interaction, b.iterations, b.free, b.spin_reversed,
                     b.converged, b.status));
  finite_state(a, b.length.precision());
}

// Unreduced equations in physical coordinates, independently of the solver's
// reflected coordinates, scaled spin roots and weak-coupling rank subtraction.
void original_equations(model::State<R> const& s)
{
  auto p = s.length.precision();
  R const zero(0, p), one(1, p), pi = uni20::pi<R>.at(p), g = s.interaction * s.length;
  R const tolerance = R(65536, p) * R(s.up + s.down) * uni20::numeric_limits<R>::epsilon(p);
  R charge = zero, spin = zero;
  auto ratio = [](R d, R w) { return C(d, w) / C(d, -w); };
  for (std::size_t j = 0; j < s.momenta.size(); ++j)
  {
    R const q = s.momenta[j] * s.length;
    R f = q - pi * R(s.quantum_numbers.charge[j].twice());
    C product(one, zero);
    for (R lambda : s.spin_rapidities)
    {
      R const d = q - lambda * s.length;
      f += 2 * uni20::math::atan(2 * d / g);
      product *= ratio(d, g / 2);
    }
    R const scale = g < one ? std::max(uni20::math::sqrt(g), uni20::math::abs(q)) : R(s.up + s.down, p);
    charge = std::max(charge, uni20::math::abs(f) / scale);
    if (s.converged)
    {
      EXPECT_REAL_NEAR(product.real(), uni20::math::cos(q), tolerance);
      EXPECT_REAL_NEAR(product.imag(), uni20::math::sin(q), tolerance);
    }
  }
  for (std::size_t a = 0; a < s.spin_rapidities.size(); ++a)
  {
    R const lambda = s.spin_rapidities[a] * s.length;
    R f = -pi * R(s.quantum_numbers.spin[a].twice());
    C lhs(one, zero), rhs(one, zero);
    for (R k : s.momenta)
    {
      R const d = lambda - k * s.length;
      f += 2 * uni20::math::atan(2 * d / g);
      lhs *= ratio(d, g / 2);
    }
    for (std::size_t b = 0; b < s.spin_rapidities.size(); ++b)
      if (a != b)
      {
        R const d = lambda - s.spin_rapidities[b] * s.length;
        f -= 2 * uni20::math::atan(d / g);
        rhs *= ratio(d, g);
      }
    spin = std::max(spin, uni20::math::abs(f) / R(s.up + s.down));
    if (s.converged)
    {
      EXPECT_REAL_NEAR(lhs.real(), rhs.real(), tolerance);
      EXPECT_REAL_NEAR(lhs.imag(), rhs.imag(), tolerance);
    }
  }
  EXPECT_REAL_NEAR(charge, s.charge_residual, tolerance);
  EXPECT_REAL_NEAR(spin, s.spin_residual, tolerance);
  EXPECT_EQ(s.residual_norm, std::max(s.charge_residual, s.spin_residual));
}

TEST_P(GaudinYangMpfr, FreeVacuumPolarizationAndFiniteZeros)
{
  model::SolverOptions<R> options(p());
  R const pi = uni20::pi<R>.at(p()), tolerance = v(256) * uni20::numeric_limits<R>::epsilon(p());
  for (std::size_t up = 0; up < 5; ++up)
    for (std::size_t down = 0; down < 5; ++down)
    {
      auto const s = model::ground_state(up, down, v(2), v(0), options);
      ASSERT_TRUE(s.converged && s.free);
      finite_state(s, p());
      R expected = v(0);
      for (auto n : {up, down})
        expected += pi * pi * R(n) * (R(n) * R(n) + (n % 2 ? -v(1) : v(2))) / 12;
      EXPECT_REAL_NEAR(s.energy, expected, tolerance * (v(1) + expected));
      EXPECT_EQ(s.momentum_index, std::int64_t((up % 2 ? 0 : up / 2) + (down % 2 ? 0 : down / 2)));
    }
  auto const a = model::ground_state(4, 0, v(2), v(9), options), b = model::ground_state(0, 4, v(2), v(9), options);
  finite_state(a, p());
  finite_state(b, p());
  EXPECT_EQ(a.energy, b.energy);
  EXPECT_TRUE(b.spin_reversed);
}

TEST_P(GaudinYangMpfr, ExplicitPrecisionSystemAndInputValidation)
{
  model::SolverOptions<R> options(p());
  model::detail::GroundSystem<R> system(6, 3, p());
  auto x = system.seed();
  auto evaluation = system.evaluate(x, v(4));
  for (auto const& value : evaluation.residual)
    finite(value, p());
  finite(evaluation.charge, p());
  finite(evaluation.spin, p());
  EXPECT_EQ(system.expand(x, v(4))[1][1], v(0));
  for (R bad : {R(1), R(1, uni20::Precision::bits(GetParam() + 1)), R()})
  {
    EXPECT_THROW((void)model::ground_state(1, 1, bad, v(1), options), std::invalid_argument);
    EXPECT_THROW((void)model::ground_state(1, 1, v(1), bad, options), std::invalid_argument);
  }
  EXPECT_THROW((void)model::ground_state(2, 2, v(1), v(1), options), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state(1, 1, v(1), -v(1), options), std::invalid_argument);
  options.residual_tolerance = R(1);
  EXPECT_THROW((void)model::ground_state(0, 0, v(1), v(0), options), std::invalid_argument);
}

TEST_P(GaudinYangMpfr, ZeroBudgetPreservesRequestedCouplingDiagnosticsWithoutLU)
{
  model::SolverOptions<R> options(p());
  options.max_iterations = 0;
  auto const state = model::ground_state(3, 3, v(6), v(1) / 100, options);
  EXPECT_FALSE(state.converged);
  EXPECT_EQ(state.iterations, 0u);
  EXPECT_EQ(state.status, model::SolveStatus::iteration_limit);
  EXPECT_GT(state.root_interaction, state.interaction);
  finite_state(state, p());
  original_equations(state);
}

TEST_P(GaudinYangMpfr, NestedEquationsSpinReversalAndScaling)
{
  model::SolverOptions<R> options(p());
  for (auto const& counts : {std::pair{1u, 1u}, std::pair{3u, 3u}, std::pair{5u, 3u}, std::pair{1u, 5u}})
    for (R c : {v(1) / 100, v(1), v(100)})
    {
      auto const a = model::ground_state(counts.first, counts.second, v(3), c, options);
      ASSERT_TRUE(a.converged);
      finite_state(a, p());
      original_equations(a);
      auto const b = model::ground_state(counts.second, counts.first, v(6), c / 2, options);
      ASSERT_TRUE(b.converged);
      EXPECT_REAL_NEAR(a.energy / 4, b.energy, v(4096) * uni20::numeric_limits<R>::epsilon(p()) * (v(1) + b.energy));
      EXPECT_EQ(a.quantum_numbers.charge, b.quantum_numbers.charge);
    }
}
TEST_P(GaudinYangMpfr, IncompleteContinuationReportsRequestedCoupling)
{
  for (std::size_t budget : {0, 1, 2})
  {
    model::SolverOptions<R> options(p());
    options.max_iterations = budget;
    auto const s = model::ground_state(3, 3, v(6), v(1) / 100, options);
    EXPECT_FALSE(s.converged);
    EXPECT_EQ(s.iterations, budget);
    EXPECT_EQ(s.status, model::SolveStatus::iteration_limit);
    EXPECT_GT(s.root_interaction, s.interaction);
    finite_state(s, p());
    original_equations(s);
  }
}

TEST(GaudinYangMpfrShared, TwoBodyAccuracyGrowthAndResolvedWeakCouplingCorrection)
{
  auto const reference_p = uni20::Precision::bits(512);
  R lo(0, reference_p), hi = uni20::pi<R>.at(reference_p), half = R(1, reference_p) / 2;
  // Independent two-body jump condition q*tan(q/2)=c/2 at c=L=1.
  for (int iteration = 0; iteration < 530; ++iteration)
  {
    R mid = (lo + hi) / 2;
    if (mid == lo || mid == hi) break;
    if (mid * uni20::math::tan(mid / 2) > half)
      hi = mid;
    else
      lo = mid;
  }
  R const expected = 2 * ((lo + hi) / 2) * ((lo + hi) / 2);
  R previous(1, reference_p);
  for (int bits : {128, 256, 333})
  {
    auto p = uni20::Precision::bits(bits);
    auto s = model::ground_state(1, 1, R(1, p), R(1, p), model::SolverOptions<R>(p));
    ASSERT_TRUE(s.converged);
    finite_state(s, p);
    R error = uni20::math::abs(s.energy.at(reference_p) - expected);
    EXPECT_LT(error, previous / R("1e15", reference_p));
    EXPECT_LT(error, (R(4096, p) * uni20::numeric_limits<R>::epsilon(p)).at(reference_p));
    previous = error;
  }
  auto p = uni20::Precision::bits(333);
  R c("1e-60", p);
  auto s = model::ground_state(1, 1, R(1, p), c, model::SolverOptions<R>(p));
  ASSERT_TRUE(s.converged);
  R correction = (s.energy / (2 * c) - 1) / c;
  EXPECT_REAL_NEAR(correction, -R(1, p) / 12, R("1e-35", p));
}

TEST(GaudinYangMpfrShared, OwningConcurrentSectorSweepAndFailureValues)
{
  auto const ambient = mpfr_get_default_prec();
  struct Request
  {
      std::size_t up, down;
      R length, c;
      model::SolverOptions<R> options;
  };
  auto submit = [](uni20::async::IAsyncScheduler& scheduler, auto& collected) {
    std::vector<Request> requests;
    for (int bits : {128, 256, 333})
    {
      auto p = uni20::Precision::bits(bits);
      model::SolverOptions<R> options(p);
      requests.push_back({5, 1, R(6, p), R(1, p), options});
      options.max_iterations = 0;
      requests.push_back({1, 5, R(6, p), R(1, p), options});
    }
    return bethe::evaluate_observable_grid_async(
        requests.size(),
        [requests](std::size_t i) {
          auto const& r = requests[i];
          return model::ground_state(r.up, r.down, r.length, r.c, r.options);
        },
        [&collected](std::size_t i, auto state) {
          EXPECT_EQ(i, collected.size());
          collected.push_back(std::move(state));
        },
        {&scheduler, 128, 2});
  };
  std::vector<model::State<R>> reference;
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  auto serial = submit(deferred, reference);
  EXPECT_EQ(serial.get_wait(deferred), 6u);
  for (int workers : {1, 4})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<model::State<R>> actual;
    auto complete = submit(scheduler, actual);
    EXPECT_EQ(complete.get_wait(scheduler), 6u);
    ASSERT_EQ(actual.size(), reference.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      same_state(actual[i], reference[i]);
      EXPECT_EQ(actual[i].converged, i % 2 == 0);
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
} // namespace
