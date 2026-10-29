// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include "xyz_finite_mpfr_references.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/xyz.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using C = uni20::complex<R>;
using P = uni20::Precision;
namespace model = bethe::xyz;
namespace math = uni20::math;
auto fields(model::GroundState<R> const& s)
{
  return std::tie(s.sites, s.iterations, s.eta, s.t, s.residual_norm, s.momentum, s.exchange.x, s.exchange.y,
                  s.exchange.z, s.quantum_numbers, s.roots, s.energy, s.converged, s.status);
}
void precision(model::GroundState<R> const& s, P p)
{
  for (auto const& x : {s.eta, s.t, s.residual_norm, s.momentum, s.exchange.x, s.exchange.y, s.exchange.z})
    EXPECT_EQ(uni20::precision_of(x), p);
  for (auto const& z : s.roots)
  {
    EXPECT_EQ(uni20::precision_of(z), p);
    EXPECT_EQ(uni20::precision_of(z.real()), p);
    EXPECT_EQ(uni20::precision_of(z.imag()), p);
  }
  if (s.energy) EXPECT_EQ(uni20::precision_of(*s.energy), p);
}
TEST(XYZMpfr, IndependentSpinHamiltonianAndFreeFermionAccuracy)
{
  auto hp = P::bits(768);
  R previous(1, hp);
  auto ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256, 333})
  {
    auto p = P::bits(bits);
    model::SolverOptions<R> controls(p);
    R worst(0, hp);
    for (auto const& f : xyz_finite_reference::fixtures)
    {
      SCOPED_TRACE(::testing::Message() << f.sites << ' ' << f.eta << ' ' << bits);
      auto s = model::ground_state(f.sites, R(f.eta, p), R(f.t, p), controls);
      ASSERT_TRUE(s.converged);
      ASSERT_TRUE(s.energy);
      precision(s, p);
      R ref(f.energy, hp);
      worst = std::max(worst, math::abs(s.energy->at(hp) - ref) / (1 + math::abs(ref)));
      EXPECT_LE(s.residual_norm, controls.residual_tolerance);
    }
    EXPECT_LT(worst, 4096 * uni20::numeric_limits<R>::epsilon(p).at(hp));
    EXPECT_LT(worst, previous * R("1e-15", hp));
    if (bits >= 256) EXPECT_LT(worst, R("1e-70", hp));
    previous = worst;
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
TEST(XYZMpfr, AnalyticTwoSitesAndFailurePrecision)
{
  for (int bits : {80, 160, 256})
  {
    auto p = P::bits(bits);
    model::SolverOptions<R> controls(p);
    controls.max_iterations = 0;
    for (char const* eta : {"0.2", "0.5", "0.8"})
    {
      auto two = model::ground_state(2, R(eta, p), R("0.7", p), controls);
      ASSERT_TRUE(two.converged);
      precision(two, p);
      EXPECT_LT(math::abs(*two.energy + (two.exchange.x + two.exchange.y + two.exchange.z) / 2),
                1024 * uni20::numeric_limits<R>::epsilon(p));
    }
    auto failed = model::ground_state(8, R("0.4", p), R("0.7", p), controls);
    precision(failed, p);
    EXPECT_EQ(failed.status, model::Status::iteration_limit);
    EXPECT_FALSE(failed.energy);
    EXPECT_FALSE(failed.converged);
    R eta("0.4", p), t("0.7", p);
    for (R bad : {R(uni20::uninitialized), R(1), R(1, P::bits(64)), R("nan", p)})
    {
      EXPECT_THROW((void)model::ground_state(4, bad, t, controls), std::invalid_argument);
      EXPECT_THROW((void)model::ground_state(4, eta, bad, controls), std::invalid_argument);
      auto invalid = controls;
      invalid.residual_tolerance = bad;
      EXPECT_THROW((void)model::ground_state(4, eta, t, invalid), std::invalid_argument);
    }
    controls.max_iterations = 100;
    auto valid = model::ground_state(4, eta, t, controls);
    ASSERT_TRUE(valid.converged);
    valid.roots[0] = valid.roots[0].at(P::bits(64));
    EXPECT_THROW((void)model::candidate_energy<R>(4, valid.roots, eta, t), std::invalid_argument);
  }
}
TEST(XYZMpfr, MultiplicativeEquationsAndJacobian)
{
  auto p = P::bits(160);
  R eta("0.4", p), t("0.7", p), zero(0, p), one(1, p), eps = uni20::numeric_limits<R>::epsilon(p);
  model::SolverOptions<R> options(p);
  auto state = model::ground_state(8, eta, t, options);
  ASSERT_TRUE(state.converged);
  auto ratio = [&](C const& a, C const& b) {
    return uni20::elliptic_theta_ratio(uni20::elliptic_theta_jet(1, a, t), 0, uni20::elliptic_theta_jet(1, b, t), 0);
  };
  for (std::size_t j = 0; j < state.roots.size(); ++j)
  {
    C product(one, zero), drive = ratio(state.roots[j] + eta / 2, state.roots[j] - eta / 2);
    for (unsigned k = 0; k < 8; ++k)
      product *= drive;
    for (std::size_t k = 0; k < state.roots.size(); ++k)
      if (k != j) product *= ratio(state.roots[j] - state.roots[k] - eta, state.roots[j] - state.roots[k] + eta);
    EXPECT_LT(math::abs(product - one), 16384 * eps);
  }
  model::detail::GroundSystem<R> system{8, eta, t};
  std::vector<R> x{state.roots[2].imag(), state.roots[3].imag()};
  uni20::DenseMatrix<R> jac(2, 2, p);
  system.evaluate(x, &jac);
  R delta = math::ldexp(one, -50);
  for (std::size_t col = 0; col < 2; ++col)
  {
    auto plus = x, minus = x;
    plus[col] += delta;
    minus[col] -= delta;
    auto a = system.evaluate(plus), b = system.evaluate(minus);
    for (std::size_t row = 0; row < 2; ++row)
      EXPECT_LT(math::abs((a.residual[row] - b.residual[row]) / (2 * delta) - jac[row, col]), 16384 * delta * delta);
  }
}
struct Request
{
    std::size_t sites;
    R eta, t;
    model::SolverOptions<R> options;
    auto operator()() const { return model::ground_state(sites, eta, t, options); }
};
TEST(XYZMpfr, OwningMixedPrecisionSweepsAndFailureIsolation)
{
  auto ambient = mpfr_get_default_prec();
  auto check = [&](uni20::async::IAsyncScheduler& scheduler) {
    std::vector<Request> requests;
    for (int bits : {80, 160, 240})
    {
      auto p = P::bits(bits);
      model::SolverOptions<R> controls(p);
      requests.push_back({8, R("0.4", p), R("0.7", p), controls});
      controls.max_iterations = 0;
      requests.push_back({8, R("0.8", p), R(1, p), controls});
    }
    std::vector<model::GroundState<R>> expected, actual;
    for (auto const& r : requests)
      expected.push_back(r());
    auto count = requests.size();
    auto done = bethe::evaluate_observable_grid_async(
        count, [requests = std::move(requests)](std::size_t i) { return requests[i](); },
        [&](std::size_t i, auto state) {
          EXPECT_EQ(i, actual.size());
          actual.push_back(std::move(state));
        },
        {.scheduler = &scheduler, .max_in_flight = 2});
    auto badp = P::bits(128);
    auto invalid = bethe::evaluate_observable_grid_async(
        1,
        [badp](std::size_t) {
          return model::ground_state(3, R("0.4", badp), R(1, badp), model::SolverOptions<R>(badp));
        },
        [](std::size_t, auto) { ADD_FAILURE(); }, {.scheduler = &scheduler, .max_in_flight = 1});
    EXPECT_THROW((void)invalid.get_wait(scheduler), std::invalid_argument);
    EXPECT_EQ(done.get_wait(scheduler), count);
    ASSERT_EQ(actual.size(), count);
    for (std::size_t i = 0; i < count; ++i)
    {
      EXPECT_TRUE(fields(actual[i]) == fields(expected[i]));
      precision(actual[i], actual[i].eta.precision());
      EXPECT_EQ(actual[i].converged, i % 2 == 0);
    }
  };
  uni20::async::DebugScheduler debug({.order = uni20::async::DebugSchedulerOrder::reverse});
  check(debug);
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    check(scheduler);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
} // namespace
