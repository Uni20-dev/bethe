// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <array>
#include <bethe/observable_grid.hpp>
#include <bethe/potts.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using C = uni20::complex<R>;
using uni20::Precision;
namespace model = bethe::potts;
namespace math = uni20::math;

void precision(model::State<R> const& s, Precision p)
{
  for (auto const* v : {&s.energy, &s.momentum, &s.residual_norm, &s.auxiliary_twist})
  {
    EXPECT_EQ(uni20::precision_of(*v), p);
    EXPECT_TRUE(uni20::isfinite(*v));
  }
  for (auto const& z : s.auxiliary_roots)
  {
    EXPECT_EQ(uni20::precision_of(z), p);
    EXPECT_TRUE(uni20::isfinite(z));
  }
}

auto fields(model::State<R> const& s)
{
  return std::tie(s.energy, s.momentum, s.residual_norm, s.sites, s.momentum_index, s.iterations, s.charge, s.converged,
                  s.auxiliary_roots, s.auxiliary_numbers, s.auxiliary_twist);
}

TEST(PottsMpfr, IndependentClockEnergiesGainAccuracyBeyondBinary128)
{
  auto const refp = Precision::bits(768);
  // Literal charge/momentum clock-Hamiltonian blocks, independent of Bethe
  // equations: scripts/reference_potts.py --sites 3 --digits 150 (mpmath 1.2.1).
  std::array<R, 3> const three_site = {R("-7."
                                         "68465843842649082473211478396111553772079883806043065159795035964243151950643"
                                         "239038179547621602644502681606140882121242295213126062380254461932221862",
                                         refp),
                                       R("-6."
                                         "92434399202024893136225131700422238789909633492884963383736869119647176861550"
                                         "31230459536457408369973802666865303219254020070193642639097800294222744",
                                         refp),
                                       R("-2."
                                         "37755164192302087519051994482424675727962147887396605970383119177753613524831"
                                         "604589633781382129967773666787234648625191317219205925550260959543697823",
                                         refp)};
  std::array<R, 2> const two_site = {-2 - 2 * math::sqrt(R(3, refp)), -(1 + math::sqrt(R(57, refp))) / 2};
  R previous(1, refp);
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256, 333})
  {
    SCOPED_TRACE(bits);
    auto const p = Precision::bits(bits);
    bethe::SolverOptions<R> options(p);
    std::array const actual{model::ground_state<R>(2, options), model::charged_one_hole<R>(2, 1, 0, options),
                            model::ground_state<R>(3, options), model::charged_one_hole<R>(3, 1, 0, options),
                            model::charged_one_hole<R>(3, -1, 1, options)};
    std::array<R, 5> const expected{two_site[0], two_site[1], three_site[0], three_site[1], three_site[2]};
    R error(0, refp);
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      ASSERT_TRUE(actual[i].converged);
      precision(actual[i], p);
      EXPECT_LE(actual[i].residual_norm, options.residual_tolerance);
      error = std::max(error, math::abs(actual[i].energy.at(refp) - expected[i]));
    }
    EXPECT_LT(error, 8192 * uni20::numeric_limits<R>::epsilon(p).at(refp));
    EXPECT_LT(error, previous * R("1e-15", refp));
    if (bits >= 256) EXPECT_LT(error, R("1e-70", refp));
    previous = error;
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(PottsMpfr, ExponentiatedTwistedEquationsAndChargePartners)
{
  auto const p = Precision::bits(160);
  bethe::SolverOptions<R> options(p);
  R const zero(0, p), one(1, p), delta = math::sqrt(R(3, p)) / 2, eps = uni20::numeric_limits<R>::epsilon(p);
  R const pi = 4 * math::atan(one);
  for (std::size_t sites : {3, 4, 7})
    for (bool vacuum : {true, false})
    {
      auto const k = sites / 2;
      auto const s = vacuum ? model::ground_state<R>(sites, options) : model::charged_one_hole<R>(sites, 1, k, options);
      ASSERT_TRUE(s.converged);
      precision(s, p);
      auto const& z = s.auxiliary_roots;
      for (std::size_t i = 0; i < sites; ++i)
      {
        C equation(math::cos(s.auxiliary_twist), math::sin(s.auxiliary_twist));
        C const drive = C(z[i], one) / C(z[i], -one);
        for (std::size_t n = 0; n < 2 * sites; ++n)
          equation *= drive;
        for (std::size_t j = 0; j < sites; ++j)
          if (i != j)
          {
            R const a = delta * (z[i] - z[j]), b = 1 + delta - (1 - delta) * z[i] * z[j];
            equation *= C(a, -b) / C(a, b);
          }
        EXPECT_LT(math::abs(equation - C(one, zero)), 512 * R(sites) * eps);
      }
      if (!vacuum)
      {
        auto const partner = model::charged_one_hole<R>(sites, -1, k, options);
        auto const reflected = model::charged_one_hole<R>(sites, 1, sites - k, options);
        EXPECT_EQ(s.energy, partner.energy);
        EXPECT_EQ(s.auxiliary_roots, partner.auxiliary_roots);
        EXPECT_LT(math::abs(s.energy - reflected.energy), 4096 * R(sites) * eps);
        EXPECT_LT(math::abs(s.momentum - 2 * pi * R(k) / R(sites)), 16 * eps);
      }
    }
}

TEST(PottsMpfr, ExplicitConstantsAndControlFailures)
{
  auto const p = Precision::bits(256), refp = Precision::bits(768);
  R const v = model::velocity<R>(p), e = model::bulk_energy_density<R>(p);
  EXPECT_EQ(uni20::precision_of(v), p);
  EXPECT_EQ(uni20::precision_of(e), p);
  R const pi = 4 * math::atan(R(1, refp)), root3 = math::sqrt(R(3, refp));
  EXPECT_LT(math::abs(v.at(refp) - 3 * root3 / 2), R("1e-75", refp));
  EXPECT_LT(math::abs(e.at(refp) + R(4, refp) / 3 + 2 * root3 / pi), R("1e-75", refp));
  EXPECT_THROW(model::velocity<R>(Precision::exact()), std::invalid_argument);
  EXPECT_THROW(model::bulk_energy_density<R>(Precision::exact()), std::invalid_argument);
  for (R const& bad : {R{}, R(1), R(1, Precision::bits(80)), R(0, p), R(-1, p), R("nan", p), R("inf", p)})
  {
    bethe::SolverOptions<R> options(p);
    options.residual_tolerance = bad;
    EXPECT_THROW(model::ground_state<R>(3, options), std::invalid_argument);
  }
  bethe::SolverOptions<R> options(p);
  options.max_iterations = 0;
  auto const initial = model::charged_one_hole<R>(4, 1, 1, options);
  EXPECT_FALSE(initial.converged);
  EXPECT_EQ(initial.iterations, 0u);
  precision(initial, p);
}

TEST(PottsMpfr, OwningMixedPrecisionGridsRetainCompleteResults)
{
  struct Request
  {
      std::size_t sites, k;
      int charge;
      bethe::SolverOptions<R> options;
  };
  auto const ambient = mpfr_get_default_prec();
  auto submit = [](uni20::async::IAsyncScheduler& scheduler, auto& output) {
    std::vector<Request> requests;
    for (int bits : {80, 160, 256})
    {
      bethe::SolverOptions<R> controls(Precision::bits(bits));
      requests.push_back({3, 0, 0, controls});
      requests.push_back({4, 1, -1, controls});
      controls.max_iterations = 0;
      requests.push_back({4, 2, 1, controls});
    }
    return bethe::evaluate_observable_grid_async(
        requests.size(),
        [owned = requests](std::size_t i) {
          auto const& r = owned[i];
          return r.charge ? model::charged_one_hole<R>(r.sites, r.charge, r.k, r.options)
                          : model::ground_state<R>(r.sites, r.options);
        },
        [&output](std::size_t i, auto result) {
          EXPECT_EQ(i, output.size());
          output.push_back(std::move(result));
        },
        {.scheduler = &scheduler, .max_in_flight = 2});
  };
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  std::vector<model::State<R>> expected;
  auto done = submit(deferred, expected);
  EXPECT_EQ(done.get_wait(deferred), 9u);
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<model::State<R>> actual;
    auto final = submit(scheduler, actual);
    EXPECT_EQ(final.get_wait(scheduler), expected.size());
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      EXPECT_TRUE(fields(actual[i]) == fields(expected[i]));
      precision(actual[i], Precision::bits(i < 3 ? 80 : i < 6 ? 160 : 256));
      EXPECT_EQ(actual[i].converged, i % 3 != 2);
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
} // namespace
