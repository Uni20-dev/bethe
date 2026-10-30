// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/detail/xxx_spinon_proposal.hpp>
#include <bethe/monte_carlo.hpp>

namespace
{
namespace d = bethe::heisenberg::detail;
using Kind = d::SpinonProposalKind;
template <class R> d::SpinonJointProposal<R> proposal(R const& one)
{
  return {8 * one,
          one / 2,
          {{Kind::background, one / 2, one, one},
           {Kind::one_pair, one / 4, 4 * one, one / 2},
           {Kind::two_pairs, one / 4, 8 * one, one}}};
}
} // namespace

TEST(SpinonGlobalProposal, FullMixtureAndEveryNativePair)
{
  auto p = proposal(1.);
  std::array<double, 6> const beta{-2.1, -1.3, -.4, .3, 1.2, 2};
  auto logs = p.global_chart_log_densities(beta);
  double density = 0, jacobian = 1;
  for (double b : beta)
    jacobian /= std::cosh(b);
  std::size_t index = 0;
  for (int i = 0; i < 6; ++i)
    for (int j = i + 1; j < 6; ++j)
    {
      std::array<double, 4> free;
      std::size_t next = 0;
      for (int k = 0; k < 6; ++k)
        if (k != i && k != j) free[next++] = beta[k];
      double direct = std::exp(p.log_density(free)) * std::exp(p.native_log_density(beta[i])) *
                      std::exp(p.native_log_density(beta[j]));
      EXPECT_NEAR(std::exp(logs[index++]), direct, 32 * std::numeric_limits<double>::epsilon() * direct);
      density += direct / 15;
    }
  EXPECT_NEAR(std::exp(p.global_log_momentum_measure(beta)), jacobian / density,
              64 * std::numeric_limits<double>::epsilon() * jacobian / density);
  std::array<int, 6> order{0, 1, 2, 3, 4, 5};
  do
  {
    std::array<double, 6> b;
    for (int j = 0; j < 6; ++j)
      b[j] = beta[order[j]];
    EXPECT_NEAR(p.global_log_momentum_measure(b), p.global_log_momentum_measure(beta), 1e-14);
  }
  while (std::next_permutation(order.begin(), order.end()));
  for (unsigned pair = 0; pair < 15; ++pair)
  {
    std::array<double, 7> u{(pair + .375) / 15, .1, .2, .7, .9, .3, .8};
    auto point = p.draw_global(u);
    EXPECT_EQ(point.native_pair, pair);
    double e = 0;
    for (auto b : point.beta)
      e += std::numbers::pi / (2 * std::cosh(b));
    EXPECT_NEAR(std::exp(point.log_energy), e, 1e-14);
    EXPECT_NEAR(point.log_momentum_measure, p.global_log_momentum_measure(point.beta), 1e-14);
  }
}

TEST(SpinonGlobalProposal, DirectLabelledVolumeAndEnergyHaveCorrectMeasures)
{
  using R = double;
  // The nearly uniform-momentum background makes the volume control very
  // sharp; the energy and nonsymmetric label test still have sampling error.
  d::SpinonJointProposal<R> p(8., 1e-6, {{Kind::background, 1., 1., 1.}});
  double const pi = std::numbers::pi, lo = 2 * std::atan(std::exp(-8.)), width = pi - 2 * lo;
  auto observe = [p, pi](bethe::UniformSample<7> const& identity, std::stop_token) {
    auto point = p.draw_global(identity.at_precision(1.));
    double measure = std::exp(point.log_momentum_measure) / (4 * pi * pi);
    double momentum0 = 2 * std::atan(std::exp(point.beta[0]));
    return bethe::MonteCarloObservation<double>{
        std::vector{measure, measure * std::exp(point.log_energy), measure * momentum0}};
  };
  auto pending = bethe::monte_carlo_async<7>({16, 256, 119}, 3, 1., observe, {.max_in_flight = 1});
  auto result = pending.get_wait();
  ASSERT_EQ(result.status, bethe::MonteCarloStatus::complete);
  ASSERT_TRUE(result.mean && result.covariance);
  std::array<double, 3> exact{std::pow(width, 6) / (4 * pi * pi),
                              6 * pi * std::cos(lo) * std::pow(width, 5) / (4 * pi * pi),
                              std::pow(width, 6) / (8 * pi)};
  for (int i = 0; i < 3; ++i)
    EXPECT_LT(std::abs((*result.mean)[i] - exact[i]),
              6 * std::sqrt((*result.covariance)[3 * i + i]) + 1e-12 * exact[i]);
}

#if UNI20_HAS_MPREAL
TEST(SpinonGlobalProposalMpfr, LogMeasureRetainsDeepSoftRange)
{
  using R = uni20::mpreal;
  for (unsigned bits : {160u, 256u})
  {
    R one(1, uni20::Precision::bits(bits));
    d::SpinonJointProposal<R> p(1040 * one, one / 2, {{Kind::background, one, one, one}});
    std::array<R, 6> beta{-1000 * one, -900 * one, -800 * one, 800 * one, 900 * one, 1000 * one};
    R log = p.global_log_momentum_measure(beta);
    EXPECT_TRUE(uni20::isfinite(log));
    EXPECT_EQ(log.precision(), one.precision());
    EXPECT_LT(log, -5000 * one);
    bethe::MonteCarloStream<7> stream(42, 0, 52);
    auto draw = p.draw_global(stream.next().at_precision(one));
    EXPECT_EQ(draw.log_momentum_measure.precision(), one.precision());
    EXPECT_TRUE(uni20::isfinite(draw.log_energy));
  }
}
#endif
