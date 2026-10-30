// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/detail/xxx_spinon_moment_proposal.hpp>
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

TEST(SpinonMomentProposal, FullLabelledMixtureAndConservation)
{
  auto const joint = proposal(1.);
  double const pi = std::numbers::pi, k = pi / 2;
  d::SpinonFixedMomentumProposal p(joint, k);
  std::array<double, 6> const beta{-2.1, -1.3, -.4, .3, 1.2, 2};
  auto const logs = p.chart_log_intensities(beta);
  double total = 0;
  std::size_t a = 0;
  for (int i = 0; i < 6; ++i)
    for (int j = i + 1; j < 6; ++j)
    {
      std::array<double, 4> free;
      std::size_t next = 0;
      for (int m = 0; m < 6; ++m)
        if (m != i && m != j) free[next++] = beta[m];
      for (auto [eliminated, native] : {std::pair{i, j}, std::pair{j, i}})
      {
        double value = std::exp(joint.log_density(free) + joint.native_log_density(beta[native]));
        for (int m = 0; m < 6; ++m)
          if (m != eliminated) value *= std::cosh(beta[m]);
        EXPECT_NEAR(std::exp(logs[a++]), value, 64 * std::numeric_limits<double>::epsilon() * value);
        total += value;
      }
    }
  EXPECT_NEAR(std::exp(p.log_momentum_measure(beta)), 30 / total,
              128 * std::numeric_limits<double>::epsilon() * 30 / total);
  std::array<int, 6> order{0, 1, 2, 3, 4, 5};
  do
  {
    std::array<double, 6> permuted;
    for (int i = 0; i < 6; ++i)
      permuted[i] = beta[order[i]];
    EXPECT_NEAR(p.log_momentum_measure(permuted), p.log_momentum_measure(beta), 2e-14);
  }
  while (std::next_permutation(order.begin(), order.end()));
  std::size_t supported = 0;
  for (unsigned chart = 0; chart < 30; ++chart)
  {
    auto const point = p.draw({(chart + .375) / 30, .1, .2, .7, .9, .3});
    EXPECT_EQ(point.eliminated, chart / 5);
    EXPECT_NE(point.eliminated, point.native);
    if (!point.inside) continue;
    ++supported;
    double momentum = 0, energy = 0;
    for (auto b : point.beta)
    {
      momentum += 2 * std::atan(std::exp(b));
      energy += pi / (2 * std::cosh(b));
    }
    EXPECT_NEAR(momentum, k + 2 * pi * point.section, 2e-14);
    EXPECT_NEAR(std::exp(point.log_energy), energy, 2e-14);
  }
  EXPECT_GT(supported, 0u);

  d::SpinonJointProposal<double> flat(3., 1e-6, {{Kind::background, 1., 1., 1.}});
  d::SpinonFixedMomentumProposal<double> sections(flat, k);
  for (int section = 0; section < 3; ++section)
  {
    double const momentum = (k + 2 * pi * section) / 6;
    double const b = std::log(std::tan(momentum / 2)), u = flat.native_cdf(b);
    auto const point = sections.draw({.375 / 30, u, u, u, u, u});
    ASSERT_TRUE(point.inside);
    EXPECT_EQ(point.section, section);
    EXPECT_NEAR(point.beta[point.eliminated], b, 2e-12);
  }
}

TEST(SpinonMomentProposal, PeriodicBoxSplineVolumeAndLabelledMomentum)
{
  // Independent analytic convolution of six uniform-momentum intervals.
  // The proposal is nearly uniform in momentum but still a complete mixture.
  double const pi = std::numbers::pi, cutoff = 3., k = pi / 2;
  d::SpinonJointProposal<double> joint(cutoff, 1e-6, {{Kind::background, 1., 1., 1.}});
  d::SpinonFixedMomentumProposal p(joint, k);
  double const lo = 2 * std::atan(std::exp(-cutoff)), width = pi - 2 * lo;
  std::array<double, 2> exact{};
  constexpr std::array<int, 7> choose{1, 6, 15, 20, 15, 6, 1};
  for (int section = 0; section < 3; ++section)
  {
    double const total = k + 2 * pi * section;
    long double sum = 0;
    for (int j = 0; j <= 6; ++j)
    {
      long double const t = std::max(0.L, (static_cast<long double>(total) - 6 * lo) / width - j);
      sum += (j % 2 ? -1 : 1) * choose[j] * std::pow(t, 5);
    }
    double const volume = std::pow(width, 5) * sum / 120;
    exact[0] += volume;
    exact[1] += volume * total / 6; // each labelled p has the same conditional mean
  }
  auto observe = [p](bethe::UniformSample<6> const& identity, std::stop_token) {
    auto const point = p.draw(identity.at_precision(1.));
    std::vector<double> values(2, 0.);
    if (point.inside)
    {
      values[0] = std::exp(point.log_momentum_measure);
      values[1] = values[0] * 2 * std::atan(std::exp(point.beta[0]));
    }
    return bethe::MonteCarloObservation<double>{values};
  };
  auto result = bethe::monte_carlo_async<6>({16, 512, 120}, 2, 1., observe).get_wait();
  ASSERT_EQ(result.status, bethe::MonteCarloStatus::complete);
  ASSERT_TRUE(result.mean && result.covariance);
  for (int i = 0; i < 2; ++i)
    EXPECT_LT(std::abs((*result.mean)[i] - exact[i]),
              6 * std::sqrt((*result.covariance)[2 * i + i]) + 1e-10 * exact[i]);
}

#if UNI20_HAS_MPREAL
TEST(SpinonMomentProposalMpfr, CutoffResolutionAndPrecisionAreExplicit)
{
  using R = uni20::mpreal;
  R const low_one(1, uni20::Precision::bits(384));
  d::SpinonJointProposal<R> low(1040 * low_one, low_one / 2, {{Kind::background, low_one, low_one, low_one}});
  EXPECT_THROW((d::SpinonFixedMomentumProposal<R>(low, 2 * uni20::math::atan(low_one))), std::runtime_error);
  std::optional<std::array<R, 6>> previous;
  for (unsigned bits : {1800u, 2048u})
  {
    R const one(1, uni20::Precision::bits(bits)), pi = 4 * uni20::math::atan(one);
    d::SpinonJointProposal<R> joint(1040 * one, one / 2, {{Kind::background, one, one, one}});
    d::SpinonFixedMomentumProposal<R> p(joint, pi / 2);
    auto const point = p.draw({one / 101, one / 10, one / 5, 7 * one / 10, 9 * one / 10, 3 * one / 10});
    ASSERT_TRUE(point.inside);
    EXPECT_EQ(point.log_momentum_measure.precision(), one.precision());
    R sum = 0 * one;
    for (auto const& b : point.beta)
      sum += 2 * uni20::math::atan(uni20::math::exp(b));
    EXPECT_LT(uni20::math::abs(sum - pi / 2 - 2 * point.section * pi), R("1e-500", one.precision()));
    if (previous)
      for (int i = 0; i < 6; ++i)
        EXPECT_LT(uni20::math::abs(point.beta[i] - (*previous)[i].at(one.precision())), R("1e-490", one.precision()));
    previous = point.beta;
  }
}
#endif
