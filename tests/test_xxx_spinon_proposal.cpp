// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "xxx_spinon_proposal_references.hpp"
#include <bethe/detail/xxx_spinon_proposal.hpp>
#include <bethe/sampling.hpp>

namespace
{
namespace math = uni20::math;
using Kind = bethe::heisenberg::detail::SpinonProposalKind;
template <class R> using Proposal = bethe::heisenberg::detail::SpinonJointProposal<R>;

template <class R> R parse(char const* text, R const& one)
{
  if constexpr (uni20::has_runtime_precision_v<R>)
    return R(text, one.precision());
  else
    return uni20::parse_real<R>(text);
}

template <class R> Proposal<R> make_proposal(R const& one, bool soft = false)
{
  auto p = [&](char const* x) { return parse(x, one); };
  return Proposal<R>(p(soft ? "1040" : "8"), p(soft ? ".5" : ".25"),
                     {{Kind::background, p(".2"), one, one},
                      {Kind::cauchy, p(".2"), p(soft ? "32" : "2"), one},
                      {Kind::one_pair, p(".3"), p(soft ? "128" : "4"), p(soft ? "2" : ".5")},
                      {Kind::two_pairs, p(".3"), p(soft ? "520" : "8"), p(soft ? "16" : "1")}});
}

template <class R> void references(R const& one)
{
  R const eps = uni20::numeric_limits<R>::epsilon(one);
  for (auto const& row : spinon_proposal_reference::points)
  {
    auto const proposal = make_proposal(one, row.soft);
    std::array<R, 6> beta;
    for (int i = 0; i < 6; ++i)
      beta[i] = parse(row.beta[i], one);
    std::array<R, 4> free{beta[0], beta[1], beta[2], beta[3]};
    R const density = proposal.log_density(free), jacobian = proposal.balanced_log_jacobian(beta);
    EXPECT_REAL_NEAR(density, parse(row.log_pdf, one), R{4096} * eps * (one + math::abs(density)));
    EXPECT_REAL_NEAR(jacobian, parse(row.log_jacobian, one), R{4096} * eps * (one + math::abs(jacobian)));
    if constexpr (uni20::has_runtime_precision_v<R>)
    {
      EXPECT_EQ(density.precision(), one.precision());
      EXPECT_EQ(jacobian.precision(), one.precision());
    }
  }
}

template <class R> void draws(R const& one)
{
  auto const proposal = make_proposal(one);
  R const eps = uni20::numeric_limits<R>::epsilon(one);
  for (auto const& row : spinon_proposal_reference::draws)
  {
    std::array<R, 5> u{parse(row.selector, one), one / R{8}, R{3} * one / R{8}, R{5} * one / R{8}, R{7} * one / R{8}};
    auto const sample = proposal.draw(u);
    for (int i = 0; i < 4; ++i)
      EXPECT_REAL_NEAR(sample.beta[i], parse(row.beta[i], one), R{512} * eps * (one + math::abs(sample.beta[i])));
    EXPECT_EQ(sample.log_density, proposal.log_density(sample.beta));
  }
}

template <class R> class SpinonProposal : public ::testing::Test {};
TYPED_TEST_SUITE(SpinonProposal, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(SpinonProposal, IndependentDensityAndCoareaReferences) { references(TypeParam{1}); }
TYPED_TEST(SpinonProposal, IndependentDrawsCoverEveryPairing) { draws(TypeParam{1}); }
TYPED_TEST(SpinonProposal, CompleteMixtureIsPermutationInvariant)
{
  using R = TypeParam;
  auto const p = make_proposal(R{1});
  std::array<R, 4> b{R{-3}, R{-1}, R{1}, R{2}};
  R const expected = p.log_density(b), eps = uni20::numeric_limits<R>::epsilon();
  do
  {
    EXPECT_REAL_NEAR(p.log_density(b), expected, R{128} * eps);
  }
  while (std::next_permutation(b.begin(), b.end()));
}
TYPED_TEST(SpinonProposal, InverseCdfAndClosedEndpointFailures)
{
  using R = TypeParam;
  auto const p = make_proposal(R{1}, true);
  R const eps = uni20::numeric_limits<R>::epsilon();
  for (int i = 1; i < 100; ++i)
  {
    R const u = R(i) / R{100};
    R const b = p.native_inverse(u);
    EXPECT_REAL_NEAR(p.native_cdf(b), u, R{64} * eps);
    EXPECT_LT(math::abs(b), p.cutoff());
  }
  EXPECT_THROW(p.native_inverse(R{0}), std::invalid_argument);
  EXPECT_THROW(p.native_inverse(R{1}), std::invalid_argument);
  EXPECT_THROW(p.native_log_density(p.cutoff()), std::invalid_argument);
  std::array<R, 5> u{R{1} / R{2}, R{1} / R{2}, R{1} / R{2}, R{1} / R{2}, R{0}};
  EXPECT_THROW(p.draw(u), std::invalid_argument);
  std::array<R, 6> degenerate{R{0}, R{0}, R{0}, R{0}, R{0}, R{0}};
  EXPECT_THROW(p.balanced_log_jacobian(degenerate), std::runtime_error);
  degenerate[0] = R{1};
  EXPECT_TRUE(uni20::isfinite(p.balanced_log_jacobian(degenerate)));
}

TEST(SpinonProposalMeasure, LiteralLabelOrbitMatchesBalancedIntensity)
{
  auto const p = make_proposal(1.0);
  std::array<double, 6> b{-2.1, -1.3, -.4, .3, 1.2, 2};
  double const expected = std::exp(-p.balanced_log_jacobian(b));
  double sum = 0;
  do
  {
    std::array<double, 4> free{b[0], b[1], b[2], b[3]};
    double h = std::exp(p.log_density(free)) * 2 * std::atan(1.) * std::abs(std::tanh(b[4]) - std::tanh(b[5]));
    for (auto const& x : free)
      h *= std::cosh(x);
    sum += h;
  }
  while (std::next_permutation(b.begin(), b.end()));
  EXPECT_NEAR(sum / (720 * expected), 1, 2e-13);
}

TEST(SpinonProposalControls, RejectsInvalidParameters)
{
  EXPECT_THROW((Proposal<double>(8, .25, {})), std::invalid_argument);
  EXPECT_THROW((Proposal<double>(8, 0, {{Kind::background, 1, 1, 1}})), std::invalid_argument);
  EXPECT_THROW((Proposal<double>(8, .25, {{Kind::background, .5, 1, 1}})), std::invalid_argument);
  EXPECT_THROW((Proposal<double>(8, .25, {{Kind::cauchy, 1, 1, 1}})), std::invalid_argument);
  EXPECT_THROW((Proposal<double>(8, .25, {{Kind::background, .5, 1, 1}, {Kind::one_pair, .5, 1, 0}})),
               std::invalid_argument);
}

#if UNI20_HAS_MPREAL
TEST(SpinonProposalMpfr, FreshMidpointsReachRepresentableQuantiles)
{
  using R = uni20::mpreal;
  for (int bits : {128, 256, 384})
  {
    SCOPED_TRACE(bits);
    R const one(1, uni20::Precision::bits(bits));
    Proposal<R> proposal(8 * one, one / 2, {{Kind::background, one, one, one}});
    bethe::MonteCarloStream<5> stream(42, 1, 52);
    for (int sample = 0; sample < 8; ++sample)
    {
      SCOPED_TRACE(sample);
      auto u = stream.next().at_precision(one);
      for (int coordinate = 1; coordinate < 5; ++coordinate)
      {
        SCOPED_TRACE(coordinate);
        SCOPED_TRACE(uni20::format_scalar(u[coordinate]));
        R beta;
        ASSERT_NO_THROW(beta = proposal.native_inverse(u[coordinate]));
        // Re-evaluate the CDF independently at a much higher precision. The
        // attainable residual includes the representable beta spacing times
        // the PDF, not just a relative error on a small uniform coordinate.
        auto hp = uni20::Precision::bits(768);
        R const b(beta, hp), input(u[coordinate], hp), h(1, hp), pi = 4 * math::atan(h);
        R const lower = 2 * math::atan(math::exp(-8 * h));
        R const norm = 2 * math::atan(4 * h);
        R const cdf =
            (2 * math::atan(math::exp(b)) - lower) / (2 * (pi - 2 * lower)) + (h / 2 + math::atan(b / 2) / norm) / 2;
        R const pdf = h / (2 * (pi - 2 * lower) * math::cosh(b)) + h / ((4 * h + b * b) * norm);
        R const bound = 256 * R(uni20::numeric_limits<R>::epsilon(one), hp) *
                        (std::min(input, h - input) + pdf * std::max(h, math::abs(b)));
        EXPECT_LE(math::abs(cdf - input), bound);
        EXPECT_LT(math::abs(beta), 8 * one);
      }
    }
  }
}

TEST(SpinonProposalMpfr, ExplicitPrecisionAndReferenceGrowth)
{
  using R = uni20::mpreal;
  auto const high = uni20::Precision::bits(768);
  R previous(1, high);
  for (int bits : {128, 256, 384})
  {
    SCOPED_TRACE(bits);
    auto const precision = uni20::Precision::bits(bits);
    R const one(1, precision);
    references(one);
    draws(one);
    auto const proposal = make_proposal(one, true);
    auto const& row = spinon_proposal_reference::points[1];
    std::array<R, 6> beta;
    for (int i = 0; i < 6; ++i)
      beta[i] = R(row.beta[i], precision);
    R const error = math::abs(R(proposal.balanced_log_jacobian(beta), high) - R(row.log_jacobian, high));
    EXPECT_LT(error, previous / 1000);
    EXPECT_LT(error, R(4096, high) * R(uni20::numeric_limits<R>::epsilon(precision), high) * R(3000, high));
    previous = error;
    EXPECT_THROW(proposal.native_inverse(R(".2", uni20::Precision::bits(bits + 32))), std::invalid_argument);
    beta.fill(R(0, precision));
    EXPECT_THROW(proposal.balanced_log_jacobian(beta), std::runtime_error);
    beta[0] = one;
    EXPECT_TRUE(uni20::isfinite(proposal.balanced_log_jacobian(beta)));
  }
  EXPECT_LT(previous, R("1e-100", high));
  EXPECT_THROW(make_proposal(R{1}), std::invalid_argument);
}
#endif
} // namespace
