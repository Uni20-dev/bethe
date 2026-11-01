// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "xxx_spinon_surface_references.hpp"
#include <bethe/detail/xxx_spinon_proposal.hpp>
#include <bethe/detail/xxx_spinon_surface.hpp>
#include <bethe/sampling.hpp>

namespace
{
namespace math = uni20::math;
template <class R> using Surface = bethe::heisenberg::detail::SpinonEnergySurface<R>;

template <class R> R parse(char const* text, R const& one)
{
  if constexpr (uni20::has_runtime_precision_v<R>)
    return R(text, one.precision());
  else
    return uni20::parse_real<R>(text);
}

template <class R> void check_reference(spinon_surface_reference::Point const& row, R const& one)
{
  SCOPED_TRACE(row.n);
  SCOPED_TRACE(row.section);
  SCOPED_TRACE(row.cutoff);
  R const pi = R{4} * math::atan(one), eps = uni20::numeric_limits<R>::epsilon(one);
  std::vector<R> beta;
  for (int i = 0; i < row.n; ++i)
    beta.push_back(parse(row.beta[i], one));
  Surface<R> const surface(row.n, R(row.cutoff) * one, pi * parse(row.q_over_pi, one), row.section,
                           parse(row.omega, one));
  auto const result = surface.solve(std::span<R const>(beta).first(row.n - 2));
  ASSERT_TRUE(result.inside);
  EXPECT_EQ(result.refinements, 0u);
  EXPECT_EQ(result.working_bits, uni20::numeric_limits<R>::digits(one));
  EXPECT_REAL_NEAR(result.log_coarea_jacobian, parse(row.log_jacobian, one), R{65536} * eps);
  for (int root = 0; root < 2; ++root)
  {
    ASSERT_EQ(result.beta[root].size(), beta.size());
    for (int i = 0; i < row.n - 2; ++i)
      EXPECT_EQ(result.beta[root][i], beta[i]);
    // All reference hard pairs are supplied in ascending rapidity order.
    EXPECT_REAL_NEAR(result.beta[root][row.n - 2], beta[row.n - 1 - root], R{65536} * eps);
    EXPECT_REAL_NEAR(result.beta[root][row.n - 1], beta[row.n - 2 + root], R{65536} * eps);
    if constexpr (uni20::has_runtime_precision_v<R>)
      for (auto const& b : result.beta[root])
        EXPECT_EQ(b.precision(), one.precision());
  }
  if constexpr (uni20::has_runtime_precision_v<R>) EXPECT_EQ(result.log_coarea_jacobian.precision(), one.precision());
}

template <class R> class SpinonSurface : public ::testing::Test {};
TYPED_TEST_SUITE(SpinonSurface, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(SpinonSurface, IndependentForwardKinematicsAndBothRoots)
{
  for (auto const& row : spinon_surface_reference::points)
    if (row.cutoff == 8 || math::exp(TypeParam{-1040}) > TypeParam{0}) check_reference(row, TypeParam{1});
}

TYPED_TEST(SpinonSurface, SupportZerosAndUnresolvedBoundariesDiffer)
{
  using R = TypeParam;
  R const pi = R{4} * math::atan(R{1}), q = pi / R{2};
  EXPECT_FALSE((Surface<R>(2, R{8}, q, 0, R{10}).solve({}).inside));
  EXPECT_FALSE((Surface<R>(2, R{8}, q, 0, R{1}).solve({}).inside));
  EXPECT_THROW((Surface<R>(2, R{8}, q, 0, pi * math::sin(q / R{2})).solve({})), std::runtime_error);
  if (math::exp(R{-1040}) == R{0})
    EXPECT_THROW((Surface<R>(2, R{1040}, q, 0, R{2})), std::runtime_error);
  else
    EXPECT_TRUE((Surface<R>(2, R{1040}, q, 0, R{2}).solve({}).inside));
  EXPECT_THROW((Surface<R>(2, R{8}, q, 1, R{2})), std::invalid_argument);
  EXPECT_THROW((Surface<R>(3, R{8}, q, 0, R{2})), std::invalid_argument);
  EXPECT_THROW((Surface<R>(6, R{8}, q, 0, R{2}).solve({})), std::invalid_argument);
  std::array<R, 4> bad{R{-8}, R{-2}, R{-1}, R{0}};
  EXPECT_THROW((Surface<R>(6, R{8}, q, 0, R{20}).solve(bad)), std::invalid_argument);
}

TYPED_TEST(SpinonSurface, FiniteIllConditionedJacobianIsNotAccepted)
{
  using R = TypeParam;
  R const one{1}, pi = R{4} * math::atan(one), q = pi / R{2};
  R const eps = uni20::numeric_limits<R>::epsilon(one);
  // Well outside the ordinary boundary guard, but conditioning amplifies
  // arithmetic error past sqrt(epsilon). A finite root is not sufficient.
  R const deficit = R{4096} * eps * pi;
  EXPECT_THROW((Surface<R>(2, R{8}, q, 0, pi * math::sin(q / R{2}) - deficit).solve({})), std::runtime_error);
}

#if UNI20_HAS_MPREAL
TEST(SpinonSurfaceMpfr, ExplicitPrecisionAndDeepSoftGeometry)
{
  using R = uni20::mpreal;
  for (int bits : {128, 256, 1800, 2048})
  {
    SCOPED_TRACE(bits);
    auto const p = uni20::Precision::bits(bits);
    R const one(1, p), pi = R{4} * math::atan(one);
    for (auto const& row : spinon_surface_reference::points)
      check_reference(row, one);
    EXPECT_THROW((Surface<R>(2, R(8, p), R(pi / R{2}, uni20::Precision::bits(bits + 32)), 0, R(2, p))),
                 std::invalid_argument);
    EXPECT_TRUE((Surface<R>(2, R(1040, p), pi / R{2}, 0, R(2, p)).solve({}).inside));
  }
  EXPECT_THROW((Surface<R>(2, R{8}, R{1}, 0, R{2})), std::invalid_argument);
}

TEST(SpinonSurfaceMpfr, ConditioningRefinesLocallyAndRetainsOutputPrecision)
{
  using R = uni20::mpreal;
  auto const p = uni20::Precision::bits(64), reference_precision = uni20::Precision::bits(1024);
  R const one(1, p), pi = R{4} * math::atan(one), q = pi / R{2};
  R const w = pi * math::sin(q / R{2}) - R("0.0000000000009094947017729282379150390625", p); // 2^-40
  auto const point = Surface<R>(2, R(1040, p), q, 0, w).solve({});
  ASSERT_TRUE(point.inside);
  EXPECT_GT(point.refinements, 0u);
  EXPECT_LE(point.refinements, 2u);
  EXPECT_GT(point.working_bits, p.bit_count());
  EXPECT_EQ(point.log_coarea_jacobian.precision(), p);
  for (auto const& branch : point.beta)
    for (auto const& b : branch)
      EXPECT_EQ(b.precision(), p);

  // Independently evaluate the two-spinon Jacobian for these EXACT represented
  // inputs; promoting cannot recreate an unrounded momentum or frequency.
  R const high_one(1, reference_precision), high_pi = R{4} * math::atan(high_one);
  R const amplitude = high_pi * math::sin(R(q, reference_precision) / R{2});
  R const energy(w, reference_precision);
  R const reference = -math::log(amplitude * amplitude - energy * energy) / R{2};
  EXPECT_REAL_NEAR(R(point.log_coarea_jacobian, reference_precision), reference, R("1e-16", reference_precision));
}

TEST(SpinonSurfaceMpfr, UnresolvedBoundaryExhaustsBoundedRefinement)
{
  using R = uni20::mpreal;
  auto const p = uni20::Precision::bits(64);
  // A boundary ambiguity cannot become a successful zero or unbounded retry.
  try
  {
    (void)Surface<R>(2, R(1040, p), R(0, p), 0, R(2, p)).solve({});
    FAIL() << "unresolved boundary accepted";
  }
  catch (std::runtime_error const& error)
  {
    EXPECT_NE(std::string(error.what()).find("refinement exhausted"), std::string::npos);
  }
}

TEST(SpinonSurfaceMpfr, RefinementCannotHideAnUnrepresentablePairGap)
{
  using R = uni20::mpreal;
  auto const p = uni20::Precision::bits(16);
  std::array<R, 4> const free{R(-500, p), R(-510, p), R(-520, p), R(-530, p)};
  // Exact dyadic inputs: the pair energy is 1.79e-9 below its maximum.
  // Refinement resolves the roots, but their 16-bit rounding changes the
  // velocity difference by more than the requested arithmetic screen.
  try
  {
    (void)Surface<R>(6, R(1040, p), R("0.49700927734375", p), 0, R("0.7726898193359375", p)).solve(free);
    FAIL() << "inaccurate returned pair gap accepted";
  }
  catch (std::runtime_error const& error)
  {
    EXPECT_NE(std::string(error.what()).find("pair gap unresolved at output precision"), std::string::npos);
  }
}

TEST(SpinonSurfaceMpfr, SavedDeepPilotGeometryNeedsNoHiddenRefinement)
{
  using R = uni20::mpreal;
  using Proposal = bethe::heisenberg::detail::SpinonJointProposal<R>;
  using Kind = bethe::heisenberg::detail::SpinonProposalKind;
  auto make_proposal = [](R const& one) {
    return Proposal(1040 * one, one / 2,
                    {{Kind::background, one / 2, one, one},
                     {Kind::one_pair, one / 4, 4 * one, one / 2},
                     {Kind::two_pairs, one / 4, 8 * one, one}});
  };
  char const* qtext = "1.5707963267948966192313216916397514420985846996875529104874722961539082031431045";
  auto const high = uni20::Precision::bits(2048);
  R const high_one(1, high);
  auto const reference_proposal = make_proposal(high_one);
  for (int bits : {128, 256})
  {
    SCOPED_TRACE(bits);
    auto const low = uni20::Precision::bits(bits);
    R const one(1, low);
    auto const proposal = make_proposal(one);
    unsigned roots = 0;
    R largest_beta(0, high), largest_weight(0, high);
    for (unsigned replicate = 0; replicate < 4; ++replicate)
    {
      bethe::MonteCarloStream<5> stream(7303, replicate, 52);
      for (unsigned sample = 0; sample < 16; ++sample)
      {
        auto const identity = stream.next();
        auto const free = proposal.draw(identity.at_precision(one)).beta;
        auto const high_free = reference_proposal.draw(identity.at_precision(high_one)).beta;
        for (int section = 0; section < 3; ++section)
        {
          auto const point = Surface<R>(6, R(1040, low), R(qtext, low), section, R(2, low)).solve(free);
          auto const reference = Surface<R>(6, R(1040, high), R(qtext, high), section, R(2, high)).solve(high_free);
          ASSERT_EQ(point.inside, reference.inside);
          EXPECT_EQ(point.refinements, 0u);
          EXPECT_EQ(point.working_bits, bits);
          if (!point.inside) continue;
          roots += 2;
          std::array<R, 6> beta, high_beta;
          std::copy(point.beta[0].begin(), point.beta[0].end(), beta.begin());
          std::copy(reference.beta[0].begin(), reference.beta[0].end(), high_beta.begin());
          for (int i = 0; i < 6; ++i)
            largest_beta = std::max(largest_beta, math::abs(R(beta[i], high) - high_beta[i]));
          R const weight(proposal.balanced_log_jacobian(beta), high);
          largest_weight =
              std::max(largest_weight, math::abs(weight - reference_proposal.balanced_log_jacobian(high_beta)));
        }
      }
    }
    EXPECT_EQ(roots, 20u);
    R const tolerance(bits == 128 ? "1e-25" : "1e-60", high);
    EXPECT_LT(largest_beta, tolerance);
    EXPECT_LT(largest_weight, tolerance);
    RecordProperty("maximum_beta_error_bits_" + std::to_string(bits),
                   uni20::format_scalar(largest_beta, {.precision = 8}));
    RecordProperty("maximum_log_weight_error_bits_" + std::to_string(bits),
                   uni20::format_scalar(largest_weight, {.precision = 8}));
  }
}
#endif
} // namespace
