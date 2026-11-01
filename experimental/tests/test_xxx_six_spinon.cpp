// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "xxx_six_spinon_fixture.hpp"
#include <bethe/experimental/xxx_six_spinon.hpp>
#include <numeric>

namespace six = bethe::experimental::xxx::six_spinon;
namespace roots = bethe::experimental::xxx::complex_roots;
namespace xxx = bethe::heisenberg;
namespace experimental_xxx = bethe::experimental::xxx;

namespace
{
six::Configuration request(six_fixture::Entry const& entry)
{
  six::Configuration result{entry.sites, entry.triple ? six::Family::triple : six::Family::two_pairs, {}, {}};
  for (int twice : entry.holes_twice)
    result.holes.push_back(uni20::from_twice(twice));
  for (int i = 0; i < (entry.triple ? 1 : 2); ++i)
    result.branches.push_back(uni20::from_twice(entry.branches_twice[i]));
  return result;
}

template <uni20::Real R> R parse(char const* text, R const& zero)
{
  if constexpr (uni20::has_runtime_precision_v<R>)
    return uni20::parse_real<R>(text, uni20::precision_of(zero));
  else
    return uni20::parse_real<R>(text);
}

template <uni20::Real R> bethe::SolverOptions<R> options(R const& zero)
{
  if constexpr (uni20::has_runtime_precision_v<R>)
    return bethe::SolverOptions<R>(uni20::precision_of(zero));
  else
    return {};
}

template <uni20::Real R> void selected_references(R const& zero)
{
  R const one = uni20::scalar_like(zero, 1), eps = uni20::numeric_limits<R>::epsilon(zero);
  for (std::size_t index = 0; index < six_fixture::entries.size(); ++index)
  {
    auto const& expected = six_fixture::entries[index];
    SCOPED_TRACE(index);
    auto const c = request(expected);
    auto const answer = six::solve<R>(c, zero);
    ASSERT_TRUE(answer.converged()) << answer.diagnostic;
    ASSERT_TRUE(answer.state);
    auto const& state = *answer.state;
    EXPECT_EQ(answer.request.holes, c.holes);
    EXPECT_EQ(state.momentum_index, six::labels(c).momentum_index);
    EXPECT_EQ(state.nonself_conjugate, index >= 6);
    EXPECT_REAL_NEAR(state.energy, parse(expected.energy, zero), 1000000 * eps);
    ASSERT_EQ(state.roots.size(), std::size_t(c.sites / 2 - 1));
    for (std::size_t i = 0; i < state.roots.size(); ++i)
    {
      EXPECT_EQ(uni20::precision_of(state.roots[i]), uni20::precision_of(zero));
      EXPECT_REAL_NEAR(state.roots[i].real(), parse(expected.roots[i][0], zero), 1000000 * eps);
      EXPECT_REAL_NEAR(state.roots[i].imag(), parse(expected.roots[i][1], zero), 1000000 * eps);
    }
    for (auto const& value : state.variables)
      EXPECT_EQ(uni20::precision_of(value), uni20::precision_of(zero));
    EXPECT_LT(state.residual, 4 * roots::residual_tolerance<R>(c.sites, state.roots, zero));
    EXPECT_LT(state.label_error, 4 * roots::residual_tolerance<R>(c.sites, state.roots, zero));
    auto const ground_state = xxx::sector_ground_state<R>(c.sites, uni20::half_int{0}, options(zero));
    ASSERT_TRUE(ground_state.converged);
    experimental_xxx::PreparedSinglet<R> const ground(c.sites, ground_state);
    ASSERT_TRUE(ground.ready());
    auto const factor = six::raising_weight(ground, answer);
    ASSERT_TRUE(factor.converged());
    ASSERT_TRUE(factor.weight);
    EXPECT_EQ(uni20::precision_of(*factor.weight), uni20::precision_of(zero));
    EXPECT_REAL_NEAR(*factor.weight / parse(expected.weight, zero), one, 1000000000 * eps);
    auto changed = answer;
    changed.state->roots.front() += one / 1000;
    EXPECT_FALSE(six::raising_weight(ground, changed).weight);
  }
}

template <uni20::Real R> void jacobians(R const& zero)
{
  R const one = uni20::scalar_like(zero, 1), h = uni20::math::cbrt(uni20::numeric_limits<R>::epsilon(zero));
  for (bool triple : {false, true})
  {
    auto const c = request(six_fixture::entries[triple ? 1 : 0]);
    auto const labels = six::labels(c);
    std::vector<R> centres{one / 10, -3 * one / 10, 4 * one / 5};
    auto const t = six::detail::takahashi<R>(12, labels, centres, zero, true);
    for (std::size_t j = 0; j < centres.size(); ++j)
    {
      auto plus = centres, minus = centres;
      plus[j] += h;
      minus[j] -= h;
      auto const fp = six::detail::takahashi<R>(12, labels, plus, zero);
      auto const fm = six::detail::takahashi<R>(12, labels, minus, zero);
      for (std::size_t i = 0; i < centres.size(); ++i)
        EXPECT_REAL_NEAR(t[i * centres.size() + j], (fp[i] - fm[i]) / (2 * h), 10000 * h * h);
    }
    std::vector<R> x = triple ? std::vector<R>{one / 10, one / 5, 3 * one / 5, 61 * one / 100, 49 * one / 50}
                              : std::vector<R>{3 * one / 10, -one / 2, 12 * one / 25, 4 * one / 5, 51 * one / 100};
    std::size_t const reals = triple ? 3 : 1;
    auto const g = six::detail::jacobian<R>(12, reals, x, zero);
    for (std::size_t j = 0; j < x.size(); ++j)
    {
      auto plus = x, minus = x;
      plus[j] += h;
      minus[j] -= h;
      auto const fp = six::detail::equations<R>(12, reals, plus, zero);
      auto const fm = six::detail::equations<R>(12, reals, minus, zero);
      for (std::size_t i = 0; i < x.size(); ++i)
        EXPECT_REAL_NEAR(g[i * x.size() + j], (fp[i] - fm[i]) / (2 * h), 10000000 * h * h);
    }
  }
}
} // namespace

TEST(ExperimentalSixLabels, RequestsAndRejections)
{
  for (auto const& entry : six_fixture::entries)
  {
    auto c = request(entry);
    auto const l = six::labels(c);
    EXPECT_EQ(std::accumulate(l.lengths.begin(), l.lengths.end(), 0), c.sites / 2 - 1);
    c.holes[1] = c.holes[0];
    EXPECT_THROW((void)six::labels(c), std::invalid_argument);
  }
  auto c = request(six_fixture::entries[0]);
  c.branches = {0, 1};
  EXPECT_THROW((void)six::solve<double>(c), std::invalid_argument);
  c = request(six_fixture::entries[0]);
  c.sites = 14;
  EXPECT_THROW((void)six::solve<double>(c), std::invalid_argument);
}

template <typename T> class ExperimentalSix : public ::testing::Test {};
TYPED_TEST_SUITE(ExperimentalSix, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(ExperimentalSix, IndependentSelectedRootsEnergiesAndWeights) { selected_references(TypeParam{0}); }
TYPED_TEST(ExperimentalSix, AnalyticJacobians) { jacobians(TypeParam{0}); }
TYPED_TEST(ExperimentalSix, SameMomentumDifferentLabelsReject)
{
  using R = TypeParam;
  six::Configuration const c{
      12, six::Family::two_pairs, {-3, -2, -1, 0, 1, 2}, {uni20::from_twice(-3), uni20::from_twice(3)}};
  auto changed = six::solve<R>(c);
  ASSERT_TRUE(changed.converged()) << changed.diagnostic;
  auto const ground_state = xxx::sector_ground_state<R>(12, uni20::half_int{0});
  experimental_xxx::PreparedSinglet<R> const ground(12, ground_state);
  ASSERT_TRUE(six::raising_weight(ground, changed).weight);
  changed.request.branches = {uni20::from_twice(-1), uni20::from_twice(1)};
  EXPECT_EQ(six::labels(changed.request).momentum_index, changed.state->momentum_index);
  EXPECT_FALSE(six::raising_weight(ground, changed).weight);
}
TYPED_TEST(ExperimentalSix, UnresolvedNativeDeviationsRetainFailure)
{
  using R = TypeParam;
  six::Configuration const c{
      20, six::Family::two_pairs, {-4, -3, -1, 0, 1, 4}, {uni20::from_twice(1), uni20::from_twice(3)}};
  auto const answer = six::solve<R>(c);
  EXPECT_FALSE(answer.converged());
  EXPECT_FALSE(answer.state);
  EXPECT_EQ(answer.status, six::SolveStatus::unresolved_geometry);
  EXPECT_EQ(answer.request.holes, c.holes);
  EXPECT_FALSE(answer.diagnostic.empty());
  auto const ground_state = xxx::sector_ground_state<R>(20, uni20::half_int{0});
  experimental_xxx::PreparedSinglet<R> const ground(20, ground_state);
  EXPECT_FALSE(six::raising_weight(ground, answer).weight);
}

#if UNI20_HAS_MPREAL
TEST(ExperimentalSixMpfr, ExplicitPrecision)
{
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {128, 256})
  {
    uni20::mpreal const zero(0, uni20::Precision::bits(bits));
    selected_references(zero);
    jacobians(zero);
    auto answer = six::solve(request(six_fixture::entries[0]), zero);
    ASSERT_TRUE(answer.converged());
    auto const ground_state = xxx::sector_ground_state<uni20::mpreal>(12, uni20::half_int{0}, options(zero));
    experimental_xxx::PreparedSinglet<uni20::mpreal> const ground(12, ground_state);
    answer.state->roots.front() = answer.state->roots.front().at(uni20::Precision::bits(bits + 32));
    EXPECT_THROW((void)six::raising_weight(ground, answer), std::invalid_argument);
    answer.state->roots.front() = uni20::complex<uni20::mpreal>(uni20::uninitialized);
    EXPECT_THROW((void)six::raising_weight(ground, answer), std::invalid_argument);
  }
  auto const c = request(six_fixture::entries[0]);
  EXPECT_THROW((void)six::solve(c, uni20::mpreal{}), std::invalid_argument);
  EXPECT_THROW((void)six::solve(c, uni20::mpreal(uni20::uninitialized)), std::invalid_argument);
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
#endif
