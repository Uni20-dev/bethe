// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/detail/eigenvalue_continuation.hpp>
#include <bethe/detail/newton_backtracking.hpp>

namespace
{
template <typename Real> class NewtonControl : public ::testing::Test {};
TYPED_TEST_SUITE(NewtonControl, test_support::RealTypes, test_support::PrecisionNames);

enum class ControlStatus
{
  iteration_limit,
  stalled,
  ill_conditioned
};
struct ControlWork
{
    std::size_t iterations = 0, stages = 0, rejected_stages = 0;
    ControlStatus status = ControlStatus::iteration_limit;
};
template <typename Real> struct ControlOptions
{
    Real residual_tolerance{};
    std::size_t max_iterations = 20, max_stages = 20;
};
template <typename Real> struct Evaluation
{
    std::vector<Real> residual;
    Real norm{};
};

TYPED_TEST(NewtonControl, ContinuationRejectsWithoutCommittingAndClampsTarget)
{
  using Real = TypeParam;
  std::vector<Real> x{Real{0}};
  Real reached{};
  ControlOptions<Real> options;
  ControlWork work;
  unsigned commits = 0;
  std::vector<Real> attempted;
  bethe::detail::continue_eigenvalues(
      x, reached, Real{0.5}, Real{1}, options, work,
      [&](auto& slope) -> std::optional<bethe::detail::ContinuationPrediction<Real>> {
        test_support::expect_exact(x[0], reached, "only committed iterates reach predictor");
        slope = {Real{1}};
        return bethe::detail::ContinuationPrediction<Real>{Real{1}, Real{1}};
      },
      [&](auto& candidate, Real next, Real, auto const&, Real) {
        attempted.push_back(next);
        ++work.iterations;
        if (attempted.size() == 1)
        {
          candidate[0] = Real{100};
          return false;
        }
        return true;
      },
      [&](Real next) {
        ++commits;
        test_support::expect_exact(next, reached, "commit after acceptance");
      });
  ASSERT_EQ(attempted.size(), 4u);
  test_support::expect_exact(attempted[0], Real{0.5}, "initial stage");
  test_support::expect_exact(attempted[1], Real{0.25}, "halved retry");
  test_support::expect_exact(attempted[2], Real{0.625}, "grown step");
  test_support::expect_exact(reached, Real{1}, "exact endpoint");
  test_support::expect_exact(x[0], Real{1}, "final iterate");
  EXPECT_EQ(work.stages, 4u);
  EXPECT_EQ(work.iterations, 4u);
  EXPECT_EQ(work.rejected_stages, 1u);
  EXPECT_EQ(commits, 3u);
}

TYPED_TEST(NewtonControl, ContinuationBudgetAndTangentFailure)
{
  using Real = TypeParam;
  for (unsigned mode = 0; mode < 3; ++mode)
  {
    std::vector<Real> x{Real{0}};
    Real reached{};
    ControlOptions<Real> options;
    if (mode == 0) options.max_iterations = 0;
    if (mode == 1) options.max_stages = 0;
    ControlWork work;
    unsigned predictions = 0;
    bethe::detail::continue_eigenvalues(
        x, reached, Real{0.5}, Real{1}, options, work,
        [&](auto&) -> std::optional<bethe::detail::ContinuationPrediction<Real>> {
          ++predictions;
          return std::nullopt;
        },
        [](auto&, Real, Real, auto const&, Real) {
          ADD_FAILURE();
          return true;
        });
    EXPECT_EQ(predictions, mode == 2 ? 1u : 0u);
    EXPECT_EQ(work.stages, 0u);
    EXPECT_EQ(work.iterations, 0u);
    EXPECT_EQ(work.status, mode == 2 ? ControlStatus::ill_conditioned : ControlStatus::iteration_limit);
    test_support::expect_exact(reached, Real{0}, "unreached target");
  }
}

TYPED_TEST(NewtonControl, ContinuationTrustLimitAndUnrepresentableStep)
{
  using Real = TypeParam;
  std::vector<Real> x{Real{0}};
  Real reached{};
  ControlOptions<Real> options;
  options.max_stages = 1;
  ControlWork work;
  auto predict = [](auto& slope) -> std::optional<bethe::detail::ContinuationPrediction<Real>> {
    slope = {Real{2}};
    return bethe::detail::ContinuationPrediction<Real>{Real{0.25}, Real{2}};
  };
  bethe::detail::continue_eigenvalues(x, reached, Real{1}, Real{2}, options, work, predict,
                                      [](auto&, Real, Real, auto const&, Real) { return true; });
  test_support::expect_exact(reached, Real{0.125}, "trust limited parameter step");
  test_support::expect_exact(x[0], Real{0.25}, "trust limited predictor");
  reached = Real{1};
  work = {};
  bethe::detail::continue_eigenvalues(x, reached, uni20::numeric_limits<Real>::epsilon() / Real{4}, Real{2}, options,
                                      work, predict, [](auto&, Real, Real, auto const&, Real) {
                                        ADD_FAILURE();
                                        return true;
                                      });
  EXPECT_EQ(work.status, ControlStatus::stalled);
  EXPECT_EQ(work.stages, 0u);
}

TYPED_TEST(NewtonControl, CorrectorCountsFailuresAndChecksFinalBudgetedUpdate)
{
  using Real = TypeParam;
  for (bool singular : {false, true})
  {
    std::vector<Real> x{Real{0}};
    ControlOptions<Real> options;
    options.max_iterations = 1;
    std::size_t iterations = 0;
    auto evaluate = [&](auto const& trial, std::vector<Real>* jac) {
      if (jac) *jac = {singular ? Real{0} : Real{1}};
      return Evaluation<Real>{{trial[0] - Real{1}}, std::abs(trial[0] - Real{1})};
    };
    EXPECT_EQ(bethe::detail::correct_eigenvalue_stage(x, Real{2}, options, iterations, evaluate,
                                                      [](auto const& trial) { return trial; }),
              !singular);
    EXPECT_EQ(iterations, 1u);
  }
}

TYPED_TEST(NewtonControl, CorrectorTrustBoxSkipsInvalidResiduals)
{
  using Real = TypeParam;
  std::vector<Real> x{Real{0}};
  ControlOptions<Real> options;
  std::size_t iterations = 0;
  unsigned evaluations = 0;
  EXPECT_FALSE(bethe::detail::correct_eigenvalue_stage(
      x, Real{0}, options, iterations,
      [&](auto const&, std::vector<Real>* jac) {
        ++evaluations;
        if (jac) *jac = {Real{1}};
        return Evaluation<Real>{{Real{-1}}, Real{1}};
      },
      [](auto const& trial) { return trial; }));
  EXPECT_EQ(evaluations, 1u);
  EXPECT_EQ(iterations, 1u);
  test_support::expect_exact(x[0], Real{0}, "failed corrector unchanged");
}

TYPED_TEST(NewtonControl, CorrectorRetainsTwelveAttemptStageLimit)
{
  using Real = TypeParam;
  std::vector<Real> x{Real{0}};
  ControlOptions<Real> options;
  options.residual_tolerance = Real{1};
  std::size_t iterations = 0;
  EXPECT_FALSE(bethe::detail::correct_eigenvalue_stage(
      x, Real{20}, options, iterations,
      [&](auto const& trial, std::vector<Real>* jac) {
        if (jac) *jac = {Real{1}};
        return Evaluation<Real>{{Real{-1}}, Real{13} - trial[0]};
      },
      [](auto const& trial) { return trial; }));
  EXPECT_EQ(iterations, 12u);
  test_support::expect_exact(x[0], Real{12}, "no thirteenth convergence check");
}

TYPED_TEST(NewtonControl, DomainRejectionHalvesFromOriginalIterate)
{
  using Real = TypeParam;
  std::vector<Real> x{Real{1}};
  unsigned trials = 0, evaluations = 0;
  EXPECT_TRUE(bethe::detail::backtrack_newton(
      x, [](std::size_t) { return Real{-2}; },
      [&](auto const& trial, Real damping) {
        ++trials;
        if (trial[0] <= Real{0}) return false;
        ++evaluations;
        return bethe::detail::newton_decreases(trial[0], Real{1}, damping, Real{0});
      }));
  EXPECT_EQ(trials, 3u);
  EXPECT_EQ(evaluations, 1u);
  test_support::expect_exact(x[0], Real{0.5}, "accepted quarter step");
}

TYPED_TEST(NewtonControl, ExhaustionKeepsIterateAndNativeHalvingBudget)
{
  using Real = TypeParam;
  std::vector<Real> x{Real{1}};
  unsigned trials = 0;
  Real last{};
  EXPECT_FALSE(bethe::detail::backtrack_newton(
      x, [](std::size_t) { return Real{1}; },
      [&](auto const&, Real damping) {
        ++trials;
        last = damping;
        return false;
      }));
  EXPECT_EQ(trials, unsigned(uni20::numeric_limits<Real>::digits + 1));
  Real expected{1};
  for (int i = 0; i < uni20::numeric_limits<Real>::digits; ++i)
    expected /= Real{2};
  test_support::expect_exact(last, expected, "last damping");
  test_support::expect_exact(x[0], Real{1}, "rejected iterate");
}

TYPED_TEST(NewtonControl, NormalizationRunsOnFreshTrials)
{
  using Real = TypeParam;
  std::vector<Real> x{Real{1}};
  unsigned trials = 0;
  EXPECT_TRUE(bethe::detail::backtrack_newton(
      x, [](std::size_t) { return Real{2}; },
      [&](auto const& trial, Real damping) {
        ++trials;
        test_support::expect_exact(trial[0], Real{11} + Real{2} * damping, "normalized trial");
        return trials == 2;
      },
      [](auto& trial) { trial[0] += Real{10}; }));
  test_support::expect_exact(x[0], Real{12}, "normalized accepted iterate");
}

TYPED_TEST(NewtonControl, ToleranceShortcutAndNativeResolution)
{
  using Real = TypeParam;
  EXPECT_TRUE(bethe::detail::newton_decreases(Real{1}, Real{1}, Real{1}, Real{1}));
  EXPECT_FALSE(bethe::detail::newton_decreases(Real{1}, Real{1}, Real{1}, Real{0}));
  EXPECT_FALSE(bethe::detail::newton_decreases(uni20::numeric_limits<Real>::quiet_NaN(), Real{1}, Real{1}, Real{1}));
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  std::vector<Real> x{Real{1}};
  EXPECT_TRUE(bethe::detail::backtrack_newton(
      x, [&](std::size_t) { return eps; }, [](auto const& trial, Real) { return trial[0] > Real{1}; }));
  test_support::expect_exact(x[0], Real{1} + eps, "native correction");
}
} // namespace
