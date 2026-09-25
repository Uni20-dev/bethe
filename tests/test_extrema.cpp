// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/detail/extrema.hpp>

namespace
{
namespace numerical = bethe::detail;
template <typename Real> class Extrema : public ::testing::Test {};
TYPED_TEST_SUITE(Extrema, test_support::RealTypes, test_support::PrecisionNames);
template <typename R, typename F> auto objective(F f)
{
  return [f](R x) -> std::optional<numerical::ObjectiveSample<R>> { return numerical::ObjectiveSample<R>{f(x), R{0}}; };
}
TYPED_TEST(Extrema, InteriorAndEndpointEdges)
{
  using R = TypeParam;
  numerical::ExtremaOptions<R> options;
  R const center = R{1} / R{3};
  auto const result =
      numerical::bounded_extrema(objective<R>([&](R x) { return (x - center) * (x - center); }), R{-1}, R{1}, options);
  ASSERT_TRUE(result.converged) << int(result.status) << " evaluations " << result.evaluations;
  EXPECT_REAL_NEAR(result.minimum->value, R{0}, options.tolerance);
  EXPECT_REAL_NEAR(result.minimum->position, center, R{4} * options.position_tolerance);
  EXPECT_REAL_NEAR(result.maximum->value, R{16} / R{9}, options.tolerance);
  EXPECT_EQ(result.maximum->position, R{-1});
  EXPECT_EQ(result.meshes, 3u);
  EXPECT_LE(result.evaluations, options.max_evaluations);
}
TYPED_TEST(Extrema, MultipleExtremaAndPeriodicSeam)
{
  using R = TypeParam;
  R const pi = R{4} * std::atan(R{1});
  numerical::ExtremaOptions<R> options;
  // Four separated minima and maxima, none generally aligned with the mesh.
  R const shift = R{1} / R{7};
  auto const result =
      numerical::bounded_extrema(objective<R>([&](R x) { return std::cos(R{4} * (x - shift)); }), -pi, pi, options);
  ASSERT_TRUE(result.converged) << int(result.status);
  EXPECT_REAL_NEAR(result.minimum->value, R{-1}, options.tolerance);
  EXPECT_REAL_NEAR(result.maximum->value, R{1}, options.tolerance);
  // A minimum close to the domain seam, within the first mesh cell.
  auto const seam =
      numerical::bounded_extrema(objective<R>([&](R x) { return std::cos(x - R{0.001}); }), -pi, pi, options);
  ASSERT_TRUE(seam.converged) << int(seam.status);
  EXPECT_REAL_NEAR(seam.minimum->value, R{-1}, options.tolerance);
  EXPECT_REAL_NEAR(seam.maximum->value, R{1}, options.tolerance);
}
TYPED_TEST(Extrema, ConstantAndCompetingWells)
{
  using R = TypeParam;
  numerical::ExtremaOptions<R> options;
  R const value = std::sqrt(R{2}); // Detect narrowing of objective values in fp128.
  auto const constant = numerical::bounded_extrema(objective<R>([&](R) { return value; }), R{0}, R{1}, options);
  ASSERT_TRUE(constant.converged);
  EXPECT_EQ(constant.minimum->value, value);
  EXPECT_EQ(constant.maximum->value, value);
  for (R tilt : {R{-0.1}, R{0.1}})
  {
    auto f = [&](R x) { return std::min((x - R{0.4}) * (x - R{0.4}) + tilt, (x + R{0.3}) * (x + R{0.3}) - tilt); };
    auto const result = numerical::bounded_extrema(objective<R>(f), R{-1}, R{1}, options);
    ASSERT_TRUE(result.converged) << int(result.status);
    EXPECT_REAL_NEAR(result.minimum->value, R{-0.1}, options.tolerance);
    EXPECT_REAL_NEAR(result.minimum->position, tilt < R{0} ? R{0.4} : R{-0.3}, R{4} * options.position_tolerance);
  }
}
TYPED_TEST(Extrema, FailuresNeverPublishCandidates)
{
  using R = TypeParam;
  auto f = objective<R>([](R x) { return x * x; });
  auto absent = [](auto const& result, numerical::ExtremaStatus status) {
    EXPECT_FALSE(result.converged);
    EXPECT_FALSE(result.minimum);
    EXPECT_FALSE(result.maximum);
    EXPECT_EQ(result.status, status);
  };
  numerical::ExtremaOptions<R> options;
  options.max_evaluations = 0;
  absent(numerical::bounded_extrema(f, R{-1}, R{1}, options), numerical::ExtremaStatus::evaluation_limit);
  options.max_evaluations = 20;
  auto const exhausted = numerical::bounded_extrema(f, R{-1}, R{1}, options);
  absent(exhausted, numerical::ExtremaStatus::evaluation_limit);
  EXPECT_EQ(exhausted.evaluations, 20u);
  options = {};
  options.max_iterations = 0;
  absent(numerical::bounded_extrema(f, R{-1}, R{1}, options), numerical::ExtremaStatus::iteration_limit);
  options = {};
  options.max_intervals = options.initial_intervals;
  absent(numerical::bounded_extrema(f, R{-1}, R{1}, options), numerical::ExtremaStatus::mesh_limit);
  auto failed = [](R) -> std::optional<numerical::ObjectiveSample<R>> { return {}; };
  absent(numerical::bounded_extrema(failed, R{-1}, R{1}), numerical::ExtremaStatus::objective_failure);
  auto nonfinite = objective<R>([](R) { return uni20::numeric_limits<R>::infinity(); });
  absent(numerical::bounded_extrema(nonfinite, R{-1}, R{1}), numerical::ExtremaStatus::objective_failure);
  auto uncertain = [](R x) -> std::optional<numerical::ObjectiveSample<R>> {
    return numerical::ObjectiveSample<R>{x * x, R{1}};
  };
  absent(numerical::bounded_extrema(uncertain, R{-1}, R{1}), numerical::ExtremaStatus::precision_limit);
  EXPECT_THROW((void)numerical::bounded_extrema(f, R{1}, R{-1}), std::invalid_argument);
  options = {};
  options.tolerance = R{0};
  EXPECT_THROW((void)numerical::bounded_extrema(f, R{-1}, R{1}, options), std::invalid_argument);
  options = {};
  options.max_intervals = 4097;
  EXPECT_THROW((void)numerical::bounded_extrema(f, R{-1}, R{1}, options), std::invalid_argument);
}
} // namespace
