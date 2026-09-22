// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/xxz_phantom_check.hpp>

#include "test_support.hpp"
#include <array>
#include <bethe/xxz_odd_continuation.hpp>

namespace
{
namespace engine = bethe::xxz::detail;
using Status = engine::PhantomLiftStatus;
template <typename Real> class XXZPhantomCheck : public ::testing::Test {};
TYPED_TEST_SUITE(XXZPhantomCheck, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXZPhantomCheck, ContinuedMixedStatesHaveResolvedWitnesses)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  for (auto sizes :
       {std::array<unsigned, 2>{7, 1}, {9, 1}, {9, 2}, {11, 1}, {11, 2}, {11, 3}, {13, 1}, {13, 2}, {13, 3}, {13, 4}})
  {
    auto const [n, p] = sizes;
    auto const m = n / 2, r = m - p;
    SCOPED_TRACE(::testing::Message() << "N=" << n << " p=" << p);
    Real const delta = p == 1 ? -Real{1} / Real{2} : -std::cos(pi / Real(2 * p + 1));
    auto const branch = engine::continue_odd_polynomial(n, delta, uni20::from_twice(std::int64_t{1}));
    ASSERT_TRUE(branch.equations_converged);
    engine::PolynomialBetheSystem<Real> const system(n, m, branch.center, branch.coordinate_scale);
    auto const result = engine::check_phantom_lift<Real>(system, branch.coefficients, delta, p, -1);
    ASSERT_EQ(result.status, Status::nonzero_witness)
        << "reduction=" << int(result.reduction.status) << " ratio=" << uni20::format_scalar(result.resolution_ratio);
    ASSERT_TRUE(result.recovered.has_value());
    EXPECT_GT(result.resolution_ratio, Real{1});
    EXPECT_GT(std::abs(result.amplitude), result.input_variation + result.roundoff_allowance);
    EXPECT_GT(result.configurations_tested, 0U);
    EXPECT_LE(result.configurations_tested, 16U);
    auto const terms = engine::phantom_dressing_terms(m, p, 1000000);
    EXPECT_EQ(result.subset_updates, result.configurations_tested * terms * (std::size_t{1} << r));
    std::vector<C> v;
    for (C x : result.recovered->roots)
    {
      C const z = branch.center + branch.coordinate_scale * x;
      v.push_back(-(Real{1} - C{0, 1} * z) / (Real{1} + C{0, 1} * z));
    }
    engine::CoordinateBetheWave<Real> const wave(n, delta, v);
    C const phase{delta, -std::sqrt(Real{1} + delta) * std::sqrt(Real{1} - delta)};
    engine::PhantomDressing<Real> const dressing(n, r, p, phase);
    C const reference =
        dressing.amplitude(result.occupied, [&](auto selected) { return wave.evaluate(selected).value; });
    EXPECT_LT(std::abs(reference - result.amplitude), Real{128} * eps * result.absolute_term_sum);
  }
}

TYPED_TEST(XXZPhantomCheck, BothChiralitiesAndNoInputMutation)
{
  using Real = TypeParam;
  Real const delta = -Real{1} / Real{2};
  auto const branch = engine::continue_odd_polynomial(9, delta, uni20::from_twice(std::int64_t{1}));
  ASSERT_TRUE(branch.equations_converged);
  auto reflected = branch.coefficients;
  for (std::size_t j = 0; j < reflected.size(); ++j)
    if ((reflected.size() - j) % 2) reflected[j] = -reflected[j];
  auto const saved = reflected;
  engine::PolynomialBetheSystem<Real> const system(9, 4, -branch.center, branch.coordinate_scale);
  auto const result = engine::check_phantom_lift<Real>(system, reflected, delta, 1, 1);
  EXPECT_EQ(result.status, Status::nonzero_witness);
  EXPECT_TRUE(reflected == saved);
  EXPECT_EQ(system.center, -branch.center);
}

TYPED_TEST(XXZPhantomCheck, BudgetsAndInconclusiveAmplitudes)
{
  using Real = TypeParam;
  Real const delta = -Real{1} / Real{2};
  auto const branch = engine::continue_odd_polynomial(7, delta, uni20::from_twice(std::int64_t{1}));
  ASSERT_TRUE(branch.equations_converged);
  engine::PolynomialBetheSystem<Real> const system(7, 3, branch.center, branch.coordinate_scale);
  engine::PhantomLiftOptions<Real> options;
  options.max_subset_updates = 11; // Each witness costs C(3,1)*2^2=12.
  auto result = engine::check_phantom_lift<Real>(system, branch.coefficients, delta, 1, -1, options);
  EXPECT_EQ(result.status, Status::work_limit);
  EXPECT_EQ(result.configurations_tested, 0U);
  EXPECT_FALSE(result.recovered.has_value());
  options.max_subset_updates = 12;
  result = engine::check_phantom_lift<Real>(system, branch.coefficients, delta, 1, -1, options);
  EXPECT_EQ(result.status, Status::nonzero_witness);
  EXPECT_EQ(result.configurations_tested, 1U);
  EXPECT_EQ(result.subset_updates, 12U);
  options.amplitude_tolerance = Real{1}; // Explicitly demand unattainable resolution.
  result = engine::check_phantom_lift<Real>(system, branch.coefficients, delta, 1, -1, options);
  EXPECT_EQ(result.status, Status::work_limit);
  EXPECT_EQ(result.configurations_tested, 1U);
  EXPECT_LE(result.resolution_ratio, Real{1});
  options.max_subset_updates = 1000000;
  options.max_configurations = 100;
  result = engine::check_phantom_lift<Real>(system, branch.coefficients, delta, 1, -1, options);
  EXPECT_EQ(result.status, Status::resolution_unresolved);
  EXPECT_TRUE(result.all_configurations_tested);
  EXPECT_EQ(result.configurations_tested, 35U); // C(7,3), not a claim the vector is zero.
  options = {};
  options.root_options.max_iterations = 0;
  result = engine::check_phantom_lift<Real>(system, branch.coefficients, delta, 1, -1, options);
  EXPECT_EQ(result.status, Status::roots_unresolved);
  EXPECT_EQ(result.configurations_tested, 0U);
  options = {};
  options.max_configurations = 0;
  EXPECT_EQ(engine::check_phantom_lift<Real>(system, branch.coefficients, delta, 1, -1, options).status,
            Status::work_limit);
}

TYPED_TEST(XXZPhantomCheck, RejectsInvalidAndOffShellCandidates)
{
  using Real = TypeParam;
  Real const delta = -Real{1} / Real{2};
  auto branch = engine::continue_odd_polynomial(7, delta, uni20::from_twice(std::int64_t{1}));
  ASSERT_TRUE(branch.equations_converged);
  engine::PolynomialBetheSystem<Real> const system(7, 3, branch.center, branch.coordinate_scale);
  auto wrong = engine::check_phantom_lift<Real>(system, branch.coefficients, delta, 1, 1);
  EXPECT_EQ(wrong.status, Status::reduction_unresolved);
  EXPECT_EQ(wrong.configurations_tested, 0U);
  branch.coefficients[0] += Real{1} / Real{100};
  auto off_shell = engine::check_phantom_lift<Real>(system, branch.coefficients, delta, 1, -1);
  EXPECT_EQ(off_shell.status, Status::reduction_unresolved);
  EXPECT_FALSE(off_shell.recovered.has_value());
  engine::PhantomLiftOptions<Real> options;
  options.amplitude_tolerance = Real{0};
  EXPECT_THROW(engine::check_phantom_lift<Real>(system, branch.coefficients, delta, 1, -1, options),
               std::invalid_argument);
  options = {};
  options.root_options.tolerance = uni20::numeric_limits<Real>::quiet_NaN();
  EXPECT_THROW(engine::check_phantom_lift<Real>(system, branch.coefficients, delta, 1, -1, options),
               std::invalid_argument);
  EXPECT_THROW(engine::check_phantom_lift<Real>(system, branch.coefficients, delta, 0, -1), std::invalid_argument);
}
} // namespace
