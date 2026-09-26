// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/sine_gordon_bethe_yang.hpp>
#include <bethe/sine_gordon_excited.hpp>

namespace
{
namespace sg = bethe::sine_gordon;
template <typename Real> class SineGordonExcited : public ::testing::Test {};
TYPED_TEST_SUITE(SineGordonExcited, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(SineGordonExcited, SourceAgainstIndependentClosedKernel)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  sg::KernelOptions<Real> options;
  options.tolerance = Real{65536} * eps;
  std::size_t evaluations = 0;
  auto const source = sg::detail::hole_table(Real{2}, Real{0.5}, Real{3}, 8, Real{0.3}, options, evaluations);
  ASSERT_TRUE(source.converged) << uni20::format_scalar(source.error);
  auto const rule = bethe::detail::gauss_legendre<Real>(96);
  auto kernel = [&](C z) { return z / (Real{2} * pi * pi * std::sinh(z)); };
  auto phase = [&](C z) {
    bethe::detail::CompensatedSum<C> sum;
    for (std::size_t i = 0; i < rule.x.size(); ++i)
      sum.add(pi * z * rule.w[i] * kernel(z * ((Real{1} + rule.x[i]) / Real{2})));
    return sum.value();
  };
  EXPECT_REAL_NEAR(source.phase, phase(C(Real{0.6})).real(), Real{131072} * eps);
  EXPECT_REAL_NEAR(source.phase_derivative, Real{2} * pi * kernel(C(Real{0.6})).real(), Real{131072} * eps);
  for (unsigned j : {0u, 4u, 8u})
  {
    Real const x = Real{0.75} * (Real(j) - Real{4});
    C const z(x, Real{0.5});
    EXPECT_REAL_NEAR(std::abs(source.source[j] + C(0, 1) * (phase(z - Real{0.3}) + phase(z + Real{0.3}))), Real{0},
                     Real{131072} * eps);
    EXPECT_REAL_NEAR(std::abs(source.counting[j] - kernel(C(Real{0.3} - x, -Real{0.5}))), Real{0}, Real{131072} * eps);
  }
}

TEST(SineGordonExcitedAudit, IndependentClosedKernelLevelsAndScaling)
{
  sg::TwoSolitonOptions<double> options;
  options.tolerance = 1e-7;
  // scripts/reference_sine_gordon_excited.py: nonuniform Gauss meshes,
  // closed p=2 kernel, 256/512 nodes, two contours and cutoffs. Rounded
  // fp64 oracle fixtures deliberately do not claim high-precision accuracy.
  for (int twice : {1, 3})
  {
    SCOPED_TRACE(twice);
    double const energy = twice == 1 ? 5.11055662981005 : 17.14609634255262;
    double const hole = twice == 1 ? 1.623489101316185 : 2.852118282211039;
    for (double mass : {1., 4.})
    {
      auto const state = sg::two_soliton_level(mass, 1 / mass, 2., uni20::from_twice(twice), options);
      ASSERT_TRUE(state.converged) << int(state.status);
      EXPECT_NEAR(*state.casimir_energy / mass, energy, 1e-7);
      EXPECT_NEAR(*state.rapidity, hole, 1e-8);
      EXPECT_NEAR(*state.scaling_function, energy, 1e-7);
      EXPECT_FALSE(state.effective_central_charge); // not a vacuum central charge
      EXPECT_LE(state.hole_error, options.tolerance / 8);
      EXPECT_LE(state.quantization_residual, options.tolerance / 8);
    }
  }
}

TEST(SineGordonExcitedAudit, UltravioletAndInfrared)
{
  sg::TwoSolitonOptions<double> options;
  options.tolerance = 1e-7;
  double const pi = 4 * std::atan(1.);
  for (int twice : {1, 3})
    for (double p : {1., 1.5, 2., 3.})
    {
      SCOPED_TRACE(twice);
      SCOPED_TRACE(p);
      double previous_error = 1.;
      for (double u : {0.1, 0.02})
      {
        auto const state = sg::two_soliton_level(1., u, p, uni20::from_twice(twice), options);
        ASSERT_TRUE(state.converged) << int(state.status);
        // FRT 1999 section 5.2.2: x=R^2+2I-1, R^2=(p+1)/(2p).
        double const error = std::abs(*state.scaling_function / (2 * pi) + 1. / 12 - ((p + 1) / (2 * p) + twice - 1));
        EXPECT_LT(error, previous_error);
        if (u == 0.02) EXPECT_LT(error, 0.008);
        previous_error = error;
      }
      auto const state = sg::two_soliton_level(1., 10., p, uni20::from_twice(twice), options);
      auto const by = sg::same_charge_pair(1., 10., p, uni20::from_twice(twice));
      ASSERT_TRUE(state.converged && by.converged);
      EXPECT_NEAR(*state.casimir_energy, *by.energy, 1e-4);
      // Finite-volume sea corrections are resolved, not an alias for BY.
      EXPECT_GT(std::abs(*state.casimir_energy - *by.energy), 1e-7);
    }
}

TYPED_TEST(SineGordonExcited, NativeFreeDiracLevel)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1});
  auto const excited = sg::two_soliton_level(Real{1}, Real{1}, Real{1});
  auto const vacuum = sg::vacuum_energy(Real{1}, Real{1}, Real{1});
  ASSERT_TRUE(excited.converged && vacuum.converged);
  EXPECT_REAL_NEAR(*excited.casimir_energy - *vacuum.casimir_energy, Real{2} * std::hypot(Real{1}, pi),
                   Real{16777216} * uni20::numeric_limits<Real>::epsilon());
}

TYPED_TEST(SineGordonExcited, WarmStartRechecksChangedSource)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  sg::VacuumOptions<Real> options;
  options.tolerance = Real{262144} * eps;
  sg::KernelOptions<Real> quadrature;
  quadrature.tolerance = Real{65536} * eps;
  std::size_t evaluations = 0;
  auto const kernel = sg::detail::kernel_table(Real{2}, Real{0.5}, Real{10} / Real{64}, 64, quadrature, evaluations);
  auto const source = sg::detail::hole_table(Real{2}, Real{0.5}, Real{5}, 64, Real{0.3}, quadrature, evaluations);
  auto const changed = sg::detail::hole_table(Real{2}, Real{0.5}, Real{5}, 64, Real{0.4}, quadrature, evaluations);
  ASSERT_TRUE(kernel.converged && source.converged && changed.converged);
  std::vector<std::complex<Real>> seed;
  sg::VacuumState<Real> first_work, cold_work, warm_work, exhausted_work;
  auto const first = sg::detail::nlie_mesh<Real>(Real{2}, Real{1}, Real{0.5}, Real{5}, 64, options, first_work, kernel,
                                                 source.source, source.counting, &seed);
  ASSERT_TRUE(first.converged);
  auto exhausted_options = options;
  exhausted_options.max_iterations = 0;
  auto const exhausted = sg::detail::nlie_mesh<Real>(Real{2}, Real{1}, Real{0.5}, Real{5}, 64, exhausted_options,
                                                     exhausted_work, kernel, changed.source, changed.counting, &seed);
  EXPECT_FALSE(exhausted.converged); // a previous converged root is not the new solution
  auto const warm = sg::detail::nlie_mesh<Real>(Real{2}, Real{1}, Real{0.5}, Real{5}, 64, options, warm_work, kernel,
                                                changed.source, changed.counting, &seed);
  auto const cold = sg::detail::nlie_mesh<Real>(Real{2}, Real{1}, Real{0.5}, Real{5}, 64, options, cold_work, kernel,
                                                changed.source, changed.counting);
  ASSERT_TRUE(warm.converged && cold.converged);
  EXPECT_REAL_NEAR(warm.value, cold.value, options.tolerance);
  EXPECT_REAL_NEAR(warm.counting_integral, cold.counting_integral, options.tolerance);
  EXPECT_GT(warm_work.iterations, 0u);
  EXPECT_LT(warm_work.iterations, cold_work.iterations);
}

TYPED_TEST(SineGordonExcited, NativeInteractingLevel)
{
  using Real = TypeParam;
  auto const state = sg::two_soliton_level(Real{1}, Real{1}, Real{2});
  ::testing::Test::RecordProperty("nonlinear_updates", std::to_string(state.iterations));
  ::testing::Test::RecordProperty("hole_trials", std::to_string(state.root_iterations));
  ::testing::Test::RecordProperty("intervals", std::to_string(state.intervals));
  ::testing::Test::RecordProperty("cutoffs", std::to_string(state.cutoffs));
  if (state.casimir_energy) ::testing::Test::RecordProperty("energy", uni20::format_scalar(*state.casimir_energy));
  ASSERT_TRUE(state.converged) << int(state.status) << " iterations " << state.iterations << " roots "
                               << state.root_iterations << " intervals " << state.intervals << " cutoffs "
                               << state.cutoffs;
  EXPECT_REAL_NEAR(*state.casimir_energy, uni20::parse_real<Real>("5.11055662981005"), Real{1e-10});
  Real const tolerance = sg::TwoSolitonOptions<Real>{}.tolerance;
  EXPECT_LE(state.hole_error, tolerance / Real{8});
  EXPECT_LE(state.quantization_residual, tolerance / Real{8});
  EXPECT_LE(state.mesh_error, tolerance / Real{4});
  EXPECT_LE(state.cutoff_error, tolerance / Real{4});
  EXPECT_LE(state.contour_error, tolerance / Real{2});
}

TEST(SineGordonExcitedAudit, CoarseHoleIsOnlyAnInitialGuess)
{
  sg::TwoSolitonOptions<double> options;
  options.tolerance = 1e-8;
  sg::TwoSolitonState<double> coarse_work, warm_work, cold_work;
  auto const coarse = sg::detail::two_soliton_mesh(2., 1., .5, 5., 16, options, coarse_work);
  ASSERT_TRUE(coarse.converged);
  auto const warm =
      sg::detail::two_soliton_mesh(2., 1., .5, 5., 64, options, warm_work, std::optional<double>(coarse.hole));
  auto const cold = sg::detail::two_soliton_mesh(2., 1., .5, 5., 64, options, cold_work);
  ASSERT_TRUE(warm.converged && cold.converged);
  EXPECT_GT(std::abs(warm.hole - coarse.hole), 1e-6);
  EXPECT_NEAR(warm.hole, cold.hole, 1e-9);
  EXPECT_NEAR(warm.value, cold.value, 1e-8);
  EXPECT_LE(warm_work.quantization_residual, options.tolerance / 8);
  EXPECT_GT(warm_work.iterations, 0u); // the sea is solved anew, not copied from the coarse grid
}

TEST(SineGordonExcitedAudit, IndependentBudgetsAndInvalidInputs)
{
  auto check = [](sg::TwoSolitonOptions<double> options, double p, sg::VacuumStatus status) {
    SCOPED_TRACE(int(status));
    auto const state = sg::two_soliton_level(1., 1., p, uni20::from_twice(1), options);
    EXPECT_EQ(state.status, status);
    EXPECT_FALSE(state.converged);
    EXPECT_FALSE(state.casimir_energy);
    EXPECT_FALSE(state.scaling_function);
    EXPECT_FALSE(state.rapidity);
  };
  sg::TwoSolitonOptions<double> options;
  options.max_root_iterations = 0;
  check(options, 2., sg::VacuumStatus::iteration_limit);
  options = {};
  options.max_root_iterations = 1;
  check(options, 1., sg::VacuumStatus::iteration_limit);
  options = {};
  options.max_iterations = 0;
  check(options, 2., sg::VacuumStatus::iteration_limit);
  options = {};
  options.max_intervals = 64;
  check(options, 1., sg::VacuumStatus::mesh_limit);
  options = {};
  options.max_cutoffs = 1;
  check(options, 1., sg::VacuumStatus::cutoff_limit);
  options = {};
  options.kernel.max_evaluations = 0;
  check(options, 2., sg::VacuumStatus::kernel_limit);
  options = {};
  options.kernel.max_levels = 0;
  check(options, 2., sg::VacuumStatus::kernel_limit);
  options = {};
  options.kernel.max_cutoffs = 0;
  check(options, 2., sg::VacuumStatus::kernel_limit);
  options = {};
  options.initial_cutoff = 1e308;
  check(options, 2., sg::VacuumStatus::precision_limit);
  options = {};
  options.initial_cutoff = .5;
  options.max_cutoffs = 2;
  options.tolerance = 1e-4;
  check(options, 1., sg::VacuumStatus::cutoff_limit);
  EXPECT_EQ(sg::two_soliton_level(1e308, 2., 1.).status, sg::VacuumStatus::precision_limit);
  EXPECT_THROW(sg::two_soliton_level(0., 1., 2.), std::invalid_argument);
  EXPECT_THROW(sg::two_soliton_level(1., 0., 2.), std::invalid_argument);
  EXPECT_THROW(sg::two_soliton_level(1., 1., .5), std::invalid_argument);
  EXPECT_THROW(sg::two_soliton_level(1., 1., 2., uni20::from_twice(2)), std::invalid_argument);
  EXPECT_THROW(sg::two_soliton_level(1., 1., 2., uni20::from_twice(5)), std::invalid_argument);
  options = {};
  options.contour_shift = 2.;
  EXPECT_THROW(sg::two_soliton_level(1., 1., 2., uni20::from_twice(1), options), std::invalid_argument);
  options = {};
  options.initial_intervals = 9;
  EXPECT_THROW(sg::two_soliton_level(1., 1., 2., uni20::from_twice(1), options), std::invalid_argument);
  options = {};
  options.max_intervals = 20000;
  EXPECT_THROW(sg::two_soliton_level(1., 1., 2., uni20::from_twice(1), options), std::length_error);
}
} // namespace
