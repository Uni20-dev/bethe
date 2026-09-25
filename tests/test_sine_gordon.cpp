// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/sine_gordon.hpp>
#include <bethe/sine_gordon_vacuum.hpp>

namespace
{
// Independent real D3 TBA at p=2 (Hegedus 2510.25344, eqs. 2.6-2.8).
// Positive-rapidity Gauss mesh, not the production shifted uniform grid.
// Integrate the constant magnon plateau analytically: integral s=1/2.
template <typename Real> Real d3_reference(Real u, std::size_t n)
{
  Real const pi = Real{4} * std::atan(Real{1}), log2 = std::log(Real{2});
  Real const eps = uni20::numeric_limits<Real>::epsilon(), b = -std::log(eps) / Real{2} + Real{4};
  auto const rule = bethe::detail::gauss_legendre<Real>(n);
  std::vector<Real> x(n), w(n), bare(n), v(n), a(n), l(n), next(n), k(n * n);
  for (std::size_t i = 0; i < n; ++i)
  {
    x[i] = b * (1 + rule.x[i]) / 2;
    w[i] = b * rule.w[i] / 2;
    bare[i] = u * std::cosh(x[i]) - log2;
  }
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      k[i * n + j] = w[j] / (2 * pi) * (1 / std::cosh(x[i] - x[j]) + 1 / std::cosh(x[i] + x[j]));
  for (std::size_t iter = 0; iter < 256; ++iter)
  {
    for (std::size_t j = 0; j < n; ++j)
      l[j] = std::log1p(std::exp(-bare[j] - v[j]));
    for (std::size_t i = 0; i < n; ++i)
    {
      bethe::detail::CompensatedSum<Real> sum;
      for (std::size_t j = 0; j < n; ++j)
        sum.add(-k[i * n + j] * l[j]);
      a[i] = std::log1p(std::exp(-sum.value())) - log2;
    }
    Real error{};
    for (std::size_t i = 0; i < n; ++i)
    {
      bethe::detail::CompensatedSum<Real> sum;
      for (std::size_t j = 0; j < n; ++j)
        sum.add(-2 * k[i * n + j] * a[j]);
      next[i] = sum.value();
      error = std::max(error, std::abs(next[i] - v[i]));
    }
    if (error < 32 * eps)
    {
      bethe::detail::CompensatedSum<Real> sum;
      for (std::size_t j = 0; j < n; ++j)
        sum.add(-w[j] * u * std::cosh(x[j]) * l[j] / pi);
      return sum.value();
    }
    v = next;
  }
  throw std::runtime_error("independent D3 reference failed to converge");
}
template <typename Real> class SineGordon : public ::testing::Test {};
TYPED_TEST_SUITE(SineGordon, test_support::RealTypes, test_support::PrecisionNames);
TEST(SineGordonVacuum, ClosedKernelReferences)
{
  for (auto const [p, expected] : {std::pair{0.5, -0.3860018985655546}, std::pair{2., -0.33346196574220655}})
  {
    auto const state = bethe::sine_gordon::vacuum_energy(1., 1., p);
    ASSERT_TRUE(state.converged) << int(state.status) << " mesh " << state.intervals << " iterations "
                                 << state.iterations << " kernel " << state.kernel_error << " mesh error "
                                 << state.mesh_error;
    EXPECT_NEAR(*state.scaling_function, expected, 1e-10);
    EXPECT_GT(state.kernel_evaluations, 0u);
  }
}
TYPED_TEST(SineGordon, VacuumFreeDiracAndBudgets)
{
  using Real = TypeParam;
  using namespace bethe::sine_gordon;
  auto const state = vacuum_energy(Real{1}, Real{1}, Real{1});
  ASSERT_TRUE(state.converged) << int(state.status) << " mesh " << state.intervals;
  EXPECT_REAL_NEAR(*state.scaling_function,
                   uni20::parse_real<Real>(
                       "-0.3456042161410258359590522072055135625147457700998595063575857703631793172742704659519"),
                   Real{1048576} * uni20::numeric_limits<Real>::epsilon());
  EXPECT_EQ(state.iterations, 0u);
  EXPECT_EQ(state.kernel_evaluations, 0u);
  auto const mesh = vacuum_energy(Real{1}, Real{1}, Real{1}, {.max_intervals = 64});
  EXPECT_EQ(mesh.status, VacuumStatus::mesh_limit);
  EXPECT_FALSE(mesh.casimir_energy);
  auto const cutoff = vacuum_energy(Real{1}, Real{1}, Real{1}, {.max_cutoffs = 1});
  EXPECT_EQ(cutoff.status, VacuumStatus::cutoff_limit);
  EXPECT_FALSE(cutoff.scaling_function);
}
TYPED_TEST(SineGordon, VacuumInteractingNative)
{
  using Real = TypeParam;
  auto const state = bethe::sine_gordon::vacuum_energy(Real{1}, Real{1}, Real{2});
  ASSERT_TRUE(state.converged) << int(state.status) << " intervals " << state.intervals << " kernel "
                               << uni20::format_scalar(state.kernel_error) << " mesh "
                               << uni20::format_scalar(state.mesh_error) << " iterations " << state.iterations;
  Real const reference = d3_reference(Real{1}, 384), refined = d3_reference(Real{1}, 512);
  Real const tolerance = Real{1048576} * uni20::numeric_limits<Real>::epsilon();
  EXPECT_REAL_NEAR(reference, refined, tolerance);
  EXPECT_REAL_NEAR(*state.scaling_function, refined, tolerance);
}
TYPED_TEST(SineGordon, DifferenceGridAgainstIndependentQuadrature)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  bethe::sine_gordon::KernelOptions<Real> options;
  options.tolerance = Real{16384} * uni20::numeric_limits<Real>::epsilon();
  std::size_t evaluations = 0;
  auto const table =
      bethe::sine_gordon::detail::kernel_table(Real{1.7}, Real{0.25}, Real{0.5}, 4, options, evaluations);
  ASSERT_TRUE(table.converged);
  for (std::size_t j : {0u, 2u, 4u})
  {
    auto const real = bethe::sine_gordon::scattering_kernel(C(Real(j) / 2, 0), Real{1.7});
    auto const shifted = bethe::sine_gordon::scattering_kernel(C(Real(j) / 2, Real{0.5}), Real{1.7});
    ASSERT_TRUE(real.converged);
    ASSERT_TRUE(shifted.converged);
    EXPECT_REAL_NEAR(table.real[j], real.value->real(), options.tolerance);
    EXPECT_REAL_NEAR(std::abs(table.shifted[j] - *shifted.value), Real{0}, options.tolerance);
  }
}
TEST(SineGordonVacuum, CouplingsScalingAndPhysicalLimits)
{
  using namespace bethe::sine_gordon;
  double const pi = 4 * std::atan(1.);
  for (double p : {.3, .7, 1.2, 3.7})
  {
    SCOPED_TRACE(p);
    VacuumOptions<double> options{.tolerance = 1e-8};
    auto const state = vacuum_energy(1., 1., p, options);
    ASSERT_TRUE(state.converged) << int(state.status);
    options.contour_shift = pi * std::min(1., p) / 5;
    auto const shifted = vacuum_energy(4., .25, p, options);
    ASSERT_TRUE(shifted.converged) << int(shifted.status);
    EXPECT_NEAR(*state.scaling_function, *shifted.scaling_function, 1e-8);
    EXPECT_NEAR(4 * *state.casimir_energy, *shifted.casimir_energy, 4e-8);
  }
  for (double p : {.5, 2.})
  {
    auto const uv = vacuum_energy(1., .0001, p, {.tolerance = 1e-8});
    ASSERT_TRUE(uv.converged) << int(uv.status);
    EXPECT_NEAR(*uv.effective_central_charge, 1., .001);
  }
  auto const ir = vacuum_energy(1., 10., .5);
  ASSERT_TRUE(ir.converged);
  double const solitons = -20 / pi * std::cyl_bessel_k(1., 10.);
  double const breather = -10 * std::sqrt(2.) / pi * std::cyl_bessel_k(1., 10 * std::sqrt(2.));
  EXPECT_NEAR(*ir.scaling_function, solitons + breather, 1e-8);
  EXPECT_GT(std::abs(*ir.scaling_function - solitons), 5e-7);
}
TEST(SineGordonVacuum, IndependentFailureBudgets)
{
  using namespace bethe::sine_gordon;
  auto const iteration = vacuum_energy(1., 1., 2., {.max_iterations = 0});
  EXPECT_EQ(iteration.status, VacuumStatus::iteration_limit);
  EXPECT_FALSE(iteration.casimir_energy);
  VacuumOptions<double> options;
  options.kernel.max_evaluations = 0;
  auto const kernel = vacuum_energy(1., 1., 2., options);
  EXPECT_EQ(kernel.status, VacuumStatus::kernel_limit);
  EXPECT_FALSE(kernel.scaling_function);
  auto const short_cutoff = vacuum_energy(1., 1., 1., {.tolerance = 1e-4, .max_cutoffs = 2, .initial_cutoff = .5});
  EXPECT_EQ(short_cutoff.status, VacuumStatus::cutoff_limit);
  EXPECT_FALSE(short_cutoff.casimir_energy);
  EXPECT_THROW(vacuum_energy(0., 1., 1.), std::invalid_argument);
  EXPECT_THROW(vacuum_energy(1., 1., 0.), std::invalid_argument);
  EXPECT_THROW(vacuum_energy(1., 1., .5, {.contour_shift = 1.}), std::invalid_argument);
  EXPECT_THROW(vacuum_energy(1., 1., 1., {.initial_intervals = 9}), std::invalid_argument);
  EXPECT_THROW(vacuum_energy(1., 1., 1., {.max_intervals = 20000}), std::length_error);
  EXPECT_THROW(vacuum_energy(1., 1., 1., {.tolerance = 0.}), std::invalid_argument);
  auto const overflow = vacuum_energy(1e308, 2., 1.);
  EXPECT_EQ(overflow.status, VacuumStatus::precision_limit);
  EXPECT_FALSE(overflow.casimir_energy);
  auto const range = vacuum_energy(1., 1., 1., {.initial_cutoff = 1e308});
  EXPECT_EQ(range.status, VacuumStatus::precision_limit);
  EXPECT_FALSE(range.scaling_function);
}
TYPED_TEST(SineGordon, KernelAgainstIndependentFourierTransforms)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const pi = Real{4} * std::atan(Real{1});
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real p : {Real{0.5}, Real{2}})
    for (C z : {C{}, C(Real{0.7}, Real{0.3}), C(Real{4}, Real{-0.6}), C(0, Real{1.2})})
    {
      auto const result = bethe::sine_gordon::scattering_kernel(z, p);
      ASSERT_TRUE(result.converged) << int(result.status);
      // Closed inverse transforms of sech(pi*k/2) and sech²(pi*k/2).
      // The implementation deliberately does not special-case either p.
      C expected;
      if (p == Real{0.5})
        expected = -Real{1} / (Real{2} * pi * std::cosh(z));
      else
        expected = z == C{} ? C(Real{1} / (Real{2} * pi * pi)) : z / (Real{2} * pi * pi * std::sinh(z));
      EXPECT_REAL_NEAR(std::abs(*result.value - expected), Real{0}, Real{4096} * eps);
      EXPECT_LE(result.quadrature_error + result.tail_bound, Real{1024} * eps);
    }
}
TYPED_TEST(SineGordon, AnalyticSymmetriesAndNearFreeCoupling)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real p : {Real{0.3}, Real{0.8}, Real{1.2}, Real{3.7}, Real{1000}})
  {
    C const z(Real{1.3}, Real{0.2});
    auto const a = bethe::sine_gordon::scattering_kernel(z, p);
    auto const b = bethe::sine_gordon::scattering_kernel(-z, p);
    auto const c = bethe::sine_gordon::scattering_kernel(std::conj(z), p);
    ASSERT_TRUE(a.converged);
    ASSERT_TRUE(b.converged);
    ASSERT_TRUE(c.converged);
    EXPECT_REAL_NEAR(std::abs(*a.value - *b.value), Real{0}, Real{32} * eps);
    EXPECT_REAL_NEAR(std::abs(*a.value - std::conj(*c.value)), Real{0}, Real{32} * eps);
  }
  // At z=0 the derivative dG/dp at p=1 is 1/8:
  // integral_0^infinity k/(2*sinh(pi*k)) dk = 1/8.
  for (Real sign : {Real{-1}, Real{1}})
  {
    Real const p = Real{1} + sign * Real{8} * eps;
    bethe::sine_gordon::KernelOptions<Real> options;
    options.tolerance *= std::abs(p - Real{1});
    auto const result = bethe::sine_gordon::scattering_kernel(C{}, p, options);
    ASSERT_TRUE(result.converged);
    EXPECT_REAL_NEAR(result.value->real() / (p - Real{1}), Real{1} / Real{8}, Real{4096} * eps);
  }
}
TYPED_TEST(SineGordon, BudgetsAndValidation)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  using namespace bethe::sine_gordon;
  auto const free = scattering_kernel(C{}, Real{1}, {.max_evaluations = 0, .max_cutoffs = 0});
  ASSERT_TRUE(free.converged);
  EXPECT_EQ(free.value->real(), Real{0});
  EXPECT_EQ(free.value->imag(), Real{0});
  EXPECT_EQ(free.evaluations, 0u);
  auto const cutoff = scattering_kernel(C{}, Real{2}, {.max_cutoffs = 0});
  EXPECT_EQ(cutoff.status, KernelStatus::cutoff_limit);
  EXPECT_FALSE(cutoff.value);
  auto const mesh = scattering_kernel(C{}, Real{2}, {.max_evaluations = 1});
  EXPECT_EQ(mesh.status, KernelStatus::quadrature_limit);
  EXPECT_FALSE(mesh.value);
  auto const levels = scattering_kernel(C{}, Real{2}, {.max_levels = 0});
  EXPECT_EQ(levels.status, KernelStatus::quadrature_limit);
  EXPECT_FALSE(levels.value);
  EXPECT_THROW(scattering_kernel(C{}, Real{0}), std::invalid_argument);
  EXPECT_THROW(scattering_kernel(C(0, 2), Real{0.5}), std::invalid_argument);
  EXPECT_THROW(scattering_kernel(C{}, Real{2}, {.tolerance = Real{0}}), std::invalid_argument);
  EXPECT_THROW(scattering_kernel(C{}, uni20::numeric_limits<Real>::infinity()), std::invalid_argument);
}
TYPED_TEST(SineGordon, NativeReferenceAndOverflowSafeContour)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  C const expected(
      uni20::parse_real<Real>(
          "0.0547765717918015314046546239425095204723847959115512148661684638805025105506472402803056409"),
      uni20::parse_real<Real>(
          "-0.00515484562906219103918414319730092148237776521629576139173090753048764604314854689936674796"));
  auto const result = bethe::sine_gordon::scattering_kernel(
      C(uni20::parse_real<Real>("0.8"), uni20::parse_real<Real>("0.4")), uni20::parse_real<Real>("2.7"));
  ASSERT_TRUE(result.converged);
  EXPECT_REAL_NEAR(std::abs(*result.value - expected), Real{0}, Real{4096} * eps);
  // Separate cosh(k*Im(z)) would overflow in fp64, although their product
  // with the exponentially small Fourier multiplier is well conditioned.
  Real const y = pi - Real{1} / Real{1000};
  auto const integrand = bethe::sine_gordon::detail::fourier_integrand(Real{1000}, C(0, y), Real{2});
  EXPECT_REAL_NEAR(integrand.real(), std::exp(-(pi - y) * Real{1000}) / (Real{2} * pi), Real{128} * eps);
  EXPECT_EQ(integrand.imag(), Real{0});
}
} // namespace
