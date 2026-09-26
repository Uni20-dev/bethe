// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/sine_gordon_bethe_yang.hpp>

namespace
{
namespace sg = bethe::sine_gordon;
template <typename Real> class SineGordonExcitations : public ::testing::Test {};
TYPED_TEST_SUITE(SineGordonExcitations, test_support::RealTypes, test_support::PrecisionNames);
constexpr sg::Particle soliton{sg::ParticleKind::soliton}, anti{sg::ParticleKind::antisoliton};
constexpr sg::Particle b1{sg::ParticleKind::breather, 1}, b2{sg::ParticleKind::breather, 2};

TEST(SineGordonExcitationLimits, UnresolvableHalfOddLabel)
{
  auto const state = sg::same_charge_pair(1., 10., 1., uni20::from_twice(std::int64_t{9007199254740993}));
  EXPECT_EQ(state.status, sg::BetheYangStatus::precision_limit);
  EXPECT_FALSE(state.energy);
}

TYPED_TEST(SineGordonExcitations, ParticlesThresholdsAndSchwingerRatios)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  sg::ParticleSpectrum<Real> const spectrum(Real{1}, Real{1} / Real{3});
  EXPECT_EQ(spectrum.particles().size(), 4u);
  EXPECT_EQ(soliton.charge(), 1);
  EXPECT_EQ(anti.charge(), -1);
  EXPECT_EQ(b1.charge(), 0);
  EXPECT_REAL_NEAR(spectrum.mass(b1), Real{1}, Real{8} * eps);
  EXPECT_REAL_NEAR(spectrum.mass(b2), std::sqrt(Real{3}), Real{8} * eps);
  EXPECT_EQ(spectrum.mass(soliton), spectrum.mass(anti));
  EXPECT_FALSE(spectrum.exists({sg::ParticleKind::breather, 3}));
  std::array<sg::Particle, 2> pair{soliton, anti};
  EXPECT_REAL_NEAR(spectrum.threshold(pair, Real{3}), std::sqrt(Real{13}), Real{8} * eps);
  EXPECT_EQ(spectrum.threshold(pair, -Real{3}), spectrum.threshold(pair, Real{3}));
  // Equal velocities minimize fixed-momentum relativistic energy, not
  // equal momenta when the masses differ.
  pair = {b1, b2};
  Real const total_mass = spectrum.mass(b1) + spectrum.mass(b2);
  Real const threshold = spectrum.threshold(pair, Real{2});
  EXPECT_REAL_NEAR(threshold,
                   spectrum.energy(b1, Real{2} * spectrum.mass(b1) / total_mass) +
                       spectrum.energy(b2, Real{2} * spectrum.mass(b2) / total_mass),
                   Real{32} * eps);
  EXPECT_EQ(sg::ParticleSpectrum<Real>(Real{1}, Real{1}).particles().size(), 2u);
  EXPECT_EQ(sg::ParticleSpectrum<Real>(Real{1}, Real{2}).particles().size(), 2u);
  EXPECT_EQ(sg::ParticleSpectrum<Real>(Real{1}, Real{0.5}).particles().size(), 3u);
  EXPECT_TRUE(sg::ParticleSpectrum<Real>(Real{1}, Real{0.5} - eps).exists(b2));
  EXPECT_FALSE(sg::ParticleSpectrum<Real>(Real{1}, Real{0.5} + eps).exists(b2));
  sg::ParticleSpectrum<Real> const native(Real{1}, uni20::parse_real<Real>("0.4"));
  EXPECT_REAL_NEAR(
      native.energy(b1, uni20::parse_real<Real>("0.7")),
      uni20::parse_real<Real>("1.368198089185226552659329247836824968488111471138021993646805247075555208134787952"),
      Real{16} * eps);
}

TYPED_TEST(SineGordonExcitations, ScatteringClosedReferencesAndDerivative)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  Real const theta = uni20::parse_real<Real>("0.7");
  auto const attractive = sg::soliton_phase(theta, Real{0.5});
  ASSERT_TRUE(attractive.converged);
  EXPECT_REAL_NEAR(*attractive.phase, -std::atan(std::sinh(theta)), Real{1024} * eps);
  auto const repulsive = sg::soliton_phase(theta, Real{2});
  ASSERT_TRUE(repulsive.converged);
  // Independent finite-interval integral of x/(pi*sinh(x)), not the
  // production Fourier quadrature. Reference generated at 90 digits.
  EXPECT_REAL_NEAR(
      *repulsive.phase,
      uni20::parse_real<Real>("0.216952007696309504212739785465189529934336381010650385652531085573082116153999453"),
      Real{1024} * eps);
  auto const negative = sg::soliton_phase(-theta, Real{2});
  ASSERT_TRUE(negative.converged);
  EXPECT_EQ(*negative.phase, -*repulsive.phase);
  EXPECT_EQ(*sg::soliton_phase(theta, Real{1}).phase, Real{0});
  EXPECT_EQ(*sg::soliton_phase(Real{0}, Real{2}).phase, Real{0});
  Real const near_free_p = Real{1} + Real{8} * eps;
  sg::KernelOptions<Real> near_free_options;
  near_free_options.tolerance = Real{4096} * eps * (near_free_p - Real{1});
  auto const near_free = sg::soliton_phase(theta, near_free_p, near_free_options);
  ASSERT_TRUE(near_free.converged);
  EXPECT_GT(*near_free.phase, Real{0});
  // d chi/dp at p=1, from the independent sech^2 Fourier transform.
  EXPECT_REAL_NEAR(*near_free.phase / (near_free_p - Real{1}), pi / Real{2} * std::tanh(theta / Real{2}),
                   Real{4096} * eps);
  Real const step = std::pow(eps, Real{1} / Real{4});
  auto const plus = sg::soliton_phase(theta + step, Real{1.7});
  auto const minus = sg::soliton_phase(theta - step, Real{1.7});
  auto const kernel = sg::scattering_kernel(std::complex<Real>(theta, Real{0}), Real{1.7});
  ASSERT_TRUE(plus.converged && minus.converged && kernel.converged);
  EXPECT_REAL_NEAR((*plus.phase - *minus.phase) / (Real{2} * step), Real{2} * pi * kernel.value->real(),
                   Real{64} * std::sqrt(eps));
}

TYPED_TEST(SineGordonExcitations, BetheYangFreeAndInteractingQuantization)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  auto const half = uni20::from_twice(1);
  auto const free = sg::same_charge_pair(Real{1}, Real{10}, Real{1}, half);
  ASSERT_TRUE(free.converged);
  EXPECT_EQ(free.evaluations, 0u);
  EXPECT_REAL_NEAR(*free.energy, Real{2} * std::hypot(Real{1}, pi / Real{10}), Real{64} * eps);
  auto const state = sg::same_charge_pair(Real{1}, Real{10}, Real{2}, half);
  ASSERT_TRUE(state.converged) << int(state.status);
  EXPECT_REAL_NEAR(
      *state.energy,
      uni20::parse_real<Real>("2.085735618785140488030977932237507957090038265403091568663327979382012603581455256"),
      Real{16384} * eps);
  EXPECT_LT(*state.energy, *free.energy);
  auto const phase = sg::soliton_phase(Real{2} * *state.rapidity, Real{2});
  ASSERT_TRUE(phase.converged);
  EXPECT_REAL_NEAR(Real{10} * std::sinh(*state.rapidity) + *phase.phase, pi, Real{16384} * eps);
  // Check the original exponential equation, including the minus sign.
  auto const exponential = -std::polar(Real{1}, Real{10} * std::sinh(*state.rapidity) + *phase.phase);
  EXPECT_REAL_NEAR(std::abs(exponential - std::complex<Real>(Real{1}, Real{0})), Real{0}, Real{32768} * eps);
  auto const scaled = sg::same_charge_pair(Real{2}, Real{5}, Real{2}, half);
  ASSERT_TRUE(scaled.converged);
  EXPECT_EQ(*scaled.rapidity, *state.rapidity);
  EXPECT_EQ(*scaled.energy, Real{2} * *state.energy);
  auto const higher = sg::same_charge_pair(Real{1}, Real{10}, Real{2}, uni20::from_twice(3));
  ASSERT_TRUE(higher.converged);
  EXPECT_GT(*higher.energy, *state.energy);
}

TYPED_TEST(SineGordonExcitations, ValidationAndDistinctNumericalBudgets)
{
  using Real = TypeParam;
  EXPECT_THROW(sg::ParticleSpectrum<Real>(Real{0}, Real{1}), std::invalid_argument);
  sg::ParticleSpectrum<Real> const spectrum(Real{1}, Real{0.5});
  EXPECT_THROW(spectrum.mass(b2), std::invalid_argument);
  EXPECT_THROW(spectrum.mass({sg::ParticleKind::soliton, 1}), std::invalid_argument);
  EXPECT_THROW(spectrum.particles(0), std::length_error);
  EXPECT_THROW(sg::soliton_phase(Real{1}, Real{0}), std::invalid_argument);
  auto const cutoff = sg::soliton_phase(Real{1}, Real{2}, {.max_cutoffs = 0});
  EXPECT_EQ(cutoff.status, sg::KernelStatus::cutoff_limit);
  EXPECT_FALSE(cutoff.phase);
  auto const quad = sg::soliton_phase(Real{1}, Real{2}, {.max_evaluations = 0});
  EXPECT_EQ(quad.status, sg::KernelStatus::quadrature_limit);
  EXPECT_FALSE(quad.phase);
  EXPECT_THROW(sg::same_charge_pair(Real{1}, Real{10}, Real{0.5}, uni20::from_twice(1)), std::invalid_argument);
  EXPECT_THROW(sg::same_charge_pair(Real{1}, Real{10}, Real{2}, uni20::half_int(1)), std::invalid_argument);
  auto const no_iterations =
      sg::same_charge_pair(Real{1}, Real{10}, Real{2}, uni20::from_twice(1), {.max_iterations = 0});
  EXPECT_EQ(no_iterations.status, sg::BetheYangStatus::iteration_limit);
  EXPECT_FALSE(no_iterations.energy);
  sg::BetheYangOptions<Real> options;
  options.phase.max_evaluations = 0;
  auto const no_phase = sg::same_charge_pair(Real{1}, Real{10}, Real{2}, uni20::from_twice(1), options);
  EXPECT_EQ(no_phase.status, sg::BetheYangStatus::phase_limit);
  EXPECT_EQ(no_phase.phase_status, sg::KernelStatus::quadrature_limit);
  EXPECT_FALSE(no_phase.energy);
  Real const tiny = uni20::numeric_limits<Real>::min(), huge = uni20::numeric_limits<Real>::max();
  EXPECT_THROW(sg::ParticleSpectrum<Real>(tiny, tiny).mass(b1), std::underflow_error);
  EXPECT_THROW(sg::ParticleSpectrum<Real>(huge, Real{0.5}).mass(b1), std::overflow_error);
  EXPECT_THROW(spectrum.energy(b1, uni20::numeric_limits<Real>::infinity()), std::invalid_argument);
  auto const underflow = sg::same_charge_pair(tiny, tiny, Real{2}, uni20::from_twice(1));
  EXPECT_EQ(underflow.status, sg::BetheYangStatus::precision_limit);
  EXPECT_FALSE(underflow.energy);
}
} // namespace
