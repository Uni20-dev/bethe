// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/hubbard_doped.hpp>

namespace
{
namespace model = bethe::hubbard::thermo;
template <typename Real> class HubbardDoped : public ::testing::Test {};
TYPED_TEST_SUITE(HubbardDoped, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(HubbardDoped, KernelReferences)
{
  using Real = TypeParam;
  struct Reference
  {
      char const *u, *x, *r, *primitive;
  };
  for (auto const& ref : {Reference{"1", "0", "0.110317800076325796698228216058998845491344874", "0"},
                          Reference{"1", "0.5", "0.101918067808414701059382271725519756019417106",
                                    "0.0537220747099031648023590169601767848974015185"},
                          Reference{"0.0625", "2", "0.00249169160345381204657485098246890837709411859",
                                    "0.24502315465044781119134400495527301838446631"},
                          Reference{"0.0001", "2", "0.00000397887359719175167859638569537560585946105257",
                                    "0.249992042252832142321128076252593915681220371"}})
  {
    Real const u = uni20::parse_real<Real>(ref.u), x = uni20::parse_real<Real>(ref.x);
    Real const expected = uni20::parse_real<Real>(ref.r), eps = uni20::numeric_limits<Real>::epsilon();
    EXPECT_REAL_NEAR(bethe::detail::hubbard_r(x, u), expected, Real{64} * eps * expected);
    EXPECT_REAL_NEAR(bethe::detail::hubbard_r_primitive(x, u), uni20::parse_real<Real>(ref.primitive), Real{64} * eps);
    EXPECT_EQ(bethe::detail::hubbard_r(-x, u), bethe::detail::hubbard_r(x, u));
    EXPECT_EQ(bethe::detail::hubbard_r_primitive(-x, u), -bethe::detail::hubbard_r_primitive(x, u));
  }
}
TYPED_TEST(HubbardDoped, NativeGaussRule)
{
  using Real = TypeParam;
  auto const rule = bethe::detail::gauss_legendre<Real>(16);
  for (int power : {0, 2, 4, 16, 30})
  {
    bethe::detail::CompensatedSum<Real> sum;
    for (std::size_t j = 0; j < rule.x.size(); ++j)
      sum.add(rule.w[j] * std::pow(rule.x[j], power));
    EXPECT_REAL_NEAR(sum.value(), Real{2} / Real(power + 1), Real{64} * uni20::numeric_limits<Real>::epsilon());
  }
}
TYPED_TEST(HubbardDoped, Backgrounds)
{
  using Real = TypeParam;
  for (Real u : {Real{1}, Real{4}, Real{20}})
    for (Real density : {Real{.25}, Real{.5}, Real{.875}})
    {
      model::DopedSolver<Real> solver(u, density);
      auto const& b = solver.background();
      ASSERT_TRUE(b.converged) << uni20::format_real(u) << ',' << uni20::format_real(density) << ',' << int(b.status)
                               << ',' << b.nodes << ',' << b.iterations << ',' << uni20::format_real(b.mesh_error);
      EXPECT_LT(*b.energy_per_site_unshifted, Real{0});
      EXPECT_REAL_NEAR(*b.mu_unshifted - *b.mu_symmetric, u / Real{2},
                       Real{32} * uni20::numeric_limits<Real>::epsilon());
    }
}

TYPED_TEST(HubbardDoped, DressedMomentumGridsAndFermiEndpoints)
{
  using Real = TypeParam;
  model::DopedSolver<Real> solver(Real{4}, Real{.5});
  ASSERT_TRUE(solver.background().converged);
  for (auto branch : {model::DopedBranch::spinon, model::DopedBranch::holon, model::DopedBranch::charge_particle})
  {
    auto const [lo, hi] = solver.momentum_range(branch);
    for (int i = 0; i <= 8; ++i)
    {
      Real const momentum = i == 8 ? hi : lo + (hi - lo) * Real(i) / Real{8};
      auto const p = solver.at_momentum(branch, momentum, model::Convention::symmetric, model::EnergyReference::fermi);
      ASSERT_TRUE(p.converged) << int(branch) << ',' << i << ',' << int(p.status) << ',' << p.iterations;
      ASSERT_TRUE(p.energy);
      EXPECT_GE(*p.energy, Real{0});
      EXPECT_EQ(p.energy, p.fermi_energy);
      if (i == 0 || i == 8)
        EXPECT_EQ(*p.energy, Real{0});
      else
        EXPECT_GT(*p.energy, Real{0});
      auto const h = solver.at_momentum(branch, momentum);
      ASSERT_TRUE(h.converged);
      EXPECT_REAL_NEAR(*h.energy - *p.energy, *solver.background().mu_symmetric * Real(p.delta_particles),
                       Real{64} * uni20::numeric_limits<Real>::epsilon());
      auto const unshifted = solver.at_momentum(branch, momentum, model::Convention::unshifted);
      ASSERT_TRUE(unshifted.converged);
      EXPECT_REAL_NEAR(*unshifted.energy - *h.energy, Real{2} * Real(p.delta_particles),
                       Real{64} * uni20::numeric_limits<Real>::epsilon());
      auto const fermi =
          solver.at_momentum(branch, momentum, model::Convention::unshifted, model::EnergyReference::fermi);
      EXPECT_EQ(fermi.energy, p.energy);
      EXPECT_EQ(unshifted.symmetric_energy, h.energy);
    }
  }
}

TYPED_TEST(HubbardDoped, DensityEnergyAdjointIdentity)
{
  using Real = TypeParam;
  auto const rule = bethe::detail::gauss_legendre<Real>(48);
  auto const mesh = model::detail::doped_mesh(Real{1}, Real{1}, rule);
  ASSERT_TRUE(mesh);
  bethe::detail::CompensatedSum<Real> sum;
  for (std::size_t j = 0; j < mesh->w.size(); ++j)
    sum.add(mesh->w[j] * mesh->epsilon[j]);
  Real const energy = sum.value() / model::detail::pi<Real>() + mesh->mu * mesh->density;
  EXPECT_REAL_NEAR(energy, mesh->e0, Real{512} * uni20::numeric_limits<Real>::epsilon());
}

TYPED_TEST(HubbardDoped, IndependentFullIntervalReferences)
{
  using Real = TypeParam;
  // mpmath, 48 decimal places, full [-Q,Q] Nyström at Q=1, U=4.
  // Digamma/log-gamma kernels, unlike the native Euler/Stirling kernels.
  // 72 and 96 nodes agree to every quoted digit; includes momentum inversion.
  auto value = [](char const* s) { return uni20::parse_real<Real>(s); };
  Real const density = value("0.38002501150130312337295694461904150571011809");
  Real const tolerance = Real{8192} * uni20::numeric_limits<Real>::epsilon();
  model::DopedSolver<Real> solver(Real{4}, density);
  auto const& b = solver.background();
  ASSERT_TRUE(b.converged);
  EXPECT_REAL_NEAR(*b.fermi_rapidity, Real{1}, tolerance);
  EXPECT_REAL_NEAR(*b.mu_unshifted, value("-1.1780462378112884478262187086713965517723693"), tolerance);
  EXPECT_REAL_NEAR(*b.energy_per_site_unshifted, value("-0.64315142876898352695897944974573494019373884"), tolerance);
  Real const offset = model::detail::pi<Real>() * density / Real{2};
  struct Reference
  {
      model::DopedBranch branch;
      Real momentum, energy, parameter;
  };
  for (auto const& ref :
       {Reference{model::DopedBranch::spinon, value("0.37467201934392443769510624564135004333897292"),
                  value("0.20027947161009807815914556412767695821058585"), Real{.5}},
        Reference{model::DopedBranch::holon, offset - value("0.61450185409484899732028805466660613845598388"),
                  value("0.68725570114182347407704028173483940897229027"), Real{.5}},
        Reference{model::DopedBranch::charge_particle, value("2.2077082495637365163877961802547245760849637") - offset,
                  value("1.9157127471381251634674330705152632804379892"), Real{2}}})
  {
    auto const point =
        solver.at_momentum(ref.branch, ref.momentum, model::Convention::symmetric, model::EnergyReference::fermi);
    ASSERT_TRUE(point.converged) << int(ref.branch) << ',' << int(point.status);
    EXPECT_REAL_NEAR(*point.energy, ref.energy, tolerance);
    EXPECT_REAL_NEAR(point.parameter, ref.parameter, tolerance);
  }
}

TYPED_TEST(HubbardDoped, StrongCouplingLimit)
{
  using Real = TypeParam;
  Real const u{1000000}, density = Real{3} / Real{8}, pi = model::detail::pi<Real>();
  model::DopedSolver<Real> solver(u, density);
  auto const& b = solver.background();
  ASSERT_TRUE(b.converged);
  EXPECT_REAL_NEAR(*b.fermi_rapidity, pi * density, Real{16} / u);
  EXPECT_REAL_NEAR(*b.energy_per_site_unshifted, -Real{2} * std::sin(pi * density) / pi, Real{16} / u);
  EXPECT_REAL_NEAR(*b.mu_unshifted, -Real{2} * std::cos(pi * density), Real{16} / u);
  for (auto branch : {model::DopedBranch::spinon, model::DopedBranch::holon, model::DopedBranch::charge_particle})
  {
    auto const [lo, hi] = solver.momentum_range(branch);
    auto const point =
        solver.at_momentum(branch, (lo + hi) / Real{2}, model::Convention::symmetric, model::EnergyReference::fermi);
    ASSERT_TRUE(point.converged) << int(point.status);
    Real const expected = branch == model::DopedBranch::spinon  ? Real{0}
                          : branch == model::DopedBranch::holon ? Real{2} * (Real{1} - std::cos(pi * density))
                                                                : Real{2} * (Real{1} + std::cos(pi * density));
    EXPECT_REAL_NEAR(*point.energy, expected, Real{32} / u);
  }
}

TYPED_TEST(HubbardDoped, HalfFillingLimitUsesLowerMottEdge)
{
  using Real = TypeParam;
  Real const doping = uni20::parse_real<Real>("0.00001"), pi = model::detail::pi<Real>();
  model::DopedSolver<Real> solver(Real{4}, Real{1} - doping);
  ASSERT_TRUE(solver.background().converged);
  auto const minimum = model::dispersion(model::Branch::holon, Real{4}, -pi / Real{2});
  ASSERT_TRUE(minimum.converged);
  EXPECT_REAL_NEAR(*solver.background().mu_symmetric, -*minimum.energy, Real{20} * doping);
  for (auto branch : {model::DopedBranch::spinon, model::DopedBranch::holon})
  {
    auto const doped = solver.at_momentum(branch, Real{1});
    auto const half = model::dispersion(
        branch == model::DopedBranch::spinon ? model::Branch::spinon : model::Branch::holon, Real{4}, Real{1});
    ASSERT_TRUE(doped.converged) << int(doped.status);
    ASSERT_TRUE(half.converged);
    EXPECT_REAL_NEAR(*doped.energy, *half.energy, Real{20} * doping);
  }
}

TYPED_TEST(HubbardDoped, FailuresNeverPublishEnergies)
{
  using Real = TypeParam;
  model::DopedOptions<Real> options;
  options.max_nodes = options.initial_nodes;
  model::DopedSolver<Real> mesh_limited(Real{4}, Real{.5}, options);
  auto const& b = mesh_limited.background();
  EXPECT_FALSE(b.converged);
  EXPECT_EQ(b.status, model::DopedStatus::mesh_limit);
  EXPECT_FALSE(b.mu_unshifted);
  EXPECT_FALSE(b.fermi_rapidity);
  auto const failed = mesh_limited.at_momentum(model::DopedBranch::holon, Real{0});
  EXPECT_FALSE(failed.energy);
  EXPECT_FALSE(failed.symmetric_energy);
  EXPECT_FALSE(failed.fermi_energy);
  options = {};
  options.max_background_iterations = 0;
  model::DopedSolver<Real> density_limited(Real{4}, Real{.5}, options);
  EXPECT_EQ(density_limited.background().status, model::DopedStatus::density_limit);
  EXPECT_FALSE(density_limited.background().energy_per_site_unshifted);
  options = {};
  options.max_iterations = 0;
  model::DopedSolver<Real> momentum_limited(Real{4}, Real{.5}, options);
  ASSERT_TRUE(momentum_limited.background().converged);
  auto const point = momentum_limited.at_momentum(model::DopedBranch::holon, Real{.3});
  EXPECT_EQ(point.status, model::DopedStatus::momentum_limit);
  EXPECT_FALSE(point.energy);
  EXPECT_FALSE(point.fermi_energy);
  // Analytic endpoints need no momentum iterations.
  EXPECT_TRUE(momentum_limited.at_momentum(model::DopedBranch::spinon, Real{0}).converged);
}

TYPED_TEST(HubbardDoped, InvalidInputs)
{
  using Real = TypeParam;
  for (Real density : {Real{0}, Real{1}, Real{-1}, uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW((model::DopedSolver<Real>(Real{4}, density)), std::invalid_argument);
  for (Real u : {Real{0}, Real{-1}, uni20::numeric_limits<Real>::infinity()})
    EXPECT_THROW((model::DopedSolver<Real>(u, Real{.5})), std::invalid_argument);
  model::DopedOptions<Real> options;
  options.initial_nodes = 3;
  EXPECT_THROW((model::DopedSolver<Real>(Real{4}, Real{.5}, options)), std::invalid_argument);
  options = {};
  options.max_nodes = 513;
  EXPECT_THROW((model::DopedSolver<Real>(Real{4}, Real{.5}, options)), std::invalid_argument);
  options = {};
  options.tolerance = Real{0};
  EXPECT_THROW((model::DopedSolver<Real>(Real{4}, Real{.5}, options)), std::invalid_argument);
  model::DopedSolver<Real> solver(Real{4}, Real{.5});
  EXPECT_THROW(solver.at_momentum(model::DopedBranch::spinon, Real{-1}), std::invalid_argument);
  EXPECT_THROW(solver.at_momentum(model::DopedBranch::holon, Real{4}), std::invalid_argument);
  EXPECT_THROW(solver.momentum_range(static_cast<model::DopedBranch>(99)), std::invalid_argument);
  EXPECT_THROW(solver.at_momentum(model::DopedBranch::holon, Real{0}, static_cast<model::Convention>(99)),
               std::invalid_argument);
  EXPECT_THROW(solver.at_momentum(model::DopedBranch::holon, Real{0}, model::Convention::symmetric,
                                  static_cast<model::EnergyReference>(99)),
               std::invalid_argument);
}

TYPED_TEST(HubbardDoped, InteriorCancellationNeverReportedAsZero)
{
  using Real = TypeParam;
  // The dilute-sea Fermi energy is quadratic in density: at n=epsilon,
  // subtracting the order-one chemical potential cannot resolve it.
  model::DopedSolver<Real> solver(Real{4}, uni20::numeric_limits<Real>::epsilon());
  ASSERT_TRUE(solver.background().converged);
  auto const [lo, hi] = solver.momentum_range(model::DopedBranch::holon);
  auto const point = solver.at_momentum(model::DopedBranch::holon, (lo + hi) / Real{2}, model::Convention::symmetric,
                                        model::EnergyReference::fermi);
  EXPECT_EQ(point.status, model::DopedStatus::precision_limit);
  EXPECT_FALSE(point.converged);
  EXPECT_FALSE(point.energy);
  EXPECT_FALSE(point.fermi_energy);
}

TEST(HubbardDopedThermodynamics, ChemicalPotentialIsGroundEnergyDerivative)
{
  double const step = 1.e-4;
  model::DopedSolver<double> center(4, .5), left(4, .5 - step), right(4, .5 + step);
  ASSERT_TRUE(center.background().converged && left.background().converged && right.background().converged);
  double const derivative =
      (*right.background().energy_per_site_unshifted - *left.background().energy_per_site_unshifted) / (2 * step);
  EXPECT_NEAR(derivative, *center.background().mu_unshifted, 1.e-7);
}
} // namespace
