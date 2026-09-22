// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/hubbard_thermo.hpp>

namespace
{
namespace model = bethe::hubbard::thermo;
template <typename Real> class HubbardThermo : public ::testing::Test {};
TYPED_TEST_SUITE(HubbardThermo, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(HubbardThermo, QuadratureAndEndpointSingularities)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  std::size_t count = 0;
  auto const integral =
      bethe::detail::tanh_sinh<Real, 3>([](Real x) { return std::array<Real, 3>{Real{1}, std::sqrt(x), std::log(x)}; },
                                        Real{0}, Real{1}, Real{64} * eps, count, 200000, 12);
  ASSERT_TRUE(integral.converged) << count;
  EXPECT_REAL_NEAR(integral.value[0], Real{1}, Real{64} * eps);
  EXPECT_REAL_NEAR(integral.value[1], Real{2} / Real{3}, Real{64} * eps);
  EXPECT_REAL_NEAR(integral.value[2], Real{-1}, Real{64} * eps);
}

TYPED_TEST(HubbardThermo, ParametricLinesAcrossCouplings)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1});
  for (Real u : {Real{1} / Real{4}, Real{1}, Real{4}, Real{20}, Real{100}})
  {
    for (Real lambda : {Real{0}, Real{.5}, Real{1}, Real{2}, Real{10}})
    {
      auto const s = model::spinon_at_rapidity(u, lambda);
      ASSERT_TRUE(s.converged) << uni20::format_real(u) << "," << uni20::format_real(lambda) << "," << s.evaluations;
      ASSERT_TRUE(s.energy);
      EXPECT_GT(*s.energy, Real{0});
    }
    for (Real k : {Real{0}, Real{.3}, pi / Real{2}, Real{2}, pi, -pi / Real{2}})
    {
      auto const h = model::charge_at_bare_momentum(model::Branch::holon, u, k);
      ASSERT_TRUE(h.converged) << uni20::format_real(u) << "," << uni20::format_real(k) << "," << h.evaluations;
      ASSERT_TRUE(h.energy);
      EXPECT_GT(*h.energy, Real{0});
    }
  }
}

TYPED_TEST(HubbardThermo, PhysicalMomentumGridAcrossCouplings)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1});
  for (Real u : {Real{.25}, Real{1}, Real{4}, Real{20}})
    for (auto branch : {model::Branch::spinon, model::Branch::holon, model::Branch::antiholon})
      for (int i = 0; i <= 8; ++i)
      {
        Real const p = branch == model::Branch::spinon ? pi * Real(i) / Real{8} : pi * (Real(i) / Real{4} - Real{1});
        auto const s = model::dispersion(branch, u, p);
        ASSERT_TRUE(s.converged) << int(branch) << "," << uni20::format_real(u) << "," << i << "," << int(s.status)
                                 << "," << s.evaluations << "," << s.iterations;
        ASSERT_TRUE(s.energy);
        EXPECT_GE(*s.energy, Real{0});
        EXPECT_LE(s.momentum_error, model::Options<Real>{}.relative_tolerance *
                                        (branch == model::Branch::spinon ? std::min(p, pi - p) : Real{1}));
      }
}

TYPED_TEST(HubbardThermo, IndependentFourierBesselReferences)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1});
  Real const tol = Real{512} * uni20::numeric_limits<Real>::epsilon();
  // Original Essler-Korepin Fourier-Bessel integrals, evaluated independently
  // with mpmath at 48 decimal digits, split into panels through omega=120/u.
  // Our implementation instead uses nonoscillatory convolution/K-series forms.
  struct Reference
  {
      Real u, lambda, k;
      char const *es, *ps, *ec, *pc;
  };
  for (auto const& ref :
       {Reference{Real{1}, Real{.5}, Real{2}, "1.62363551638715649341503641522580392287950694",
                  "1.01875209564840617568054248964328257088350928", "0.407368939609261750044336500912450251711767674",
                  "-1.38101056299565894567309850812931709644086994"},
        Reference{Real{4}, Real{0}, pi, "1.24228148941956239984591553194945497655659639",
                  "1.5707963267948966192313216916397514420985847", "0.643363511006452197366294884644872250891658422",
                  "-1.5707963267948966192313216916397514420985847"},
        Reference{Real{4}, Real{2}, Real{1}, "0.180195419367165283255675569344877169275747983",
                  "0.147435530616030589362559903539606986354219786", "3.62482643594807997383907104208212004914561354",
                  "0.0869244707432127794663054310826567415930967806"},
        Reference{Real{8}, Real{1}, pi / Real{2}, "0.583988260387869747204901399931155905248506344",
                  "0.925124904640531102000028987528110039741546307", "4.31565065357205341928655078233197306290634754",
                  "-0.325690058692956251671241842752280188869145233"}})
  {
    SCOPED_TRACE(uni20::format_real(ref.u));
    auto const s = model::spinon_at_rapidity(ref.u, ref.lambda);
    auto const c = model::charge_at_bare_momentum(model::Branch::holon, ref.u, ref.k);
    ASSERT_TRUE(s.converged);
    ASSERT_TRUE(c.converged);
    EXPECT_REAL_NEAR(*s.energy, uni20::parse_real<Real>(ref.es), tol);
    EXPECT_REAL_NEAR(s.momentum, uni20::parse_real<Real>(ref.ps), tol);
    EXPECT_REAL_NEAR(*c.energy, uni20::parse_real<Real>(ref.ec), tol);
    EXPECT_REAL_NEAR(c.momentum, uni20::parse_real<Real>(ref.pc), tol);
    auto const inverted_spin = model::dispersion(model::Branch::spinon, ref.u, uni20::parse_real<Real>(ref.ps));
    auto const inverted_charge = model::dispersion(model::Branch::holon, ref.u, uni20::parse_real<Real>(ref.pc));
    ASSERT_TRUE(inverted_spin.converged);
    ASSERT_TRUE(inverted_charge.converged);
    EXPECT_REAL_NEAR(*inverted_spin.energy, uni20::parse_real<Real>(ref.es), Real{4} * tol);
    EXPECT_REAL_NEAR(*inverted_charge.energy, uni20::parse_real<Real>(ref.ec), Real{4} * tol);
  }
}

TYPED_TEST(HubbardThermo, QuantumNumbersAndEnergyConventions)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1});
  Real const tol = Real{512} * uni20::numeric_limits<Real>::epsilon();
  auto const h = model::charge_at_bare_momentum(model::Branch::holon, Real{4}, Real{1});
  auto const a = model::charge_at_bare_momentum(model::Branch::antiholon, Real{4}, Real{1});
  ASSERT_TRUE(h.converged);
  ASSERT_TRUE(a.converged);
  EXPECT_EQ(h.delta_particles, -1);
  EXPECT_EQ(a.delta_particles, 1);
  EXPECT_EQ(h.spin, uni20::half_int{0});
  EXPECT_EQ(a.spin, uni20::half_int{0});
  EXPECT_EQ(*h.energy, *a.energy);
  EXPECT_REAL_NEAR(h.momentum - a.momentum, pi, tol);
  for (auto branch : {model::Branch::spinon, model::Branch::holon, model::Branch::antiholon})
  {
    auto const sym = model::dispersion(branch, Real{4}, Real{1});
    auto const raw = model::dispersion(branch, Real{4}, Real{1}, model::Convention::unshifted);
    ASSERT_TRUE(sym.converged);
    ASSERT_TRUE(raw.converged);
    EXPECT_EQ(sym.symmetric_energy, raw.symmetric_energy);
    EXPECT_EQ(raw.energy_offset, Real{2} * Real(raw.delta_particles));
    EXPECT_REAL_NEAR(*raw.energy - *sym.energy, raw.energy_offset, tol);
    if (branch == model::Branch::spinon) EXPECT_EQ(sym.spin.twice(), 1);
  }
  // In the unshifted convention, removal need not cost positive energy.
  auto const removal = model::dispersion(model::Branch::holon, Real{4}, -pi / Real{2}, model::Convention::unshifted);
  ASSERT_TRUE(removal.converged);
  EXPECT_LT(*removal.energy, Real{0});
}

TYPED_TEST(HubbardThermo, WeakCouplingGapRetainsRelativePrecision)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1});
  auto const gap = model::charge_at_bare_momentum(model::Branch::holon, Real{.25}, pi);
  ASSERT_TRUE(gap.converged);
  // Independently summed K_1 series at 60 decimal digits (four terms suffice).
  auto const expected = uni20::parse_real<Real>("7.856419420831164720513491069438241188549399224515880816e-12");
  EXPECT_REAL_NEAR(*gap.energy / expected, Real{1}, Real{512} * uni20::numeric_limits<Real>::epsilon());
  EXPECT_LT(gap.energy_error, model::Options<Real>{}.relative_tolerance * expected);
}

TYPED_TEST(HubbardThermo, SpinonEndpointsReflectionAndStrongCoupling)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1});
  for (Real p : {Real{0}, pi})
  {
    auto const s = model::dispersion(model::Branch::spinon, Real{4}, p);
    ASSERT_TRUE(s.converged);
    EXPECT_EQ(*s.energy, Real{0});
    EXPECT_EQ(s.evaluations, 0u);
  }
  auto const a = model::spinon_at_rapidity(Real{4}, Real{1});
  auto const b = model::spinon_at_rapidity(Real{4}, Real{-1});
  ASSERT_TRUE(a.converged);
  ASSERT_TRUE(b.converged);
  EXPECT_EQ(*a.energy, *b.energy);
  EXPECT_REAL_NEAR(a.momentum + b.momentum, pi, Real{16} * uni20::numeric_limits<Real>::epsilon());
  Real const u = Real{10000}, p = pi / Real{3};
  auto const s = model::dispersion(model::Branch::spinon, u, p);
  ASSERT_TRUE(s.converged);
  EXPECT_REAL_NEAR(*s.energy * u, Real{2} * pi * std::sin(p), Real{1} / Real{100000});
}

TYPED_TEST(HubbardThermo, BudgetFailuresDoNotPublishEnergies)
{
  using Real = TypeParam;
  model::Options<Real> options;
  options.max_evaluations = 10;
  auto const limited = model::dispersion(model::Branch::holon, Real{4}, Real{0}, model::Convention::symmetric, options);
  EXPECT_FALSE(limited.converged);
  EXPECT_FALSE(limited.energy);
  EXPECT_FALSE(limited.symmetric_energy);
  EXPECT_EQ(limited.evaluations, 10u);
  EXPECT_EQ(limited.status, model::Status::quadrature_limit);
  options = {};
  options.max_iterations = 0;
  auto const iteration =
      model::dispersion(model::Branch::spinon, Real{4}, Real{1}, model::Convention::symmetric, options);
  EXPECT_FALSE(iteration.energy);
  EXPECT_EQ(iteration.status, model::Status::momentum_limit);
  EXPECT_THROW((void)model::dispersion(model::Branch::spinon, Real{0}, Real{1}), std::invalid_argument);
  EXPECT_THROW((void)model::dispersion(model::Branch::spinon, Real{-1}, Real{1}), std::invalid_argument);
  EXPECT_THROW((void)model::dispersion(model::Branch::spinon, Real{4}, Real{-1}), std::invalid_argument);
  EXPECT_THROW((void)model::dispersion(model::Branch::holon, Real{4}, Real{4}), std::invalid_argument);
  EXPECT_THROW((void)model::charge_at_bare_momentum(model::Branch::spinon, Real{4}, Real{1}), std::invalid_argument);
}
} // namespace
