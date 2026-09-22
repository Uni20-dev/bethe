// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "su_fermions_ed.hpp"
#include "test_support.hpp"
#include <bethe/gaudin_yang.hpp>
#include <bethe/lieb_liniger.hpp>
#include <bethe/su3.hpp>
#include <bethe/su_fermions.hpp>
#include <complex>

namespace
{
namespace model = bethe::su_fermions;
template <typename Real> class SUFermions : public ::testing::Test {};
TYPED_TEST_SUITE(SUFermions, test_support::RealTypes, test_support::PrecisionNames);

template <uni20::Real Real> void check_equations(model::State<Real> const& s)
{
  using C = std::complex<Real>;
  auto const& r = s.rapidities;
  Real const c = *s.reached_interaction, eps = uni20::numeric_limits<Real>::epsilon();
  auto ratio = [](Real d, Real w) { return C{d, w} / C{d, -w}; };
  for (std::size_t a = 0; a < r.size(); ++a)
    for (std::size_t j = 0; j < r[a].size(); ++j)
    {
      EXPECT_EQ(r[a][j], -r[a][r[a].size() - 1 - j]);
      if (j) EXPECT_GT(r[a][j], r[a][j - 1]);
      C lhs{Real{1}, Real{0}}, rhs{Real{1}, Real{0}};
      if (a == 0)
      {
        lhs = std::exp(C{Real{0}, r[0][j] * s.length});
        for (Real v : r[1])
          rhs *= ratio(r[0][j] - v, c / Real{2});
      }
      else
      {
        for (Real v : r[a - 1])
          lhs *= ratio(r[a][j] - v, c / Real{2});
        if (a + 1 < r.size())
          for (Real v : r[a + 1])
            lhs *= ratio(r[a][j] - v, c / Real{2});
        for (std::size_t k = 0; k < r[a].size(); ++k)
          if (k != j) rhs *= ratio(r[a][j] - r[a][k], c);
      }
      EXPECT_REAL_NEAR(lhs.real(), rhs.real(), Real{2048} * Real(s.particles) * eps);
      EXPECT_REAL_NEAR(lhs.imag(), rhs.imag(), Real{2048} * Real(s.particles) * eps);
    }
}

TYPED_TEST(SUFermions, NestedCouplingSweepAndOriginalEquations)
{
  using Real = TypeParam;
  for (auto const& pop :
       {std::vector<std::size_t>{1, 1, 1}, std::vector<std::size_t>{3, 1, 1}, std::vector<std::size_t>{5, 3, 1},
        std::vector<std::size_t>{3, 3, 3}, std::vector<std::size_t>{1, 1, 1, 1}, std::vector<std::size_t>{3, 1, 1, 1},
        std::vector<std::size_t>{3, 3, 3, 3}, std::vector<std::size_t>{1, 1, 1, 1, 1, 1}})
    for (Real c : {Real{1} / Real{1000}, Real{1}, Real{100}, Real{10000}})
    {
      SCOPED_TRACE(::testing::Message() << pop.size() << "," << pop.front() << "," << double(c));
      auto const s = model::ground_state<Real>(pop, Real{1}, c);
      ASSERT_TRUE(s.converged) << int(s.status) << " iterations=" << s.iterations
                               << " residual=" << double(s.target_residual_norm);
      EXPECT_EQ(*s.reached_interaction, c);
      EXPECT_EQ(s.momentum, Real{0});
      ASSERT_NO_FATAL_FAILURE(check_equations(s));
    }
}

TYPED_TEST(SUFermions, TwoComponentReduction)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (auto const& pop :
       {std::vector<std::size_t>{1, 1}, std::vector<std::size_t>{3, 3}, std::vector<std::size_t>{3, 5}})
    for (Real c : {Real{1} / Real{100}, Real{1}, Real{100}})
    {
      auto const s = model::ground_state<Real>(pop, Real{2}, c);
      auto const old = bethe::gaudin_yang::ground_state(pop[0], pop[1], Real{2}, c);
      ASSERT_TRUE(s.converged && old.converged);
      EXPECT_REAL_NEAR(*s.energy, old.energy, Real{8192} * eps * (Real{1} + old.energy));
      for (std::size_t j = 0; j < old.momenta.size(); ++j)
        EXPECT_REAL_NEAR(s.rapidities[0][j], old.momenta[j], Real{8192} * eps);
    }
  // Independent two-body jump condition: c=pi, ell=1 gives k=+/-pi/2.
  Real const pi = Real{4} * std::atan(Real{1}), exact = pi * pi / Real{2};
  auto const analytic = model::ground_state<Real>(std::vector<std::size_t>{1, 1, 0}, Real{1}, pi);
  ASSERT_TRUE(analytic.converged);
  EXPECT_REAL_NEAR(*analytic.energy, exact, Real{512} * eps);
  if constexpr (uni20::numeric_limits<Real>::digits > 53)
    EXPECT_GT(std::abs(Real(double(exact)) - exact), Real{512} * eps);
}

TYPED_TEST(SUFermions, AntisymmetricSpinSingletEqualsBosons)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t colors : {3, 4, 5, 6})
    for (Real c : {Real{1} / Real{1000}, Real{1}, Real{1000}})
    {
      auto const s = model::ground_state<Real>(std::vector<std::size_t>(colors, 1), Real{1}, c);
      auto const bosons = bethe::lieb_liniger::ground_state(colors, Real{1}, c);
      ASSERT_TRUE(s.converged && bosons.converged);
      EXPECT_REAL_NEAR(*s.energy, bosons.energy, Real{8192} * eps * bosons.energy);
    }
  Real const tiny = uni20::parse_real<Real>("1e-40");
  auto const s = model::ground_state<Real>(std::vector<std::size_t>{1, 1, 1}, Real{1}, tiny);
  ASSERT_TRUE(s.converged) << int(s.status) << " iterations=" << s.iterations;
  EXPECT_REAL_NEAR(*s.energy / (Real{6} * tiny), Real{1}, Real{8192} * eps);
}

TYPED_TEST(SUFermions, FreeLimitsAndComponentPermutation)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (auto const& pop : {std::vector<std::size_t>{0}, std::vector<std::size_t>{0, 0, 0},
                          std::vector<std::size_t>{2, 1, 4}, std::vector<std::size_t>{0, 6, 0}})
  {
    auto const s = model::ground_state<Real>(pop, Real{2}, Real{0});
    ASSERT_TRUE(s.converged && s.free);
    Real exact{};
    for (auto n : pop)
      exact += pi * pi * Real(n) * (Real(n) * Real(n) + (n % 2 ? Real{-1} : Real{2})) / Real{12};
    EXPECT_REAL_NEAR(*s.energy, exact, Real{256} * eps * (Real{1} + exact));
    EXPECT_TRUE(s.quantum_numbers.empty());
    EXPECT_EQ(s.free_modes.size(), pop.size());
  }
  auto const polarized = model::ground_state<Real>(std::vector<std::size_t>{0, 6, 0}, Real{2}, Real{10});
  EXPECT_TRUE(polarized.free && polarized.converged);
  auto const s = model::ground_state<Real>(std::vector<std::size_t>{3, 1, 1}, Real{1}, Real{1});
  auto const perm = model::ground_state<Real>(std::vector<std::size_t>{1, 0, 3, 1}, Real{1}, Real{1});
  ASSERT_TRUE(s.converged && perm.converged);
  EXPECT_EQ(*s.energy, *perm.energy);
  EXPECT_EQ(perm.component_order, (std::vector<std::size_t>{2, 0, 3}));
  auto const scaled = model::ground_state<Real>(std::vector<std::size_t>{3, 1, 1}, Real{2}, Real{0.5});
  ASSERT_TRUE(scaled.converged);
  EXPECT_REAL_NEAR(*s.energy, Real{4} * (*scaled.energy), Real{256} * eps * (*s.energy));
}

TYPED_TEST(SUFermions, WeakPerturbationAndStrongSpinChain)
{
  using Real = TypeParam;
  std::vector<std::size_t> pop{5, 3, 1};
  Real const c = Real{1} / Real{100000};
  auto const weak = model::ground_state<Real>(pop, Real{1}, c);
  auto const free = model::ground_state<Real>(pop, Real{1}, Real{0});
  ASSERT_TRUE(weak.converged);
  EXPECT_REAL_NEAR((*weak.energy - *free.energy) / (Real{2} * c), Real{23}, Real{1} / Real{1000});
  Real const strong_c = Real{1000000}, pi = Real{4} * std::atan(Real{1});
  for (std::size_t n : {3, 9})
  {
    auto const s = model::ground_state<Real>(std::vector<std::size_t>(3, n / 3), Real{1}, strong_c);
    auto const spin = bethe::su3::ground_state<Real>(n);
    ASSERT_TRUE(s.converged && spin.converged);
    Real const infinity = pi * pi * Real(n) * Real(n * n - 1) / Real{3};
    Real const omega = Real(n) - spin.energy;
    EXPECT_REAL_NEAR((Real{1} - *s.energy / infinity) * strong_c / Real{2}, omega, Real{1} / Real{500});
    for (std::size_t a = 0; a < 2; ++a)
      for (std::size_t j = 0; j < spin.rapidities[a].size(); ++j)
        EXPECT_REAL_NEAR(s.rapidities[a + 1][j] / strong_c, spin.rapidities[a][j], Real{1} / Real{100000});
  }
}

TYPED_TEST(SUFermions, AnalyticJacobian)
{
  using Real = TypeParam;
  model::detail::GroundSystem<Real> const system({9, 6, 3});
  auto const x = system.seed();
  Real const eps = uni20::numeric_limits<Real>::epsilon(), h = std::cbrt(eps);
  for (Real g : {Real{1} / Real{5}, Real{1}, Real{20}})
  {
    std::vector<Real> jac;
    (void)system.evaluate(x, g, &jac);
    auto const roots = system.expand(x, g);
    auto raw = [&](std::vector<Real> const& values, std::size_t a, std::size_t j) {
      auto const r = system.expand(values, g);
      auto const full = system.counts[a] - system.counts[a] / 2 + j;
      Real const value = r[a][full];
      Real f = (a ? Real{0} : value) - system.pi * Real(system.labels[a][full].twice());
      for (std::size_t b = a ? a - 1 : 1; b < std::min(r.size(), a + 2); ++b)
        for (std::size_t k = 0; k < r[b].size(); ++k)
          if (!(a == b && k == full))
            f += (a == b ? Real{-2} : Real{2}) * std::atan((value - r[b][k]) / (a == b ? g : g / Real{2}));
      return f;
    };
    for (std::size_t a = 0; a < 3; ++a)
      for (std::size_t j = 0; j < system.counts[a] / 2; ++j)
      {
        auto const row = system.offsets[a] + j, full = system.counts[a] - system.counts[a] / 2 + j;
        Real const normalization = g < Real{1} ? std::max(std::sqrt(g), std::abs(roots[a][full])) : Real{9};
        for (std::size_t col = 0; col < x.size(); ++col)
        {
          auto plus = x, minus = x;
          plus[col] += h;
          minus[col] -= h;
          EXPECT_REAL_NEAR((raw(plus, a, j) - raw(minus, a, j)) / (Real{2} * h * normalization),
                           jac[row * x.size() + col], Real{16384} * h * h);
        }
      }
  }
}

TYPED_TEST(SUFermions, BudgetsAndInputValidation)
{
  using Real = TypeParam;
  std::vector<std::size_t> pop{3, 1, 1};
  model::SolverOptions<Real> options;
  options.max_iterations = 0;
  auto const seed = model::ground_state<Real>(pop, Real{1}, Real{1}, options);
  EXPECT_FALSE(seed.converged);
  EXPECT_FALSE(seed.energy);
  EXPECT_FALSE(seed.reached_interaction);
  EXPECT_TRUE(seed.rapidities.empty());
  EXPECT_EQ(seed.iterations, 0U);
  EXPECT_EQ(seed.stages, 0U);
  auto const free = model::ground_state<Real>(pop, Real{1}, Real{0}, options);
  EXPECT_TRUE(free.free && free.converged);
  options.max_iterations = 1;
  auto const exhausted = model::ground_state<Real>(pop, Real{1}, Real{1}, options);
  EXPECT_FALSE(exhausted.converged);
  EXPECT_EQ(exhausted.iterations, 1U);
  options.max_iterations = 10000;
  options.max_stages = 1;
  auto const partial = model::ground_state<Real>(pop, Real{1}, Real{1}, options);
  EXPECT_FALSE(partial.converged);
  ASSERT_TRUE(partial.energy);
  ASSERT_TRUE(partial.reached_interaction);
  EXPECT_EQ(partial.status, model::SolveStatus::stage_limit);
  EXPECT_GT(*partial.reached_interaction, Real{1});
  EXPECT_GT(partial.target_residual_norm, Real{1} / Real{100});
  auto const reference = model::ground_state<Real>(pop, Real{1}, *partial.reached_interaction);
  ASSERT_TRUE(reference.converged);
  EXPECT_EQ(*partial.energy, *reference.energy);
  ASSERT_NO_FATAL_FAILURE(check_equations(partial));
  options.max_stages = 0;
  auto const no_stage = model::ground_state<Real>(pop, Real{1}, Real{1}, options);
  EXPECT_EQ(no_stage.status, model::SolveStatus::stage_limit);
  EXPECT_FALSE(no_stage.energy);
  EXPECT_THROW((void)model::ground_state<Real>({}, Real{1}, Real{1}), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state<Real>(std::vector<std::size_t>{2, 1, 1}, Real{1}, Real{1}),
               std::invalid_argument);
  EXPECT_THROW((void)model::ground_state<Real>(std::vector<std::size_t>{2, 2, 2}, Real{1}, Real{1}),
               std::invalid_argument);
  EXPECT_THROW((void)model::ground_state<Real>(pop, Real{0}, Real{1}), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state<Real>(pop, Real{1}, Real{-1}), std::invalid_argument);
  Real const nan = uni20::numeric_limits<Real>::quiet_NaN();
  EXPECT_THROW((void)model::ground_state<Real>(pop, Real{1}, nan), std::invalid_argument);
  options.residual_tolerance = Real{0};
  EXPECT_THROW((void)model::ground_state<Real>(pop, Real{1}, Real{1}, options), std::invalid_argument);
  std::vector<std::size_t> huge{std::size_t(std::numeric_limits<std::int64_t>::max()), 1, 1};
  EXPECT_THROW((void)model::ground_state<Real>(huge, Real{1}, Real{1}), std::invalid_argument);
  EXPECT_THROW((void)model::ground_quantum_numbers(std::vector<std::size_t>{0, 1, 0}), std::invalid_argument);
}

TYPED_TEST(SUFermions, LabelsAndLargerSixComponentSystem)
{
  using Real = TypeParam;
  auto const labels = model::ground_quantum_numbers(std::vector<std::size_t>{1, 0, 3, 1});
  ASSERT_EQ(labels.size(), 3U);
  EXPECT_EQ(labels[0].size(), 5U);
  EXPECT_EQ(labels[1].size(), 2U);
  EXPECT_EQ(labels[2].size(), 1U);
  EXPECT_EQ(labels[0].front().twice(), -4);
  EXPECT_EQ(labels[1].front().twice(), -1);
  EXPECT_EQ(labels[2].front().twice(), 0);
  for (Real c : {Real{1} / Real{1000}, Real{1}, Real{1000}})
  {
    auto const s = model::ground_state<Real>(std::vector<std::size_t>(6, 5), Real{30}, c);
    ASSERT_TRUE(s.converged) << int(s.status) << " iterations=" << s.iterations;
    EXPECT_EQ(s.rapidities.size(), 6U);
    EXPECT_EQ(s.particles, 30U);
    ASSERT_NO_FATAL_FAILURE(check_equations(s));
  }
}

TEST(SUFermionsHamiltonian, PlaneWaveVariationalConvergence)
{
  std::vector<std::size_t> pop{3, 1, 1};
  auto const s = model::ground_state<double>(pop, 1., 1.);
  ASSERT_TRUE(s.converged);
  std::vector<unsigned> counts{3, 1, 1};
  EXPECT_NEAR(bethe::test::su_fermions_cutoff_ground(counts, 2, 1., 0.), 8 * std::numbers::pi * std::numbers::pi,
              1e-10);
  double previous = 1000, first = 0;
  for (int cutoff : {2, 3, 4})
  {
    double const e = bethe::test::su_fermions_cutoff_ground(counts, cutoff, 1., 1.);
    EXPECT_GT(e, *s.energy);
    EXPECT_LT(e, previous);
    if (cutoff == 2) first = e - *s.energy;
    if (cutoff == 4)
    {
      EXPECT_LT(e - *s.energy, first * .6);
      EXPECT_LT(e - *s.energy, .4);
    }
    previous = e;
  }
}
} // namespace
