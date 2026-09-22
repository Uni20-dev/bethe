// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "hubbard_ed.hpp"
#include "test_support.hpp"
#include <array>
#include <bethe/heisenberg_open.hpp>
#include <bethe/hubbard_open.hpp>
#include <string>
#include <string_view>

namespace
{
namespace model = bethe::hubbard::open;
template <typename T>
concept HasMomentum = requires(T value) { value.momentum; };
template <typename T>
concept HasMomentumIndex = requires(T value) { value.momentum_index; };
static_assert(!HasMomentum<model::State<double>> && !HasMomentumIndex<model::State<double>>);

// Independently written full equations, using returned conventional Lambda.
template <uni20::Real Real> void check_state(model::State<Real> const& state)
{
  using std::abs;
  using std::atan;
  using std::cos;
  using std::sin;
  Real const pi = Real{4} * atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  auto const n = state.root_particles, m = state.root_down_spins;
  auto const& k = state.charge_momenta;
  auto const& lambda = state.spin_rapidities;
  ASSERT_TRUE(k.size() == n && state.converged == (state.status == model::SolveStatus::converged))
      << "root count/status";
  Real const allowance = Real{128} * Real(state.sites) * eps;
  Real energy = Real{0};
  for (Real value : k)
    energy -= Real{2} * cos(value);
  ASSERT_TRUE(abs(energy + state.energy_offset - state.energy) <= allowance * (Real{1} + abs(state.energy_offset)))
      << "OBC energy from returned roots and exact mapping offset";
  if (state.free_fermion)
  {
    ASSERT_TRUE((state.root_interaction == Real{0} || m == 0 || m == n) && lambda.empty() &&
                state.quantum_numbers.charge.empty() && state.quantum_numbers.spin.empty())
        << "free metadata";
    ASSERT_TRUE(state.converged && state.residual_norm == Real{0} && state.iterations == 0 &&
                state.continuation_steps == 0)
        << "exact free path does not iterate";
    return;
  }
  ASSERT_TRUE(lambda.size() == m && state.quantum_numbers.charge.size() == n && state.quantum_numbers.spin.size() == m)
      << "OBC nested labels and rapidities";
  Real const u = state.root_interaction / Real{4}, normalization = Real{2} * Real(state.sites + 1);
  Real charge = Real{0}, spin = Real{0};
  for (std::size_t j = 0; j < n; ++j)
  {
    SCOPED_TRACE(::testing::Message() << " j=" << j);
    ASSERT_TRUE(k[j] > Real{0} && k[j] < pi && (j == 0 || k[j] > k[j - 1])) << "positive ordered charge roots";
    ASSERT_TRUE(state.quantum_numbers.charge[j] == uni20::half_int(static_cast<std::int64_t>(j + 1)))
        << "charge labels";
    Real f = normalization * k[j] - Real{2} * pi * Real(j + 1);
    for (Real l : lambda)
      f += Real{2} * (atan((sin(k[j]) - l) / u) + atan((sin(k[j]) + l) / u));
    charge = std::max(charge, abs(f) / normalization);
  }
  for (std::size_t a = 0; a < m; ++a)
  {
    SCOPED_TRACE(::testing::Message() << " a=" << a);
    ASSERT_TRUE(lambda[a] > Real{0} && (a == 0 || lambda[a] > lambda[a - 1])) << "positive ordered spin roots";
    ASSERT_TRUE(state.quantum_numbers.spin[a] == uni20::half_int(static_cast<std::int64_t>(a + 1))) << "spin labels";
    Real f = -Real{2} * pi * Real(a + 1);
    for (Real value : k)
      f += Real{2} * (atan((lambda[a] - sin(value)) / u) + atan((lambda[a] + sin(value)) / u));
    for (std::size_t b = 0; b < m; ++b)
      if (b != a)
        f -= Real{2} * (atan((lambda[a] - lambda[b]) / (Real{2} * u)) + atan((lambda[a] + lambda[b]) / (Real{2} * u)));
    spin = std::max(spin, abs(f) / normalization);
  }

  ASSERT_REAL_NEAR(charge, state.charge_residual, allowance) << "OBC full-equation residuals at final auxiliary U";
  ASSERT_REAL_NEAR(spin, state.spin_residual, allowance) << "OBC full-equation residuals at final auxiliary U";

  ASSERT_TRUE(state.residual_norm == std::max(state.charge_residual, state.spin_residual)) << "residual max norm";
}

template <typename Real> class HubbardOpen : public ::testing::Test {};
TYPED_TEST_SUITE(HubbardOpen, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(HubbardOpen, Jacobian)
{
  using Real = TypeParam;

  using std::abs;
  using std::cbrt;
  Real const step = cbrt(uni20::numeric_limits<Real>::epsilon());
  for (auto [sites, n, m] : {std::array<std::size_t, 3>{2, 2, 1}, {3, 3, 1}, {4, 4, 2}, {5, 4, 1}, {6, 5, 2}})
    for (Real u : {Real{1} / Real{4}, Real{1}, Real{4}})
    {
      SCOPED_TRACE(::testing::Message() << " sites=" << sites << " n=" << n << " m=" << m);
      SCOPED_TRACE(::testing::Message() << " u=" << uni20::format_scalar(u));
      model::detail::GroundSystem<Real> system(sites, n, m);
      auto x = system.seed();
      for (Real& value : x)
        value *= Real{9} / Real{10};
      uni20::DenseMatrix<Real> a(system.order, system.order);
      (void)system.evaluate(x, u, &a);
      for (std::size_t col = 0; col < system.order; ++col)
      {
        SCOPED_TRACE(::testing::Message() << " col=" << col);
        auto plus = x, minus = x;
        plus[col] += step;
        minus[col] -= step;
        auto const fp = system.evaluate(plus, u), fm = system.evaluate(minus, u);
        for (std::size_t row = 0; row < system.order; ++row)
        {
          SCOPED_TRACE(::testing::Message() << " row=" << row);
          Real const numerical = (fp.residual[row] - fm.residual[row]) / (Real{2} * step);
          EXPECT_REAL_NEAR((a[row, col]), numerical, Real{1000} * step * step * (Real{1} + abs(numerical)))
              << "open Hubbard analytic Jacobian vs central differences";
        }
      }
    }
}

TYPED_TEST(HubbardOpen, SectorEnergiesAgainstED)
{
  using Real = TypeParam;

  using std::abs;
  // Every physical spin population, not just balanced or half-filled sectors.
  for (unsigned sites : {1, 2, 3, 4, 5, 6})
    for (unsigned up = 0; up <= sites; ++up)
      for (unsigned down = 0; down <= sites; ++down)
        for (Real interaction :
             {Real{0}, Real{1} / Real{10}, Real{4}, Real{100}, Real{-1} / Real{10}, Real{-4}, Real{-100}})
        {
          SCOPED_TRACE(::testing::Message() << " sites=" << sites);
          SCOPED_TRACE(::testing::Message() << " up=" << up);
          SCOPED_TRACE(::testing::Message() << " down=" << down);
          SCOPED_TRACE(::testing::Message() << " interaction=" << uni20::format_scalar(interaction));
          auto const sz = uni20::from_twice(static_cast<std::int64_t>(up) - down);
          auto const state = model::sector_ground_state<Real>(sites, up + down, sz, interaction);
          auto const ed = bethe::test::exact_ground(sites, static_cast<double>(interaction), up, down, false);
          ASSERT_TRUE(state.converged) << "open Hubbard convergence vs independent fermionic ED";
          EXPECT_LE(std::abs(static_cast<double>(state.energy) - ed.energy), 3e-10)
              << "open Hubbard energy vs independent fermionic ED";
          ASSERT_TRUE(state.particles == up + down && state.down_spins == down && state.interaction == interaction)
              << "OBC physical metadata survives symmetry mappings";
          ASSERT_NO_FATAL_FAILURE(check_state(state));
        }
}

TYPED_TEST(HubbardOpen, DimerPrecision)
{
  using Real = TypeParam;

  using std::sqrt;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real u :
       {Real{0}, Real{1} / Real{100}, Real{4}, Real{100}, uni20::parse_real<Real>("4.12345678901234567890123456789")})
  {
    SCOPED_TRACE(::testing::Message() << " u=" << uni20::format_scalar(u));
    auto const positive = model::ground_state<Real>(2, u, {.residual_tolerance = Real{2} * eps});
    auto const negative = model::ground_state<Real>(2, -u, {.residual_tolerance = Real{2} * eps});
    Real const exact = -Real{8} / (u + sqrt(u * u + Real{16}));

    ASSERT_TRUE(positive.converged) << "native-precision ordinary dimer";
    ASSERT_TRUE(negative.converged) << "native-precision ordinary dimer";
    ASSERT_REAL_NEAR(positive.energy, exact, Real{64} * eps) << "native-precision ordinary dimer";
    ASSERT_REAL_NEAR(negative.energy, (exact - u), Real{64} * eps * (Real{1} + u)) << "native-precision ordinary dimer";

    test_support::expect_exact(uni20::parse_real<Real>(uni20::format_real(positive.energy)), positive.energy,
                               "energy scalar I/O");
  }
}

TYPED_TEST(HubbardOpen, LargeSystems)
{
  using Real = TypeParam;

  using std::abs;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  // Larger systems exercise all filling parities, both signs and asymmetric spin populations.
  for (auto [sites, particles, down] : {std::array<std::size_t, 3>{7, 7, 3},
                                        {8, 4, 2},
                                        {15, 8, 3},
                                        {16, 9, 4},
                                        {16, 24, 12},
                                        {32, 32, 16},
                                        {31, 18, 7}})
    for (Real interaction : {Real{0}, Real{1} / Real{1000000}, Real{4}, Real{100}, Real{-4}})
    {
      SCOPED_TRACE(::testing::Message() << " sites=" << sites << " particles=" << particles << " down=" << down);
      SCOPED_TRACE(::testing::Message() << " interaction=" << uni20::format_scalar(interaction));
      auto const sz = uni20::from_twice(static_cast<std::int64_t>(particles) - 2 * static_cast<std::int64_t>(down));
      auto const state = model::sector_ground_state<Real>(sites, particles, sz, interaction);
      ASSERT_TRUE(state.converged && state.residual_norm <= Real{32} * eps)
          << "larger OBC native-precision convergence";
      ASSERT_NO_FATAL_FAILURE(check_state(state));
      if (interaction > Real{0} && interaction < Real{1})
      {
        auto const free = model::sector_ground_state<Real>(sites, particles, sz, Real{0});
        ASSERT_TRUE(state.energy > free.energy && state.energy - free.energy < Real(sites) * interaction)
            << "OBC weak repulsion approaches free-fermion energy from above";
      }
      std::size_t stages = interaction == Real{0} || state.free_fermion ? 0 : 1;
      for (Real current = Real{8}; interaction != Real{0} && !state.free_fermion && current > abs(interaction);
           current /= Real{2})
        ++stages;
      ASSERT_TRUE(state.continuation_steps == stages) << "all continuation stages reach the requested magnitude of U";
    }
}

TYPED_TEST(HubbardOpen, StrongCoupling)
{
  using Real = TypeParam;

  for (std::size_t sites : {2, 3, 4, 5, 8, 15})
    for (std::size_t down = 0; down <= sites / 2; ++down)
    {
      SCOPED_TRACE(::testing::Message() << " sites=" << sites);
      SCOPED_TRACE(::testing::Message() << " down=" << down);
      auto const sz = uni20::from_twice(static_cast<std::int64_t>(sites - 2 * down));
      auto const xxx = bethe::heisenberg::open::sector_ground_state<Real>(sites, sz);
      auto const state = model::sector_ground_state<Real>(sites, sites, sz, Real{10000});

      ASSERT_TRUE(xxx.converged) << "open Hubbard strong coupling has L-1 bonds and the correct XXX offset";
      ASSERT_TRUE(state.converged) << "open Hubbard strong coupling has L-1 bonds and the correct XXX offset";
      ASSERT_REAL_NEAR(Real{2500} * state.energy, (xxx.energy - Real(sites - 1) / Real{4}), Real(sites) / Real{1000000})
          << "open Hubbard strong coupling has L-1 bonds and the correct XXX offset";
    }
}

TYPED_TEST(HubbardOpen, FreeSpectrumAndParity)
{
  using Real = TypeParam;

  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t sites : {1, 3, 4, 9})
  {
    SCOPED_TRACE(::testing::Message() << " sites=" << sites);
    auto const ground = model::ground_state<Real>(sites, Real{4});
    ASSERT_TRUE(ground.particles == sites && ground.down_spins == sites / 2)
        << "half-filled default chooses Sz by parity";
    for (std::size_t count = 0; count <= sites; ++count)
    {
      SCOPED_TRACE(::testing::Message() << " count=" << count);
      auto const free =
          model::sector_ground_state<Real>(sites, 2 * count, uni20::half_int{0}, Real{0}, {.max_iterations = 0});
      Real reference = Real{0}, pi = Real{4} * std::atan(Real{1});
      for (std::size_t j = 1; j <= count; ++j)
        reference -= Real{4} * std::cos(pi * Real(j) / Real(sites + 1));
      EXPECT_REAL_NEAR(free.energy, reference, Real{32} * Real(sites) * eps) << "native free standing-wave spectrum";
      ASSERT_NO_FATAL_FAILURE(check_state(free));
    }
  }
}

TYPED_TEST(HubbardOpen, FailureBudgets)
{
  using Real = TypeParam;

  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real interaction : {Real{1} / Real{100}, -Real{1} / Real{100}})
    for (std::size_t budget : {0, 1, 5})
    {
      SCOPED_TRACE(::testing::Message() << " interaction=" << uni20::format_scalar(interaction));
      SCOPED_TRACE(::testing::Message() << " budget=" << budget);
      auto const state = model::sector_ground_state<Real>(7, 5, uni20::half_int::parse("-1/2"), interaction,
                                                          {.max_iterations = budget});
      ASSERT_TRUE(!state.converged && state.status == model::SolveStatus::iteration_limit && state.iterations == budget)
          << "global OBC budget is honored including mapped estimates";
      ASSERT_NO_FATAL_FAILURE(check_state(state));
    }
  auto const weak = model::ground_state<Real>(6, eps * eps * eps * eps);
  ASSERT_NO_FATAL_FAILURE(check_state(weak));
  if (weak.converged)
  {
    EXPECT_REAL_NEAR(weak.energy, model::ground_state<Real>(6, Real{0}).energy, Real{4096} * eps)
        << "unresolved weak coupling retains the free energy to native precision";
  }
  auto const stalled = model::ground_state<Real>(6, Real{4}, {.residual_tolerance = eps * eps});
  ASSERT_TRUE(!stalled.converged && stalled.status == model::SolveStatus::stalled)
      << "unattainable OBC tolerance stalls explicitly";
  ASSERT_NO_FATAL_FAILURE(check_state(stalled));
}

TYPED_TEST(HubbardOpen, InvalidInputs)
{
  using Real = TypeParam;

  model::detail::GroundSystem<Real> overflow_system(4, 4, 2);
  auto overflow_roots = overflow_system.seed();
  overflow_roots[overflow_system.nk] = (uni20::numeric_limits<Real>::max() / Real{4}) * Real{3};
  EXPECT_THROW(([&] { (void)overflow_system.evaluate(overflow_roots, Real{1}); })(), std::overflow_error);
  EXPECT_THROW(([&] { (void)model::ground_state<Real>(0, Real{4}); })(), std::invalid_argument);
  EXPECT_THROW(([&] { (void)model::ground_quantum_numbers(std::numeric_limits<std::size_t>::max()); })(),
               std::invalid_argument);
  EXPECT_THROW(([&] { (void)model::ground_quantum_numbers(4, 5, 1); })(), std::invalid_argument);
  EXPECT_THROW(([&] { (void)model::ground_quantum_numbers(4, 4, 3); })(), std::invalid_argument);
  EXPECT_THROW(([&] {
                 (void)model::ground_state<Real>(
                     static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() / 4) - 1, Real{4});
               })(),
               std::invalid_argument);
  for (Real value : {uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW(([&] { (void)model::ground_state<Real>(2, value); })(), std::invalid_argument);
  for (Real value :
       {Real{0}, Real{-1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW(([&] { (void)model::ground_state<Real>(2, Real{0}, {.residual_tolerance = value}); })(),
                 std::invalid_argument);
  EXPECT_THROW(([&] { (void)model::ground_state<Real>(2, uni20::numeric_limits<Real>::min() / Real{16}); })(),
               std::invalid_argument);
  EXPECT_THROW(([&] { (void)model::sector_ground_state<Real>(4, 9, uni20::half_int{0}, Real{4}); })(),
               std::invalid_argument);
  EXPECT_THROW(([&] { (void)model::sector_ground_state<Real>(4, 3, uni20::half_int{0}, Real{4}); })(),
               std::invalid_argument);
  EXPECT_THROW(([&] { (void)model::sector_ground_state<Real>(4, 4, uni20::half_int{3}, Real{4}); })(),
               std::invalid_argument);
  EXPECT_THROW(([&] {
                 (void)model::sector_ground_state<Real>(
                     4, 4, uni20::from_twice(std::numeric_limits<std::int64_t>::min()), Real{4});
               })(),
               std::invalid_argument);
  EXPECT_THROW(([&] {
                 (void)model::sector_ground_state<Real>(3, 6, uni20::half_int{0}, -uni20::numeric_limits<Real>::max());
               })(),
               std::overflow_error);
}

} // namespace
