// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include "test_support.hpp"
#include <bethe/heisenberg_excitations.hpp>
#include <bethe/xxz_excitations.hpp>

#include "exact_spectrum.hpp"

#include <set>
#include <string>
#include <string_view>

namespace
{
namespace model = bethe::xxz;
using uni20::half_int;

// Independent check in the conventional hyperbolic rapidity, not the rational
// scattering formula used in the implementation. Energy from magnon momenta.
template <uni20::Real Real> void check_state(std::size_t n, model::RealState<Real> const& state)
{
  using std::abs;
  using std::acos;
  using std::atan;
  using std::atanh;
  using std::cos;
  using std::sin;
  using std::tan;
  using std::tanh;
  Real const pi = Real{4} * atan(Real{1});
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const d = state.delta;
  Real const gamma = d > Real{1} ? Real{0} : acos(d);
  Real const scale = d > Real{1} ? std::sqrt((d - Real{1}) / (d + Real{1})) : tan(gamma / Real{2});
  Real residual = Real{0};
  Real energy = Real(n) * d / Real{4};
  Real momentum = Real{0};
  auto const& z = state.rapidities;
  ASSERT_TRUE(z.size() == state.quantum_numbers.size()) << "XXZ root/label count";
  for (std::size_t i = 0; i < z.size(); ++i)
  {
    SCOPED_TRACE(::testing::Message() << " i=" << i);
    Real const theta = Real{2} * atan(z[i]);
    Real f = Real(n) * theta - pi * Real(state.quantum_numbers[i].twice());
    if (d != Real{0})
      for (std::size_t j = 0; j < z.size(); ++j)
        if (i != j)
        {
          if (d > Real{1})
          {
            Real const difference = atan(scale * z[i]) - atan(scale * z[j]);
            Real const t = Real{2} * scale / (Real{1} + scale * scale);
            f -= Real{2} * std::atan2(sin(difference), t * cos(difference));
          }
          else if (d == Real{1})
            f -= Real{2} * atan((z[i] - z[j]) / Real{2});
          else
          {
            Real const difference = atanh(scale * z[i]) - atanh(scale * z[j]);
            f -= Real{2} * atan(tanh(difference) / tan(gamma));
          }
        }
    residual = std::max(residual, abs(f) / Real(n));
    energy += cos(pi - theta) - d;
    momentum += pi - theta;
    ASSERT_TRUE(uni20::isfinite(z[i])) << "finite XXZ rapidities";
    if (d < Real{1})
    {
      EXPECT_REAL_NEAR(scale * z[i], Real{0}, Real{1}) << "XXZ roots stay on finite real branch";
    }
    if (state.converged && i > 0)
    {
      ASSERT_TRUE(z[i - 1] < z[i]) << "ordered XXZ roots";
    }
  }
  Real const allowance = Real{128} * Real(n) * eps;
  EXPECT_REAL_NEAR(energy, state.energy, allowance) << "XXZ returned energy vs roots";
  EXPECT_REAL_NEAR(residual, state.residual_norm, allowance) << "XXZ independent logarithmic residual";
  ASSERT_TRUE(state.momentum_index < n) << "XXZ momentum index range";
  if (state.converged)
  {
    Real const phase_allowance = allowance + Real(n) * state.residual_norm;

    ASSERT_REAL_NEAR(cos(momentum), cos(state.momentum), phase_allowance) << "XXZ momentum from roots";
    ASSERT_REAL_NEAR(sin(momentum), sin(state.momentum), phase_allowance) << "XXZ momentum from roots";
  }
}

template <typename Real> class XXZExcitations : public ::testing::Test {};
TYPED_TEST_SUITE(XXZExcitations, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXZExcitations, EnergiesAgainstED)
{
  using Real = TypeParam;

  using std::abs;
  using std::atan;
  using std::cos;
  using Numbers = model::QuantumNumbers;
  auto half = [](std::int64_t twice) { return uni20::from_twice(twice); };
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const pi = Real{4} * atan(Real{1});
  auto const all_count = std::numeric_limits<std::size_t>::max();
  for (Real d : {Real{0}, Real{1} / Real{10}, Real{1} / Real{2}, Real{9} / Real{10}, Real{1}})
    for (unsigned n = 2; n <= 9; ++n)
      for (unsigned m = 0; m <= n / 2; ++m)
      {
        SCOPED_TRACE(::testing::Message() << " d=" << uni20::format_scalar(d));
        SCOPED_TRACE(::testing::Message() << " n=" << n);
        SCOPED_TRACE(::testing::Message() << " m=" << m);
        auto const sz = half(n - 2 * m);
        auto const scan = model::real_excitations<Real>(n, d, sz, {.count = all_count});
        ASSERT_TRUE(scan.converged() && !scan.first_unconverged) << "XXZ excitation family convergence";
        ASSERT_TRUE(scan.sz == sz && scan.delta == d && scan.candidate_count == scan.levels.size() &&
                    scan.candidate_count == model::real_excitation_count(n, d, sz))
            << "XXZ scan accounting";
        auto ed = test_support::exact_spectrum(n, m, 0.371, true, static_cast<double>(d));
        std::set<Numbers> identities;
        for (std::size_t i = 0; i < scan.levels.size(); ++i)
        {
          SCOPED_TRACE(::testing::Message() << " i=" << i);
          auto const& level = scan.levels[i];
          auto const& state = level.state;
          ASSERT_TRUE(state.sz == sz && !state.spin_reversed && state.rapidities.size() == m && state.delta == d)
              << "XXZ excitation metadata";
          ASSERT_TRUE(identities.insert(state.quantum_numbers).second) << "unique XXZ excitation identity";
          ASSERT_TRUE(level.gap && *level.gap == state.energy - scan.ground_state.energy &&
                      *level.gap >= -Real{256} * Real(n) * eps)
              << "XXZ global excitation gap";
          if (i > 0)
          {
            auto const& previous = scan.levels[i - 1].state;
            ASSERT_TRUE(previous.energy <= state.energy) << "XXZ energy ordering";
            if (previous.energy == state.energy)
            {
              ASSERT_TRUE(previous.quantum_numbers < state.quantum_numbers) << "XXZ deterministic exact ties";
            }
          }
          ASSERT_NO_FATAL_FAILURE(check_state(n, state));
          double const target =
              static_cast<double>(state.energy) + 0.371 * std::cos(static_cast<double>(state.momentum));
          auto found = std::min_element(
              ed.begin(), ed.end(), [&](double a, double b) { return std::abs(a - target) < std::abs(b - target); });

          ASSERT_TRUE(found != ed.end()) << "XXZ excitation E/P and multiplicities vs ED";
          ASSERT_REAL_NEAR(*found, target, 3e-11) << "XXZ excitation E/P and multiplicities vs ED";

          ed.erase(found);
          auto const restart =
              model::solve_real<Real>(n, d, state.quantum_numbers, {.max_iterations = 0}, state.rapidities);
          ASSERT_TRUE(restart.converged && restart.iterations == 0 && restart.energy == state.energy)
              << "XXZ restart residual and energy";
          if (d == Real{0})
          {
            Real exact = Real{0};
            for (auto number : state.quantum_numbers)
              exact += cos(pi - pi * Real(number.twice()) / Real(n));

            ASSERT_REAL_NEAR(state.energy, exact, Real{128} * Real(n) * eps)
                << "XX excitation free-particle energy and one-update convergence";
            ASSERT_TRUE(state.iterations <= 1) << "XX excitation free-particle energy and one-update convergence";
          }
        }
        // Enumerate subsets of the full XXX grid independently, applying exact
        // integer infinity bounds at Delta=0,1/2,1. No production window helper.
        if (d == Real{0} || d == Real{1} / Real{2} || d == Real{1})
        {
          std::size_t expected = 0;
          auto const slots = n - m;
          for (unsigned mask = 0; mask < (1U << slots); ++mask)
            if (std::popcount(mask) == static_cast<int>(m))
            {
              Numbers numbers;
              bool allowed = true;
              for (unsigned j = 0; j < slots; ++j)
                if (mask & (1U << j))
                {
                  int const q = -static_cast<int>(slots - 1) + 2 * static_cast<int>(j);
                  if (d == Real{0}) allowed = allowed && 2 * std::abs(q) < static_cast<int>(n);
                  if (d == Real{1} / Real{2}) allowed = allowed && 3 * std::abs(q) < static_cast<int>(2 * n - m + 1);
                  numbers.push_back(half(q));
                }
              ASSERT_TRUE(identities.contains(numbers) == allowed) << "independent XXZ excitation window enumeration";
              expected += allowed;
            }
          ASSERT_TRUE(expected == scan.candidate_count) << "independent XXZ candidate count";
        }
        auto const few = model::real_excitations<Real>(n, d, sz, {.count = 2});
        auto const reversed = model::real_excitations<Real>(n, d, -sz, {.count = 2});
        ASSERT_TRUE(reversed.converged());
        ASSERT_EQ(reversed.levels.size(), few.levels.size());
        ASSERT_TRUE(few.candidate_count == scan.candidate_count &&
                    few.levels.size() == std::min(std::size_t{2}, scan.levels.size()))
            << "XXZ bounded heap scans entire family";
        for (std::size_t i = 0; i < few.levels.size(); ++i)
        {
          SCOPED_TRACE(::testing::Message() << " i=" << i);
          ASSERT_TRUE(few.levels[i].state.quantum_numbers == scan.levels[i].state.quantum_numbers &&
                      few.levels[i].state.energy == scan.levels[i].state.energy)
              << "XXZ lowest prefix";
          auto const& rev = reversed.levels[i].state;
          ASSERT_TRUE(rev.sz == -sz && rev.spin_reversed == (sz.twice() != 0) &&
                      rev.energy == few.levels[i].state.energy && rev.rapidities == few.levels[i].state.rapidities &&
                      rev.momentum_index == few.levels[i].state.momentum_index &&
                      reversed.levels[i].gap == few.levels[i].gap)
              << "XXZ spin-reversed excitations";
        }
        if (d == Real{1})
        {
          auto const xxx = bethe::heisenberg::real_excitations<Real>(n, sz, {.count = all_count});
          ASSERT_EQ(xxx.levels.size(), scan.levels.size());
          ASSERT_TRUE(xxx.candidate_count == scan.candidate_count) << "XXX endpoint candidate count";
          for (std::size_t i = 0; i < scan.levels.size(); ++i)
          {
            SCOPED_TRACE(::testing::Message() << " i=" << i);
            auto const& x = xxx.levels[i].state;
            auto const& z = scan.levels[i].state;
            ASSERT_TRUE(x.rapidities == z.rapidities && x.energy == z.energy && x.iterations == z.iterations &&
                        x.residual_norm == z.residual_norm && x.momentum_index == z.momentum_index &&
                        xxx.levels[i].gap == scan.levels[i].gap)
                << "XXX endpoint exact excitation delegation";
          }
        }
      }
}

TYPED_TEST(XXZExcitations, WindowThreshold)
{
  using Real = TypeParam;

  auto const all_count = std::numeric_limits<std::size_t>::max();
  ASSERT_TRUE(model::real_excitation_count(8, Real{0}, half_int{2}) == 6 &&
              model::real_excitation_count(8, Real{1} / Real{2}, half_int{2}) == 6 &&
              model::real_excitation_count(8, Real{1}, half_int{2}) == 15)
      << "anisotropy-dependent family size";
  // At Delta=1/2 the outer M=2, N=8 labels reach infinity. Exclude the exact
  // boundary in every precision, include them safely above it, never round down.
  Real const offset = Real{1} / Real{1000000};
  for (Real d : {Real{1} / Real{2} - offset, Real{1} / Real{2}, Real{1} / Real{2} + offset})
  {
    SCOPED_TRACE(::testing::Message() << " d=" << uni20::format_scalar(d));
    auto const scan = model::real_excitations<Real>(8, d, half_int{2}, {.count = all_count});
    ASSERT_TRUE(scan.converged() && scan.candidate_count == (d > Real{1} / Real{2} ? 15 : 6))
        << "infinity threshold crossing";
    for (auto const& level : scan.levels)
      ASSERT_NO_FATAL_FAILURE(check_state(8, level.state));
  }
}

TYPED_TEST(XXZExcitations, NativeGap)
{
  using Real = TypeParam;

  using std::abs;
  using std::sqrt;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  auto const scan = model::real_excitations<Real>(4, Real{3} / Real{4}, half_int{1}, {.count = 1});
  ASSERT_TRUE(scan.converged());
  ASSERT_EQ(scan.levels.size(), 1u);
  auto const gap = scan.levels[0].gap;
  Real const exact = (Real{3} / Real{4} + sqrt(Real{9} / Real{16} + Real{8})) / Real{2} - Real{1};

  ASSERT_TRUE(gap) << "XXZ excitation gap native precision";
  ASSERT_REAL_NEAR(*gap, exact, Real{128} * eps) << "XXZ excitation gap native precision";

  test_support::expect_exact(uni20::parse_real<Real>(uni20::format_real(*gap)), *gap, "XXZ gap I/O round trip");
  if constexpr (uni20::numeric_limits<Real>::digits > uni20::numeric_limits<double>::digits)
  {
    ASSERT_TRUE(abs(Real(static_cast<double>(exact)) - exact) > Real{128} * eps)
        << "XXZ gap oracle rejects double narrowing";
  }
}

TYPED_TEST(XXZExcitations, FailureSemantics)
{
  using Real = TypeParam;

  auto const failed =
      model::real_excitations<Real>(6, Real{1} / Real{2}, half_int{1}, {.count = 1}, {.max_iterations = 0});
  ASSERT_TRUE(!failed.converged() && failed.first_unconverged && failed.levels.empty() && failed.converged_count == 0 &&
              failed.candidate_count == 6)
      << "XXZ failed candidates excluded, not ranked";
  auto const vacuum = model::real_excitations<Real>(6, Real{1} / Real{2}, half_int{-3}, {}, {.max_iterations = 0});
  ASSERT_TRUE(vacuum.family_converged() && !vacuum.converged() && vacuum.levels.size() == 1 && !vacuum.levels[0].gap &&
              vacuum.levels[0].state.spin_reversed)
      << "XXZ failed ground does not invalidate vacuum energy";
  auto const half_filling = model::real_excitations<Real>(8, Real{1} / Real{2}, half_int{0});
  ASSERT_EQ(half_filling.levels.size(), 1u);
  ASSERT_TRUE(half_filling.candidate_count == 1 && half_filling.levels[0].gap == Real{0})
      << "no complete Sz=0 spectrum claim";
}

TYPED_TEST(XXZExcitations, LargeChains)
{
  using Real = TypeParam;

  for (auto sz : {half_int{0}, half_int{1}})
    ASSERT_TRUE(model::real_excitations<Real>(24, Real{3} / Real{4}, sz, {.count = 2}).converged())
        << "larger XXZ excitation scan";
}

TYPED_TEST(XXZExcitations, InvalidInputs)
{
  using Real = TypeParam;

  using Numbers = model::QuantumNumbers;
  auto half = [](std::int64_t twice) { return uni20::from_twice(twice); };
  Real const offset = Real{1} / Real{1000000};
  EXPECT_THROW(([&] { (void)model::real_excitations<Real>(8, Real{1}, half_int{2}, {.max_candidates = 14}); })(),
               std::length_error);
  ASSERT_TRUE(model::real_excitations<Real>(8, Real{1}, half_int{2}, {.max_candidates = 15}).candidate_count == 15)
      << "exact limit accepted";
  EXPECT_THROW(([&] { (void)model::real_excitation_count(1000, Real{1}, half_int{100}); })(), std::length_error);
  EXPECT_THROW(([&] { (void)model::real_excitations<Real>(1000000000, Real{1} / Real{2}, half_int{1}); })(),
               std::length_error);
  EXPECT_THROW(([&] { (void)model::real_excitations<Real>(4, Real{1} / Real{2}, half_int{1}, {.count = 0}); })(),
               std::invalid_argument);
  EXPECT_THROW(
      ([&] { (void)model::real_excitations<Real>(4, Real{1} / Real{2}, half_int{1}, {.max_candidates = 0}); })(),
      std::invalid_argument);
  for (auto const& numbers :
       {Numbers{half(1)}, Numbers{half(2), half(2)}, Numbers{half(1), half(-1)}, Numbers{half(-3), half(-1), half(1)},
        Numbers{half(6)}, Numbers{half(std::numeric_limits<std::int64_t>::min())}})
    EXPECT_THROW(([&] { (void)model::solve_real<Real>(4, Real{1} / Real{2}, numbers); })(), std::invalid_argument);
  Numbers const edge{half(-5), half(5)};
  EXPECT_THROW(([&] { (void)model::solve_real<Real>(8, Real{1} / Real{2}, edge); })(), std::invalid_argument);
  ASSERT_TRUE(model::solve_real<Real>(8, Real{1} / Real{2} + offset, edge).converged)
      << "finite labels above threshold accepted";
  Numbers const one{half_int{0}};
  for (Real bad : {Real{1}, Real{2}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    SCOPED_TRACE(::testing::Message() << " bad=" << uni20::format_scalar(bad));
    std::vector<Real> roots{bad};
    EXPECT_THROW(([&] { (void)model::solve_real<Real>(4, Real{0}, one, {}, roots); })(), std::invalid_argument);
  }
  std::vector<Real> wrong_size(2);
  EXPECT_THROW(([&] { (void)model::solve_real<Real>(4, Real{0}, one, {}, wrong_size); })(), std::invalid_argument);
}

TYPED_TEST(XXZExcitations, SpecifiedBudget)
{
  using Real = TypeParam;

  using Numbers = model::QuantumNumbers;
  auto half = [](std::int64_t twice) { return uni20::from_twice(twice); };
  Numbers const excited{half(-3), half(1)};
  auto const zero = model::solve_real<Real>(8, Real{1} / Real{2}, excited, {.max_iterations = 0});
  auto const first = model::solve_real<Real>(8, Real{1} / Real{2}, excited, {.max_iterations = 1});
  ASSERT_TRUE(!zero.converged && zero.iterations == 0 && !first.converged && first.iterations == 1)
      << "XXZ excitation budget semantics";
  ASSERT_NO_FATAL_FAILURE(check_state(8, zero));
  ASSERT_NO_FATAL_FAILURE(check_state(8, first));
}

template <typename Real> class XXZ : public ::testing::Test {};
TYPED_TEST_SUITE(XXZ, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXZ, SectorsAgainstED)
{
  using Real = TypeParam;

  using std::abs;
  using std::atan;
  using std::cos;
  using Options = model::SolverOptions<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const pi = Real{4} * atan(Real{1});

  for (Real d : {Real{0}, Real{1} / Real{10}, Real{1} / Real{2}, Real{9} / Real{10}, Real{1}, Real{101} / Real{100},
                 Real{2}, Real{10}})
    for (unsigned n = 2; n <= 9; ++n)
    {
      SCOPED_TRACE(::testing::Message() << " d=" << uni20::format_scalar(d));
      SCOPED_TRACE(::testing::Message() << " n=" << n);
      auto const states = model::sector_ground_states<Real>(n, d);
      ASSERT_TRUE(states.size() == n + 1) << "XXZ sector count";
      for (unsigned m = 0; m <= n / 2; ++m)
      {
        SCOPED_TRACE(::testing::Message() << " m=" << m);
        auto const& state = states[n - m];
        auto const spin = uni20::from_twice(static_cast<std::int64_t>(n - 2 * m));
        ASSERT_TRUE(state.converged && state.residual_norm <= Options{}.residual_tolerance) << "XXZ sector convergence";
        ASSERT_TRUE(state.sz == spin && state.delta == d && state.rapidities.size() == m) << "XXZ sector metadata";
        auto const ed = test_support::exact_spectrum(n, m, 0, true, static_cast<double>(d));
        EXPECT_REAL_NEAR(static_cast<double>(state.energy), ed.front(), 3e-11) << "XXZ sector minimum vs ED";
        auto const ed_p = test_support::exact_spectrum(n, m, 0.371, true, static_cast<double>(d));
        double const target = static_cast<double>(state.energy) + 0.371 * std::cos(static_cast<double>(state.momentum));
        ASSERT_TRUE(std::any_of(ed_p.begin(), ed_p.end(), [&](double e) { return std::abs(e - target) < 3e-11; }))
            << "XXZ joint energy/momentum vs ED";
        ASSERT_TRUE(states[m].energy == state.energy && states[m].rapidities == state.rapidities &&
                    states[m].sz == -spin && states[m].spin_reversed == (m != n - m) &&
                    states[m].momentum_index == state.momentum_index)
            << "XXZ spin reversal";
        auto const direct = model::sector_ground_state<Real>(n, d, -spin);
        ASSERT_TRUE(direct.energy == state.energy && direct.rapidities == state.rapidities) << "XXZ direct negative Sz";
        ASSERT_NO_FATAL_FAILURE(check_state(n, state));
        if (d == Real{1})
        {
          auto const xxx = bethe::heisenberg::sector_ground_state<Real>(n, spin);
          ASSERT_TRUE(state.energy == xxx.energy && state.rapidities == xxx.rapidities &&
                      state.quantum_numbers == xxx.quantum_numbers && state.momentum_index == xxx.momentum_index &&
                      state.momentum == xxx.momentum && state.iterations == xxx.iterations &&
                      state.residual_norm == xxx.residual_norm)
              << "Delta=1 delegates exactly to XXX";
        }
        if (d == Real{0})
        {
          // Jordan-Wigner: periodic fermion momenta for odd M, antiperiodic
          // for even M. Fill the M lowest cos(k), independently of Bethe labels.
          std::vector<Real> levels;
          for (unsigned j = 0; j < n; ++j)
            levels.push_back(cos(pi * Real(2 * j + (m % 2 == 0 ? 1 : 0)) / Real(n)));
          std::sort(levels.begin(), levels.end());
          Real exact = Real{0};
          for (unsigned j = 0; j < m; ++j)
            exact += levels[j];
          EXPECT_REAL_NEAR(state.energy, exact, Real{128} * Real(n) * eps) << "XX free-fermion spectrum";
          ASSERT_TRUE(state.iterations <= 1) << "XX endpoint needs at most one update";
        }
      }
      ASSERT_TRUE(model::ground_state<Real>(n, d).energy == states[(n + 1) / 2].energy) << "XXZ ground wrapper";
    }
}

TYPED_TEST(XXZ, NativePrecision)
{
  using Real = TypeParam;

  using std::abs;
  using std::sqrt;
  Real const eps = uni20::numeric_limits<Real>::epsilon();

  // N=4 analytic ground energy -(Delta+sqrt(Delta^2+8))/2 tests native
  // precision, parsing beyond double, and both nonsingular endpoint limits.
  for (Real d :
       {Real{0}, Real{1} / Real{2}, Real{3} / Real{4}, uni20::parse_real<Real>("0.1234567890123456789012345678901234"),
        Real{1024} * eps, Real{1} - Real{1024} * eps, Real{1}, Real{1} + Real{1024} * eps, Real{2}, Real{10}})
  {
    SCOPED_TRACE(::testing::Message() << " d=" << uni20::format_scalar(d));
    auto const state = model::ground_state<Real>(4, d);
    Real const exact = -(d + sqrt(d * d + Real{8})) / Real{2};

    ASSERT_TRUE(state.converged) << "irrational XXZ ground energy precision";
    ASSERT_REAL_NEAR(state.energy, exact, Real{128} * eps) << "irrational XXZ ground energy precision";

    test_support::expect_exact(uni20::parse_real<Real>(uni20::format_real(state.energy)), state.energy,
                               "XXZ energy I/O round trip");
    ASSERT_NO_FATAL_FAILURE(check_state(4, state));
    EXPECT_REAL_NEAR(model::ground_state<Real>(2, d).energy + Real{1} + d / Real{2}, Real{0}, Real{32} * eps)
        << "XXZ N=2 double-bond convention";
    EXPECT_REAL_NEAR(model::ground_state<Real>(3, d).energy + Real{1} / Real{2} + d / Real{4}, Real{0}, Real{64} * eps)
        << "XXZ N=3 frustrated odd-chain minimum";
    if constexpr (uni20::numeric_limits<Real>::digits > uni20::numeric_limits<double>::digits)
      if (d == Real{3} / Real{4} || d == Real{2})
      {
        ASSERT_TRUE(abs(static_cast<Real>(static_cast<double>(exact)) - exact) > Real{128} * eps)
            << "XXZ oracle must discriminate double narrowing";
      }
  }
}

TYPED_TEST(XXZ, IterationBudgets)
{
  using Real = TypeParam;

  using std::sqrt;
  using Options = model::SolverOptions<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();

  // Exactly at Delta=1, budget semantics and roots must also match XXX.
  for (Real d : {Real{0}, Real{1} / Real{2}, Real{1}, Real{2}, Real{10}})
  {
    SCOPED_TRACE(::testing::Message() << " d=" << uni20::format_scalar(d));
    auto const zero = model::ground_state<Real>(4, d, Options{.max_iterations = 0});
    ASSERT_TRUE(!zero.converged && zero.iterations == 0 && zero.energy == -Real{2} - d) << "XXZ zero budget";
    ASSERT_NO_FATAL_FAILURE(check_state(4, zero));
    auto const one = model::ground_state<Real>(4, d, Options{.max_iterations = 1});
    ASSERT_TRUE(one.iterations == 1 && one.converged == (d == Real{0})) << "XXZ one-update budget";
    EXPECT_REAL_NEAR(one.energy + d + sqrt(Real{2}), Real{0}, Real{32} * eps) << "XXZ simultaneous update";
    ASSERT_NO_FATAL_FAILURE(check_state(4, one));
    auto const vacuum = model::sector_ground_state<Real>(8, d, half_int{-4}, Options{.max_iterations = 0});
    ASSERT_TRUE(vacuum.converged && vacuum.rapidities.empty() && vacuum.energy == Real{2} * d && vacuum.spin_reversed)
        << "XXZ polarized vacuum";
  }
}

TYPED_TEST(XXZ, LargeChains)
{
  using Real = TypeParam;

  for (std::size_t n : {16, 65, 128})
    for (Real d : {Real{0}, Real{1} / Real{2}, Real{99} / Real{100}, Real{101} / Real{100}, Real{2}, Real{10}})
    {
      SCOPED_TRACE(::testing::Message() << " n=" << n);
      SCOPED_TRACE(::testing::Message() << " d=" << uni20::format_scalar(d));
      auto const state = model::ground_state<Real>(n, d);
      ASSERT_TRUE(state.converged) << "larger XXZ chain convergence";
      ASSERT_NO_FATAL_FAILURE(check_state(n, state));
    }
}

TYPED_TEST(XXZ, InvalidInputs)
{
  using Real = TypeParam;

  using Options = model::SolverOptions<Real>;

  for (Real bad : {-Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    SCOPED_TRACE(::testing::Message() << " bad=" << uni20::format_scalar(bad));
    EXPECT_THROW(([&] { (void)model::ground_state<Real>(4, bad); })(), std::invalid_argument);
    EXPECT_THROW(([&] { (void)model::sector_ground_states<Real>(4, bad); })(), std::invalid_argument);
  }
  for (Real bad :
       {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    EXPECT_THROW(([&] { (void)model::ground_state<Real>(4, Real{0}, Options{.residual_tolerance = bad}); })(),
                 std::invalid_argument);
  for (std::size_t n : {std::size_t{0}, std::size_t{1}, std::numeric_limits<std::size_t>::max()})
    EXPECT_THROW(([&] { (void)model::ground_state<Real>(n, Real{0}); })(), std::invalid_argument);
  for (auto spin :
       {half_int{3}, half_int{-3}, half_int::parse("1/2"), uni20::from_twice(std::numeric_limits<std::int64_t>::min())})
    EXPECT_THROW(([&] { (void)model::sector_ground_state<Real>(4, Real{0}, spin); })(), std::invalid_argument);
  EXPECT_THROW(([&] { (void)model::sector_ground_state<Real>(5, Real{0}, half_int{0}); })(), std::invalid_argument);
}

TYPED_TEST(XXZ, MassiveRationalEquations)
{
  using Real = TypeParam;
  using Complex = uni20::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {7, 8, 16})
    for (std::size_t m = 1; m <= n / 2; ++m)
      for (Real d : {Real{101} / Real{100}, Real{2}, Real{10}})
      {
        auto const state = model::sector_ground_state(n, d, uni20::from_twice(std::int64_t(n - 2 * m)));
        ASSERT_TRUE(state.converged);
        Real const t1 = std::sqrt((d - Real{1}) / (d + Real{1})), t2 = Real{2} * t1 / (Real{1} + t1 * t1);
        // sin(lambda+i*a)/sin(lambda-i*a), divided by cosh(a).
        auto const ratio = [](Real lambda, Real t) {
          return Complex{std::sin(lambda), t * std::cos(lambda)} / Complex{std::sin(lambda), -t * std::cos(lambda)};
        };
        std::vector<Real> lambda;
        for (Real z : state.rapidities)
          lambda.push_back(std::atan(t1 * z));
        bool winding = false;
        for (std::size_t j = 0; j < m; ++j)
        {
          Complex lhs{Real{1}, Real{0}}, rhs = lhs;
          auto const bare = ratio(lambda[j], t1);
          for (std::size_t i = 0; i < n; ++i)
            lhs *= bare;
          for (std::size_t a = 0; a < m; ++a)
            if (a != j)
            {
              Real const difference = lambda[j] - lambda[a];
              rhs *= ratio(difference, t2);
              winding = winding || std::cos(difference) < Real{0};
            }
          EXPECT_REAL_NEAR(lhs.real(), rhs.real(), Real{256} * Real(n) * eps);
          EXPECT_REAL_NEAR(lhs.imag(), rhs.imag(), Real{256} * Real(n) * eps);
        }
        if (n == 16 && m == 8 && d == Real{10}) EXPECT_TRUE(winding);
      }
}

TYPED_TEST(XXZ, IsingLimitAndLargeAnisotropy)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n : {7, 8, 16})
    for (std::size_t m = 1; m <= n / 2; ++m)
    {
      Real previous = Real{1};
      for (Real d : {Real{1000}, Real{10000}, Real{100000}})
      {
        auto const state = model::sector_ground_state(n, d, uni20::from_twice(std::int64_t(n - 2 * m)));
        ASSERT_TRUE(state.converged);
        Real sum = Real{0};
        for (auto q : state.quantum_numbers)
          sum += Real(q.twice()) / Real{2};
        Real error = Real{0};
        for (std::size_t j = 0; j < m; ++j)
        {
          // At Delta=infinity theta_1=2*lambda and theta_2=2*(lambda-lambda').
          Real const exact = pi * (Real(state.quantum_numbers[j].twice()) / Real{2} - sum / Real(n)) / Real(n - m);
          Real const lambda = std::atan(std::sqrt((d - Real{1}) / (d + Real{1})) * state.rapidities[j]);
          error = std::max(error, std::abs(lambda - exact));
        }
        EXPECT_LT(error, previous / Real{4} + Real{128} * eps);
        previous = error;
        EXPECT_LT(std::abs(state.energy / d - (Real(n) / Real{4} - Real(m))), Real(m) / d + Real{128} * eps);
      }
    }
  // Avoid overflow in scattering and in N*Delta/4 when the answer is finite.
  Real const huge = uni20::numeric_limits<Real>::max();
  auto const ground = model::ground_state(4, huge);
  ASSERT_TRUE(ground.converged);
  EXPECT_EQ(ground.energy, -huge);
  EXPECT_EQ(model::sector_ground_state(4, huge, half_int{1}).energy, -Real{1});
  EXPECT_EQ(model::sector_ground_state(4, huge, half_int{2}).energy, huge);
  EXPECT_THROW((void)model::sector_ground_state(5, huge, uni20::from_twice(std::int64_t{5})), std::runtime_error);
  EXPECT_THROW((void)model::solve_real<Real>(4, Real{2}, model::sector_ground_quantum_numbers(4, half_int{0})),
               std::invalid_argument);
  EXPECT_THROW((void)model::real_excitation_count(4, Real{2}, half_int{0}), std::invalid_argument);
  EXPECT_TRUE(model::open::ground_state(4, Real{2}).converged);
}

} // namespace
