// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include <bethe/heisenberg_excitations.hpp>
#include <bethe/xxz_excitations.hpp>
#include <uni20/core/scalar_io.hpp>

#include "exact_spectrum.hpp"

#include <iostream>
#include <set>
#include <string>
#include <string_view>

namespace
{
namespace model = bethe::xxz;
using uni20::half_int;
void require(bool condition, std::string_view message)
{
  if (!condition) throw std::runtime_error(std::string(message));
}
template <typename Exception = std::invalid_argument, typename Function> void rejects(Function&& function)
{
  try
  {
    function();
  }
  catch (Exception const&)
  {
    return;
  }
  throw std::runtime_error("expected exception");
}

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
  Real const gamma = acos(d);
  Real const scale = tan(gamma / Real{2});
  Real residual = Real{0};
  Real energy = Real(n) * d / Real{4};
  Real momentum = Real{0};
  auto const& z = state.rapidities;
  require(z.size() == state.quantum_numbers.size(), "XXZ root/label count");
  for (std::size_t i = 0; i < z.size(); ++i)
  {
    Real const theta = Real{2} * atan(z[i]);
    Real f = Real(n) * theta - pi * Real(state.quantum_numbers[i].twice());
    if (d != Real{0})
      for (std::size_t j = 0; j < z.size(); ++j)
        if (i != j)
        {
          if (d == Real{1})
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
    require(uni20::isfinite(z[i]), "finite XXZ rapidities");
    if (d < Real{1}) require(abs(scale * z[i]) < Real{1}, "XXZ roots stay on finite real branch");
    if (state.converged && i > 0) require(z[i - 1] < z[i], "ordered XXZ roots");
  }
  Real const allowance = Real{128} * Real(n) * eps;
  require(abs(energy - state.energy) < allowance, "XXZ returned energy vs roots");
  require(abs(residual - state.residual_norm) < allowance, "XXZ independent logarithmic residual");
  require(state.momentum_index < n, "XXZ momentum index range");
  if (state.converged)
  {
    Real const phase_allowance = allowance + Real(n) * state.residual_norm;
    require(abs(cos(momentum) - cos(state.momentum)) < phase_allowance &&
                abs(sin(momentum) - sin(state.momentum)) < phase_allowance,
            "XXZ momentum from roots");
  }
}

template <uni20::Real Real> void check_excitations()
{
  using std::abs;
  using std::atan;
  using std::cos;
  using std::sqrt;
  using Numbers = model::QuantumNumbers;
  auto half = [](std::int64_t twice) { return uni20::from_twice(twice); };
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const pi = Real{4} * atan(Real{1});
  auto const all_count = std::numeric_limits<std::size_t>::max();
  for (Real d : {Real{0}, Real{1} / Real{10}, Real{1} / Real{2}, Real{9} / Real{10}, Real{1}})
    for (unsigned n = 2; n <= 9; ++n)
      for (unsigned m = 0; m <= n / 2; ++m)
      {
        auto const sz = half(n - 2 * m);
        auto const scan = model::real_excitations<Real>(n, d, sz, {.count = all_count});
        require(scan.converged() && !scan.first_unconverged, "XXZ excitation family convergence");
        require(scan.sz == sz && scan.delta == d && scan.candidate_count == scan.levels.size() &&
                    scan.candidate_count == model::real_excitation_count(n, d, sz),
                "XXZ scan accounting");
        auto ed = test_support::exact_spectrum(n, m, 0.371, true, static_cast<double>(d));
        std::set<Numbers> identities;
        for (std::size_t i = 0; i < scan.levels.size(); ++i)
        {
          auto const& level = scan.levels[i];
          auto const& state = level.state;
          require(state.sz == sz && !state.spin_reversed && state.rapidities.size() == m && state.delta == d,
                  "XXZ excitation metadata");
          require(identities.insert(state.quantum_numbers).second, "unique XXZ excitation identity");
          require(level.gap && *level.gap == state.energy - scan.ground_state.energy &&
                      *level.gap >= -Real{256} * Real(n) * eps,
                  "XXZ global excitation gap");
          if (i > 0)
          {
            auto const& previous = scan.levels[i - 1].state;
            require(previous.energy <= state.energy, "XXZ energy ordering");
            if (previous.energy == state.energy)
              require(previous.quantum_numbers < state.quantum_numbers, "XXZ deterministic exact ties");
          }
          check_state(n, state);
          double const target =
              static_cast<double>(state.energy) + 0.371 * std::cos(static_cast<double>(state.momentum));
          auto found = std::min_element(
              ed.begin(), ed.end(), [&](double a, double b) { return std::abs(a - target) < std::abs(b - target); });
          require(found != ed.end() && std::abs(*found - target) < 3e-11,
                  "XXZ excitation E/P and multiplicities vs ED");
          ed.erase(found);
          auto const restart =
              model::solve_real<Real>(n, d, state.quantum_numbers, {.max_iterations = 0}, state.rapidities);
          require(restart.converged && restart.iterations == 0 && restart.energy == state.energy,
                  "XXZ restart residual and energy");
          if (d == Real{0})
          {
            Real exact = Real{0};
            for (auto number : state.quantum_numbers)
              exact += cos(pi - pi * Real(number.twice()) / Real(n));
            require(abs(state.energy - exact) < Real{128} * Real(n) * eps && state.iterations <= 1,
                    "XX excitation free-particle energy and one-update convergence");
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
              require(identities.contains(numbers) == allowed, "independent XXZ excitation window enumeration");
              expected += allowed;
            }
          require(expected == scan.candidate_count, "independent XXZ candidate count");
        }
        auto const few = model::real_excitations<Real>(n, d, sz, {.count = 2});
        auto const reversed = model::real_excitations<Real>(n, d, -sz, {.count = 2});
        require(few.candidate_count == scan.candidate_count &&
                    few.levels.size() == std::min(std::size_t{2}, scan.levels.size()),
                "XXZ bounded heap scans entire family");
        for (std::size_t i = 0; i < few.levels.size(); ++i)
        {
          require(few.levels[i].state.quantum_numbers == scan.levels[i].state.quantum_numbers &&
                      few.levels[i].state.energy == scan.levels[i].state.energy,
                  "XXZ lowest prefix");
          auto const& rev = reversed.levels[i].state;
          require(rev.sz == -sz && rev.spin_reversed == (sz.twice() != 0) && rev.energy == few.levels[i].state.energy &&
                      rev.rapidities == few.levels[i].state.rapidities &&
                      rev.momentum_index == few.levels[i].state.momentum_index &&
                      reversed.levels[i].gap == few.levels[i].gap,
                  "XXZ spin-reversed excitations");
        }
        if (d == Real{1})
        {
          auto const xxx = bethe::heisenberg::real_excitations<Real>(n, sz, {.count = all_count});
          require(xxx.candidate_count == scan.candidate_count, "XXX endpoint candidate count");
          for (std::size_t i = 0; i < scan.levels.size(); ++i)
          {
            auto const& x = xxx.levels[i].state;
            auto const& z = scan.levels[i].state;
            require(x.rapidities == z.rapidities && x.energy == z.energy && x.iterations == z.iterations &&
                        x.residual_norm == z.residual_norm && x.momentum_index == z.momentum_index &&
                        xxx.levels[i].gap == scan.levels[i].gap,
                    "XXX endpoint exact excitation delegation");
          }
        }
      }
  require(model::real_excitation_count(8, Real{0}, half_int{2}) == 6 &&
              model::real_excitation_count(8, Real{1} / Real{2}, half_int{2}) == 6 &&
              model::real_excitation_count(8, Real{1}, half_int{2}) == 15,
          "anisotropy-dependent family size");
  // At Delta=1/2 the outer M=2, N=8 labels reach infinity. Exclude the exact
  // boundary in every precision, include them safely above it, never round down.
  Real const offset = Real{1} / Real{1000000};
  for (Real d : {Real{1} / Real{2} - offset, Real{1} / Real{2}, Real{1} / Real{2} + offset})
  {
    auto const scan = model::real_excitations<Real>(8, d, half_int{2}, {.count = all_count});
    require(scan.converged() && scan.candidate_count == (d > Real{1} / Real{2} ? 15 : 6),
            "infinity threshold crossing");
    for (auto const& level : scan.levels)
      check_state(8, level.state);
  }
  auto const gap = model::real_excitations<Real>(4, Real{3} / Real{4}, half_int{1}, {.count = 1}).levels[0].gap;
  Real const exact = (Real{3} / Real{4} + sqrt(Real{9} / Real{16} + Real{8})) / Real{2} - Real{1};
  require(gap && abs(*gap - exact) < Real{128} * eps, "XXZ excitation gap native precision");
  require(uni20::parse_real<Real>(uni20::format_real(*gap)) == *gap, "XXZ gap I/O round trip");
  if constexpr (uni20::numeric_limits<Real>::digits > uni20::numeric_limits<double>::digits)
    require(abs(Real(static_cast<double>(exact)) - exact) > Real{128} * eps, "XXZ gap oracle rejects double narrowing");
  auto const failed =
      model::real_excitations<Real>(6, Real{1} / Real{2}, half_int{1}, {.count = 1}, {.max_iterations = 0});
  require(!failed.converged() && failed.first_unconverged && failed.levels.empty() && failed.converged_count == 0 &&
              failed.candidate_count == 6,
          "XXZ failed candidates excluded, not ranked");
  auto const vacuum = model::real_excitations<Real>(6, Real{1} / Real{2}, half_int{-3}, {}, {.max_iterations = 0});
  require(vacuum.family_converged() && !vacuum.converged() && vacuum.levels.size() == 1 && !vacuum.levels[0].gap &&
              vacuum.levels[0].state.spin_reversed,
          "XXZ failed ground does not invalidate vacuum energy");
  auto const half_filling = model::real_excitations<Real>(8, Real{1} / Real{2}, half_int{0});
  require(half_filling.candidate_count == 1 && half_filling.levels[0].gap == Real{0},
          "no complete Sz=0 spectrum claim");
  for (auto sz : {half_int{0}, half_int{1}})
    require(model::real_excitations<Real>(24, Real{3} / Real{4}, sz, {.count = 2}).converged(),
            "larger XXZ excitation scan");
  rejects<std::length_error>(
      [&] { (void)model::real_excitations<Real>(8, Real{1}, half_int{2}, {.max_candidates = 14}); });
  require(model::real_excitations<Real>(8, Real{1}, half_int{2}, {.max_candidates = 15}).candidate_count == 15,
          "exact limit accepted");
  rejects<std::length_error>([&] { (void)model::real_excitation_count(1000, Real{1}, half_int{100}); });
  rejects<std::length_error>([&] { (void)model::real_excitations<Real>(1000000000, Real{1} / Real{2}, half_int{1}); });
  rejects([&] { (void)model::real_excitations<Real>(4, Real{1} / Real{2}, half_int{1}, {.count = 0}); });
  rejects([&] { (void)model::real_excitations<Real>(4, Real{1} / Real{2}, half_int{1}, {.max_candidates = 0}); });
  for (auto const& numbers :
       {Numbers{half(1)}, Numbers{half(2), half(2)}, Numbers{half(1), half(-1)}, Numbers{half(-3), half(-1), half(1)},
        Numbers{half(6)}, Numbers{half(std::numeric_limits<std::int64_t>::min())}})
    rejects([&] { (void)model::solve_real<Real>(4, Real{1} / Real{2}, numbers); });
  Numbers const edge{half(-5), half(5)};
  rejects([&] { (void)model::solve_real<Real>(8, Real{1} / Real{2}, edge); });
  require(model::solve_real<Real>(8, Real{1} / Real{2} + offset, edge).converged,
          "finite labels above threshold accepted");
  Numbers const one{half_int{0}};
  for (Real bad : {Real{1}, Real{2}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    std::vector<Real> roots{bad};
    rejects([&] { (void)model::solve_real<Real>(4, Real{0}, one, {}, roots); });
  }
  std::vector<Real> wrong_size(2);
  rejects([&] { (void)model::solve_real<Real>(4, Real{0}, one, {}, wrong_size); });
  Numbers const excited{half(-3), half(1)};
  auto const zero = model::solve_real<Real>(8, Real{1} / Real{2}, excited, {.max_iterations = 0});
  auto const first = model::solve_real<Real>(8, Real{1} / Real{2}, excited, {.max_iterations = 1});
  require(!zero.converged && zero.iterations == 0 && !first.converged && first.iterations == 1,
          "XXZ excitation budget semantics");
  check_state(8, zero);
  check_state(8, first);
}

template <uni20::Real Real> void tests()
{
  check_excitations<Real>();
  using std::abs;
  using std::atan;
  using std::cos;
  using std::sqrt;
  using Options = model::SolverOptions<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const pi = Real{4} * atan(Real{1});

  for (Real d : {Real{0}, Real{1} / Real{10}, Real{1} / Real{2}, Real{9} / Real{10}, Real{1}})
    for (unsigned n = 2; n <= 9; ++n)
    {
      auto const states = model::sector_ground_states<Real>(n, d);
      require(states.size() == n + 1, "XXZ sector count");
      for (unsigned m = 0; m <= n / 2; ++m)
      {
        auto const& state = states[n - m];
        auto const spin = uni20::from_twice(static_cast<std::int64_t>(n - 2 * m));
        require(state.converged && state.residual_norm <= Options{}.residual_tolerance, "XXZ sector convergence");
        require(state.sz == spin && state.delta == d && state.rapidities.size() == m, "XXZ sector metadata");
        auto const ed = test_support::exact_spectrum(n, m, 0, true, static_cast<double>(d));
        require(std::abs(static_cast<double>(state.energy) - ed.front()) < 3e-11, "XXZ sector minimum vs ED");
        auto const ed_p = test_support::exact_spectrum(n, m, 0.371, true, static_cast<double>(d));
        double const target = static_cast<double>(state.energy) + 0.371 * std::cos(static_cast<double>(state.momentum));
        require(std::any_of(ed_p.begin(), ed_p.end(), [&](double e) { return std::abs(e - target) < 3e-11; }),
                "XXZ joint energy/momentum vs ED");
        require(states[m].energy == state.energy && states[m].rapidities == state.rapidities && states[m].sz == -spin &&
                    states[m].spin_reversed == (m != n - m) && states[m].momentum_index == state.momentum_index,
                "XXZ spin reversal");
        auto const direct = model::sector_ground_state<Real>(n, d, -spin);
        require(direct.energy == state.energy && direct.rapidities == state.rapidities, "XXZ direct negative Sz");
        check_state(n, state);
        if (d == Real{1})
        {
          auto const xxx = bethe::heisenberg::sector_ground_state<Real>(n, spin);
          require(state.energy == xxx.energy && state.rapidities == xxx.rapidities &&
                      state.quantum_numbers == xxx.quantum_numbers && state.momentum_index == xxx.momentum_index &&
                      state.momentum == xxx.momentum && state.iterations == xxx.iterations &&
                      state.residual_norm == xxx.residual_norm,
                  "Delta=1 delegates exactly to XXX");
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
          require(abs(state.energy - exact) < Real{128} * Real(n) * eps, "XX free-fermion spectrum");
          require(state.iterations <= 1, "XX endpoint needs at most one update");
        }
      }
      require(model::ground_state<Real>(n, d).energy == states[(n + 1) / 2].energy, "XXZ ground wrapper");
    }

  // N=4 analytic ground energy -(Delta+sqrt(Delta^2+8))/2 tests native
  // precision, parsing beyond double, and both nonsingular endpoint limits.
  for (Real d :
       {Real{0}, Real{1} / Real{2}, Real{3} / Real{4}, uni20::parse_real<Real>("0.1234567890123456789012345678901234"),
        Real{1024} * eps, Real{1} - Real{1024} * eps, Real{1}})
  {
    auto const state = model::ground_state<Real>(4, d);
    Real const exact = -(d + sqrt(d * d + Real{8})) / Real{2};
    require(state.converged && abs(state.energy - exact) < Real{128} * eps, "irrational XXZ ground energy precision");
    require(uni20::parse_real<Real>(uni20::format_real(state.energy)) == state.energy, "XXZ energy I/O round trip");
    check_state(4, state);
    require(abs(model::ground_state<Real>(2, d).energy + Real{1} + d / Real{2}) < Real{32} * eps,
            "XXZ N=2 double-bond convention");
    require(abs(model::ground_state<Real>(3, d).energy + Real{1} / Real{2} + d / Real{4}) < Real{64} * eps,
            "XXZ N=3 frustrated odd-chain minimum");
    if constexpr (uni20::numeric_limits<Real>::digits > uni20::numeric_limits<double>::digits)
      if (d == Real{3} / Real{4})
        require(abs(static_cast<Real>(static_cast<double>(exact)) - exact) > Real{128} * eps,
                "XXZ oracle must discriminate double narrowing");
  }
  // Exactly at Delta=1, budget semantics and roots must also match XXX.
  for (Real d : {Real{0}, Real{1} / Real{2}, Real{1}})
  {
    auto const zero = model::ground_state<Real>(4, d, Options{.max_iterations = 0});
    require(!zero.converged && zero.iterations == 0 && zero.energy == -Real{2} - d, "XXZ zero budget");
    check_state(4, zero);
    auto const one = model::ground_state<Real>(4, d, Options{.max_iterations = 1});
    require(one.iterations == 1 && one.converged == (d == Real{0}), "XXZ one-update budget");
    require(abs(one.energy + d + sqrt(Real{2})) < Real{32} * eps, "XXZ simultaneous update");
    check_state(4, one);
    auto const vacuum = model::sector_ground_state<Real>(8, d, half_int{-4}, Options{.max_iterations = 0});
    require(vacuum.converged && vacuum.rapidities.empty() && vacuum.energy == Real{2} * d && vacuum.spin_reversed,
            "XXZ polarized vacuum");
  }
  for (std::size_t n : {16, 65, 128})
    for (Real d : {Real{0}, Real{1} / Real{2}, Real{99} / Real{100}})
    {
      auto const state = model::ground_state<Real>(n, d);
      require(state.converged, "larger XXZ chain convergence");
      check_state(n, state);
    }
  for (Real bad :
       {-Real{1}, Real{2}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    rejects([&] { (void)model::ground_state<Real>(4, bad); });
    rejects([&] { (void)model::sector_ground_states<Real>(4, bad); });
  }
  for (Real bad :
       {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    rejects([&] { (void)model::ground_state<Real>(4, Real{0}, Options{.residual_tolerance = bad}); });
  for (std::size_t n : {std::size_t{0}, std::size_t{1}, std::numeric_limits<std::size_t>::max()})
    rejects([&] { (void)model::ground_state<Real>(n, Real{0}); });
  for (auto spin :
       {half_int{3}, half_int{-3}, half_int::parse("1/2"), uni20::from_twice(std::numeric_limits<std::int64_t>::min())})
    rejects([&] { (void)model::sector_ground_state<Real>(4, Real{0}, spin); });
  rejects([&] { (void)model::sector_ground_state<Real>(5, Real{0}, half_int{0}); });
}
} // namespace

int main(int argc, char** argv)
{
  try
  {
    require(argc == 2, "pass a precision");
    std::string_view const precision = argv[1];
    if (precision == "fp64")
      tests<double>();
    else if (precision == "long-double")
      tests<long double>();
#if UNI20_HAS_FLOAT128
    else if (precision == "fp128")
      tests<uni20::float128>();
#endif
    else
      throw std::invalid_argument("unavailable test precision");
    std::cout << "XXZ checks passed for " << precision << '\n';
    return 0;
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
