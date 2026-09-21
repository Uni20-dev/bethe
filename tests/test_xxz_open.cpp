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
namespace model = bethe::xxz::open;
using uni20::half_int;
using Numbers = model::QuantumNumbers;
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
template <typename State> constexpr bool has_momentum = requires(State s) { s.momentum; };
static_assert(!has_momentum<model::RealState<double>>);

// Independently evaluate the conventional hyperbolic equations, including
// the boundary reflection phase, and energy from one-magnon wave numbers.
template <uni20::Real Real> void check_state(std::size_t n, model::RealState<Real> const& state)
{
  using std::abs;
  using std::acos;
  using std::atan;
  using std::atanh;
  using std::cos;
  using std::tan;
  using std::tanh;
  Real const pi = Real{4} * atan(Real{1});
  Real const d = state.delta, gamma = acos(d), scale = tan(gamma / Real{2});
  Real residual = Real{0}, energy = Real(n - 1) * d / Real{4};
  auto const& z = state.rapidities;
  require(z.size() == state.quantum_numbers.size(), "open XXZ root/label count");
  for (std::size_t i = 0; i < z.size(); ++i)
  {
    Real const theta = Real{2} * atan(z[i]);
    Real const lambda = d == Real{1} ? z[i] / Real{2} : atanh(scale * z[i]);
    Real f = Real(2 * n) * theta - pi * Real(state.quantum_numbers[i].twice());
    if (d != Real{1}) f += Real{4} * atan(tanh(lambda) * scale);
    if (d != Real{0})
      for (std::size_t j = 0; j < z.size(); ++j)
        if (i != j)
        {
          Real const other = d == Real{1} ? z[j] / Real{2} : atanh(scale * z[j]);
          for (Real difference : {lambda - other, lambda + other})
            f -= Real{2} * atan(d == Real{1} ? difference : tanh(difference) / tan(gamma));
        }
    residual = std::max(residual, abs(f) / Real(2 * n));
    energy += cos(pi - theta) - d;
    require(uni20::isfinite(z[i]) && z[i] >= Real{0}, "finite nonnegative open XXZ roots");
    if (d < Real{1}) require(scale * z[i] < Real{1}, "finite hyperbolic rapidity");
    if (state.converged) require(z[i] > Real{0} && (i == 0 || z[i - 1] < z[i]), "positive ordered physical roots");
  }
  Real const allowance = Real{128} * Real(n) * uni20::numeric_limits<Real>::epsilon();
  require(abs(energy - state.energy) < allowance, "open XXZ returned energy vs roots");
  require(abs(residual - state.residual_norm) < allowance, "open XXZ independent boundary equations");
}

template <uni20::Real Real> void spectra()
{
  using std::abs;
  using std::atan;
  using std::cos;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * atan(Real{1});
  auto const all = std::numeric_limits<std::size_t>::max();
  for (Real d : {Real{0}, Real{1} / Real{10}, Real{1} / Real{2}, Real{9} / Real{10}, Real{1}})
    for (unsigned n = 2; n <= 9; ++n)
    {
      auto const sectors = model::sector_ground_states<Real>(n, d);
      require(sectors.size() == n + 1, "open XXZ sector count");
      for (unsigned m = 0; m <= n / 2; ++m)
      {
        auto const sz = uni20::from_twice(static_cast<std::int64_t>(n - 2 * m));
        auto ed = test_support::exact_spectrum(n, m, 0, false, static_cast<double>(d));
        auto const& ground = sectors[n - m];
        require(ground.converged && ground.sz == sz && ground.delta == d && ground.rapidities.size() == m &&
                    std::abs(static_cast<double>(ground.energy) - ed.front()) < 3e-11,
                "open XXZ sector minimum vs ED");
        require(sectors[m].sz == -sz && sectors[m].spin_reversed == (m != n - m) &&
                    sectors[m].rapidities == ground.rapidities && sectors[m].energy == ground.energy,
                "open XXZ sector spin reversal");
        auto const negative = model::sector_ground_state<Real>(n, d, -sz);
        require(negative.rapidities == ground.rapidities && negative.energy == ground.energy, "direct negative sector");
        check_state(n, ground);
        auto const scan = model::real_excitations<Real>(n, d, sz, {.count = all});
        require(scan.converged() && !scan.first_unconverged && scan.sz == sz && scan.delta == d &&
                    scan.candidate_count == scan.levels.size() &&
                    scan.candidate_count == model::real_excitation_count(n, d, sz),
                "open XXZ scan accounting");
        require(scan.levels.front().state.energy == ground.energy, "scan contains sector minimum");
        std::set<Numbers> identities;
        for (std::size_t i = 0; i < scan.levels.size(); ++i)
        {
          auto const& level = scan.levels[i];
          auto const& state = level.state;
          require(state.sz == sz && !state.spin_reversed && state.rapidities.size() == m && state.delta == d,
                  "open XXZ excitation metadata");
          require(identities.insert(state.quantum_numbers).second, "unique excitation labels");
          require(level.gap && *level.gap == state.energy - scan.ground_state.energy &&
                      *level.gap >= -Real{256} * Real(n) * eps,
                  "open XXZ global gaps");
          if (i != 0)
          {
            auto const& previous = scan.levels[i - 1].state;
            require(previous.energy <= state.energy, "energy ordering");
            if (previous.energy == state.energy)
              require(previous.quantum_numbers < state.quantum_numbers, "deterministic exact ties");
          }
          check_state(n, state);
          double const target = static_cast<double>(state.energy);
          auto found = std::min_element(
              ed.begin(), ed.end(), [&](double a, double b) { return std::abs(a - target) < std::abs(b - target); });
          require(found != ed.end() && std::abs(*found - target) < 3e-11, "open XXZ excitation multiplicities vs ED");
          ed.erase(found);
          auto const restart =
              model::solve_real<Real>(n, d, state.quantum_numbers, {.max_iterations = 0}, state.rapidities);
          require(restart.converged && restart.iterations == 0 && restart.energy == state.energy,
                  "restart uses returned-root residual");
          if (d == Real{0})
          {
            Real exact = Real{0};
            for (auto number : state.quantum_numbers)
              exact -= cos(pi * Real(number.twice()) / Real(2 * (n + 1)));
            require(abs(state.energy - exact) < Real{128} * Real(n) * eps && state.iterations <= 1,
                    "open XX standing waves use N+1");
          }
        }
        // Independent integer thresholds, without the production window helper.
        if (d == Real{0} || d == Real{1} / Real{2} || d == Real{1})
        {
          std::size_t expected = 0;
          for (unsigned mask = 0; mask < (1U << (n - m)); ++mask)
            if (std::popcount(mask) == static_cast<int>(m))
            {
              Numbers numbers;
              bool allowed = true;
              for (unsigned j = 1; j <= n - m; ++j)
                if (mask & (1U << (j - 1)))
                {
                  if (d == Real{0}) allowed = allowed && 2 * j < n + 1;
                  if (d == Real{1} / Real{2}) allowed = allowed && 3 * j < 2 * n - m + 2;
                  numbers.emplace_back(j);
                }
              require(identities.contains(numbers) == allowed, "independent positive label window");
              expected += allowed;
            }
          require(expected == scan.candidate_count, "independent open XXZ candidate count");
        }
        auto const few = model::real_excitations<Real>(n, d, sz, {.count = 2});
        auto const reversed = model::real_excitations<Real>(n, d, -sz, {.count = 2});
        require(few.candidate_count == scan.candidate_count &&
                    few.levels.size() == std::min(std::size_t{2}, scan.levels.size()),
                "bounded heap accounting");
        for (std::size_t i = 0; i < few.levels.size(); ++i)
        {
          auto const& state = few.levels[i].state;
          auto const& rev = reversed.levels[i].state;
          require(state.energy == scan.levels[i].state.energy &&
                      state.quantum_numbers == scan.levels[i].state.quantum_numbers,
                  "lowest prefix");
          require(rev.sz == -sz && rev.spin_reversed == (sz.twice() != 0) && rev.energy == state.energy &&
                      rev.rapidities == state.rapidities && reversed.levels[i].gap == few.levels[i].gap,
                  "spin reversed excitations");
        }
        if (d == Real{1})
        {
          auto const xxx = bethe::heisenberg::open::real_excitations<Real>(n, sz, {.count = all});
          require(xxx.candidate_count == scan.candidate_count, "open XXX endpoint count");
          for (std::size_t i = 0; i < scan.levels.size(); ++i)
          {
            auto const& x = xxx.levels[i].state;
            auto const& z = scan.levels[i].state;
            require(x.rapidities == z.rapidities && x.quantum_numbers == z.quantum_numbers && x.energy == z.energy &&
                        x.iterations == z.iterations && x.residual_norm == z.residual_norm &&
                        xxx.levels[i].gap == scan.levels[i].gap,
                    "exact XXX delegation");
          }
        }
      }
      require(model::ground_state<Real>(n, d).energy == sectors[(n + 1) / 2].energy, "ground wrapper");
    }
}

template <uni20::Real Real> void tests()
{
  spectra<Real>();
  using std::abs;
  using std::sqrt;
  using Options = model::SolverOptions<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (Real d :
       {Real{0}, Real{1} / Real{2}, Real{3} / Real{4}, uni20::parse_real<Real>("0.1234567890123456789012345678901234"),
        Real{1024} * eps, Real{1} - Real{1024} * eps, Real{1}})
  {
    auto const two = model::ground_state<Real>(2, d);
    require(two.converged && abs(two.energy + Real{1} / Real{2} + d / Real{4}) < Real{128} * eps &&
                abs(two.rapidities[0] - sqrt((Real{1} + d) / (Real{3} - d))) < Real{128} * eps,
            "two-site energy and boundary-phase root");
    auto const three = model::ground_state<Real>(3, d);
    Real const exact = -(d + sqrt(d * d + Real{8})) / Real{4};
    require(three.converged && abs(three.energy - exact) < Real{128} * eps, "native precision irrational ground");
    require(uni20::parse_real<Real>(uni20::format_real(three.energy)) == three.energy, "energy round trip");
    check_state(2, two);
    check_state(3, three);
    if (d == Real{3} / Real{4})
    {
      auto const scan = model::real_excitations<Real>(3, d, uni20::from_twice(std::int64_t{1}), {.count = 2});
      auto const gap = scan.levels.at(1).gap;
      require(gap && abs(*gap + exact) < Real{128} * eps, "native precision irrational gap");
      require(uni20::parse_real<Real>(uni20::format_real(*gap)) == *gap, "gap round trip");
      if constexpr (uni20::numeric_limits<Real>::digits > uni20::numeric_limits<double>::digits)
        require(abs(Real(static_cast<double>(exact)) - exact) > Real{128} * eps,
                "analytic oracle rejects double narrowing");
    }
  }
  for (Real d : {Real{0}, Real{1} / Real{2}, Real{1}})
  {
    auto const zero = model::ground_state<Real>(4, d, {.max_iterations = 0});
    require(!zero.converged && zero.iterations == 0 && zero.energy == -Real{2} - Real{5} * d / Real{4}, "zero budget");
    auto const first = model::ground_state<Real>(4, d, {.max_iterations = 1});
    require(first.iterations == 1 && first.converged == (d == Real{0}), "one-update budget");
    if (d < Real{1})
      require(abs(first.energy + sqrt(Real{5}) / Real{2} + Real{5} * d / Real{4}) < Real{128} * eps,
              "simultaneous update and N+1 denominator");
    else
    {
      auto const xxx = bethe::heisenberg::open::ground_state<Real>(4, {.max_iterations = 1});
      require(first.energy == xxx.energy && first.rapidities == xxx.rapidities &&
                  first.residual_norm == xxx.residual_norm,
              "XXX budget delegation");
    }
    check_state(4, zero);
    check_state(4, first);
  }
  Real const offset = Real{1} / Real{1000000};
  for (Real d : {Real{1} / Real{2} - offset, Real{1} / Real{2}, Real{1} / Real{2} + offset})
  {
    auto const scan = model::real_excitations<Real>(8, d, half_int{1}, {.count = 100});
    require(scan.converged() && scan.candidate_count == (d > Real{1} / Real{2} ? 10 : 4), "infinity threshold");
    for (auto const& level : scan.levels)
      check_state(8, level.state);
    require(model::real_excitation_count(4, d, half_int{1}) == (d > Real{1} / Real{2} ? 3 : 2),
            "one-magnon infinity threshold");
  }
  for (std::size_t n : {16, 65, 128})
    for (Real d : {Real{0}, Real{1} / Real{2}, Real{99} / Real{100}})
    {
      auto const state = model::ground_state<Real>(n, d);
      require(state.converged, "larger open XXZ chain convergence");
      check_state(n, state);
    }
  require(model::real_excitations<Real>(24, Real{3} / Real{4}, half_int{1}, {.count = 2}).converged(),
          "larger excitation scan");
  auto const failed = model::real_excitations<Real>(6, Real{1} / Real{2}, half_int{1}, {}, {.max_iterations = 0});
  require(!failed.converged() && failed.first_unconverged && failed.levels.empty() && failed.converged_count == 0,
          "failed candidates excluded");
  auto const vacuum = model::real_excitations<Real>(6, Real{1} / Real{2}, half_int{-3}, {}, {.max_iterations = 0});
  require(vacuum.family_converged() && !vacuum.converged() && vacuum.levels.size() == 1 && !vacuum.levels[0].gap &&
              vacuum.levels[0].state.spin_reversed && vacuum.levels[0].state.energy == Real{5} / Real{8},
          "vacuum valid with unavailable ground gap");
  require(model::real_excitation_count(8, Real{0}, half_int{2}) == 6 &&
              model::real_excitation_count(8, Real{1} / Real{2}, half_int{2}) == 10 &&
              model::real_excitation_count(8, Real{1}, half_int{2}) == 15,
          "anisotropy-dependent window");
  rejects<std::length_error>(
      [&] { (void)model::real_excitations<Real>(8, Real{1}, half_int{2}, {.max_candidates = 14}); });
  require(model::real_excitations<Real>(8, Real{1}, half_int{2}, {.max_candidates = 15}).candidate_count == 15,
          "exact cap accepted");
  rejects<std::length_error>([&] { (void)model::real_excitation_count(1000, Real{1}, half_int{100}); });
  rejects<std::length_error>([&] { (void)model::real_excitations<Real>(1000000000, Real{1} / Real{2}, half_int{1}); });
  rejects([&] { (void)model::real_excitations<Real>(4, Real{0}, half_int{1}, {.count = 0}); });
  rejects([&] { (void)model::real_excitations<Real>(4, Real{0}, half_int{1}, {.max_candidates = 0}); });
  for (auto const& numbers : {Numbers{half_int{0}}, Numbers{half_int{-1}}, Numbers{uni20::from_twice(std::int64_t{1})},
                              Numbers{half_int{1}, half_int{1}}, Numbers{half_int{2}, half_int{1}},
                              Numbers{half_int{1}, half_int{2}, half_int{3}}, Numbers{half_int{3}}})
    rejects([&] { (void)model::solve_real<Real>(4, Real{1} / Real{2}, numbers); });
  Numbers const one{half_int{1}};
  for (Real bad :
       {-Real{1}, Real{1}, Real{2}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    std::vector<Real> roots{bad};
    rejects([&] { (void)model::solve_real<Real>(4, Real{0}, one, {}, roots); });
  }
  std::vector<Real> wrong_size(2);
  rejects([&] { (void)model::solve_real<Real>(4, Real{0}, one, {}, wrong_size); });
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
  for (auto sz : {half_int{3}, half_int{-3}, uni20::from_twice(std::int64_t{1}),
                  uni20::from_twice(std::numeric_limits<std::int64_t>::min())})
    rejects([&] { (void)model::sector_ground_state<Real>(4, Real{0}, sz); });
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
    std::cout << "Open XXZ checks passed for " << precision << '\n';
    return 0;
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
