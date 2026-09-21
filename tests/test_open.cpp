// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include <bethe/heisenberg_open.hpp>
#include <uni20/core/scalar_io.hpp>

#include "exact_spectrum.hpp"

#include <iostream>
#include <string>
#include <string_view>

namespace
{
namespace model = bethe::heisenberg::open;
using bethe::heisenberg::QuantumNumbers;

void require(bool condition, std::string_view message)
{
  if (!condition) throw std::runtime_error(std::string(message));
}

template <typename Function> void invalid_argument(Function&& function)
{
  try
  {
    function();
  }
  catch (std::invalid_argument const&)
  {
    return;
  }
  throw std::runtime_error("expected invalid_argument");
}

template <typename State>
concept HasMomentum = requires(State state) {
  state.momentum;
  state.momentum_index;
};
static_assert(!HasMomentum<model::RealState<double>>);

template <uni20::Real Real> void check_state(std::size_t n, model::RealState<Real> const& state)
{
  using std::abs;
  using std::atan;
  Real const pi = Real{4} * atan(Real{1});
  Real residual = Real{0};
  Real energy = static_cast<Real>(n - 1) / Real{4};
  require(state.quantum_numbers.size() == state.rapidities.size(), "open root/quantum-number count");
  for (std::size_t i = 0; i < state.rapidities.size(); ++i)
  {
    Real const z = state.rapidities[i];
    // Direct logarithmic equation, independently of the update angle and
    // compensated sum. Both direct and reflected self scattering are absent.
    Real f = Real{4} * Real(n) * atan(z) - pi * Real(state.quantum_numbers[i].twice());
    for (std::size_t j = 0; j < state.rapidities.size(); ++j)
      if (i != j)
        f -= Real{2} * (atan((z - state.rapidities[j]) / Real{2}) + atan((z + state.rapidities[j]) / Real{2}));
    residual = std::max(residual, abs(f) / (Real{2} * Real(n)));
    energy -= Real{2} / (Real{1} + z * z);
    if (state.converged)
    {
      require(z > Real{0}, "open physical roots must be positive");
      if (i > 0) require(state.rapidities[i - 1] < z, "open roots strictly ordered");
    }
  }
  Real const allowance = Real{64} * Real(n) * uni20::numeric_limits<Real>::epsilon();
  require(abs(energy - state.energy) < allowance, "open energy matches returned roots");
  require(abs(residual - state.residual_norm) < allowance, "open residual matches returned roots");
}

template <uni20::Real Real> void tests()
{
  using std::abs;
  using std::atan;
  using std::cos;
  using std::sqrt;
  using uni20::half_int;
  using Options = bethe::heisenberg::SolverOptions<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const pi = Real{4} * atan(Real{1});

  auto const two = model::ground_state<Real>(2);
  require(two.converged && abs(two.energy + Real{3} / Real{4}) < Real{16} * eps, "N=2 one-bond normalization");
  require(abs(two.rapidities[0] - Real{1}) < Real{8} * eps, "N=2 exact open rapidity");
  auto const three = model::ground_state<Real>(3);
  require(three.converged && abs(three.energy + Real{1}) < Real{16} * eps, "N=3 open ground energy");
  require(abs(three.rapidities[0] - Real{1} / sqrt(Real{3})) < Real{8} * eps, "N=3 exact open rapidity");
  auto const four = model::ground_state<Real>(4);
  Real const exact_four = -Real{3} / Real{4} - sqrt(Real{3}) / Real{2};
  Real const precision_tolerance = Real{128} * eps;
  require(four.converged && abs(four.energy - exact_four) < precision_tolerance, "N=4 irrational open energy");
  if constexpr (uni20::numeric_limits<Real>::digits > uni20::numeric_limits<double>::digits)
    require(abs(static_cast<Real>(static_cast<double>(exact_four)) - exact_four) > precision_tolerance,
            "open precision oracle must reject double narrowing");
  require(uni20::parse_real<Real>(uni20::format_real(four.energy)) == four.energy, "open energy text round trip");

  // Every sector, both parities of N, with independently assembled N-1 bonds.
  for (unsigned n = 2; n <= 9; ++n)
  {
    auto const states = model::sector_ground_states<Real>(n);
    require(states.size() == n + 1, "open sector count");
    for (unsigned down = 0; down <= n; ++down)
    {
      auto const& state = states[n - down];
      auto const spectrum = test_support::exact_spectrum(n, down, 0, false);
      require(state.converged && state.residual_norm <= Options{}.residual_tolerance, "open sector convergence");
      require(std::abs(static_cast<double>(state.energy) - spectrum.front()) < 2e-11, "open sector energy vs ED");
      require(state.sz == uni20::from_twice(static_cast<std::int64_t>(n) - 2 * down), "open sector label");
      require(state.spin_reversed == (2 * down > n), "open spin-reversed reference");
      require(state.rapidities == states[down].rapidities && state.energy == states[down].energy,
              "open spin reversal reuses roots and energy");
      for (std::size_t i = 0; i < state.quantum_numbers.size(); ++i)
        require(state.quantum_numbers[i] == half_int(static_cast<std::int64_t>(i + 1)), "open consecutive integers");
      check_state(n, state);
    }
    auto const ground = model::ground_state<Real>(n);
    require(ground.energy == states[(n + 1) / 2].energy, "open ground wrapper selects minimal |Sz|");
  }

  // Every supported real-root configuration through N=8, not only ground states.
  for (unsigned n = 2; n <= 8; ++n)
    for (unsigned m = 1; m <= n / 2; ++m)
    {
      auto const spectrum = test_support::exact_spectrum(n, m, 0, false);
      QuantumNumbers numbers;
      auto visit = [&](auto&& self, unsigned first) -> void {
        if (numbers.size() == m)
        {
          auto const state = model::solve_real<Real>(n, numbers);
          require(state.converged, "open excited-state convergence");
          double nearest = std::numeric_limits<double>::infinity();
          for (double e : spectrum)
            nearest = std::min(nearest, std::abs(e - static_cast<double>(state.energy)));
          require(nearest < 2e-11, "open excited real-root energy vs ED");
          require(state.quantum_numbers == numbers, "open state identity");
          check_state(n, state);
          return;
        }
        for (unsigned i = first; i <= n - m; ++i)
        {
          numbers.emplace_back(i);
          self(self, i + 1);
          numbers.pop_back();
        }
      };
      visit(visit, 1);
    }

  // The M=1 standing waves have q=(N-I)*pi/N; no fictitious self reflection.
  for (std::size_t n = 2; n <= 12; ++n)
    for (std::size_t i = 1; i < n; ++i)
    {
      auto const state = model::solve_real<Real>(n, QuantumNumbers{half_int(static_cast<std::int64_t>(i))});
      Real const exact = Real(n - 1) / Real{4} - Real{1} - cos(pi * Real(i) / Real(n));
      require(state.converged && abs(state.energy - exact) < Real{64} * Real(n) * eps, "open one-magnon spectrum");
    }
  for (std::size_t n : {16, 65, 128})
  {
    auto const state = model::ground_state<Real>(n);
    require(state.converged, "larger open-chain convergence");
    check_state(n, state);
  }

  auto const zero = model::ground_state<Real>(4, Options{.max_iterations = 0});
  require(!zero.converged && zero.iterations == 0 && zero.energy == -Real{13} / Real{4}, "open zero budget");
  check_state(4, zero);
  auto const one = model::ground_state<Real>(4, Options{.max_iterations = 1});
  require(!one.converged && one.iterations == 1, "open one-update budget");
  require(abs(one.energy + Real{5} / Real{4} + sqrt(Real{2}) / Real{2}) < Real{16} * eps, "open simultaneous update");
  check_state(4, one);
  auto const restart = model::solve_real<Real>(4, four.quantum_numbers, Options{.max_iterations = 0}, four.rapidities);
  require(restart.converged && restart.iterations == 0 && restart.energy == four.energy, "open restart guess");
  auto const vacuum = model::solve_real<Real>(5, QuantumNumbers{}, Options{.max_iterations = 0});
  require(vacuum.converged && vacuum.energy == Real{1} && vacuum.sz == half_int::parse("5/2"), "open empty state");

  for (std::size_t n : {std::size_t{0}, std::size_t{1}, std::numeric_limits<std::size_t>::max()})
    invalid_argument([&] { (void)model::ground_state<Real>(n); });
  for (auto sz :
       {half_int::parse("1/2"), half_int{3}, half_int{-3}, uni20::from_twice(std::numeric_limits<std::int64_t>::min())})
    invalid_argument([&] { (void)model::sector_ground_state<Real>(4, sz); });
  invalid_argument([&] { (void)model::sector_ground_state<Real>(5, half_int{0}); });
  for (QuantumNumbers const& numbers :
       {QuantumNumbers{half_int{0}}, QuantumNumbers{half_int{-1}}, QuantumNumbers{half_int::parse("1/2")},
        QuantumNumbers{half_int{4}}, QuantumNumbers{half_int{1}, half_int{1}}, QuantumNumbers{half_int{2}, half_int{1}},
        QuantumNumbers{half_int{1}, half_int{3}}, QuantumNumbers{half_int{1}, half_int{2}, half_int{3}}})
    invalid_argument([&] { (void)model::solve_real<Real>(4, numbers); });
  invalid_argument([&] { (void)model::solve_real<Real>(4, four.quantum_numbers, {}, std::vector<Real>{Real{1}}); });
  for (Real bad : {-Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    invalid_argument(
        [&] { (void)model::solve_real<Real>(4, four.quantum_numbers, {}, std::vector<Real>{Real{1}, bad}); });
  for (Real bad :
       {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    invalid_argument([&] { (void)model::ground_state<Real>(4, Options{.residual_tolerance = bad}); });
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
    std::cout << "Open-chain checks passed for " << precision << '\n';
    return 0;
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
