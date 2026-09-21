// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include <bethe/heisenberg_excitations.hpp>
#include <uni20/core/scalar_io.hpp>

#include "exact_spectrum.hpp"

#include <iostream>
#include <set>
#include <string>
#include <string_view>

namespace
{
using namespace bethe::heisenberg;
using uni20::half_int;

half_int half(std::int64_t twice) { return uni20::from_twice(twice); }

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

void remove_energy(std::vector<double>& spectrum, double target)
{
  auto const found = std::min_element(spectrum.begin(), spectrum.end(),
                                      [&](double a, double b) { return std::abs(a - target) < std::abs(b - target); });
  require(found != spectrum.end() && std::abs(*found - target) < 3e-11,
          "multiplet energy/momentum vs ED multiplicities");
  spectrum.erase(found);
}

// Subtract the Sz=S+1 spectrum as a multiset, leaving one entry per spin-S
// multiplet. Translation commutes with SU(2), so the shifted oracle works too.
auto spin_spectrum(unsigned n, unsigned m, bool periodic)
{
  auto result = test_support::exact_spectrum(n, m, periodic ? 0.371 : 0, periodic);
  if (m > 0)
    for (double energy : test_support::exact_spectrum(n, m - 1, periodic ? 0.371 : 0, periodic))
      remove_energy(result, energy);
  return result;
}

template <uni20::Real Real, bool Periodic> void check_boundary()
{
  auto scan = [](std::size_t n, half_int spin, RealExcitationOptions const& enumeration = {},
                 SolverOptions<Real> const& solver = {}) {
    if constexpr (Periodic)
      return real_excitations<Real>(n, spin, enumeration, solver);
    else
      return open::real_excitations<Real>(n, spin, enumeration, solver);
  };
  using std::abs;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (unsigned n = 2; n <= 8; ++n)
    for (unsigned m = 0; m <= n / 2; ++m)
    {
      auto const spin = half(n - 2 * m);
      auto const all = scan(n, spin, {.count = 1000});
      require(all.converged() && all.family_converged() && !all.first_unconverged, "exhaustive scan converged");
      require(all.spin == spin && all.candidate_count == all.converged_count, "scan spin and candidate accounting");
      require(all.levels.size() == real_excitation_count(n, spin), "full family returned");
      auto ed = spin_spectrum(n, m, Periodic);
      std::set<QuantumNumbers> identities;
      for (std::size_t i = 0; i < all.levels.size(); ++i)
      {
        auto const& level = all.levels[i];
        auto const& state = level.state;
        require(state.converged && state.sz == spin && !state.spin_reversed, "one highest-weight state per multiplet");
        require(state.quantum_numbers.size() == m && state.rapidities.size() == m, "excitation root count");
        require(level.gap && *level.gap == state.energy - all.ground_state.energy, "gap uses the global ground energy");
        require(*level.gap >= -Real{256} * Real(n) * eps, "excitation energy not below global ground beyond roundoff");
        require(identities.insert(state.quantum_numbers).second, "unique quantum-number configurations");
        if (i > 0)
        {
          auto const& previous = all.levels[i - 1].state;
          require(previous.energy <= state.energy, "energy ordering");
          if (previous.energy == state.energy)
            require(previous.quantum_numbers < state.quantum_numbers, "deterministic exact-tie ordering");
        }
        double target = static_cast<double>(state.energy);
        if constexpr (Periodic) target += 0.371 * std::cos(static_cast<double>(state.momentum));
        remove_energy(ed, target);
      }
      // Independent bit-subset enumeration tests the actual sets, not just a
      // binomial formula shared with the implementation.
      std::size_t expected = 0;
      auto const slots = n - m;
      for (unsigned mask = 0; mask < (1U << slots); ++mask)
        if (std::popcount(mask) == static_cast<int>(m))
        {
          QuantumNumbers numbers;
          for (unsigned i = 0; i < slots; ++i)
            if (mask & (1U << i))
              numbers.push_back(half(Periodic ? -static_cast<int>(slots - 1) + 2 * static_cast<int>(i)
                                              : 2 * static_cast<int>(i + 1)));
          require(identities.contains(numbers), "every independently generated configuration was visited");
          ++expected;
        }
      require(expected == all.candidate_count, "independent candidate count");
      auto const few = scan(n, spin, {.count = 2});
      require(few.levels.size() == std::min(std::size_t{2}, expected), "bounded retained count");
      require(few.candidate_count == expected && few.converged_count == expected, "count never truncates the scan");
      for (std::size_t i = 0; i < few.levels.size(); ++i)
        require(few.levels[i].state.quantum_numbers == all.levels[i].state.quantum_numbers &&
                    few.levels[i].state.energy == all.levels[i].state.energy && few.levels[i].gap == all.levels[i].gap,
                "bounded heap returns exactly the lowest prefix");
    }

  auto const singlet = scan(4, half_int{0});
  require(singlet.candidate_count == 1 && singlet.levels.size() == 1 && *singlet.levels[0].gap == Real{0},
          "even singlet real-root family contains only the ground configuration");
  require(spin_spectrum(4, 2, Periodic).size() == 2, "ED detects the omitted excited singlet");
  auto const vacuum = scan(4, half_int{2}, {}, {.max_iterations = 0});
  require(vacuum.family_converged() && !vacuum.converged() && vacuum.levels.size() == 1 && !vacuum.levels[0].gap,
          "failed ground reference leaves the vacuum energy valid but gap unavailable");
  auto const failed = scan(6, half_int{1}, {.count = 1}, {.max_iterations = 0});
  require(!failed.converged() && !failed.family_converged() && failed.first_unconverged && failed.levels.empty() &&
              failed.candidate_count == 6 && failed.converged_count == 0,
          "failed candidates are counted, never ranked");

  // Precision-sensitive gap, not just precision-sensitive absolute energies.
  if constexpr (!Periodic)
  {
    using std::sqrt;
    auto const triplet = scan(4, half_int{1}, {.count = 1});
    Real const exact_gap = (Real{1} + sqrt(Real{3}) - sqrt(Real{2})) / Real{2};
    require(abs(*triplet.levels[0].gap - exact_gap) < Real{64} * eps, "irrational gap at selected precision");
    if constexpr (uni20::numeric_limits<Real>::digits > uni20::numeric_limits<double>::digits)
      require(abs(static_cast<Real>(static_cast<double>(exact_gap)) - exact_gap) > Real{64} * eps,
              "gap oracle rejects double narrowing");
    require(uni20::parse_real<Real>(uni20::format_real(*triplet.levels[0].gap)) == *triplet.levels[0].gap,
            "excitation gap text round trip");
  }
  else
  {
    auto const triplets = scan(4, half_int{1});
    require(triplets.levels.size() == 3 &&
                abs(triplets.levels[1].state.energy - triplets.levels[2].state.energy) < Real{64} * eps,
            "degenerate triplets remain distinct multiplets");
    require(triplets.levels[1].state.momentum_index != triplets.levels[2].state.momentum_index,
            "degenerate momentum partners retained");
  }
  for (auto spin : {half_int{-1}, half_int{3}, half(1), uni20::from_twice(std::numeric_limits<std::int64_t>::min())})
    rejects([&] { (void)scan(4, spin); });
  rejects([&] { (void)scan(1, half_int{0}); });
  rejects([&] { (void)scan(4, half_int{1}, {.count = 0}); });
  rejects([&] { (void)scan(4, half_int{1}, {.max_candidates = 0}); });
  rejects([&] { (void)scan(4, half_int{1}, {}, {.residual_tolerance = Real{0}}); });
  rejects<std::length_error>([&] { (void)scan(4, half_int{1}, {.max_candidates = 2}); });
  require(scan(4, half_int{1}, {.max_candidates = 3}).candidate_count == 3, "exact candidate limit accepted");
}

template <uni20::Real Real> void tests()
{
  check_boundary<Real, true>();
  check_boundary<Real, false>();
  require(real_excitation_count(64, half_int{1}) == 528, "two-hole triplet count");
  require(real_excitation_count(65, half(1)) == 33, "odd one-spinon count");
  require(real_excitation_count(64, half_int{0}) == 1 && real_excitation_count(64, half_int{32}) == 1,
          "empty combination edge cases");
  rejects<std::length_error>([] { (void)real_excitation_count(1000, half_int{100}); });
  rejects<std::length_error>([] { (void)real_excitation_count(64, half_int{1}, 527); });
  // Preflight must reject before calling either solver, even with a huge N.
  rejects<std::length_error>([] {
    (void)detail::scan_real_excitations<RealState<Real>>(
        1000000000, half_int{1}, {}, true,
        [](QuantumNumbers const&) -> RealState<Real> { throw std::runtime_error("unexpected candidate solve"); },
        []() -> RealState<Real> { throw std::runtime_error("unexpected ground solve"); });
  });

  // Inject a failed high-energy candidate: the retained lowest state converges,
  // but that must not hide an incomplete scan. Also test a failed ground alone.
  auto const partial = detail::scan_real_excitations<open::RealState<Real>>(
      4, half_int{1}, {.count = 1}, false,
      [](QuantumNumbers const& numbers) {
        auto state = open::solve_real<Real>(4, numbers);
        if (numbers == QuantumNumbers{half_int{3}}) state.converged = false;
        return state;
      },
      [] { return open::ground_state<Real>(4); });
  require(partial.ground_state.converged && !partial.converged() && partial.converged_count == 2 &&
              partial.candidate_count == 3 && partial.levels.size() == 1 && partial.levels[0].state.converged &&
              partial.levels[0].gap && partial.first_unconverged->quantum_numbers == QuantumNumbers{half_int{3}},
          "non-retained failure invalidates family ordering, not valid gaps");
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
    std::cout << "Excitation checks passed for " << precision << '\n';
    return 0;
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
