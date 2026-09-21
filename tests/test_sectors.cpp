// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include <bethe/heisenberg.hpp>
#include <uni20/core/scalar_io.hpp>

#include <bit>
#include <iostream>
#include <string>
#include <string_view>

namespace
{
void require(bool condition, std::string_view message)
{
  if (!condition)
    throw std::runtime_error(std::string(message));
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

// Independent small-chain oracle: construct H directly in the Sz bit basis.
// Optionally add a*(T+T^-1)/2, which shifts an eigenvalue by a*cos(P), to
// check momentum as well as energy without a momentum-space Bethe formula.
// Double is intentional only in this ED oracle; separate analytic regressions
// below check the solver and dispersions in the selected precision.
std::vector<double> exact_spectrum(unsigned n, unsigned down, double translation_weight = 0)
{
  std::vector<unsigned> basis;
  std::vector<std::size_t> index(1U << n);
  for (unsigned bits = 0; bits < (1U << n); ++bits)
    if (std::popcount(bits) == static_cast<int>(down))
    {
      index[bits] = basis.size();
      basis.push_back(bits);
    }
  auto const dim = basis.size();
  std::vector<double> matrix(dim * dim, 0);
  auto at = [&](std::size_t i, std::size_t j) -> double& { return matrix[i * dim + j]; };
  for (std::size_t i = 0; i < dim; ++i)
  {
    auto const bits = basis[i];
    for (unsigned site = 0; site < n; ++site)
    {
      unsigned const next = (site + 1) % n;
      bool const opposite = ((bits >> site) & 1U) != ((bits >> next) & 1U);
      at(i, i) += opposite ? -0.25 : 0.25;
      if (opposite)
        at(index[bits ^ (1U << site) ^ (1U << next)], i) += 0.5;
    }
    unsigned const translated = ((bits << 1) & ((1U << n) - 1)) | (bits >> (n - 1));
    at(index[translated], i) += translation_weight / 2;
    at(i, index[translated]) += translation_weight / 2;
  }
  // Cyclic Jacobi diagonalization of a real symmetric matrix.
  bool diagonal = false;
  for (int sweep = 0; sweep < 80; ++sweep)
  {
    double largest = 0;
    for (std::size_t p = 0; p < dim; ++p)
      for (std::size_t q = p + 1; q < dim; ++q)
      {
        double const off = at(p, q);
        largest = std::max(largest, std::abs(off));
        if (std::abs(off) < 1e-14)
          continue;
        double const tau = (at(q, q) - at(p, p)) / (2 * off);
        double const t = std::copysign(1.0, tau) / (std::abs(tau) + std::hypot(1.0, tau));
        double const c = 1 / std::sqrt(1 + t * t);
        double const s = t * c;
        at(p, p) -= t * off;
        at(q, q) += t * off;
        at(p, q) = at(q, p) = 0;
        for (std::size_t r = 0; r < dim; ++r)
          if (r != p && r != q)
          {
            double const rp = at(r, p), rq = at(r, q);
            at(r, p) = at(p, r) = c * rp - s * rq;
            at(r, q) = at(q, r) = s * rp + c * rq;
          }
      }
    if (largest < 1e-13)
    {
      diagonal = true;
      break;
    }
  }
  require(diagonal, "ED oracle diagonalization failed");
  std::vector<double> energies(dim);
  for (std::size_t i = 0; i < dim; ++i)
    energies[i] = at(i, i);
  std::sort(energies.begin(), energies.end());
  return energies;
}

template <uni20::Real Real> void check_state(std::size_t n, bethe::heisenberg::RealState<Real> const& state)
{
  using std::abs;
  using std::atan;
  using std::cos;
  using std::sin;
  Real const pi = Real{4} * atan(Real{1});
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real residual = Real{0};
  Real energy = static_cast<Real>(n) / Real{4};
  Real momentum = Real{0};
  require(state.quantum_numbers.size() == state.rapidities.size(), "root/quantum-number count");
  for (std::size_t i = 0; i < state.rapidities.size(); ++i)
  {
    auto const z = state.rapidities[i];
    Real f = Real{2} * static_cast<Real>(n) * atan(z) - pi * static_cast<Real>(state.quantum_numbers[i].twice());
    for (std::size_t j = 0; j < state.rapidities.size(); ++j)
      if (i != j)
        f -= Real{2} * atan((z - state.rapidities[j]) / Real{2});
    residual = std::max(residual, abs(f) / static_cast<Real>(n));
    energy -= Real{2} / (Real{1} + z * z);
    momentum += pi - Real{2} * atan(z);
    if (i > 0 && state.converged)
      require(state.rapidities[i - 1] < z, "ordered converged roots");
  }
  Real const allowance = Real{64} * static_cast<Real>(n) * eps;
  require(abs(energy - state.energy) < allowance, "returned energy matches roots");
  require(abs(residual - state.residual_norm) < allowance, "returned residual matches roots");
  require(state.momentum_index < n, "momentum index in range");
  if (state.converged)
  {
    // Summing individual momentum errors costs O(N) times the phase residual.
    Real const momentum_allowance =
        Real{4} * static_cast<Real>(n * n) * eps + static_cast<Real>(n) * state.residual_norm;
    require(abs(cos(momentum) - cos(state.momentum)) < momentum_allowance, "cos(P) from rapidities");
    require(abs(sin(momentum) - sin(state.momentum)) < momentum_allowance, "sin(P) from rapidities");
  }
}

template <uni20::Real Real> void tests()
{
  using namespace bethe::heisenberg;
  using std::abs;
  using std::atan;
  using std::cos;
  using std::sin;
  using std::sqrt;
  using uni20::half_int;
  using Options = SolverOptions<Real>;
  auto half = [](std::int64_t twice) { return uni20::from_twice(twice); };
  Real const pi = Real{4} * atan(Real{1});
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const tolerance = Real{4096} * eps;

  for (unsigned n = 2; n <= 9; ++n)
  {
    auto const sectors = sector_ground_states<Real>(n);
    require(sectors.size() == n + 1, "all Sz sectors present");
    for (unsigned m = 0; m <= n / 2; ++m)
    {
      auto const ed = exact_spectrum(n, m);
      auto const& state = sectors[n - m];
      require(state.converged && abs(state.energy - static_cast<Real>(ed.front())) < Real{1e-11},
              "sector minimum disagrees with exact diagonalization");
      require(state.sz == half(n - 2 * m) && state.rapidities.size() == m, "sector labels");
      auto const& reversed = sectors[m];
      require(reversed.energy == state.energy && reversed.sz == -state.sz, "spin-reversed sectors");
      require(reversed.momentum_index == state.momentum_index, "spin reversal preserves momentum");
      require(reversed.spin_reversed == (reversed.sz.twice() < 0), "spin-reversed reference label");
      check_state(n, state);
      auto const direct = sector_ground_state<Real>(n, -state.sz);
      require(direct.energy == state.energy, "direct negative Sz solve");
      if (n % 2 != 0 && m != 0)
      {
        auto reflected_numbers = state.quantum_numbers;
        std::reverse(reflected_numbers.begin(), reflected_numbers.end());
        for (auto& number : reflected_numbers)
          number = -number;
        auto const reflected = solve_real<Real>(n, reflected_numbers);
        require(reflected.converged && abs(reflected.energy - state.energy) < tolerance,
                "odd-N reflection-related sector minima");
        require((reflected.momentum_index + state.momentum_index) % n == 0, "reflected momentum");
      }
    }
    auto const gs = ground_state<Real>(n);
    require(gs.converged && gs.energy == sectors[(n + 1) / 2].energy, "ground-state wrapper");
  }

  // Precision-sensitive one-magnon sector oracle: E=N/4-1-cos(pi/N), odd N.
  auto const magnon = sector_ground_state<Real>(5, half(3));
  Real const exact_magnon = -sqrt(Real{5}) / Real{4};
  require(abs(magnon.energy - exact_magnon) < Real{64} * eps, "one-magnon irrational energy precision");
  if constexpr (uni20::numeric_limits<Real>::digits > uni20::numeric_limits<double>::digits)
    require(abs(static_cast<Real>(static_cast<double>(exact_magnon)) - exact_magnon) > Real{64} * eps,
            "one-magnon oracle discriminates double narrowing");

  // Exercise general specified states, not only the sector-ground or spinon
  // generators: every supported two-root configuration for these small chains.
  for (unsigned n : {6, 7, 8})
  {
    auto const ed = exact_spectrum(n, 2, 0.371);
    auto const bound = static_cast<std::int64_t>(n) - 3;
    for (std::int64_t first = -bound; first <= bound; first += 2)
      for (std::int64_t second = first + 2; second <= bound; second += 2)
      {
        QuantumNumbers const numbers{half(first), half(second)};
        auto const state = solve_real<Real>(n, numbers);
        require(state.converged && state.quantum_numbers == numbers, "specified two-root state");
        check_state(n, state);
        double const target = static_cast<double>(state.energy) + 0.371 * std::cos(static_cast<double>(state.momentum));
        require(std::any_of(ed.begin(), ed.end(), [&](double e) { return std::abs(e - target) < 1e-11; }),
                "specified real-root energy/momentum disagrees with ED");
      }
  }

  for (unsigned n : {3, 5, 7, 9})
  {
    auto const branch = one_spinon_branch<Real>(n);
    auto const ed = exact_spectrum(n, (n - 1) / 2, 0.371);
    require(branch.size() == (n + 1) / 2, "one-spinon state count");
    std::vector<std::size_t> momenta;
    for (std::size_t i = 0; i < branch.size(); ++i)
    {
      auto const& point = branch[i];
      auto const& state = point.state;
      require(state.converged && state.sz == half(1), "one-spinon convergence and Sz");
      require(point.hole.twice() == static_cast<std::int64_t>((n - 1) / 2) - 2 * static_cast<std::int64_t>(i),
              "hole enumeration");
      require(std::find(state.quantum_numbers.begin(), state.quantum_numbers.end(), point.hole) ==
                  state.quantum_numbers.end(),
              "specified hole is absent");
      require(point.spinon_momentum > Real{0} && point.spinon_momentum < pi, "spinon k range");
      if (i > 0)
        require(branch[i - 1].spinon_momentum < point.spinon_momentum, "ascending spinon k");
      auto const& mirror = branch[branch.size() - 1 - i];
      require(abs(point.spinon_momentum + mirror.spinon_momentum - pi) < tolerance, "spinon k reflection");
      require(abs(state.energy - mirror.state.energy) < tolerance, "spinon energy reflection");
      require(abs(point.bulk_subtracted_energy - (state.energy - Real(n) * bulk_energy_density<Real>())) < tolerance,
              "bulk energy subtraction");
      Real const p = pi * static_cast<Real>((n - 1) / 2) + pi / Real{2} - point.spinon_momentum;
      require(abs(cos(p) - cos(state.momentum)) < tolerance && abs(sin(p) - sin(state.momentum)) < tolerance,
              "spinon/lattice momentum offset");
      double const target = static_cast<double>(state.energy) + 0.371 * std::cos(static_cast<double>(state.momentum));
      require(std::any_of(ed.begin(), ed.end(), [&](double e) { return std::abs(e - target) < 1e-11; }),
              "one-spinon energy/momentum disagrees with ED translation oracle");
      momenta.push_back(state.momentum_index);
      check_state(n, state);
    }
    std::sort(momenta.begin(), momenta.end());
    require(std::adjacent_find(momenta.begin(), momenta.end()) == momenta.end(), "unique lattice momenta");
    require(abs(branch.front().state.energy - ground_state<Real>(n).energy) < tolerance, "branch endpoint is GS");
  }

  // N=5 central-hole state has roots +/-1, E=-3/4.
  auto const center = one_spinon_state<Real>(5, half(0));
  require(abs(center.state.energy + Real{3} / Real{4}) < tolerance, "N=5 central spinon exact energy");
  require(abs(center.state.rapidities[1] - Real{1}) < tolerance, "N=5 central spinon exact root");
  require(center.state.momentum_index == 0, "N=5 central spinon momentum");
  // At fixed k=pi/2, bulk-subtracted finite-size energies approach pi/2.
  auto const c17 = one_spinon_state<Real>(17, half(0));
  auto const c65 = one_spinon_state<Real>(65, half(0));
  require(c17.state.converged && c65.state.converged, "larger one-spinon convergence");
  require(abs(c65.bulk_subtracted_energy - pi / Real{2}) < abs(c17.bulk_subtracted_energy - pi / Real{2}),
          "one-spinon thermodynamic convergence");
  require(abs(c65.bulk_subtracted_energy - pi / Real{2}) < Real{0.02}, "one-spinon thermodynamic normalization");

  auto const no_updates = one_spinon_state<Real>(7, half(1), Options{.max_iterations = 0});
  require(!no_updates.state.converged && no_updates.state.iterations == 0, "spinon zero budget");
  check_state(7, no_updates.state);
  auto const one_update = sector_ground_state<Real>(7, half(1), Options{.max_iterations = 1});
  require(!one_update.converged && one_update.iterations == 1, "sector one-update budget");
  check_state(7, one_update);
  auto const vacuum = solve_real<Real>(8, QuantumNumbers{}, Options{.max_iterations = 0});
  require(vacuum.converged && vacuum.energy == Real{2} && vacuum.momentum_index == 0, "empty-root vacuum");
  auto const restarted = solve_real<Real>(5, center.state.quantum_numbers, {}, center.state.rapidities);
  require(restarted.converged && restarted.iterations == 0 && restarted.energy == center.state.energy,
          "initial roots remain separate from state identity");

  require(spinon_energy(Real{0}) == Real{0} && spinon_energy(pi) == Real{0}, "exact dispersion endpoints");
  require(abs(spinon_energy(pi / Real{2}) - pi / Real{2}) < Real{8} * eps, "XXX dispersion maximum");
  require(abs(spinon_energy(pi / Real{6}) - pi / Real{4}) < Real{8} * eps, "XXX dispersion interior");
  require(spinon_energy(pi / Real{3}, Real{2}) == Real{2} * spinon_energy(pi / Real{3}), "exchange scaling");
  for (Real const k : {Real{0}, pi / Real{7}, pi / Real{2}, pi})
  {
    require(bethe::xxz::spinon_energy(k, Real{1}) == spinon_energy(k), "XXZ isotropic limit");
    Real const xx = k == pi ? Real{0} : sin(k);
    require(abs(bethe::xxz::spinon_energy(k, Real{0}) - xx) < Real{8} * eps, "XX free-fermion dispersion");
    require(abs(bethe::xxz::spinon_energy(k, Real{1} - Real{8} * eps) - spinon_energy(k)) < Real{32} * eps,
            "XXZ stable near-isotropic limit");
    require(abs(bethe::xxz::spinon_energy(k, Real{1} / Real{2}) - Real{3} * sqrt(Real{3}) / Real{4} * xx) <
                Real{16} * eps,
            "XXZ Delta=1/2 analytic coefficient");
    require(abs(bethe::xxz::spinon_energy(k, -Real{1} / Real{2}) - Real{3} * sqrt(Real{3}) / Real{8} * xx) <
                Real{16} * eps,
            "XXZ Delta=-1/2 analytic coefficient");
  }

  for (auto const sz : {half(1), half(10), half(-10), half(std::numeric_limits<std::int64_t>::min())})
    invalid_argument([&] { (void)sector_ground_state<Real>(4, sz); });
  invalid_argument([&] { (void)sector_ground_state<Real>(5, half(0)); });
  for (QuantumNumbers const& numbers :
       {QuantumNumbers{half(-1), half(-1)}, QuantumNumbers{half(1), half(-1)}, QuantumNumbers{half(-2), half(2)},
        QuantumNumbers{half(-3), half(3)}, QuantumNumbers{half(-2), half(0), half(2)}})
    invalid_argument([&] { (void)solve_real<Real>(4, numbers); });
  invalid_argument([&] { (void)one_spinon_branch<Real>(4); });
  invalid_argument([&] { (void)one_spinon_state<Real>(5, half(1)); });
  invalid_argument([&] { (void)one_spinon_state<Real>(5, half(4)); });
  invalid_argument([&] { (void)solve_real<Real>(5, center.state.quantum_numbers, {}, std::vector<Real>{Real{0}}); });
  invalid_argument(
      [&]
      {
        (void)solve_real<Real>(5, center.state.quantum_numbers, {},
                               std::vector<Real>{Real{0}, uni20::numeric_limits<Real>::infinity()});
      });
  for (Real const bad : {-Real{1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
  {
    invalid_argument([&] { (void)spinon_energy(bad); });
    invalid_argument([&] { (void)spinon_energy(Real{1}, bad); });
    invalid_argument([&] { (void)bethe::xxz::spinon_energy(Real{1}, bad); });
  }
  invalid_argument([&] { (void)spinon_energy(pi + Real{1}); });
  invalid_argument([&] { (void)spinon_energy(Real{1}, Real{0}); });
  invalid_argument([&] { (void)bethe::xxz::spinon_energy(Real{1}, Real{2}); });
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
    std::cout << "Sector and spinon checks passed for " << precision << '\n';
    return 0;
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
