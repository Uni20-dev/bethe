// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "hubbard_ed.hpp"
#include <array>
#include <bethe/heisenberg.hpp>
#include <bethe/hubbard.hpp>
#include <bit>
#include <iostream>
#include <string>
#include <string_view>
#include <uni20/core/scalar_io.hpp>
#include <uni20/linalg/ops/self_adjoint_eigh.hpp>

namespace
{
namespace model = bethe::hubbard;
void require(bool condition, std::string_view message)
{
  if (!condition) throw std::runtime_error(std::string(message));
}
template <typename Exception = std::invalid_argument, typename Function> void rejects(Function&& f)
{
  try
  {
    f();
  }
  catch (Exception const&)
  {
    return;
  }
  throw std::runtime_error("expected exception");
}

using bethe::test::exact_ground;

// All original Lieb-Wu equations, not the reduced implementation or Jacobian.
template <uni20::Real Real> void check_state(model::State<Real> const& s)
{
  using std::abs;
  using std::atan;
  using std::cos;
  using std::sin;
  auto const sites = s.sites, n = s.root_particles, m = s.root_down_spins;
  Real const pi = Real{4} * atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  auto const& k = s.charge_momenta;
  auto const& l = s.spin_rapidities;
  require(k.size() == n, "Hubbard root particle metadata");
  Real energy = s.energy_offset, momentum = Real{2} * pi * Real(s.momentum_offset) / Real(sites);
  for (Real value : k)
  {
    energy -= Real{2} * cos(value);
    momentum += value;
  }
  Real const allowance = Real{64} * Real(sites) * eps * (Real{1} + abs(s.energy_offset));
  require(abs(energy - s.energy) < allowance, "Hubbard energy from returned momenta");
  require(abs(cos(momentum) - cos(s.momentum)) < allowance && abs(sin(momentum) - sin(s.momentum)) < allowance,
          "Hubbard momentum from returned roots");
  require(s.converged == (s.status == model::SolveStatus::converged), "consistent solve status");
  if (s.free_fermion)
  {
    require((s.root_interaction == Real{0} || m == 0 || m == n) && l.empty() && s.quantum_numbers.charge.empty() &&
                s.quantum_numbers.spin.empty(),
            "free fermions do not have interacting Bethe labels or Lambda");
    require(s.converged && s.iterations == 0 && s.continuation_steps == 0 && s.residual_norm == Real{0},
            "exact free-fermion path does not iterate");
    return;
  }
  require(l.size() == m && s.quantum_numbers.charge.size() == n && s.quantum_numbers.spin.size() == m,
          "nested root and label counts");
  Real const u = s.root_interaction / Real{4};
  std::vector<Real> sine(n);
  for (std::size_t j = 0; j < n; ++j)
    sine[j] = k[j] == Real{0} || k[j] == pi ? Real{0} : sin(k[j]);
  Real rc = Real{0}, rs = Real{0};
  std::int64_t twice_labels = 0;
  for (std::size_t j = 0; j < n; ++j)
  {
    auto const label = s.quantum_numbers.charge[j].twice();
    require((label % 2 == 0) == (m % 2 == 0), "charge integer/half-integer parity");
    if (j)
      require(k[j - 1] < k[j] && s.quantum_numbers.charge[j] - s.quantum_numbers.charge[j - 1] == uni20::half_int{1},
              "ordered distinct charge roots and consecutive labels");
    Real f = Real(sites) * k[j] - pi * Real(label);
    for (Real lambda : l)
      f += Real{2} * atan((sine[j] - lambda) / u);
    rc = std::max(rc, abs(f) / Real(sites));
    twice_labels += label;
  }
  for (std::size_t a = 0; a < m; ++a)
  {
    auto const label = s.quantum_numbers.spin[a].twice();
    require((label % 2 == 0) == ((n - m) % 2 != 0), "spin integer/half-integer parity");
    if (a)
      require(l[a - 1] < l[a] && s.quantum_numbers.spin[a] - s.quantum_numbers.spin[a - 1] == uni20::half_int{1},
              "ordered distinct spin roots and consecutive labels");
    Real f = -pi * Real(label);
    for (Real value : sine)
      f += Real{2} * atan((l[a] - value) / u);
    for (Real lambda : l)
      f -= Real{2} * atan((l[a] - lambda) / (Real{2} * u));
    rs = std::max(rs, abs(f) / Real(sites));
    twice_labels += label;
  }
  require(abs(rc - s.charge_residual) < allowance && abs(rs - s.spin_residual) < allowance,
          "independent full charge/spin residuals at requested U");
  require(s.residual_norm == std::max(s.charge_residual, s.spin_residual), "nested max norm");
  require((static_cast<std::size_t>(twice_labels / 2) + s.momentum_offset) % sites == s.momentum_index,
          "exact label momentum");
}

template <uni20::Real Real> void jacobian()
{
  using std::abs;
  using std::cbrt;
  Real const step = cbrt(uni20::numeric_limits<Real>::epsilon());
  for (auto [sites, particles, down] : {std::array<std::size_t, 3>{2, 2, 1},
                                        {4, 4, 2},
                                        {6, 6, 3},
                                        {8, 8, 4},
                                        {10, 10, 5},
                                        {8, 8, 2},
                                        {8, 8, 3},
                                        {8, 2, 1},
                                        {10, 6, 1},
                                        {10, 6, 3}})
    for (Real u : {Real{1} / Real{4}, Real{1}, Real{4}})
    {
      model::detail::GroundSystem<Real> system(sites, particles, down);
      auto x = system.seed();
      for (std::size_t j = 0; j < x.size(); ++j)
        x[j] *= Real{9} / Real{10};
      uni20::DenseMatrix<Real> a(system.order, system.order);
      (void)system.evaluate(x, u, &a);
      for (std::size_t col = 0; col < system.order; ++col)
      {
        auto plus = x, minus = x;
        plus[col] += step;
        minus[col] -= step;
        auto fp = system.evaluate(plus, u), fm = system.evaluate(minus, u);
        for (std::size_t row = 0; row < system.order; ++row)
        {
          Real const numerical = (fp.residual[row] - fm.residual[row]) / (Real{2} * step);
          require(abs(a[row, col] - numerical) < Real{1000} * step * step * (Real{1} + abs(numerical)),
                  "analytic reduced Hubbard Jacobian vs central differences");
        }
      }
    }
}

template <uni20::Real Real> void tests()
{
  using std::abs;
  using std::atan;
  using std::cos;
  using std::sqrt;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * atan(Real{1});
  jacobian<Real>();
  model::detail::GroundSystem<Real> overflow_system(4);
  auto overflow_roots = overflow_system.seed();
  overflow_roots[overflow_system.nk] = (uni20::numeric_limits<Real>::max() / Real{4}) * Real{3};
  rejects<std::overflow_error>([&] { (void)overflow_system.evaluate(overflow_roots, Real{1}); });
  // ED is an independent fp64 oracle; high-precision accuracy is tested below
  // against native-precision analytic values, not by promoting a double result.
  for (unsigned n : {2, 4, 6})
    for (Real u : {Real{0}, Real{1} / Real{10}, Real{1}, Real{4}, Real{8}, Real{100}})
    {
      auto const state = model::ground_state<Real>(n, u);
      auto const ed = exact_ground(n, static_cast<double>(u));
      require(state.converged && std::abs(static_cast<double>(state.energy) - ed.energy) < 2e-11,
              "Hubbard half-filled ground energy vs fermionic ED");
      if (u != Real{0})
        require(std::abs(ed.translation - (n % 4 == 0 ? -1.0 : 1.0)) < 2e-11,
                "Hubbard ground momentum vs fermionic ED translation");
      check_state(state);
    }
  for (std::size_t n : {2, 4, 6, 8, 16, 32})
    for (Real u : {Real{0}, Real{1} / Real{100}, Real{1}, Real{4}, Real{8}, Real{100}})
    {
      auto const s = model::ground_state<Real>(n, u);
      require(s.converged && s.residual_norm <= Real{32} * eps,
              "Hubbard default convergence across both parity classes");
      std::size_t stages = u == Real{0} ? 0 : 1;
      for (Real current = Real{8}; u != Real{0} && current > u; current /= Real{2})
        ++stages;
      require(s.continuation_steps == stages, "continuation reaches requested U");
      check_state(s);
    }
  for (Real u : {Real{0}, Real{1} / Real{100}, Real{4}, Real{100}})
  {
    auto const state = model::ground_state<Real>(2, u, {.residual_tolerance = Real{2} * eps});
    Real const exact = -Real{32} / (u + sqrt(u * u + Real{64}));
    require(state.converged && abs(state.energy - exact) < Real{32} * eps,
            "native-precision L=2 irrational energy (doubled periodic hopping)");
    require(uni20::parse_real<Real>(uni20::format_real(state.energy)) == state.energy, "energy round-trip scalar I/O");
  }
  Real const parsed = uni20::parse_real<Real>("4.12345678901234567890123456789");
  auto const parsed_state = model::ground_state<Real>(2, parsed, {.residual_tolerance = Real{2} * eps});
  Real const exact = -Real{32} / (parsed + sqrt(parsed * parsed + Real{64}));
  require(parsed_state.interaction == parsed && abs(parsed_state.energy - exact) < Real{32} * eps,
          "U parsing and solve keep selected arithmetic");
  if constexpr (uni20::numeric_limits<Real>::digits > 53)
    require(parsed != Real(static_cast<double>(parsed)), "high precision U has digits beyond double");

  for (std::size_t n : {2, 4, 6, 8, 10, 32})
  {
    std::vector<Real> levels;
    for (std::size_t q = 0; q < n; ++q)
      levels.push_back(-Real{2} * cos(Real{2} * pi * Real(q) / Real(n)));
    std::sort(levels.begin(), levels.end());
    Real free = Real{0};
    for (std::size_t j = 0; j < n / 2; ++j)
      free += Real{2} * levels[j];
    auto const zero = model::ground_state<Real>(n, Real{0}, {.max_iterations = 0});
    require(abs(zero.energy - free) < Real{32} * Real(n) * eps, "U=0 filled free-fermion spectrum");
    check_state(zero);
    auto const weak = model::ground_state<Real>(n, Real{1} / Real{1000000});
    require(weak.converged && weak.energy > free && weak.energy - free < Real(n) / Real{1000000},
            "weak repulsion approaches free-fermion energy from above");
    check_state(weak);
  }
  for (std::size_t n : {4, 6, 8, 16})
  {
    auto const xxx = bethe::heisenberg::ground_state<Real>(n);
    require(xxx.converged, "XXX strong-coupling reference converges");
    Real const u = Real{10000};
    auto const strong = model::ground_state<Real>(n, u);
    require(strong.converged &&
                abs(u * strong.energy / Real{4} - (xxx.energy - Real(n) / Real{4})) < Real(n) / Real{1000000},
            "Hubbard strong coupling includes XXX -L/4 shift");
    check_state(strong);
  }
  for (std::size_t n : {4, 6, 16})
    for (std::size_t budget : {0, 1, 5, 10})
    {
      auto const state = model::ground_state<Real>(n, Real{1} / Real{100}, {.max_iterations = budget});
      require(!state.converged && state.status == model::SolveStatus::iteration_limit && state.iterations == budget,
              "global continuation budget and explicit nonconvergence");
      check_state(state);
    }
  auto const stalled = model::ground_state<Real>(8, eps * eps * eps * eps);
  require(!stalled.converged && stalled.status == model::SolveStatus::stalled,
          "unresolved weak-coupling roots return a stalled state");
  check_state(stalled);
  for (std::size_t n : {0, 1, 3, 5})
    rejects([&] { (void)model::ground_state<Real>(n, Real{4}); });
  rejects([&] { (void)model::ground_quantum_numbers(std::numeric_limits<std::size_t>::max()); });
  rejects([&] {
    (void)model::ground_state<Real>(static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() / 4) - 1,
                                    Real{4});
  });
  for (Real u : {uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    rejects([&] { (void)model::ground_state<Real>(4, u); });
  for (Real t : {Real{0}, Real{-1}, uni20::numeric_limits<Real>::infinity(), uni20::numeric_limits<Real>::quiet_NaN()})
    rejects([&] { (void)model::ground_state<Real>(4, Real{4}, {.residual_tolerance = t}); });
  rejects([&] { (void)model::ground_state<Real>(4, uni20::numeric_limits<Real>::min() / Real{16}); });
  // All physical spin populations on small even rings, including both signs
  // of U, above half filling, empty/filled bands and degenerate free shells.
  for (unsigned sites : {2, 4, 6})
    for (unsigned up = 0; up <= sites; ++up)
      for (unsigned down = 0; down <= sites; ++down)
        for (Real interaction :
             {Real{0}, Real{1} / Real{10}, Real{4}, Real{100}, Real{-1} / Real{10}, Real{-4}, Real{-100}})
        {
          auto root_up = up, root_down = interaction < Real{0} ? sites - down : down;
          if (root_up + root_down > sites)
          {
            root_up = sites - root_up;
            root_down = sites - root_down;
          }
          auto const m = std::min(root_up, root_down), n = root_up + root_down;
          bool const supported =
              interaction == Real{0} || up == 0 || down == 0 || m == 0 || n == sites || (n % 2 == 0 && m % 2 == 1);
          auto solve = [&] {
            return model::sector_ground_state<Real>(
                sites, up + down, uni20::from_twice(static_cast<std::int64_t>(up) - down), interaction);
          };
          if (!supported)
          {
            rejects(solve);
            continue;
          }
          auto const state = solve();
          auto const ed = exact_ground(sites, static_cast<double>(interaction), up, down);
          if (!state.converged || std::abs(static_cast<double>(state.energy) - ed.energy) > 3e-10 ||
              ed.momentum_weights[state.momentum_index] < 0.9)
          {
            std::cerr << "sector L=" << sites << " up=" << up << " down=" << down
                      << " U=" << uni20::format_real(interaction) << " E=" << uni20::format_real(state.energy)
                      << " ED=" << ed.energy << " P=" << state.momentum_index
                      << " weight=" << ed.momentum_weights[state.momentum_index] << '\n';
            require(false, "Hubbard sector energy and momentum vs independent ED");
          }
          require(state.particles == up + down && state.down_spins == down && state.interaction == interaction,
                  "physical sector metadata survives mappings");
          check_state(state);
        }

  // Native-precision mapping identities, including balanced attraction away
  // from half filling and both signs of physical Sz.
  for (std::size_t pairs = 1; pairs < 16; ++pairs)
  {
    Real const u = Real{4};
    auto const attractive = model::sector_ground_state<Real>(16, 2 * pairs, uni20::half_int{0}, -u);
    auto const repulsive =
        model::sector_ground_state<Real>(16, 16, uni20::from_twice(2 * static_cast<std::int64_t>(pairs) - 16), u);
    require(attractive.converged && repulsive.converged && attractive.shiba_transformed &&
                attractive.root_particles == 16 && attractive.root_down_spins == std::min(pairs, 16 - pairs) &&
                attractive.root_interaction == u && attractive.energy_offset == -u * Real(pairs) &&
                abs(attractive.energy - (repulsive.energy - u * Real(pairs))) < Real{256} * eps,
            "balanced attraction maps to the correct half-filled spin sector in native precision");
    check_state(attractive);
    check_state(repulsive);
  }
  for (Real u : {Real{1} / Real{100}, Real{4}, Real{100}})
  {
    auto const state = model::ground_state<Real>(2, -u, {.residual_tolerance = Real{2} * eps});
    Real const exact = -(u + sqrt(u * u + Real{64})) / Real{2};
    require(state.converged && abs(state.energy - exact) < Real{32} * eps * (Real{1} + u),
            "native-precision attractive dimer energy");
  }
  for (auto [sites, particles, down] : {std::array<std::size_t, 3>{16, 6, 3}, {32, 10, 1}, {32, 14, 5}, {16, 16, 2}})
    for (Real interaction : {Real{1} / Real{1000000}, Real{4}, Real{100}})
    {
      auto const sz = uni20::from_twice(static_cast<std::int64_t>(particles - 2 * down));
      auto const state = model::sector_ground_state<Real>(sites, particles, sz, interaction);
      require(state.converged, "doped/polarized convergence beyond ED sizes");
      check_state(state);
      if (interaction < Real{1})
      {
        auto const free = model::sector_ground_state<Real>(sites, particles, sz, Real{0});
        require(state.energy > free.energy && state.energy - free.energy < Real(sites) * interaction,
                "doped/polarized weak coupling approaches independently filled free bands");
      }
      for (std::size_t budget : {0, 1})
      {
        auto const unfinished =
            model::sector_ground_state<Real>(sites, particles, sz, interaction, {.max_iterations = budget});
        require(!unfinished.converged && unfinished.iterations == budget, "doped/polarized global update budget");
        check_state(unfinished);
      }
    }
  auto const mapped_unfinished =
      model::sector_ground_state<Real>(16, 10, uni20::half_int{0}, -Real{1} / Real{100}, {.max_iterations = 1});
  require(!mapped_unfinished.converged && mapped_unfinished.auxiliary_roots() &&
              mapped_unfinished.root_interaction == Real{1} / Real{100},
          "mapped failure retains final root U");
  check_state(mapped_unfinished);
  rejects([&] { (void)model::sector_ground_state<Real>(4, 9, uni20::half_int{0}, Real{4}); });
  rejects([&] { (void)model::sector_ground_state<Real>(4, 3, uni20::half_int{0}, Real{4}); });
  rejects([&] { (void)model::sector_ground_state<Real>(4, 4, uni20::half_int{3}, Real{4}); });
  rejects([&] { (void)model::sector_ground_state<Real>(4, 4, uni20::half_int{-3}, Real{4}); });
  rejects<std::overflow_error>(
      [&] { (void)model::sector_ground_state<Real>(4, 8, uni20::half_int{0}, -uni20::numeric_limits<Real>::max()); });
  std::cout << "Hubbard ED, momentum, Jacobian, native precision, limits and failure checks passed\n";
}
} // namespace
int main(int argc, char** argv)
{
  try
  {
    if (argc != 2) throw std::invalid_argument("pass a precision");
    std::string_view const p = argv[1];
    if (p == "fp64")
      tests<double>();
    else if (p == "long-double")
      tests<long double>();
#if UNI20_HAS_FLOAT128
    else if (p == "fp128")
      tests<uni20::float128>();
#endif
    else
      throw std::invalid_argument("unavailable precision");
    return 0;
  }
  catch (std::exception const& e)
  {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
