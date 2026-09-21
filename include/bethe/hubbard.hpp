// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/hubbard_common.hpp>

namespace bethe::hubbard
{
/// Periodic sector ground state, t=1, even L>=2 (supported sectors below).
/// H=-sum_(j,sigma)(c^dagger_j c_(j+1)+h.c.)+U sum_j n_up n_down.
template <uni20::Real Real> struct State : StateData<Real>
{
    std::size_t momentum_index = 0, momentum_offset = 0;
    Real momentum = Real{0};
};

namespace detail
{
inline std::int64_t checked_sites(std::size_t sites)
{
  if (sites < 2 || sites % 2 != 0 || sites > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() / 4))
    throw std::invalid_argument("Hubbard ground state requires even 2 <= sites <= INT64_MAX/4");
  return static_cast<std::int64_t>(sites);
}

inline void checked_root_sector(std::size_t sites, std::size_t particles, std::size_t down)
{
  checked_sites(sites);
  if (particles > sites || down > particles / 2)
    throw std::invalid_argument("repulsive root sector requires N <= L and 0 <= M <= N/2");
  if (down != 0 && particles != sites && (particles % 2 != 0 || down % 2 == 0))
    throw std::invalid_argument(
        "unsupported Hubbard sector: after symmetry mapping, interacting doped sectors require odd "
        "N_up and N_down; other shell parities need competing spin branches (not implemented)");
}
} // namespace detail

/// Labels in a supported normalized repulsive sector. At half filling, charge
/// labels fill a Brillouin zone (the +pi endpoint for even M). Doped interacting
/// sectors have even N, odd M and centered half-odd I / integer J. Spin labels
/// are centered consecutively. These are auxiliary labels if a mapping is used.
[[nodiscard]] inline NestedQuantumNumbers ground_quantum_numbers(std::size_t sites, std::size_t particles,
                                                                 std::size_t down)
{
  detail::checked_root_sector(sites, particles, down);
  auto const n = static_cast<std::int64_t>(particles), m = static_cast<std::int64_t>(down);
  NestedQuantumNumbers result;
  result.charge.reserve(particles);
  result.spin.reserve(down);
  for (std::int64_t j = 0; j < n; ++j)
    result.charge.push_back(uni20::from_twice(2 * j - n + (m % 2 == 0 && n % 2 == 0 ? 2 : 1)));
  for (std::int64_t j = 0; j < m; ++j)
    result.spin.push_back(uni20::from_twice(2 * j - m + 1));
  return result;
}

[[nodiscard]] inline NestedQuantumNumbers ground_quantum_numbers(std::size_t sites)
{
  return ground_quantum_numbers(sites, sites, sites / 2);
}

namespace detail
{
/// Reflection-symmetric ground-state equations. Only positive k and Lambda
/// vary; k=0,pi (half filling, even M) or Lambda=0 (odd M) are held exactly fixed.
/// Spin variables are Lambda/max(1,U/4), to condition the large-U solve.
template <uni20::Real Real> class GroundSystem {
  public:
    explicit GroundSystem(std::size_t sites) : GroundSystem(sites, sites, sites / 2) {}
    GroundSystem(std::size_t sites, std::size_t particles, std::size_t down)
        : sites(sites), n(particles), m(down), nk(n / 2 - (m % 2 == 0 ? 1 : 0)), ns(m / 2), order(nk + ns),
          pi(Real{4} * std::atan(Real{1})), labels(ground_quantum_numbers(sites, particles, down))
    {}

    std::vector<Real> seed() const
    {
      using std::tan;
      std::vector<Real> x(order);
      for (std::size_t j = 0; j < nk; ++j)
        x[j] = Real{2} * pi * charge_label(j) / Real(sites);
      for (std::size_t a = 0; a < ns; ++a)
        x[nk + a] = tan(pi * spin_label(a) / Real(n));
      return x;
    }
    bool physical(std::vector<Real> const& x) const
    {
      for (std::size_t i = 0; i < order; ++i)
        if (!uni20::isfinite(x[i]) || x[i] <= Real{0} || (i < nk && x[i] >= pi) ||
            (i != 0 && i != nk && x[i - 1] >= x[i]))
          return false;
      return true;
    }
    void expand(std::vector<Real> const& x, Real u, std::vector<Real>& k, std::vector<Real>& lambda) const
    {
      k.clear();
      lambda.clear();
      for (std::size_t j = nk; j-- > 0;)
        k.push_back(-x[j]);
      if (m % 2 == 0) k.push_back(Real{0});
      for (std::size_t j = 0; j < nk; ++j)
        k.push_back(x[j]);
      if (m % 2 == 0) k.push_back(pi);
      Real const scale = std::max(Real{1}, u);
      for (std::size_t a = ns; a-- > 0;)
        lambda.push_back(-scale * x[nk + a]);
      if (m % 2 != 0) lambda.push_back(Real{0});
      for (std::size_t a = 0; a < ns; ++a)
        lambda.push_back(scale * x[nk + a]);
    }
    struct Evaluation
    {
        std::vector<Real> residual;
        Real charge = Real{0}, spin = Real{0};
        Real norm() const { return std::max(charge, spin); }
    };
    Evaluation evaluate(std::vector<Real> const& x, Real u, uni20::DenseMatrix<Real>* jacobian = nullptr) const
    {
      using std::abs;
      using std::atan;
      using std::cos;
      using std::sin;
      std::vector<Real> k, lambda;
      expand(x, u, k, lambda);
      // Reflection creates both +/-Lambda: keep their differences representable
      // instead of letting overflow masquerade as a saturated scattering phase.
      for (Real value : lambda)
        if (!uni20::isfinite(value) || abs(value) > uni20::numeric_limits<Real>::max() / Real{2})
          throw std::overflow_error("Hubbard rapidities exceed the representable scattering range");
      std::vector<Real> sine(n), cosine(n);
      for (std::size_t j = 0; j < n; ++j)
      {
        // The fixed points are exact, not floating sin(pi) divided by small U.
        sine[j] = k_variable(j) < 0 ? Real{0} : sin(k[j]);
        cosine[j] = cos(k[j]);
      }
      if (jacobian)
        for (std::size_t i = 0; i < order; ++i)
          for (std::size_t j = 0; j < order; ++j)
            (*jacobian)[i, j] = Real{0};
      Real const scale = std::max(Real{1}, u), weight = Real{2} / Real(sites);
      Evaluation result{.residual = std::vector<Real>(order)};
      for (std::size_t i = 0; i < nk; ++i)
      {
        Real const s = sin(x[i]), c = cos(x[i]);
        bethe::detail::CompensatedSum<Real> phases;
        if (jacobian) (*jacobian)[i, i] = Real{1};
        for (std::size_t b = 0; b < m; ++b)
        {
          Real const d = s - lambda[b];
          phases.add(atan(d / u));
          if (jacobian)
          {
            Real const derivative = weight * kernel(d, u);
            (*jacobian)[i, i] += c * derivative;
            auto const variable = spin_variable(b);
            if (variable >= 0) (*jacobian)[i, nk + variable] -= spin_sign(b) * scale * derivative;
          }
        }
        result.residual[i] = x[i] - Real{2} * pi * charge_label(i) / Real(sites) + weight * phases.value();
        result.charge = std::max(result.charge, abs(result.residual[i]));
      }
      for (std::size_t a = 0; a < ns; ++a)
      {
        auto const row = nk + a;
        Real const l = scale * x[row];
        bethe::detail::CompensatedSum<Real> phases;
        for (std::size_t j = 0; j < n; ++j)
        {
          Real const d = l - sine[j];
          phases.add(atan(d / u));
          if (jacobian)
          {
            Real const derivative = weight * kernel(d, u);
            (*jacobian)[row, row] += scale * derivative;
            auto const variable = k_variable(j);
            if (variable >= 0) (*jacobian)[row, variable] -= k_sign(j) * cosine[j] * derivative;
          }
        }
        for (std::size_t b = 0; b < m; ++b)
        {
          Real const d = l - lambda[b];
          phases.add(-atan(d / (Real{2} * u)));
          if (jacobian)
          {
            Real const derivative = weight * scale * kernel(d, Real{2} * u);
            (*jacobian)[row, row] -= derivative;
            auto const variable = spin_variable(b);
            if (variable >= 0) (*jacobian)[row, nk + variable] += spin_sign(b) * derivative;
          }
        }
        result.residual[row] = weight * (phases.value() - pi * spin_label(a));
        result.spin = std::max(result.spin, abs(result.residual[row]));
      }
      for (Real f : result.residual)
        if (!uni20::isfinite(f)) throw std::runtime_error("nonfinite Hubbard equation residual");
      if (jacobian)
        for (std::size_t i = 0; i < order; ++i)
          for (std::size_t j = 0; j < order; ++j)
            if (!uni20::isfinite((*jacobian)[i, j])) throw std::runtime_error("nonfinite Hubbard Jacobian");
      return result;
    }
    std::size_t const sites, n, m, nk, ns, order;
    Real const pi;
    NestedQuantumNumbers const labels;

  private:
    Real charge_label(std::size_t i) const { return Real(i) + (m % 2 == 0 ? Real{1} : Real{1} / Real{2}); }
    Real spin_label(std::size_t i) const { return Real(i) + (m % 2 == 0 ? Real{1} / Real{2} : Real{1}); }
    std::int64_t k_variable(std::size_t j) const
    {
      if (j < nk) return static_cast<std::int64_t>(nk - 1 - j);
      auto const first = nk + (m % 2 == 0 ? 1 : 0);
      if (j >= first && j < first + nk) return static_cast<std::int64_t>(j - first);
      return -1;
    }
    Real k_sign(std::size_t j) const { return j < nk ? Real{-1} : Real{1}; }
    std::int64_t spin_variable(std::size_t b) const
    {
      if (b < ns) return static_cast<std::int64_t>(ns - 1 - b);
      auto const first = ns + m % 2;
      if (b >= first) return static_cast<std::int64_t>(b - first);
      return -1;
    }
    Real spin_sign(std::size_t b) const { return b < ns ? Real{-1} : Real{1}; }
};

template <uni20::Real Real>
State<Real> free_state(std::size_t sites, std::size_t particles, std::size_t down, Real interaction)
{
  using std::atan;
  using std::cos;
  using std::sin;
  State<Real> result;
  result.sites = sites;
  result.particles = result.root_particles = particles;
  result.down_spins = result.root_down_spins = down;
  result.interaction = result.root_interaction = interaction;
  result.free_fermion = result.converged = true;
  result.status = SolveStatus::converged;
  Real const pi = Real{4} * atan(Real{1}), angle = pi / Real(sites);
  std::int64_t momentum = 0, length = static_cast<std::int64_t>(sites);
  for (auto count : {particles - down, down})
  {
    if (count == 0) continue;
    auto const first = -(static_cast<std::int64_t>(count) - 1) / 2;
    for (std::size_t j = 0; j < count; ++j)
    {
      auto const q = first + static_cast<std::int64_t>(j);
      result.charge_momenta.push_back(Real{2} * pi * Real(q) / Real(sites));
      momentum = (momentum + q) % length;
    }
    // A filled band has exactly zero kinetic energy; do not evaluate sin(pi).
    if (count != sites)
      result.energy -= Real{2} * sin(Real(count) * angle) / sin(angle) * (count % 2 == 0 ? cos(angle) : Real{1});
  }
  result.momentum_index = static_cast<std::size_t>((momentum + length) % length);
  result.momentum = Real{2} * pi * Real(result.momentum_index) / Real(sites);
  return result;
}

template <uni20::Real Real>
State<Real> repulsive_ground_state(std::size_t sites, std::size_t particles, std::size_t down, Real interaction,
                                   SolverOptions<Real> const& options)
{
  checked_root_sector(sites, particles, down);
  if (down == 0) return free_state(sites, particles, down, interaction);
  auto const order = particles / 2 - (down % 2 == 0 ? 1 : 0) + down / 2;
  check_newton_size<Real>(order);
  GroundSystem<Real> const system(sites, particles, down);
  auto result = solve_ground_system<Real, State<Real>>(system, interaction, options);
  result.momentum_index = down % 2 == 0 ? sites / 2 : 0;
  result.momentum = result.momentum_index == 0 ? Real{0} : Real{4} * std::atan(Real{1});
  return result;
}
} // namespace detail

/// Lowest energy in a supported (N,Sz) sector on an even periodic ring.
/// Normalize U<0 by down-spin particle-hole (Shiba), N>L by full particle-hole,
/// and negative Sz by spin reversal. The repulsive solver covers half filling,
/// doped odd/odd spin populations, and fully polarized sectors. U=0 is unrestricted.
/// Roots/residuals belong to root_* metadata; energy/momentum to the requested sector.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> sector_ground_state(std::size_t sites, std::size_t particles, uni20::half_int sz,
                                              Real interaction, SolverOptions<Real> const& options = {})
{
  detail::checked_sites(sites);
  auto const down = detail::validate_sector(sites, particles, sz, interaction, options);
  if (interaction == Real{0} || down == 0 || down == particles)
    return detail::free_state(sites, particles, down, interaction);
  auto const mapping = detail::map_sector(sites, particles, down, interaction);
  auto result =
      detail::repulsive_ground_state(sites, mapping.up + mapping.down, mapping.down, mapping.interaction, options);
  detail::apply_mapping(result, particles, down, interaction, mapping);
  std::size_t momentum_shift = mapping.shiba && down % 2 == 0 ? sites / 2 : 0;
  auto const before_particle_hole = mapping.shiba ? particles + sites - 2 * down : particles;
  if (mapping.particle_hole && before_particle_hole % 2 != 0) momentum_shift = (momentum_shift + sites / 2) % sites;
  result.momentum_offset = momentum_shift;
  result.momentum_index = (result.momentum_index + momentum_shift) % sites;
  result.momentum = Real{8} * std::atan(Real{1}) * Real(result.momentum_index) / Real(sites);
  return result;
}

/// Convenience wrapper: half filling, Sz=0. Both signs of U are supported.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> ground_state(std::size_t sites, Real interaction, SolverOptions<Real> const& options = {})
{
  return sector_ground_state(sites, sites, uni20::half_int{0}, interaction, options);
}
} // namespace bethe::hubbard
