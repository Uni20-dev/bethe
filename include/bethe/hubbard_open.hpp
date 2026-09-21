// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/hubbard_common.hpp>

namespace bethe::hubbard::open
{
using hubbard::NestedQuantumNumbers;
using hubbard::QuantumNumbers;
using hubbard::SolveStatus;
template <uni20::Real Real> using SolverOptions = hubbard::SolverOptions<Real>;

/// Free-end Hubbard sector ground state, t=1, L>=1, either sign of U.
/// charge_momenta contains standing-wave k, not conserved lattice momentum.
/// Residuals are normalized by 2*(L+1); there are no momentum members.
template <uni20::Real Real> struct State : hubbard::StateData<Real>
{};

/// Positive consecutive integer labels in a normalized repulsive root sector.
[[nodiscard]] inline NestedQuantumNumbers ground_quantum_numbers(std::size_t sites, std::size_t particles,
                                                                 std::size_t down)
{
  hubbard::detail::checked_length(sites);
  if (particles > sites || down > particles / 2)
    throw std::invalid_argument("open Hubbard root sector requires N <= L and 0 <= M <= N/2");
  NestedQuantumNumbers result;
  for (std::size_t j = 0; j < particles; ++j)
    result.charge.emplace_back(static_cast<std::int64_t>(j + 1));
  for (std::size_t a = 0; a < down; ++a)
    result.spin.emplace_back(static_cast<std::int64_t>(a + 1));
  return result;
}

[[nodiscard]] inline NestedQuantumNumbers ground_quantum_numbers(std::size_t sites)
{
  return ground_quantum_numbers(sites, sites, sites / 2);
}

namespace detail
{
/// Free-end logarithmic equations, with both direct and reflected scattering.
/// All k and Lambda are positive; variables are k and Lambda/max(1,U/4).
template <uni20::Real Real> class GroundSystem {
  public:
    GroundSystem(std::size_t sites, std::size_t particles, std::size_t down)
        : sites(sites), n(particles), m(down), nk(n), ns(m), order(n + m), pi(Real{4} * std::atan(Real{1})),
          labels(open::ground_quantum_numbers(sites, particles, down))
    {}

    std::vector<Real> seed() const
    {
      using std::tan;
      std::vector<Real> x(order);
      for (std::size_t j = 0; j < n; ++j)
        x[j] = pi * Real(j + 1) / Real(sites + 1);
      for (std::size_t a = 0; a < m; ++a)
        x[n + a] = tan(pi * Real(a + 1) / (Real{2} * Real(n)));
      return x;
    }

    bool physical(std::vector<Real> const& x) const
    {
      for (std::size_t i = 0; i < order; ++i)
        if (!uni20::isfinite(x[i]) || x[i] <= Real{0} || (i < n && x[i] >= pi) ||
            (i != 0 && i != n && x[i - 1] >= x[i]))
          return false;
      return true;
    }

    void expand(std::vector<Real> const& x, Real u, std::vector<Real>& k, std::vector<Real>& lambda) const
    {
      k.assign(x.begin(), x.begin() + n);
      lambda.resize(m);
      Real const scale = std::max(Real{1}, u);
      for (std::size_t a = 0; a < m; ++a)
        lambda[a] = scale * x[n + a];
    }

    struct Evaluation
    {
        std::vector<Real> residual;
        Real charge = Real{0}, spin = Real{0};
        Real norm() const { return std::max(charge, spin); }
    };

    Evaluation evaluate(std::vector<Real> const& x, Real u, uni20::DenseMatrix<Real>* jacobian = nullptr) const
    {
      using hubbard::detail::kernel;
      using std::abs;
      using std::atan;
      using std::cos;
      using std::sin;
      Real const scale = std::max(Real{1}, u), weight = Real{1} / Real(sites + 1);
      std::vector<Real> sine(n), cosine(n), lambda(m);
      for (std::size_t j = 0; j < n; ++j)
      {
        sine[j] = sin(x[j]);
        cosine[j] = cos(x[j]);
      }
      for (std::size_t a = 0; a < m; ++a)
      {
        lambda[a] = scale * x[n + a];
        if (!uni20::isfinite(lambda[a]) || abs(lambda[a]) > uni20::numeric_limits<Real>::max() / Real{2})
          throw std::overflow_error("Hubbard rapidities exceed the representable scattering range");
      }
      if (jacobian)
        for (std::size_t i = 0; i < order; ++i)
          for (std::size_t j = 0; j < order; ++j)
            (*jacobian)[i, j] = Real{0};
      Evaluation result{.residual = std::vector<Real>(order)};
      for (std::size_t j = 0; j < n; ++j)
      {
        bethe::detail::CompensatedSum<Real> phases;
        if (jacobian) (*jacobian)[j, j] = Real{1};
        for (std::size_t a = 0; a < m; ++a)
        {
          Real const minus = sine[j] - lambda[a], plus = sine[j] + lambda[a];
          phases.add(atan(minus / u));
          phases.add(atan(plus / u));
          if (jacobian)
          {
            Real const km = weight * kernel(minus, u), kp = weight * kernel(plus, u);
            (*jacobian)[j, j] += cosine[j] * (km + kp);
            (*jacobian)[j, n + a] = scale * (kp - km);
          }
        }
        result.residual[j] = x[j] - pi * Real(j + 1) * weight + weight * phases.value();
        result.charge = std::max(result.charge, abs(result.residual[j]));
      }
      for (std::size_t a = 0; a < m; ++a)
      {
        auto const row = n + a;
        bethe::detail::CompensatedSum<Real> phases;
        for (std::size_t j = 0; j < n; ++j)
        {
          Real const minus = lambda[a] - sine[j], plus = lambda[a] + sine[j];
          phases.add(atan(minus / u));
          phases.add(atan(plus / u));
          if (jacobian)
          {
            Real const km = weight * kernel(minus, u), kp = weight * kernel(plus, u);
            (*jacobian)[row, row] += scale * (km + kp);
            (*jacobian)[row, j] = cosine[j] * (kp - km);
          }
        }
        for (std::size_t b = 0; b < m; ++b)
        {
          // Both factors exclude b=a: reflected self-scattering is absent too.
          if (b == a) continue;
          Real const minus = lambda[a] - lambda[b], plus = lambda[a] + lambda[b];
          phases.add(-atan(minus / (Real{2} * u)));
          phases.add(-atan(plus / (Real{2} * u)));
          if (jacobian)
          {
            Real const km = weight * scale * kernel(minus, Real{2} * u);
            Real const kp = weight * scale * kernel(plus, Real{2} * u);
            (*jacobian)[row, row] -= km + kp;
            (*jacobian)[row, n + b] = km - kp;
          }
        }
        result.residual[row] = weight * (phases.value() - pi * Real(a + 1));
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
};

template <uni20::Real Real>
State<Real> free_state(std::size_t sites, std::size_t particles, std::size_t down, Real interaction)
{
  using std::cos;
  State<Real> result;
  result.sites = sites;
  result.particles = result.root_particles = particles;
  result.down_spins = result.root_down_spins = down;
  result.interaction = result.root_interaction = interaction;
  result.free_fermion = result.converged = true;
  result.status = SolveStatus::converged;
  Real const pi = Real{4} * std::atan(Real{1});
  bethe::detail::CompensatedSum<Real> energy;
  for (auto count : {particles - down, down})
    for (std::size_t j = 0; j < count; ++j)
    {
      Real const k = pi * Real(j + 1) / Real(sites + 1);
      result.charge_momenta.push_back(k);
      // A filled band has exactly zero kinetic energy, including L=1.
      if (count != sites) energy.add(-Real{2} * cos(k));
    }
  result.energy = energy.value();
  return result;
}
} // namespace detail

/// Ground state in any physical (N,Sz) sector, with free ends and no boundary fields.
/// Odd and even L are bipartite. Mapped roots/residuals are explicitly auxiliary.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> sector_ground_state(std::size_t sites, std::size_t particles, uni20::half_int sz,
                                              Real interaction, SolverOptions<Real> const& options = {})
{
  auto const down = hubbard::detail::validate_sector(sites, particles, sz, interaction, options);
  if (interaction == Real{0} || down == 0 || down == particles)
    return detail::free_state(sites, particles, down, interaction);
  auto const mapping = hubbard::detail::map_sector(sites, particles, down, interaction);
  auto const n = mapping.up + mapping.down, m = mapping.down;
  State<Real> result;
  if (m == 0)
    result = detail::free_state(sites, n, m, mapping.interaction);
  else
  {
    hubbard::detail::check_newton_size<Real>(n + m);
    detail::GroundSystem<Real> const system(sites, n, m);
    result = hubbard::detail::solve_ground_system<Real, State<Real>>(system, mapping.interaction, options);
  }
  hubbard::detail::apply_mapping(result, particles, down, interaction, mapping);
  return result;
}

/// Half filling; choose Sz=0 for even L and Sz=1/2 for odd L.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> ground_state(std::size_t sites, Real interaction, SolverOptions<Real> const& options = {})
{
  return sector_ground_state(sites, sites, uni20::from_twice(static_cast<std::int64_t>(sites % 2)), interaction,
                             options);
}
} // namespace bethe::hubbard::open
