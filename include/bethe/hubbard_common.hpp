// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/newton_backtracking.hpp>

#include <algorithm>
#include <bethe/solver.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <uni20/common/half_int.hpp>
#include <uni20/linalg/ops/linear_solve.hpp>
#include <utility>
#include <vector>

namespace bethe::hubbard
{
template <uni20::Real Real> using SolverOptions = bethe::SolverOptions<Real>;
using QuantumNumbers = std::vector<uni20::half_int>;
struct NestedQuantumNumbers
{
    QuantumNumbers charge;
    QuantumNumbers spin;
};
enum class SolveStatus
{
  converged,
  iteration_limit,
  stalled
};

/// Boundary-independent observables and root-sector diagnostics. No lattice momentum.
template <uni20::Real Real> struct StateData
{
    std::size_t sites = 0, particles = 0, down_spins = 0;
    Real interaction = Real{0};
    /// Roots and residuals describe this sector, which may be symmetry-mapped.
    std::size_t root_particles = 0, root_down_spins = 0;
    Real root_interaction = Real{0};
    /// E_physical = E_roots + energy_offset.
    Real energy_offset = Real{0};
    bool shiba_transformed = false, particle_hole_transformed = false, spin_reversed = false;
    bool auxiliary_roots() const { return shiba_transformed || particle_hole_transformed || spin_reversed; }
    std::vector<Real> charge_momenta;
    /// Conventional Lieb-Wu Lambda (not the internally scaled Newton variable).
    std::vector<Real> spin_rapidities;
    NestedQuantumNumbers quantum_numbers;
    Real energy = Real{0};
    /// Boundary-normalized equation maxima at root_interaction (not a continuation stage):
    /// divide by L for PBC, by 2*(L+1) for free ends.
    Real charge_residual = Real{0}, spin_residual = Real{0}, residual_norm = Real{0};
    std::size_t iterations = 0, continuation_steps = 0;
    bool converged = false;
    SolveStatus status = SolveStatus::iteration_limit;
    /// Exact free occupations at U=0 or with only one spin species in the root
    /// sector. Bethe labels and spin rapidities are absent in either case.
    bool free_fermion = false;
};

namespace detail
{
inline void checked_length(std::size_t sites)
{
  if (sites < 1 || sites > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() / 4))
    throw std::invalid_argument("Hubbard requires 1 <= sites <= INT64_MAX/4");
}

template <uni20::Real Real>
std::size_t validate_sector(std::size_t sites, std::size_t particles, uni20::half_int sz, Real interaction,
                            SolverOptions<Real> const& options)
{
  checked_length(sites);
  if (!uni20::isfinite(interaction)) throw std::invalid_argument("Hubbard ground state requires finite U");
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");
  if (particles > 2 * sites) throw std::invalid_argument("Hubbard particles must satisfy 0 <= N <= 2L");
  auto const bound = static_cast<std::int64_t>(std::min(particles, 2 * sites - particles)), twice_sz = sz.twice();
  if (twice_sz < -bound || twice_sz > bound || (static_cast<std::int64_t>(particles) - twice_sz) % 2 != 0)
    throw std::invalid_argument("Hubbard N and Sz must give integer N_up and N_down between 0 and L");
  return static_cast<std::size_t>((static_cast<std::int64_t>(particles) - twice_sz) / 2);
}

/// Exact bipartite-lattice mappings. Boundary-specific momentum is not part of this type.
template <uni20::Real Real> struct SectorMapping
{
    std::size_t up, down;
    Real interaction, energy_offset;
    bool shiba, particle_hole, spin_reversed;
};

template <uni20::Real Real>
SectorMapping<Real> map_sector(std::size_t sites, std::size_t particles, std::size_t physical_down, Real interaction)
{
  auto up = particles - physical_down, down = physical_down;
  bool const shiba = interaction < Real{0};
  Real const root_u = shiba ? -interaction : interaction;
  // Combine exact integer coefficients before multiplying by U.
  std::int64_t shift = 0;
  if (shiba)
  {
    shift = static_cast<std::int64_t>(up);
    down = sites - down;
  }
  bool const particle_hole = up + down > sites;
  if (particle_hole)
  {
    shift += (shiba ? -1 : 1) * static_cast<std::int64_t>(up + down - sites);
    up = sites - up;
    down = sites - down;
  }
  bool const spin_reversed = up < down;
  if (spin_reversed) std::swap(up, down);
  Real const offset = interaction * Real(shift);
  if (!uni20::isfinite(offset)) throw std::overflow_error("Hubbard symmetry energy offset exceeds this precision");
  return {up, down, root_u, offset, shiba, particle_hole, spin_reversed};
}

template <typename StateType, uni20::Real Real>
void apply_mapping(StateType& result, std::size_t particles, std::size_t down, Real interaction,
                   SectorMapping<Real> const& mapping)
{
  result.particles = particles;
  result.down_spins = down;
  result.interaction = interaction;
  result.energy_offset = mapping.energy_offset;
  result.shiba_transformed = mapping.shiba;
  result.particle_hole_transformed = mapping.particle_hole;
  result.spin_reversed = mapping.spin_reversed;
  result.energy += mapping.energy_offset;
  if (!uni20::isfinite(result.energy)) throw std::overflow_error("Hubbard energy exceeds this precision");
}

template <uni20::Real Real> void check_newton_size(std::size_t order)
{
  auto const max_elements = static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (order != 0 && order > max_elements / order) throw std::invalid_argument("Hubbard Newton matrix is too large");
}

// Width/(width^2+d^2), without squaring a potentially huge value.
template <uni20::Real Real> Real kernel(Real d, Real width)
{
  using std::abs;
  if (abs(d) > width)
  {
    Real const ratio = width / d;
    return (ratio / d) / (Real{1} + ratio * ratio);
  }
  Real const ratio = d / width;
  return (Real{1} / width) / (Real{1} + ratio * ratio);
}

/// Shared damped Newton/continuation driver for boundary-specific root systems.
/// Analytic Newton with U=8 -> requested U continuation when U<8.
/// The update budget is shared by all stages; final residuals use the target U.
template <uni20::Real Real, typename StateType, typename System>
StateType solve_ground_system(System const& system, Real interaction, SolverOptions<Real> const& options)
{
  using std::cos;
  StateType result;
  result.sites = system.sites;
  result.particles = result.root_particles = system.n;
  result.down_spins = result.root_down_spins = system.m;
  result.interaction = result.root_interaction = interaction;
  Real const target_u = interaction / Real{4};
  if (target_u == Real{0} || !uni20::isfinite(Real{1} / target_u))
    throw std::invalid_argument("U is too small for this precision; U/4 and its reciprocal must be representable");
  auto x = system.seed();
  Real u = std::max(Real{2}, target_u);
  uni20::DenseMatrix<Real> jacobian(system.order, system.order), step(system.order, 1);
  for (;;)
  {
    auto evaluation = system.evaluate(x, u, &jacobian);
    if (evaluation.norm() <= options.residual_tolerance)
    {
      ++result.continuation_steps;
      if (u == target_u) break;
      Real const next = std::max(target_u, u / Real{2});
      Real const factor = std::max(Real{1}, u) / std::max(Real{1}, next);
      for (std::size_t a = 0; a < system.ns; ++a)
        x[system.nk + a] *= factor;
      u = next;
      continue;
    }
    if (result.iterations == options.max_iterations) break;
    for (std::size_t i = 0; i < system.order; ++i)
      step[i, 0] = -evaluation.residual[i];
    uni20::linalg::solve_inplace(jacobian, step);
    bool const accepted = bethe::detail::backtrack_newton(
        x, [&](std::size_t j) { return step[j, 0]; },
        [&](auto const& trial, Real damping) {
          return system.physical(trial) &&
                 bethe::detail::newton_decreases(system.evaluate(trial, u).norm(), evaluation.norm(), damping,
                                                 options.residual_tolerance);
        });
    if (!accepted)
    {
      result.status = SolveStatus::stalled;
      break;
    }
    ++result.iterations;
  }
  // Preserve conventional Lambda while reevaluating at the requested U.
  Real const factor = std::max(Real{1}, u) / std::max(Real{1}, target_u);
  for (std::size_t a = 0; a < system.ns; ++a)
    x[system.nk + a] *= factor;
  auto const evaluation = system.evaluate(x, target_u);
  result.charge_residual = evaluation.charge;
  result.spin_residual = evaluation.spin;
  result.residual_norm = evaluation.norm();
  result.converged = result.residual_norm <= options.residual_tolerance;
  if (result.converged) result.status = SolveStatus::converged;
  system.expand(x, target_u, result.charge_momenta, result.spin_rapidities);
  result.quantum_numbers = system.labels;
  bethe::detail::CompensatedSum<Real> energy;
  for (Real k : result.charge_momenta)
    energy.add(-Real{2} * cos(k));
  result.energy = energy.value();
  if (!uni20::isfinite(result.energy)) throw std::runtime_error("nonfinite Hubbard energy");
  for (Real lambda : result.spin_rapidities)
    if (!uni20::isfinite(lambda)) throw std::runtime_error("nonfinite Hubbard spin rapidity");
  return result;
}
} // namespace detail
} // namespace bethe::hubbard
