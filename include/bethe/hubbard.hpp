// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/solver.hpp>
#include <uni20/common/half_int.hpp>
#include <uni20/linalg/ops/linear_solve.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
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

/// Periodic sector ground state, t=1, even L>=2 (supported sectors below).
/// H=-sum_(j,sigma)(c^dagger_j c_(j+1)+h.c.)+U sum_j n_up n_down.
template <uni20::Real Real> struct State
{
    std::size_t sites = 0, particles = 0, down_spins = 0;
    Real interaction = Real{0};
    /// Roots and residuals describe this sector, which may be symmetry-mapped.
    std::size_t root_particles = 0, root_down_spins = 0;
    Real root_interaction = Real{0};
    /// E_physical = E_roots + energy_offset; likewise for momentum modulo L.
    Real energy_offset = Real{0};
    std::size_t momentum_offset = 0;
    bool shiba_transformed = false, particle_hole_transformed = false, spin_reversed = false;
    bool auxiliary_roots() const { return shiba_transformed || particle_hole_transformed || spin_reversed; }
    std::vector<Real> charge_momenta;
    /// Conventional Lieb-Wu Lambda (not the internally scaled Newton variable).
    std::vector<Real> spin_rapidities;
    NestedQuantumNumbers quantum_numbers;
    Real energy = Real{0};
    std::size_t momentum_index = 0;
    Real momentum = Real{0};
    /// max|F_charge|/L and max|F_spin|/L, at root_interaction (not a continuation stage).
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

/// Normalized repulsive sector. Caller validates the public parameters.
/// Analytic reduced Newton with U=8 -> requested U continuation when U<8.
/// The update budget is shared by all stages; final residuals use the target U.
template <uni20::Real Real>
State<Real> repulsive_ground_state(std::size_t sites, std::size_t particles, std::size_t down, Real interaction,
                                   SolverOptions<Real> const& options)
{
  using std::abs;
  using std::atan;
  using std::cos;
  using std::sin;
  checked_root_sector(sites, particles, down);
  if (down == 0) return free_state(sites, particles, down, interaction);
  State<Real> result;
  result.sites = sites;
  result.particles = result.root_particles = particles;
  result.down_spins = result.root_down_spins = down;
  result.interaction = result.root_interaction = interaction;
  Real const pi = Real{4} * atan(Real{1});
  result.momentum_index = down % 2 == 0 ? sites / 2 : 0;
  result.momentum = result.momentum_index == 0 ? Real{0} : pi;
  Real const target_u = interaction / Real{4};
  if (target_u == Real{0} || !uni20::isfinite(Real{1} / target_u))
    throw std::invalid_argument("U is too small for this precision; U/4 and its reciprocal must be representable");
  auto const order = particles / 2 - (down % 2 == 0 ? 1 : 0) + down / 2;
  auto const max_elements = static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (order > max_elements / order) throw std::invalid_argument("Hubbard Newton matrix is too large");
  GroundSystem<Real> const system(sites, particles, down);
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
    bool accepted = false;
    auto trial = x;
    Real damping = Real{1};
    for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits; ++backtrack)
    {
      for (std::size_t i = 0; i < system.order; ++i)
        trial[i] = x[i] + damping * step[i, 0];
      if (system.physical(trial))
      {
        auto const norm = system.evaluate(trial, u).norm();
        if (norm <= options.residual_tolerance || norm < (Real{1} - damping / Real{10000}) * evaluation.norm())
        {
          accepted = true;
          break;
        }
      }
      damping /= Real{2};
    }
    if (!accepted)
    {
      result.status = SolveStatus::stalled;
      break;
    }
    x = std::move(trial);
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
  if (!uni20::isfinite(interaction)) throw std::invalid_argument("Hubbard ground state requires finite U");
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");
  if (particles > 2 * sites) throw std::invalid_argument("Hubbard particles must satisfy 0 <= N <= 2L");
  auto const bound = static_cast<std::int64_t>(std::min(particles, 2 * sites - particles)), twice_sz = sz.twice();
  if (twice_sz < -bound || twice_sz > bound || (static_cast<std::int64_t>(particles) - twice_sz) % 2 != 0)
    throw std::invalid_argument("Hubbard N and Sz must give integer N_up and N_down between 0 and L");
  auto const physical_down = static_cast<std::size_t>((static_cast<std::int64_t>(particles) - twice_sz) / 2);
  auto up = particles - physical_down, down = physical_down;
  if (interaction == Real{0} || up == 0 || down == 0) return detail::free_state(sites, particles, down, interaction);

  bool const shiba = interaction < Real{0};
  Real const root_u = shiba ? -interaction : interaction;
  // Combine integer coefficients BEFORE multiplying by U, avoiding cancellation
  // between potentially huge Shiba and particle-hole energy offsets.
  std::int64_t shift = 0;
  std::size_t momentum_shift = 0;
  if (shiba)
  {
    shift = static_cast<std::int64_t>(up);
    if (down % 2 == 0) momentum_shift = sites / 2;
    down = sites - down;
  }
  bool const particle_hole = up + down > sites;
  if (particle_hole)
  {
    shift += (shiba ? -1 : 1) * static_cast<std::int64_t>(up + down - sites);
    if ((up + down) % 2 != 0) momentum_shift = (momentum_shift + sites / 2) % sites;
    up = sites - up;
    down = sites - down;
  }
  bool const spin_reversed = up < down;
  if (spin_reversed) std::swap(up, down);
  Real const offset = interaction * Real(shift);
  if (!uni20::isfinite(offset)) throw std::overflow_error("Hubbard symmetry energy offset exceeds this precision");
  auto result = detail::repulsive_ground_state(sites, up + down, down, root_u, options);
  result.particles = particles;
  result.down_spins = physical_down;
  result.interaction = interaction;
  result.energy_offset = offset;
  result.momentum_offset = momentum_shift;
  result.shiba_transformed = shiba;
  result.particle_hole_transformed = particle_hole;
  result.spin_reversed = spin_reversed;
  result.energy += offset;
  result.momentum_index = (result.momentum_index + momentum_shift) % sites;
  result.momentum = Real{8} * std::atan(Real{1}) * Real(result.momentum_index) / Real(sites);
  if (!uni20::isfinite(result.energy)) throw std::overflow_error("Hubbard energy exceeds this precision");
  return result;
}

/// Convenience wrapper: half filling, Sz=0. Both signs of U are supported.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> ground_state(std::size_t sites, Real interaction, SolverOptions<Real> const& options = {})
{
  return sector_ground_state(sites, sites, uni20::half_int{0}, interaction, options);
}
} // namespace bethe::hubbard
