// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/detail/continuum_newton.hpp>
#include <bethe/real_excitations.hpp>
#include <bethe/solver.hpp>
#include <cmath>
#include <span>
#include <uni20/linalg/ops/linear_solve.hpp>

namespace bethe::lieb_liniger
{
using QuantumNumbers = std::vector<uni20::half_int>;
template <uni20::Real Real> using SolverOptions = bethe::SolverOptions<Real>;
enum class SolveStatus
{
  converged,
  iteration_limit,
  stalled
};

/// Repulsive continuum bosons on a ring. Units: hbar^2/(2m)=1,
/// H=-sum d_j^2 + 2c sum_(i<j) delta(x_i-x_j).
template <uni20::Real Real> struct State
{
    std::size_t particles = 0;
    Real length = Real{1}, interaction = Real{1};
    QuantumNumbers quantum_numbers;
    /// Physical k, not dimensionless k*length. Energy is sum(k*k).
    std::vector<Real> momenta;
    Real energy = Real{0}, momentum = Real{0};
    /// Exact sum of I, signed and NOT reduced modulo N.
    std::int64_t momentum_index = 0;
    /// Dimensionless, component-scaled equation residual; see the model guide.
    Real residual_norm = Real{0};
    std::size_t iterations = 0;
    bool converged = false;
    SolveStatus status = SolveStatus::iteration_limit;
};

namespace detail
{
inline void check_particles(std::size_t n)
{
  if (n > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() / 4))
    throw std::invalid_argument("Lieb-Liniger particle count exceeds the quantum-number range");
}

template <uni20::Real Real>
Real validate_parameters(std::size_t n, Real length, Real c, SolverOptions<Real> const& options)
{
  check_particles(n);
  if (!uni20::isfinite(length) || length <= Real{0} || !uni20::isfinite(c) || c <= Real{0})
    throw std::invalid_argument("repulsive Lieb-Liniger requires finite length > 0 and c > 0");
  Real const g = c * length;
  if (!uni20::isfinite(g) || g <= Real{0})
    throw std::invalid_argument("c*length must be finite and nonzero in the selected precision");
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");
  auto const elements = static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (n && n > elements / n) throw std::length_error("Lieb-Liniger Newton matrix is too large");
  return g;
}

inline std::int64_t validate_numbers(QuantumNumbers const& numbers)
{
  auto const n = numbers.size();
  check_particles(n);
  // A wider accumulator permits cancellation without signed intermediate overflow.
  __int128 twice_sum = 0;
  for (std::size_t j = 0; j < n; ++j)
  {
    auto const value = numbers[j].twice();
    if ((value % 2 != 0) != (n % 2 == 0) || (j && numbers[j] <= numbers[j - 1]))
      throw std::invalid_argument("I must be strictly increasing integers for odd N, half-odd integers for even N");
    twice_sum += value;
  }
  auto const sum = twice_sum / 2;
  if (sum < std::numeric_limits<std::int64_t>::min() || sum > std::numeric_limits<std::int64_t>::max())
    throw std::overflow_error("total Lieb-Liniger momentum index exceeds int64");
  return static_cast<std::int64_t>(sum);
}

template <uni20::Real Real> struct System
{
    QuantumNumbers const& numbers;
    Real g;
    Real const pi = Real{4} * std::atan(Real{1});

    std::int64_t ground_twice(std::size_t j) const
    {
      return 2 * static_cast<std::int64_t>(j) - static_cast<std::int64_t>(numbers.size() - 1);
    }
    Real reduced_twice(std::size_t j) const
    {
      return static_cast<Real>(static_cast<__int128>(numbers[j].twice()) - ground_twice(j));
    }
    auto seed() const
    {
      using std::sqrt;
      Real const factor = std::min(Real{1}, sqrt(g) / (pi * sqrt(Real(std::max(std::size_t{1}, numbers.size())))));
      std::vector<Real> q(numbers.size());
      for (std::size_t j = 0; j < q.size(); ++j)
        q[j] = pi * (reduced_twice(j) + factor * Real(ground_twice(j)));
      return q;
    }
    bool physical(std::span<Real const> q) const
    {
      for (std::size_t j = 0; j < q.size(); ++j)
        if (!uni20::isfinite(q[j]) || (j && q[j] <= q[j - 1])) return false;
      return true;
    }
    // g/(g*g+d*d), avoiding overflow/underflow from the squares.
    Real kernel(Real d) const { return bethe::detail::rational_scattering_kernel(d, g); }
    struct Evaluation
    {
        std::vector<Real> residual;
        Real norm = Real{0};
    };
    Evaluation evaluate(std::span<Real const> q, uni20::DenseMatrix<Real>* jacobian = nullptr) const
    {
      using std::abs;
      using std::atan;
      using std::sqrt;
      auto const n = q.size();
      Evaluation result{std::vector<Real>(n)};
      for (std::size_t j = 0; j < n; ++j)
      {
        bethe::detail::CompensatedSum<Real> sum, diagonal;
        Real const target = pi * (g < Real{1} ? reduced_twice(j) : Real(numbers[j].twice()));
        sum.add(q[j] - target);
        diagonal.add(Real{1});
        for (std::size_t k = 0; k < n; ++k)
          if (j != k)
          {
            Real const difference = q[j] - q[k];
            if (!uni20::isfinite(difference)) throw std::overflow_error("Lieb-Liniger root difference overflow");
            // For g<1 subtract the exact rank contribution from I first.
            // This avoids cancelling pi-sized phases to recover sqrt(g) roots.
            Real phase;
            if (g < Real{1})
              phase = -Real{2} * (abs(difference) > g
                                      ? atan(g / difference)
                                      : (difference > Real{0} ? pi / Real{2} : -pi / Real{2}) - atan(difference / g));
            else
              phase = Real{2} * (abs(difference) > g
                                     ? (difference > Real{0} ? pi / Real{2} : -pi / Real{2}) - atan(g / difference)
                                     : atan(difference / g));
            sum.add(phase);
            if (jacobian)
            {
              Real const value = Real{2} * kernel(difference);
              if (!uni20::isfinite(value)) throw std::overflow_error("Lieb-Liniger Jacobian exceeds this precision");
              (*jacobian)[j, k] = -value;
              diagonal.add(value);
            }
          }
        if (jacobian)
        {
          (*jacobian)[j, j] = diagonal.value();
          if (!uni20::isfinite(diagonal.value())) throw std::overflow_error("Lieb-Liniger Jacobian overflow");
        }
        result.residual[j] = sum.value();
        if (!uni20::isfinite(result.residual[j])) throw std::overflow_error("Lieb-Liniger residual overflow");
        Real const scale = std::max({std::min(Real{1}, sqrt(g)), abs(q[j]), abs(target)});
        result.norm = std::max(result.norm, abs(sum.value()) / scale);
      }
      return result;
    }
};
} // namespace detail

/// Consecutive, reflection-symmetric labels I_j=j-(N-1)/2, including N=0.
inline QuantumNumbers ground_quantum_numbers(std::size_t particles)
{
  detail::check_particles(particles);
  QuantumNumbers result(particles);
  for (std::size_t j = 0; j < particles; ++j)
    result[j] = uni20::from_twice(2 * static_cast<std::int64_t>(j) - static_cast<std::int64_t>(particles - 1));
  return result;
}

/// Solve a specified repulsive state. The label count specifies N. Optional
/// initial momenta must be finite, strictly ordered, and have the same size.
/// A zero budget evaluates only the seed. Failed solves remain unranked estimates.
template <uni20::Real Real>
State<Real> solve_real(Real length, Real c, QuantumNumbers const& numbers, SolverOptions<Real> const& options = {},
                       std::span<Real const> initial_momenta = {})
{
  Real const g = detail::validate_parameters(numbers.size(), length, c, options);
  auto const momentum_index = detail::validate_numbers(numbers);
  detail::System<Real> const system{numbers, g};
  auto q = system.seed();
  if (!initial_momenta.empty())
  {
    if (initial_momenta.size() != q.size()) throw std::invalid_argument("initial momentum count differs from N");
    for (std::size_t j = 0; j < q.size(); ++j)
      q[j] = initial_momenta[j] * length;
  }
  if (!system.physical(q))
    throw std::invalid_argument("initial momenta must be finite and strictly ordered at the selected precision");
  State<Real> result;
  result.particles = numbers.size();
  result.length = length;
  result.interaction = c;
  result.quantum_numbers = numbers;
  result.momentum_index = momentum_index;
  result.momentum = (Real{2} * system.pi * Real(momentum_index)) / length;
  if (!uni20::isfinite(result.momentum)) throw std::overflow_error("Lieb-Liniger total momentum overflow");
  auto const n = q.size();
  result.iterations = bethe::detail::continuum_newton(system, q, options);
  result.momenta.resize(n);
  bethe::detail::CompensatedSum<Real> energy;
  bool nonzero = false;
  for (std::size_t j = 0; j < n; ++j)
  {
    Real const k = q[j] / length, term = k * k;
    if (!uni20::isfinite(k) || !uni20::isfinite(term)) throw std::overflow_error("Lieb-Liniger energy overflow");
    nonzero = nonzero || q[j] != Real{0};
    result.momenta[j] = k;
    energy.add(term);
    q[j] = k * length;
  }
  result.energy = energy.value();
  if (!uni20::isfinite(result.energy)) throw std::overflow_error("Lieb-Liniger energy sum overflow");
  if (nonzero && result.energy == Real{0}) throw std::underflow_error("Lieb-Liniger energy underflow");
  result.residual_norm = system.physical(q) ? system.evaluate(q).norm : uni20::numeric_limits<Real>::infinity();
  result.converged = system.physical(q) && result.residual_norm <= options.residual_tolerance;
  if (result.converged)
    result.status = SolveStatus::converged;
  else if (result.iterations < options.max_iterations)
    result.status = SolveStatus::stalled;
  return result;
}

template <uni20::Real Real>
State<Real> ground_state(std::size_t particles, Real length, Real c, SolverOptions<Real> const& options = {})
{
  (void)detail::validate_parameters(particles, length, c, options);
  return solve_real(length, c, ground_quantum_numbers(particles), options);
}

/// Add padding label slots at EACH edge of the ground-state interval.
/// Choose N distinct labels from N+2*padding slots. There are infinitely many
/// states without this explicit cutoff. Padding 0 contains only the ground state.
inline std::size_t excitation_count(std::size_t particles, std::size_t padding,
                                    std::size_t limit = std::numeric_limits<std::size_t>::max())
{
  detail::check_particles(particles);
  auto const bound = static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() / 4);
  if (padding > (bound - particles) / 2) throw std::length_error("Lieb-Liniger label window exceeds integer range");
  return bethe::detail::bounded_binomial(particles + 2 * padding, particles, limit);
}

template <uni20::Real Real> struct ExcitationScan : bethe::RealExcitationScan<State<Real>>
{
    std::size_t particles = 0, padding = 0;
    Real length = Real{1}, interaction = Real{1};
};

template <uni20::Real Real>
ExcitationScan<Real> real_excitations(std::size_t particles, Real length, Real c, std::size_t padding,
                                      RealExcitationOptions const& enumeration = {},
                                      SolverOptions<Real> const& options = {})
{
  (void)detail::validate_parameters(particles, length, c, options);
  if (enumeration.count == 0 || enumeration.max_candidates == 0)
    throw std::invalid_argument("excitation count and max_candidates must be positive");
  (void)excitation_count(particles, padding, enumeration.max_candidates);
  auto const first = particles ? -static_cast<std::int64_t>(particles - 1) - 2 * static_cast<std::int64_t>(padding) : 0;
  auto result = bethe::detail::scan_real_combinations<ExcitationScan<Real>>(
      particles + 2 * padding, particles, first, enumeration,
      [&](QuantumNumbers const& numbers) { return solve_real(length, c, numbers, options); },
      [&] { return ground_state(particles, length, c, options); });
  result.particles = particles;
  result.padding = padding;
  result.length = length;
  result.interaction = c;
  return result;
}
} // namespace bethe::lieb_liniger
