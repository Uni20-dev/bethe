// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/lieb_liniger.hpp>

namespace bethe::lieb_liniger::open
{
/// Dirichlet walls at x=0,L. Same units/coupling as the ring, but NO conserved
/// total momentum: positive Bethe wave numbers are standing-wave coordinates.
template <uni20::Real Real> struct State
{
    std::size_t particles = 0;
    Real length{1}, interaction{1};
    QuantumNumbers quantum_numbers;
    std::vector<Real> momenta{};
    Real energy{}, residual_norm{};
    std::size_t iterations = 0;
    bool converged = false;
    SolveStatus status = SolveStatus::iteration_limit;
};

namespace detail
{
inline void validate_numbers(QuantumNumbers const& numbers)
{
  lieb_liniger::detail::check_particles(numbers.size());
  for (std::size_t j = 0; j < numbers.size(); ++j)
    if (numbers[j].twice() <= 0 || numbers[j].twice() % 2 || (j && numbers[j] <= numbers[j - 1]))
      throw std::invalid_argument("hard-wall Lieb-Liniger labels must be strictly increasing positive integers");
}

template <uni20::Real Real> struct System
{
    QuantumNumbers const& numbers;
    Real g;
    Real const pi = Real{4} * std::atan(Real{1});
    auto seed() const
    {
      Real const factor =
          std::min(Real{1}, std::sqrt(g) / (pi * std::sqrt(Real(std::max(std::size_t{1}, numbers.size())))));
      std::vector<Real> q(numbers.size());
      for (std::size_t j = 0; j < q.size(); ++j)
        q[j] = pi * (Real(numbers[j].twice() / 2 - std::int64_t(j)) + factor * Real(j));
      return q;
    }
    bool physical(std::span<Real const> q) const
    {
      for (std::size_t j = 0; j < q.size(); ++j)
        if (!uni20::isfinite(q[j]) || q[j] <= Real{0} || (j && q[j] <= q[j - 1])) return false;
      return true;
    }
    typename lieb_liniger::detail::System<Real>::Evaluation evaluate(std::span<Real const> q,
                                                                     uni20::DenseMatrix<Real>* jacobian = nullptr) const
    {
      typename lieb_liniger::detail::System<Real>::Evaluation result{std::vector<Real>(q.size())};
      auto phase = [&](Real x) {
        if (g >= Real{1}) return std::atan2(x, g);
        Real const complement = std::atan2(g, std::abs(x));
        return x > Real{0} ? -complement : complement;
      };
      for (std::size_t j = 0; j < q.size(); ++j)
      {
        auto const label = numbers[j].twice() / 2;
        Real const target = pi * Real(g < Real{1} ? label - std::int64_t(j) : label);
        bethe::detail::CompensatedSum<Real> sum, diagonal;
        sum.add(q[j] - target);
        diagonal.add(Real{1});
        for (std::size_t k = 0; k < q.size(); ++k)
          if (j != k) // Neither direct self-scattering nor the self-image at 2q_j.
          {
            Real const difference = q[j] - q[k], image = q[j] + q[k];
            if (!uni20::isfinite(image)) throw std::overflow_error("hard-wall root sum overflow");
            sum.add(phase(difference));
            sum.add(phase(image));
            if (jacobian)
            {
              Real const a = bethe::detail::rational_scattering_kernel(difference, g);
              Real const b = bethe::detail::rational_scattering_kernel(image, g);
              if (!uni20::isfinite(a) || !uni20::isfinite(b))
                throw std::overflow_error("hard-wall Jacobian exceeds this precision");
              diagonal.add(a);
              diagonal.add(b);
              (*jacobian)[j, k] = b - a;
            }
          }
        if (jacobian)
        {
          (*jacobian)[j, j] = diagonal.value();
          if (!uni20::isfinite(diagonal.value())) throw std::overflow_error("hard-wall Jacobian overflow");
        }
        result.residual[j] = sum.value();
        if (!uni20::isfinite(sum.value())) throw std::overflow_error("hard-wall residual overflow");
        result.norm = std::max(result.norm, std::abs(sum.value()) / std::max({Real{1}, std::abs(target), q[j]}));
      }
      return result;
    }
};
} // namespace detail

/// I_j=j+1, j=0,...,N-1; vacuum allowed. No parity shift as on the ring.
inline QuantumNumbers ground_quantum_numbers(std::size_t particles)
{
  lieb_liniger::detail::check_particles(particles);
  QuantumNumbers result(particles);
  for (std::size_t j = 0; j < particles; ++j)
    result[j] = uni20::half_int(std::int64_t(j + 1));
  return result;
}

/// Positive real roots of q_j + sum_(k!=j)[atan((q_j-q_k)/g)+
/// atan((q_j+q_k)/g)] = pi I_j, q=kL, g=cL. Repulsion c>0 only.
template <uni20::Real Real>
State<Real> solve_real(Real length, Real c, QuantumNumbers const& numbers, SolverOptions<Real> const& options = {},
                       std::span<Real const> initial_momenta = {})
{
  Real const g = lieb_liniger::detail::validate_parameters(numbers.size(), length, c, options);
  detail::validate_numbers(numbers);
  detail::System<Real> const system{numbers, g};
  auto q = system.seed();
  if (!initial_momenta.empty())
  {
    if (initial_momenta.size() != q.size()) throw std::invalid_argument("initial wave-number count differs from N");
    for (std::size_t j = 0; j < q.size(); ++j)
      q[j] = initial_momenta[j] * length;
  }
  if (!system.physical(q))
    throw std::invalid_argument("hard-wall seed must be finite, positive and distinctly ordered at this precision");
  State<Real> result{.particles = numbers.size(), .length = length, .interaction = c, .quantum_numbers = numbers};
  result.iterations = bethe::detail::continuum_newton(system, q, options);
  result.momenta.resize(q.size());
  bethe::detail::CompensatedSum<Real> energy;
  for (std::size_t j = 0; j < q.size(); ++j)
  {
    Real const k = q[j] / length, term = k * k;
    if (!uni20::isfinite(k) || !uni20::isfinite(term)) throw std::overflow_error("hard-wall energy overflow");
    result.momenta[j] = k;
    energy.add(term);
    q[j] = k * length;
  }
  result.energy = energy.value();
  if (!uni20::isfinite(result.energy)) throw std::overflow_error("hard-wall energy sum overflow");
  if (!q.empty() && result.energy == Real{0}) throw std::underflow_error("hard-wall energy underflow");
  result.residual_norm = system.physical(q) ? system.evaluate(q).norm : uni20::numeric_limits<Real>::infinity();
  result.converged = result.residual_norm <= options.residual_tolerance;
  result.status = result.converged                              ? SolveStatus::converged
                  : result.iterations == options.max_iterations ? SolveStatus::iteration_limit
                                                                : SolveStatus::stalled;
  return result;
}

template <uni20::Real Real>
State<Real> ground_state(std::size_t particles, Real length, Real c, SolverOptions<Real> const& options = {})
{
  (void)lieb_liniger::detail::validate_parameters(particles, length, c, options);
  return solve_real(length, c, ground_quantum_numbers(particles), options);
}

/// Choose N distinct labels from 1,...,N+padding (one upper edge, unlike PBC).
inline std::size_t excitation_count(std::size_t particles, std::size_t padding,
                                    std::size_t limit = std::numeric_limits<std::size_t>::max())
{
  lieb_liniger::detail::check_particles(particles);
  auto const bound = static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() / 4);
  if (padding > bound - particles) throw std::length_error("hard-wall label window exceeds integer range");
  return bethe::detail::bounded_binomial(particles + padding, particles, limit);
}

template <uni20::Real Real> using ExcitationScan = bethe::RealExcitationScan<State<Real>>;

template <uni20::Real Real>
ExcitationScan<Real> real_excitations(std::size_t particles, Real length, Real c, std::size_t padding,
                                      RealExcitationOptions const& enumeration = {},
                                      SolverOptions<Real> const& options = {})
{
  (void)lieb_liniger::detail::validate_parameters(particles, length, c, options);
  if (enumeration.count == 0 || enumeration.max_candidates == 0)
    throw std::invalid_argument("excitation count and max_candidates must be positive");
  (void)excitation_count(particles, padding, enumeration.max_candidates);
  return bethe::detail::scan_real_combinations<ExcitationScan<Real>>(
      particles + padding, particles, 2, enumeration,
      [&](QuantumNumbers const& numbers) { return solve_real(length, c, numbers, options); },
      [&] { return ground_state(particles, length, c, options); });
}
} // namespace bethe::lieb_liniger::open
