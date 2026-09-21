// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include <bethe/heisenberg.hpp>

namespace bethe::heisenberg::open
{
/// Free-end XXX result for H=sum_{i=0}^{N-2} S_i.S_{i+1}, J=1, h=0.
/// Unlike a periodic state, this has no translation momentum or momentum index.
template <uni20::Real Real> struct RealState
{
    std::vector<Real> rapidities;
    Real energy = Real{0};
    /// max|F_i|/(2*N), where F_i is the open-chain logarithmic equation.
    Real residual_norm = Real{0};
    std::size_t iterations = 0;
    bool converged = false;
    QuantumNumbers quantum_numbers;
    uni20::half_int sz;
    bool spin_reversed = false;
};

namespace detail
{
inline void validate_quantum_numbers(std::size_t sites, std::span<uni20::half_int const> numbers)
{
  auto const n = heisenberg::detail::checked_sites(sites);
  if (numbers.size() > sites / 2) throw std::invalid_argument("open finite real-root states require M <= N/2");
  auto const bound = n - static_cast<std::int64_t>(numbers.size());
  for (std::size_t i = 0; i < numbers.size(); ++i)
  {
    auto const twice = numbers[i].twice();
    if (twice < 2 || twice > 2 * bound || twice % 2 != 0)
      throw std::invalid_argument("open Bethe quantum numbers must be integers in [1,N-M]");
    if (i != 0 && numbers[i - 1] >= numbers[i])
      throw std::invalid_argument("open Bethe quantum numbers must be strictly increasing");
  }
}
} // namespace detail

/// Solve a free-end, highest-weight XXX state with M positive, finite real roots.
/// Supported quantum numbers are sorted distinct integers 1 <= I_i <= N-M,
/// with M <= N/2. Zero roots and their spurious Bethe vectors are excluded;
/// complex strings and infinite-root SU(2) descendants are not implemented.
///
/// With phi(z)=2*atan(z), solve
/// F_i=2*N*phi(z_i)-2*pi*I_i-sum_{j!=i}[phi((z_i-z_j)/2)+phi((z_i+z_j)/2)].
/// Neither the direct self term nor the reflected self term is included.
/// E=(N-1)/4-sum_i 2/(1+z_i^2). Residual tolerance applies to max|F|/(2*N),
/// not to the energy error. O(M^2) per simultaneous update and O(M) storage.
/// The initial zero guess is allowed; a zero budget evaluates it without updates.
template <uni20::Real Real = double>
[[nodiscard]] RealState<Real> solve_real(std::size_t sites, std::span<uni20::half_int const> numbers,
                                         SolverOptions<Real> const& options = {},
                                         std::span<Real const> initial_roots = {})
{
  using std::abs;
  using std::atan;
  using std::tan;
  detail::validate_quantum_numbers(sites, numbers);
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");

  Real const n = static_cast<Real>(sites);
  Real const pi = bethe::detail::pi<Real>();
  auto const roots = numbers.size();
  RealState<Real> result;
  result.quantum_numbers.assign(numbers.begin(), numbers.end());
  result.sz = uni20::from_twice(static_cast<std::int64_t>(sites - 2 * roots));
  result.rapidities.resize(roots, Real{0});
  if (!initial_roots.empty())
  {
    if (initial_roots.size() != roots) throw std::invalid_argument("initial root count must match the quantum numbers");
    for (Real const root : initial_roots)
      if (!uni20::isfinite(root) || root < Real{0})
        throw std::invalid_argument("open initial roots must be finite and nonnegative");
    result.rapidities.assign(initial_roots.begin(), initial_roots.end());
  }
  std::vector<Real> angles(roots);
  for (;;)
  {
    result.residual_norm = Real{0};
    for (std::size_t i = 0; i < roots; ++i)
    {
      heisenberg::detail::CompensatedSum<Real> phase;
      for (std::size_t j = 0; j < roots; ++j)
        if (i != j)
        {
          phase.add(atan((result.rapidities[i] - result.rapidities[j]) / Real{2}));
          // Divide before adding, to avoid overflow for large finite guesses.
          phase.add(atan(result.rapidities[i] / Real{2} + result.rapidities[j] / Real{2}));
        }
      Real const quantum_number = static_cast<Real>(numbers[i].twice()) / Real{2};
      angles[i] = (pi * quantum_number + phase.value()) / (Real{2} * n);
      Real const residual = Real{2} * abs(atan(result.rapidities[i]) - angles[i]);
      if (!uni20::isfinite(residual)) throw std::runtime_error("nonfinite open Heisenberg equation residual");
      result.residual_norm = std::max(result.residual_norm, residual);
    }
    result.converged = result.residual_norm <= options.residual_tolerance;
    if (result.converged || result.iterations == options.max_iterations) break;
    for (std::size_t i = 0; i < roots; ++i)
    {
      result.rapidities[i] = tan(angles[i]);
      if (!uni20::isfinite(result.rapidities[i])) throw std::runtime_error("nonfinite open Heisenberg rapidity");
    }
    ++result.iterations;
  }
  heisenberg::detail::CompensatedSum<Real> magnons;
  for (Real const z : result.rapidities)
    magnons.add(-Real{2} / (Real{1} + z * z));
  result.energy = (n - Real{1}) / Real{4} + magnons.value();
  if (!uni20::isfinite(result.energy)) throw std::runtime_error("nonfinite open Heisenberg energy");
  return result;
}

/// The sector minimum uses I=1,...,M with M=N/2-|Sz|, for either parity of N.
[[nodiscard]] inline QuantumNumbers sector_ground_quantum_numbers(std::size_t sites, uni20::half_int sz)
{
  auto const m = heisenberg::detail::sector_roots(sites, sz);
  QuantumNumbers numbers;
  numbers.reserve(static_cast<std::size_t>(m));
  for (std::int64_t i = 1; i <= m; ++i)
    numbers.emplace_back(i);
  return numbers;
}

template <uni20::Real Real = double>
[[nodiscard]] RealState<Real> sector_ground_state(std::size_t sites, uni20::half_int sz,
                                                  SolverOptions<Real> const& options = {})
{
  auto const numbers = open::sector_ground_quantum_numbers(sites, sz);
  auto result = open::solve_real<Real>(sites, numbers, options);
  result.sz = sz;
  result.spin_reversed = sz.twice() < 0;
  return result;
}

/// All N+1 sector minima in ascending Sz. Reuses spin-reversed partners.
/// Retaining the roots of every sector takes O(N^2) storage.
template <uni20::Real Real = double>
[[nodiscard]] std::vector<RealState<Real>> sector_ground_states(std::size_t sites,
                                                                SolverOptions<Real> const& options = {})
{
  auto const n = heisenberg::detail::checked_sites(sites);
  std::vector<RealState<Real>> states(sites + 1);
  for (std::size_t m = 0; m <= sites / 2; ++m)
  {
    auto const sz = uni20::from_twice(n - 2 * static_cast<std::int64_t>(m));
    states[sites - m] = open::sector_ground_state<Real>(sites, sz, options);
    if (m != sites - m)
    {
      states[m] = states[sites - m];
      states[m].sz = -sz;
      states[m].spin_reversed = true;
    }
  }
  return states;
}

/// Free-end ground state, N>=2: Sz=0 for even N, Sz=1/2 for odd N.
template <uni20::Real Real = double>
[[nodiscard]] RealState<Real> ground_state(std::size_t sites, SolverOptions<Real> const& options = {})
{
  return open::sector_ground_state<Real>(sites, uni20::from_twice(static_cast<std::int64_t>(sites % 2)), options);
}
} // namespace bethe::heisenberg::open
