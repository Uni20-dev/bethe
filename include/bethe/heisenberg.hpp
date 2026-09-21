// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2015-2023, 2026 Ian McCulloch
// Successor to MPToolkit's misc/heisenberg-energy.cpp. See CITATIONS.md.

#pragma once

#include <bethe/spinon.hpp>
#include <uni20/common/half_int.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace bethe::heisenberg
{

/// Controls the real-root, zero-field fixed-point iteration.
template <uni20::Real Real> struct SolverOptions
{
  Real residual_tolerance = Real{32} * uni20::numeric_limits<Real>::epsilon();
  std::size_t max_iterations = 10000;
};

/// The last iterate, including on iteration-limit exhaustion.
/// Energy includes the ferromagnetic reference N/4 for J=1.
template <uni20::Real Real> struct RealState
{
  std::vector<Real> rapidities;
  Real energy = Real{0};
  /// Infinity norm of the logarithmic Bethe equations divided by N (radians).
  Real residual_norm = Real{0};
  /// Number of simultaneous updates; evaluating the initial guess is not an update.
  std::size_t iterations = 0;
  bool converged = false;
  std::vector<uni20::half_int> quantum_numbers;
  uni20::half_int sz;
  /// Exact lattice momentum index: P = 2*pi*momentum_index/sites in [0,2*pi).
  std::size_t momentum_index = 0;
  Real momentum = Real{0};
  /// Roots describe overturned spins relative to the all-down vacuum if true.
  bool spin_reversed = false;
};

// Preserve the original public result name.
template <uni20::Real Real> using GroundState = RealState<Real>;
using QuantumNumbers = std::vector<uni20::half_int>;

namespace detail
{
inline std::int64_t checked_sites(std::size_t sites)
{
  // Leave room for modular momentum arithmetic, without signed overflow.
  if (sites < 2 || sites > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() / 4))
    throw std::invalid_argument("Heisenberg chain requires 2 <= sites <= INT64_MAX/4");
  return static_cast<std::int64_t>(sites);
}

inline std::int64_t sector_roots(std::size_t sites, uni20::half_int sz)
{
  auto const n = checked_sites(sites);
  auto const twice_sz = sz.twice();
  if (twice_sz < -n || twice_sz > n || (n - twice_sz) % 2 != 0)
    throw std::invalid_argument("Sz must lie in [-N/2,N/2] with the same half-integer parity as N/2");
  return (n - (twice_sz < 0 ? -twice_sz : twice_sz)) / 2;
}

inline void validate_quantum_numbers(std::size_t sites, std::span<uni20::half_int const> numbers)
{
  auto const n = checked_sites(sites);
  if (numbers.size() > sites / 2)
    throw std::invalid_argument("finite real-root states require M <= N/2");
  auto const m = static_cast<std::int64_t>(numbers.size());
  auto const bound = n - m - 1; // Conventional all-1-string window, in doubled units.
  for (std::size_t i = 0; i < numbers.size(); ++i)
  {
    auto const twice = numbers[i].twice();
    if (twice < -bound || twice > bound || (twice + bound) % 2 != 0)
      throw std::invalid_argument("Bethe quantum number outside the supported finite-real window or wrong parity");
    if (i != 0 && numbers[i - 1] >= numbers[i])
      throw std::invalid_argument("Bethe quantum numbers must be strictly increasing");
  }
}

inline std::size_t momentum_index(std::size_t sites, std::span<uni20::half_int const> numbers)
{
  auto const n = static_cast<std::int64_t>(sites);
  auto const modulus = 2 * n;
  // 2*mom_index = N*M - sum(2*I), modulo 2*N. Reduce at each step
  // rather than forming the potentially overflowing product or sum.
  std::int64_t index = numbers.size() % 2 == 0 ? 0 : n;
  for (auto const number : numbers)
  {
    index = (index - number.twice()) % modulus;
    if (index < 0)
      index += modulus;
  }
  return static_cast<std::size_t>(index / 2);
}

// Keep roundoff in the phase and energy sums from growing with the root count.
template <typename Real> class CompensatedSum
{
public:
  void add(Real value)
  {
    Real const corrected = value - correction_;
    Real const next = sum_ + corrected;
    correction_ = (next - sum_) - corrected;
    sum_ = next;
  }

  Real value() const { return sum_; }

private:
  Real sum_ = Real{0};
  Real correction_ = Real{0};
};
} // namespace detail

/// Solve a specified periodic XXX highest-weight state, J=1, h=0.
/// M <= N/2 sorted quantum numbers, |2I| <= N-M-1 and 2I == N-M-1 mod 2.
/// This is the conventional all-1-string window, not a complete classification
/// of every finite real solution. Complex roots and infinite-root descendants
/// are not supported. The empty set gives the fully polarized state.
///
/// The optional initial roots are a numerical guess, not part of state identity.
/// O(M^2) work per simultaneous sweep and O(M) storage. Convergence tests
/// max|F|/N; it does not bound the energy error. Budget exhaustion returns the
/// last iterate with converged=false. A zero budget only evaluates the guess.
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

  std::size_t const roots = numbers.size();
  Real const n = static_cast<Real>(sites);
  Real const pi = Real{4} * atan(Real{1});
  RealState<Real> result;
  result.quantum_numbers.assign(numbers.begin(), numbers.end());
  result.sz = uni20::from_twice(static_cast<std::int64_t>(sites - 2 * roots));
  result.momentum_index = detail::momentum_index(sites, numbers);
  result.momentum = Real{2} * pi * (static_cast<Real>(result.momentum_index) / n);
  result.rapidities.resize(roots, Real{0});
  if (!initial_roots.empty())
  {
    if (initial_roots.size() != roots)
      throw std::invalid_argument("initial root count must match the quantum numbers");
    for (Real const root : initial_roots)
      if (!uni20::isfinite(root))
        throw std::invalid_argument("initial roots must be finite");
    result.rapidities.assign(initial_roots.begin(), initial_roots.end());
  }
  std::vector<Real> angles(roots);

  for (;;)
  {
    // Evaluate both the residual and the next fixed-point angles at the same
    // current roots. Never report a residual belonging to the previous iterate.
    result.residual_norm = Real{0};
    for (std::size_t i = 0; i < roots; ++i)
    {
      detail::CompensatedSum<Real> phase;
      for (std::size_t j = 0; j < roots; ++j)
      {
        if (i != j)
          phase.add(atan((result.rapidities[i] - result.rapidities[j]) / Real{2}));
      }
      Real const quantum_number = static_cast<Real>(numbers[i].twice()) / Real{2};
      angles[i] = (pi * quantum_number + phase.value()) / n;
      Real const residual = Real{2} * abs(atan(result.rapidities[i]) - angles[i]);
      if (!uni20::isfinite(residual))
        throw std::runtime_error("nonfinite Heisenberg equation residual");
      result.residual_norm = std::max(result.residual_norm, residual);
    }

    result.converged = result.residual_norm <= options.residual_tolerance;
    if (result.converged || result.iterations == options.max_iterations)
      break;

    // All angles were computed before any roots change: this is a Jacobi update.
    for (std::size_t i = 0; i < roots; ++i)
    {
      result.rapidities[i] = tan(angles[i]);
      if (!uni20::isfinite(result.rapidities[i]))
        throw std::runtime_error("nonfinite Heisenberg rapidity");
    }
    ++result.iterations;
  }

  detail::CompensatedSum<Real> magnons;
  for (Real const z : result.rapidities)
    magnons.add(-Real{2} / (Real{1} + z * z));
  result.energy = n / Real{4} + magnons.value();
  if (!uni20::isfinite(result.energy))
    throw std::runtime_error("nonfinite Heisenberg energy");
  return result;
}

/// Consecutive quantum numbers of the lowest state in an Sz sector.
/// Negative Sz uses global spin reversal, not infinite-root descendants.
/// For odd N and M>0 there are two reflection-related minima; choose the
/// sequence centered at -1/2. Negate and reverse it for the other momentum.
[[nodiscard]] inline QuantumNumbers sector_ground_quantum_numbers(std::size_t sites, uni20::half_int sz)
{
  auto const m = detail::sector_roots(sites, sz);
  QuantumNumbers numbers;
  numbers.reserve(static_cast<std::size_t>(m));
  auto const first = -(m - 1) - static_cast<std::int64_t>(sites % 2);
  for (std::int64_t i = 0; i < m; ++i)
    numbers.push_back(uni20::from_twice(first + 2 * i));
  return numbers;
}

template <uni20::Real Real = double>
[[nodiscard]] RealState<Real> sector_ground_state(std::size_t sites, uni20::half_int sz,
                                                  SolverOptions<Real> const& options = {})
{
  auto const numbers = sector_ground_quantum_numbers(sites, sz);
  auto result = solve_real<Real>(sites, numbers, options);
  result.sz = sz;
  result.spin_reversed = sz.twice() < 0;
  return result;
}

/// One representative per Sz sector, ordered from -N/2 to N/2. Spin-reversed
/// partners reuse the same solve. Retains roots for every result: O(N^2) storage.
template <uni20::Real Real = double>
[[nodiscard]] std::vector<RealState<Real>> sector_ground_states(std::size_t sites,
                                                                SolverOptions<Real> const& options = {})
{
  auto const n = detail::checked_sites(sites);
  std::vector<RealState<Real>> states(sites + 1);
  for (std::size_t m = 0; m <= sites / 2; ++m)
  {
    auto const sz = uni20::from_twice(n - 2 * static_cast<std::int64_t>(m));
    states[sites - m] = sector_ground_state<Real>(sites, sz, options);
    if (m != sites - m)
    {
      states[m] = states[sites - m];
      states[m].sz = -sz;
      states[m].spin_reversed = true;
    }
  }
  return states;
}

/// Periodic XXX zero-field ground state, J=1, N>=2. Even N uses Sz=0;
/// odd N returns one Sz=1/2 representative. For N=2 the bond is counted twice.
template <uni20::Real Real = double>
[[nodiscard]] GroundState<Real> ground_state(std::size_t sites, SolverOptions<Real> const& options = {})
{
  return sector_ground_state<Real>(sites, uni20::from_twice(static_cast<std::int64_t>(sites % 2)), options);
}

/// Odd N, Sz=1/2 one-spinon state: occupy M=(N-1)/2 of the M+1 slots
/// -M/2, -M/2+1, ..., M/2, leaving exactly the specified hole.
[[nodiscard]] inline QuantumNumbers one_spinon_quantum_numbers(std::size_t sites, uni20::half_int hole)
{
  auto const n = detail::checked_sites(sites);
  if (sites % 2 == 0)
    throw std::invalid_argument("one-spinon states require odd N >= 3");
  auto const m = (n - 1) / 2;
  if (hole.twice() < -m || hole.twice() > m || (hole.twice() + m) % 2 != 0)
    throw std::invalid_argument("spinon hole outside the allowed window or wrong parity");
  QuantumNumbers numbers;
  numbers.reserve(static_cast<std::size_t>(m));
  for (std::int64_t i = 0; i <= m; ++i)
  {
    auto const number = uni20::from_twice(-m + 2 * i);
    if (number != hole)
      numbers.push_back(number);
  }
  return numbers;
}

template <uni20::Real Real> struct SpinonState
{
  RealState<Real> state;
  uni20::half_int hole;
  /// k=pi/2-2*pi*I_h/N in (0,pi). P=pi*M+pi/2-k modulo 2*pi.
  /// A finite-size labeling convention; not momentum relative to the odd-N GS.
  Real spinon_momentum = Real{0};
  /// E_N - N*(1/4-log(2)), J=1. Includes finite-size corrections.
  Real bulk_subtracted_energy = Real{0};
};

template <uni20::Real Real = double>
[[nodiscard]] SpinonState<Real> one_spinon_state(std::size_t sites, uni20::half_int hole,
                                                 SolverOptions<Real> const& options = {})
{
  auto const numbers = one_spinon_quantum_numbers(sites, hole);
  SpinonState<Real> result;
  result.state = solve_real<Real>(sites, numbers, options);
  result.hole = hole;
  auto const n = static_cast<std::int64_t>(sites);
  result.spinon_momentum =
      bethe::detail::pi<Real>() * (static_cast<Real>(n - 2 * hole.twice()) / (Real{2} * static_cast<Real>(n)));
  result.bulk_subtracted_energy = result.state.energy - static_cast<Real>(sites) * bulk_energy_density<Real>();
  return result;
}

/// All (N+1)/2 one-spinon states with Sz=1/2, in ascending spinon momentum.
/// Not the complete Sz=1/2 spectrum. Retains all roots: O(N^2) storage.
template <uni20::Real Real = double>
[[nodiscard]] std::vector<SpinonState<Real>> one_spinon_branch(std::size_t sites,
                                                               SolverOptions<Real> const& options = {})
{
  auto const n = detail::checked_sites(sites);
  if (sites % 2 == 0)
    throw std::invalid_argument("one-spinon states require odd N >= 3");
  auto const m = (n - 1) / 2;
  std::vector<SpinonState<Real>> branch;
  branch.reserve(static_cast<std::size_t>(m + 1));
  for (std::int64_t i = 0; i <= m; ++i)
    branch.push_back(one_spinon_state<Real>(sites, uni20::from_twice(m - 2 * i), options));
  return branch;
}

} // namespace bethe::heisenberg
