// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include <bethe/heisenberg.hpp>
#include <bethe/solver.hpp>

#include <utility>

namespace bethe::xxz
{
template <uni20::Real Real> using SolverOptions = bethe::SolverOptions<Real>;
using QuantumNumbers = std::vector<uni20::half_int>;

/// Periodic XXZ finite-real-root state, J=1, h=0, Delta >= 0.
/// Delta>1 currently supports sector ground states, not general excitations.
/// Not an SU(2) highest-weight/multiplet classification away from Delta=1.
template <uni20::Real Real> struct RealState
{
    /// z=tanh(lambda)/tan(gamma/2), Delta=cos(gamma); z=2*lambda_XXX at Delta=1.
    /// For Delta>1: z=tan(lambda)/tanh(eta/2), Delta=cosh(eta).
    /// This scaled coordinate, not lambda itself, is used at every anisotropy.
    std::vector<Real> rapidities;
    QuantumNumbers quantum_numbers;
    Real delta = Real{0};
    Real energy = Real{0};
    /// max|F|/N in radians, evaluated at the returned roots.
    Real residual_norm = Real{0};
    std::size_t iterations = 0;
    bool converged = false;
    uni20::half_int sz;
    std::size_t momentum_index = 0;
    Real momentum = Real{0};
    bool spin_reversed = false;
};

namespace detail
{
inline std::int64_t checked_sites(std::size_t sites)
{
  if (sites < 2 || sites > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() / 4))
    throw std::invalid_argument("XXZ chain requires 2 <= sites <= INT64_MAX/4");
  return static_cast<std::int64_t>(sites);
}

inline std::size_t sector_roots(std::size_t sites, uni20::half_int sz)
{
  auto const n = checked_sites(sites);
  auto const twice = sz.twice();
  if (twice < -n || twice > n || (n - twice) % 2 != 0)
    throw std::invalid_argument("Sz must lie in [-N/2,N/2] with the same half-integer parity as N/2");
  return static_cast<std::size_t>((n - (twice < 0 ? -twice : twice)) / 2);
}

template <uni20::Real Real> void validate_delta(Real delta)
{
  if (!uni20::isfinite(delta) || delta < Real{0} || delta > Real{1})
    throw std::invalid_argument("finite-chain XXZ solver requires finite Delta in [0,1]");
}

template <uni20::Real Real> void validate_ground_delta(Real delta)
{
  if (!uni20::isfinite(delta) || delta < Real{0})
    throw std::invalid_argument("periodic XXZ ground states require finite Delta >= 0");
}
} // namespace detail

/// Supported symmetric window of unit-spaced quantum numbers I. M=0 uses
/// an empty window. Restrict to the conventional XXX window AND |I|<I_infinity,
/// where 2*I_infinity=N/2+(N-2*M+2)*asin(Delta)/pi. This is a supported family,
/// not a classification of all real or complex XXZ solutions.
/// At 0<Delta<1 exclude a 32-epsilon relative band at the infinity threshold;
/// an exactly marginal label must not enter through transcendental roundoff.
struct RealQuantumNumberWindow
{
    uni20::half_int first;
    uni20::half_int last;
    std::size_t slots = 0;
};

template <uni20::Real Real = double>
[[nodiscard]] RealQuantumNumberWindow real_quantum_number_window(std::size_t sites, Real delta, uni20::half_int sz)
{
  using std::asin;
  using std::atan;
  detail::validate_delta(delta);
  auto const m = detail::sector_roots(sites, sz);
  if (m == 0) return {};
  auto const n = static_cast<std::int64_t>(sites);
  auto const bound = n - static_cast<std::int64_t>(m) - 1;
  auto const parity = bound % 2;
  Real const infinity = Real(n) / Real{2} + Real(sites - 2 * m + 2) * asin(delta) / (Real{4} * atan(Real{1}));
  Real const margin = Real{32} * uni20::numeric_limits<Real>::epsilon() * std::max(Real{1}, infinity);
  auto allowed = [&](std::int64_t q) {
    if (delta == Real{1}) return true;
    if (delta == Real{0}) return 2 * q < n; // Exact integer endpoint, |I|<N/4.
    return Real(q) < infinity - margin;
  };
  // Binary search needs O(log N) work and no root/slot allocation.
  std::int64_t lo = 0, hi = (bound - parity) / 2 + 1;
  while (lo < hi)
  {
    auto const mid = lo + (hi - lo) / 2;
    if (allowed(parity + 2 * mid))
      lo = mid + 1;
    else
      hi = mid;
  }
  if (lo == 0) return {};
  auto const maximum = parity + 2 * (lo - 1);
  return {uni20::from_twice(-maximum), uni20::from_twice(maximum), static_cast<std::size_t>(maximum + 1)};
}

namespace detail
{
template <uni20::Real Real>
void validate_quantum_numbers(std::size_t sites, Real delta, std::span<uni20::half_int const> numbers)
{
  auto const n = checked_sites(sites);
  if (numbers.size() > sites / 2) throw std::invalid_argument("finite real-root states require M <= N/2");
  auto const sz = uni20::from_twice(n - 2 * static_cast<std::int64_t>(numbers.size()));
  auto const window = xxz::real_quantum_number_window(sites, delta, sz);
  for (std::size_t i = 0; i < numbers.size(); ++i)
  {
    auto const q = numbers[i].twice();
    if (window.slots == 0 || numbers[i] < window.first || numbers[i] > window.last ||
        (q - window.first.twice()) % 2 != 0)
      throw std::invalid_argument("Bethe quantum number outside the supported XXZ finite-real window or wrong parity");
    if (i != 0 && numbers[i - 1] >= numbers[i])
      throw std::invalid_argument("Bethe quantum numbers must be strictly increasing");
  }
}

inline std::size_t momentum_index(std::size_t sites, std::span<uni20::half_int const> numbers)
{
  auto const n = static_cast<std::int64_t>(sites);
  auto const modulus = 2 * n;
  std::int64_t index = numbers.size() % 2 == 0 ? 0 : n;
  for (auto const number : numbers)
  {
    index = (index - number.twice()) % modulus;
    if (index < 0) index += modulus;
  }
  return static_cast<std::size_t>(index / 2);
}
} // namespace detail

/// Consecutive labels for the supported sector minimum. For odd N, select
/// the sequence centered at -1/2, one of two reflection-related minima.
/// Excited-state labels are validated separately by solve_real.
[[nodiscard]] inline QuantumNumbers sector_ground_quantum_numbers(std::size_t sites, uni20::half_int sz)
{
  auto const m = static_cast<std::int64_t>(detail::sector_roots(sites, sz));
  QuantumNumbers numbers;
  numbers.reserve(static_cast<std::size_t>(m));
  auto const first = -(m - 1) - static_cast<std::int64_t>(sites % 2);
  for (std::int64_t i = 0; i < m; ++i)
    numbers.push_back(uni20::from_twice(first + 2 * i));
  return numbers;
}

namespace detail
{
// The caller validates the label family. Massive ground states use the same
// fixed-point map, but the scattering phase must retain its winding branch.
template <uni20::Real Real>
[[nodiscard]] RealState<Real> solve_validated(std::size_t sites, Real delta, std::span<uni20::half_int const> numbers,
                                              SolverOptions<Real> const& options,
                                              std::span<Real const> initial_roots = {})
{
  using std::abs;
  using std::atan;
  using std::tan;
  auto const m = numbers.size();
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");

  RealState<Real> result;
  result.delta = delta;
  result.sz = uni20::from_twice(static_cast<std::int64_t>(sites - 2 * m));
  if (delta == Real{1})
  {
    // Exact endpoint only: never snap nearby anisotropies to the XXX model.
    auto state = heisenberg::solve_real<Real>(sites, numbers, options, initial_roots);
    result.rapidities = std::move(state.rapidities);
    result.quantum_numbers = std::move(state.quantum_numbers);
    result.energy = state.energy;
    result.residual_norm = state.residual_norm;
    result.iterations = state.iterations;
    result.converged = state.converged;
    result.momentum_index = state.momentum_index;
    result.momentum = state.momentum;
    return result;
  }

  result.quantum_numbers.assign(numbers.begin(), numbers.end());
  result.momentum_index = detail::momentum_index(sites, numbers);
  Real const n = static_cast<Real>(sites);
  Real const pi = Real{4} * atan(Real{1});
  result.momentum = Real{2} * pi * (static_cast<Real>(result.momentum_index) / n);
  result.rapidities.resize(m, Real{0});
  std::vector<Real> angles(m);
  Real const plus = Real{1} + delta;
  Real const minus = Real{1} - delta;
  bool const massive = delta > Real{1};
  Real const inverse = massive ? Real{1} / delta : Real{0};
  if (!initial_roots.empty())
  {
    if (initial_roots.size() != m) throw std::invalid_argument("initial root count must match the quantum numbers");
    for (Real const z : initial_roots)
      if (!uni20::isfinite(z) || (!massive && minus * z * z >= plus))
        throw std::invalid_argument("initial roots must be finite and inside the XXZ real-rapidity branch");
    result.rapidities.assign(initial_roots.begin(), initial_roots.end());
  }
  for (;;)
  {
    result.residual_norm = Real{0};
    for (std::size_t i = 0; i < m; ++i)
    {
      bethe::detail::CompensatedSum<Real> phase;
      if (delta != Real{0})
        for (std::size_t j = 0; j < m; ++j)
          if (i != j)
          {
            if (massive)
            {
              // Divide both atan2 arguments by Delta to avoid overflow at
              // large anisotropy. A negative denominator is NOT a pole exit.
              phase.add(
                  std::atan2(result.rapidities[i] - result.rapidities[j],
                             Real{1} + inverse + (Real{1} - inverse) * result.rapidities[i] * result.rapidities[j]));
            }
            else
            {
              Real const denominator = plus - minus * result.rapidities[i] * result.rapidities[j];
              if (denominator <= Real{0}) throw std::runtime_error("XXZ iterate left the finite real-rapidity branch");
              phase.add(atan(delta * (result.rapidities[i] - result.rapidities[j]) / denominator));
            }
          }
      Real const number = static_cast<Real>(result.quantum_numbers[i].twice()) / Real{2};
      angles[i] = (pi * number + phase.value()) / n;
      Real const residual = Real{2} * abs(atan(result.rapidities[i]) - angles[i]);
      if (!uni20::isfinite(residual)) throw std::runtime_error("nonfinite XXZ equation residual");
      result.residual_norm = std::max(result.residual_norm, residual);
    }
    result.converged = result.residual_norm <= options.residual_tolerance;
    if (result.converged || result.iterations == options.max_iterations) break;
    for (std::size_t i = 0; i < m; ++i)
    {
      if (massive && abs(angles[i]) >= pi / Real{2})
        throw std::runtime_error("XXZ iterate reached the massive rapidity endpoint");
      result.rapidities[i] = tan(angles[i]);
      Real const z = result.rapidities[i];
      if (!uni20::isfinite(z) || (!massive && minus * z * z >= plus))
        throw std::runtime_error("XXZ iterate left the finite real-rapidity branch");
    }
    ++result.iterations;
  }
  bethe::detail::CompensatedSum<Real> magnons;
  for (Real const z : result.rapidities)
    if (massive)
      magnons.add((z * z - Real{1}) / (Real{1} + z * z));
    else
      magnons.add(-(plus - minus * z * z) / (Real{1} + z * z));
  result.energy = massive ? (n / Real{4} - Real(m)) * delta + magnons.value() : n * delta / Real{4} + magnons.value();
  if (!uni20::isfinite(result.energy)) throw std::runtime_error("nonfinite XXZ energy");
  return result;
}
} // namespace detail

/// H=sum_i (Sx_i Sx_{i+1}+Sy_i Sy_{i+1}+Delta Sz_i Sz_{i+1}), periodic.
/// N=2 counts the bond twice, consistently with the XXX solver.
/// Explicit-label solves remain restricted to 0<=Delta<=1; use the ground-state
/// entry points for Delta>1, whose excitation labels are not classified here.
/// The normalized logarithmic equation in the stored coordinate is
/// F_i/N = 2 atan(z_i) - (2*pi*I_i + sum_{j!=i} 2 atan(A_ij))/N,
/// A_ij = Delta*(z_i-z_j)/[1+Delta-(1-Delta)*z_i*z_j].
/// The iteration needs no inverse hyperbolic functions. Starts with z=0 unless
/// supplied with finite, strictly in-branch initial roots (a numerical guess).
/// Jacobi updates, O(M^2) work per sweep and O(M) storage.
/// Exhaustion returns a consistent last iterate with converged=false.
template <uni20::Real Real = double>
[[nodiscard]] RealState<Real> solve_real(std::size_t sites, Real delta, std::span<uni20::half_int const> numbers,
                                         SolverOptions<Real> const& options = {},
                                         std::span<Real const> initial_roots = {})
{
  detail::validate_quantum_numbers(sites, delta, numbers);
  return detail::solve_validated<Real>(sites, delta, numbers, options, initial_roots);
}

/// Lowest state in an Sz sector for finite Delta>=0. Negative Sz uses spin reversal.
/// At Delta>1, theta_2/2=atan2(z_i-z_j,1+1/Delta+(1-1/Delta)*z_i*z_j).
template <uni20::Real Real = double>
[[nodiscard]] RealState<Real> sector_ground_state(std::size_t sites, Real delta, uni20::half_int sz,
                                                  SolverOptions<Real> const& options = {})
{
  detail::validate_ground_delta(delta);
  auto const numbers = xxz::sector_ground_quantum_numbers(sites, sz);
  auto result = delta > Real{1} ? detail::solve_validated<Real>(sites, delta, numbers, options)
                                : xxz::solve_real<Real>(sites, delta, numbers, options);
  result.sz = sz;
  result.spin_reversed = sz.twice() < 0;
  return result;
}

/// One representative per Sz, ordered -N/2,...,N/2. Reuses spin reversal;
/// these are not SU(2) multiplets. Retaining all roots uses O(N^2) storage.
template <uni20::Real Real = double>
[[nodiscard]] std::vector<RealState<Real>> sector_ground_states(std::size_t sites, Real delta,
                                                                SolverOptions<Real> const& options = {})
{
  auto const n = detail::checked_sites(sites);
  detail::validate_ground_delta(delta);
  std::vector<RealState<Real>> states(sites + 1);
  for (std::size_t m = 0; m <= sites / 2; ++m)
  {
    auto const sz = uni20::from_twice(n - 2 * static_cast<std::int64_t>(m));
    states[sites - m] = xxz::sector_ground_state<Real>(sites, delta, sz, options);
    if (m != sites - m)
    {
      states[m] = states[sites - m];
      states[m].sz = -sz;
      states[m].spin_reversed = true;
    }
  }
  return states;
}

/// Even N uses Sz=0; odd N returns one Sz=1/2 ground-state representative.
template <uni20::Real Real = double>
[[nodiscard]] RealState<Real> ground_state(std::size_t sites, Real delta, SolverOptions<Real> const& options = {})
{
  return xxz::sector_ground_state<Real>(sites, delta, uni20::from_twice(static_cast<std::int64_t>(sites % 2)), options);
}
} // namespace bethe::xxz
