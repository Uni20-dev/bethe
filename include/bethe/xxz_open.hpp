// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include <bethe/heisenberg_open.hpp>
#include <bethe/xxz.hpp>
#include <bethe/xxz_negative.hpp>
#include <bethe/xxz_open_massive.hpp>

namespace bethe::xxz::open
{
using xxz::QuantumNumbers;
using xxz::RealQuantumNumberWindow;
using xxz::SolverOptions;

/// Free-end XXZ state, J=1, h=0, no boundary fields, 0<=Delta<=1.
/// No translation momentum and no SU(2) multiplet interpretation.
template <uni20::Real Real> struct RealState
{
    /// Same scaled coordinate as periodic XXZ: z=tanh(lambda)/tan(gamma/2).
    /// At Delta=1, z=2*lambda_XXX. Physical roots are strictly positive.
    std::vector<Real> rapidities;
    QuantumNumbers quantum_numbers;
    Real delta = Real{0};
    Real energy = Real{0};
    /// max|F|/(2*N), evaluated at the returned roots, not an energy-error bound.
    Real residual_norm = Real{0};
    std::size_t iterations = 0;
    bool converged = false;
    uni20::half_int sz;
    bool spin_reversed = false;
};

/// Integer labels 1<=I<=N-M intersected with I<I_infinity, where
/// I_infinity=(N+1)/2+(N-2*M+1)*asin(Delta)/pi. Excludes a conservative
/// 32-epsilon relative band at the infinity threshold for 0<Delta<1.
/// The endpoints use exact integer bounds. This is NOT the complete spectrum.
template <uni20::Real Real = double>
[[nodiscard]] RealQuantumNumberWindow real_quantum_number_window(std::size_t sites, Real delta, uni20::half_int sz)
{
  using std::asin;
  using std::atan;
  xxz::detail::validate_delta(delta);
  auto const m = xxz::detail::sector_roots(sites, sz);
  if (m == 0) return {};
  auto const n = static_cast<std::int64_t>(sites);
  auto const bound = sites - m;
  if (delta == Real{1}) return {uni20::half_int{1}, uni20::half_int(static_cast<std::int64_t>(bound)), bound};
  if (delta == Real{0}) return {uni20::half_int{1}, uni20::half_int(n / 2), sites / 2};
  Real const infinity = Real(sites + 1) / Real{2} + Real(sites - 2 * m + 1) * asin(delta) / (Real{4} * atan(Real{1}));
  Real const margin = Real{32} * uni20::numeric_limits<Real>::epsilon() * std::max(Real{1}, infinity);
  std::size_t lo = 1, hi = bound + 1;
  while (lo < hi)
  {
    auto const mid = lo + (hi - lo) / 2;
    if (Real(mid) < infinity - margin)
      lo = mid + 1;
    else
      hi = mid;
  }
  auto const slots = lo - 1;
  if (slots == 0) return {};
  return {uni20::half_int{1}, uni20::half_int(static_cast<std::int64_t>(slots)), slots};
}

namespace detail
{
template <uni20::Real Real>
void validate_quantum_numbers(std::size_t sites, Real delta, std::span<uni20::half_int const> numbers)
{
  auto const n = xxz::detail::checked_sites(sites);
  if (numbers.size() > sites / 2) throw std::invalid_argument("open XXZ finite real-root states require M <= N/2");
  auto const sz = uni20::from_twice(n - 2 * static_cast<std::int64_t>(numbers.size()));
  auto const window = open::real_quantum_number_window(sites, delta, sz);
  for (std::size_t i = 0; i < numbers.size(); ++i)
  {
    if (window.slots == 0 || numbers[i] < window.first || numbers[i] > window.last || numbers[i].twice() % 2 != 0)
      throw std::invalid_argument(
          "open XXZ Bethe quantum numbers must be integers in the supported positive finite-real window");
    if (i != 0 && numbers[i - 1] >= numbers[i])
      throw std::invalid_argument("open XXZ Bethe quantum numbers must be strictly increasing");
  }
}
} // namespace detail

/// H=sum_{i=0}^{N-2}(Sx_i*Sx_{i+1}+Sy_i*Sy_{i+1}+Delta*Sz_i*Sz_{i+1}).
/// With p=1+Delta, q=1-Delta and c=q/p, the logarithmic equation is
/// F_i=4*N*atan(z_i)+4*atan(c*z_i)-2*pi*I_i
///     -2*sum_{j!=i}[atan(Delta*(z_i-z_j)/(p-q*z_i*z_j))
///                  +atan(Delta*(z_i+z_j)/(p+q*z_i*z_j))].
/// Both direct and reflected self scattering are excluded; the boundary term
/// is essential even with zero boundary fields. Energy has reference (N-1)*Delta/4.
///
/// A simultaneous fixed-point update uses atan(z)-atan(c*z)
/// =atan(2*Delta*z/(p+q*z*z)), putting the driving denominator at 2*(N+1).
/// This keeps the update in the finite positive branch for supported labels.
/// Starts at zero unless given finite, nonnegative in-branch initial roots.
/// O(M^2) per update (O(M) at Delta=0), O(M) storage. A zero budget evaluates
/// the guess; exhaustion returns a consistent iterate with converged=false.
template <uni20::Real Real = double>
[[nodiscard]] RealState<Real> solve_real(std::size_t sites, Real delta, std::span<uni20::half_int const> numbers,
                                         SolverOptions<Real> const& options = {},
                                         std::span<Real const> initial_roots = {})
{
  using std::abs;
  using std::atan;
  using std::tan;
  detail::validate_quantum_numbers(sites, delta, numbers);
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");
  RealState<Real> result;
  result.delta = delta;
  result.sz = uni20::from_twice(static_cast<std::int64_t>(sites - 2 * numbers.size()));
  if (delta == Real{1})
  {
    // Exact endpoint only: preserve the existing open XXX solver and diagnostics.
    auto state = heisenberg::open::solve_real<Real>(sites, numbers, options, initial_roots);
    result.rapidities = std::move(state.rapidities);
    result.quantum_numbers = std::move(state.quantum_numbers);
    result.energy = state.energy;
    result.residual_norm = state.residual_norm;
    result.iterations = state.iterations;
    result.converged = state.converged;
    return result;
  }
  Real const n = Real(sites), pi = Real{4} * atan(Real{1});
  Real const plus = Real{1} + delta, minus = Real{1} - delta;
  auto const m = numbers.size();
  result.quantum_numbers.assign(numbers.begin(), numbers.end());
  result.rapidities.resize(m, Real{0});
  if (!initial_roots.empty())
  {
    if (initial_roots.size() != m) throw std::invalid_argument("initial root count must match the quantum numbers");
    for (Real const z : initial_roots)
      if (!uni20::isfinite(z) || z < Real{0} || minus * z * z >= plus)
        throw std::invalid_argument(
            "open XXZ initial roots must be finite, nonnegative, and inside the real-rapidity branch");
    result.rapidities.assign(initial_roots.begin(), initial_roots.end());
  }
  std::vector<Real> angles(m);
  for (;;)
  {
    result.residual_norm = Real{0};
    for (std::size_t i = 0; i < m; ++i)
    {
      Real const z = result.rapidities[i];
      bethe::detail::CompensatedSum<Real> phase;
      if (delta != Real{0})
      {
        for (std::size_t j = 0; j < m; ++j)
          if (i != j)
          {
            Real const w = result.rapidities[j];
            Real const product = minus * z * w;
            if (plus - product <= Real{0})
              throw std::runtime_error("open XXZ iterate left the finite real-rapidity branch");
            phase.add(atan(delta * (z - w) / (plus - product)));
            phase.add(atan(delta * (z + w) / (plus + product)));
          }
        phase.add(Real{2} * atan(Real{2} * delta * z / (plus + minus * z * z)));
      }
      Real const number = Real(numbers[i].twice()) / Real{2};
      angles[i] = (pi * number + phase.value()) / (Real{2} * (n + Real{1}));
      Real const residual = (Real{2} * (n + Real{1}) / n) * abs(atan(z) - angles[i]);
      if (!uni20::isfinite(residual)) throw std::runtime_error("nonfinite open XXZ equation residual");
      result.residual_norm = std::max(result.residual_norm, residual);
    }
    result.converged = result.residual_norm <= options.residual_tolerance;
    if (result.converged || result.iterations == options.max_iterations) break;
    for (std::size_t i = 0; i < m; ++i)
    {
      Real const z = tan(angles[i]);
      if (!uni20::isfinite(z) || z <= Real{0} || minus * z * z >= plus)
        throw std::runtime_error("open XXZ iterate left the finite positive real-rapidity branch");
      result.rapidities[i] = z;
    }
    ++result.iterations;
  }
  bethe::detail::CompensatedSum<Real> magnons;
  for (Real const z : result.rapidities)
    magnons.add(-(plus - minus * z * z) / (Real{1} + z * z));
  result.energy = Real(sites - 1) * delta / Real{4} + magnons.value();
  if (!uni20::isfinite(result.energy)) throw std::runtime_error("nonfinite open XXZ energy");
  return result;
}

/// The sector minimum fills I=1,...,M, M=N/2-|Sz|, for either parity of N.
/// This returns the full label sequence. A massive even zero-Sz GroundState
/// stores its final label separately in boundary_root, not quantum_numbers.
[[nodiscard]] inline QuantumNumbers sector_ground_quantum_numbers(std::size_t sites, uni20::half_int sz)
{
  auto const m = xxz::detail::sector_roots(sites, sz);
  QuantumNumbers numbers;
  numbers.reserve(m);
  for (std::size_t i = 1; i <= m; ++i)
    numbers.emplace_back(static_cast<std::int64_t>(i));
  return numbers;
}

/// Ground states cover Delta>-1; explicit all-real states and scans remain
/// restricted to 0<=Delta<=1. GroundState carries the representation's
/// coordinates and residual convention, including any massive boundary root.
template <uni20::Real Real = double>
[[nodiscard]] GroundState<Real> sector_ground_state(std::size_t sites, Real delta, uni20::half_int sz,
                                                    SolverOptions<Real> const& options = {})
{
  if (!uni20::isfinite(delta) || delta <= -Real{1})
    throw std::invalid_argument("open XXZ ground states require finite Delta > -1");
  if (delta < Real{0})
  {
    auto roots = xxz::detail::negative_ground_roots<Real>(sites, delta, sz, true, options);
    GroundState<Real> result;
    result.rapidities = std::move(roots.rapidities);
    result.log_rapidities = std::move(roots.log_rapidities);
    result.quantum_numbers = std::move(roots.quantum_numbers);
    result.delta = result.root_delta = delta;
    result.energy = roots.energy;
    result.residual_norm = roots.residual_norm;
    result.residual_convention = GroundResidualConvention::negative_rank_scaled;
    result.iterations = roots.iterations;
    result.converged = roots.converged;
    switch (roots.status)
    {
      case xxz::detail::NegativeSolveStatus::converged:
        result.status = GroundSolveStatus::converged;
        break;
      case xxz::detail::NegativeSolveStatus::iteration_limit:
        result.status = GroundSolveStatus::iteration_limit;
        break;
      case xxz::detail::NegativeSolveStatus::stalled:
        result.status = GroundSolveStatus::stalled;
        break;
    }
    result.sz = sz;
    result.spin_reversed = sz.twice() < 0;
    return result;
  }
  if (delta > Real{1}) return massive::sector_ground_state<Real>(sites, delta, sz, options);
  auto const numbers = open::sector_ground_quantum_numbers(sites, sz);
  auto real = open::solve_real<Real>(sites, delta, numbers, options);
  GroundState<Real> result;
  result.rapidities = std::move(real.rapidities);
  result.quantum_numbers = std::move(real.quantum_numbers);
  result.delta = result.root_delta = delta;
  result.energy = real.energy;
  result.residual_norm = real.residual_norm;
  result.iterations = real.iterations;
  result.converged = real.converged;
  result.status = real.converged ? GroundSolveStatus::converged : GroundSolveStatus::iteration_limit;
  result.sz = sz;
  result.spin_reversed = sz.twice() < 0;
  return result;
}

/// One representative per Sz, ordered -N/2,...,N/2. Reuses spin reversal;
/// retains all roots, O(N^2) storage. These are not SU(2) multiplets.
template <uni20::Real Real = double>
[[nodiscard]] std::vector<GroundState<Real>> sector_ground_states(std::size_t sites, Real delta,
                                                                  SolverOptions<Real> const& options = {})
{
  auto const n = xxz::detail::checked_sites(sites);
  if (!uni20::isfinite(delta) || delta <= -Real{1})
    throw std::invalid_argument("open XXZ ground states require finite Delta > -1");
  std::vector<GroundState<Real>> states(sites + 1);
  for (std::size_t m = 0; m <= sites / 2; ++m)
  {
    auto const sz = uni20::from_twice(n - 2 * static_cast<std::int64_t>(m));
    states[sites - m] = open::sector_ground_state<Real>(sites, delta, sz, options);
    if (m != sites - m)
    {
      states[m] = states[sites - m];
      states[m].sz = -sz;
      states[m].spin_reversed = true;
    }
  }
  return states;
}

template <uni20::Real Real = double>
[[nodiscard]] GroundState<Real> ground_state(std::size_t sites, Real delta, SolverOptions<Real> const& options = {})
{
  return open::sector_ground_state<Real>(sites, delta, uni20::from_twice(static_cast<std::int64_t>(sites % 2)),
                                         options);
}
} // namespace bethe::xxz::open
