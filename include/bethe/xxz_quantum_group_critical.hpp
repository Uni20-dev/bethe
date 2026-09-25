// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/xxz_quantum_group.hpp>
#include <optional>

namespace bethe::xxz::quantum_group::critical
{
enum class Status
{
  converged,
  iteration_limit,
  precision_limit
};
template <uni20::Real Real> struct RealState
{
    std::size_t sites = 0, iterations = 0;
    Real delta{}, residual_norm{};
    QuantumNumbers quantum_numbers;
    std::vector<Real> scaled_roots; // z=tanh(lambda)/tan(gamma/2).
    std::vector<Real> rapidities;   // Positive finite lambda, Delta=cos(gamma).
    std::optional<Real> energy, energy_shift;
    bool converged = false;
    Status status = Status::iteration_limit;
};

/// Positive finite-real-root family of the non-Hermitian quantum-group chain:
/// H=sum(SxSx+SySy+Delta SzSz)+i*sqrt(1-Delta^2)/2*(Sz_1-Sz_N).
/// 0<Delta<1; no ordinary SU(2) multiplicity or Jordan structure is inferred.
/// Log residual: max|2N theta_1-sum(theta_2^-+theta_2^+)-2pi I|/(2N).
/// Necessary label-window bounds do not guarantee a finite regular solution.
template <uni20::Real Real = double>
RealState<Real> solve_real(std::size_t sites, Real delta, std::span<uni20::half_int const> numbers,
                           SolverOptions<Real> const& options = {})
{
  xxz::detail::checked_sites(sites);
  if (!uni20::isfinite(delta) || !(delta > Real{0}) || !(delta < Real{1}))
    throw std::invalid_argument("critical quantum-group real roots require 0<Delta<1");
  if (!uni20::isfinite(options.residual_tolerance) || !(options.residual_tolerance > Real{0}))
    throw std::invalid_argument("residual tolerance must be finite and positive");
  auto const m = numbers.size();
  if (m > sites / 2) throw std::invalid_argument("critical quantum-group real roots require M<=N/2");
  Real const pi = Real{4} * std::atan(Real{1}), gamma = std::acos(delta);
  Real const bound = Real(sites - m + 1) - Real(sites - 2 * m + 2) * gamma / pi;
  Real const margin = Real{32} * uni20::numeric_limits<Real>::epsilon() * std::max(Real{1}, bound);
  for (std::size_t i = 0; i < m; ++i)
  {
    auto const twice = numbers[i].twice();
    if (twice <= 0 || twice % 2 || (i && numbers[i] <= numbers[i - 1]) || Real(twice) / Real{2} > Real(sites - m) ||
        Real(twice) / Real{2} >= bound - margin)
      throw std::invalid_argument("labels must be increasing positive integers below the finite-root threshold");
  }
  RealState<Real> out;
  out.sites = sites;
  out.delta = delta;
  out.quantum_numbers.assign(numbers.begin(), numbers.end());
  Real const plus = Real{1} + delta, minus = Real{1} - delta;
  std::vector<Real> next(m), angles(m);
  out.scaled_roots.resize(m);
  for (std::size_t i = 0; i < m; ++i)
    out.scaled_roots[i] = std::tan(pi * Real(numbers[i].twice()) / (Real{4} * Real(sites)));
  auto physical = [&](std::vector<Real> const& z) {
    for (std::size_t i = 0; i < m; ++i)
      if (!uni20::isfinite(z[i]) || z[i] <= Real{0} || minus * z[i] * z[i] >= plus || (i && z[i] <= z[i - 1]))
        return false;
    return true;
  };
  if (!physical(out.scaled_roots)) throw std::invalid_argument("seed lies outside the finite-real branch");
  for (;;)
  {
    out.residual_norm = Real{0};
    for (std::size_t i = 0; i < m; ++i)
    {
      Real const z = out.scaled_roots[i];
      bethe::detail::CompensatedSum<Real> phase;
      for (std::size_t j = 0; j < m; ++j)
        if (j != i)
        {
          Real const w = out.scaled_roots[j];
          phase.add(std::atan(delta * (z - w) / (plus - minus * z * w)));
          phase.add(std::atan(delta * (z + w) / (plus + minus * z * w)));
        }
      angles[i] = (pi * Real(numbers[i].twice()) / Real{2} + phase.value()) / (Real{2} * Real(sites));
      out.residual_norm = std::max(out.residual_norm, Real{2} * std::abs(std::atan(z) - angles[i]));
    }
    if (out.residual_norm <= options.residual_tolerance)
    {
      out.converged = true;
      out.status = Status::converged;
      break;
    }
    if (out.iterations == options.max_iterations) break;
    for (std::size_t i = 0; i < m; ++i)
      next[i] = std::tan(angles[i]);
    if (!physical(next) || next == out.scaled_roots)
    {
      out.status = Status::precision_limit;
      break;
    }
    out.scaled_roots.swap(next);
    ++out.iterations;
  }
  Real const t = std::sqrt(minus / plus);
  for (Real z : out.scaled_roots)
  {
    auto const lambda = std::atanh(t * z);
    out.rapidities.push_back(lambda);
    if (!uni20::isfinite(lambda))
    {
      out.converged = false;
      out.status = Status::precision_limit;
    }
  }
  if (out.converged)
  {
    bethe::detail::CompensatedSum<Real> shift;
    for (Real z : out.scaled_roots)
      shift.add(-(plus - minus * z * z) / (Real{1} + z * z));
    out.energy_shift = shift.value();
    out.energy = Real(sites - 1) * delta / Real{4} + *out.energy_shift;
  }
  return out;
}

/// Consecutive-label sea for M=(N-ell)/2. Not a claim of spectrum completeness
/// or generic-q multiplicities at roots of unity.
template <uni20::Real Real = double>
RealState<Real> sea_state(std::size_t sites, Real delta, std::size_t through_lines = 0,
                          SolverOptions<Real> const& options = {})
{
  auto const m = quantum_group::detail::real_sector_roots(sites, through_lines);
  return solve_real<Real>(sites, delta, quantum_group::detail::consecutive(m), options);
}
} // namespace bethe::xxz::quantum_group::critical
