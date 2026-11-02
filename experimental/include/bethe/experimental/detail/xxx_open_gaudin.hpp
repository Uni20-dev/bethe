// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/heisenberg_open.hpp>

namespace bethe::experimental::xxx::open::boundary_detail
{
namespace model = bethe::heisenberg::open;

template <uni20::Real R> void validate_singlet_shape(std::size_t sites, model::RealState<R> const& ground)
{
  (void)bethe::heisenberg::detail::checked_sites(sites);
  if (sites % 2 || ground.rapidities.size() != sites / 2 ||
      ground.quantum_numbers != model::sector_ground_quantum_numbers(sites, uni20::half_int{0}) ||
      ground.sz != uni20::half_int{0} || ground.spin_reversed)
    throw std::invalid_argument("boundary singlet observables require an even free-end XXX ground state");
}

/// Original equations in solve_real's max|F|/(2*N) units; no cached residual.
/// Empty excited roots use the prepared ground's explicit working precision.
template <uni20::Real R>
std::optional<R> measured_residual(std::size_t sites, model::RealState<R> const& state, R const& zero)
{
  using uni20::math::abs;
  using uni20::math::atan;
  auto const& z = state.rapidities;
  if (state.quantum_numbers.size() != z.size()) throw std::invalid_argument("root and quantum-number counts differ");
  model::detail::validate_quantum_numbers(sites, state.quantum_numbers);
  for (auto const& root : z)
    bethe::detail::require_same_precision(root, zero);
  for (std::size_t i = 0; i < z.size(); ++i)
    if (!uni20::isfinite(z[i]) || z[i] <= zero || (i && z[i] <= z[i - 1])) return {};
  R result = zero;
  R const pi = bethe::detail::pi_like(zero), n = uni20::scalar_like(zero, sites);
  for (std::size_t i = 0; i < z.size(); ++i)
  {
    bethe::detail::CompensatedSum<R> phase(zero);
    for (std::size_t j = 0; j < z.size(); ++j)
      if (i != j)
      {
        phase.add(atan((z[i] - z[j]) / 2));
        phase.add(atan(z[i] / 2 + z[j] / 2));
      }
    R const angle = (pi * uni20::scalar_like(zero, state.quantum_numbers[i].twice()) / 2 + phase.value()) / (2 * n);
    R const residual = 2 * abs(atan(z[i]) - angle);
    if (!uni20::isfinite(residual)) return {};
    result = std::max(result, residual);
  }
  return result;
}

/// Row-major Jacobian with respect to standard ABA u=z/2, for norms/bond response.
template <uni20::Real R> std::vector<R> gaudin_matrix(std::size_t sites, std::span<R const> u, R const& zero)
{
  auto const m = u.size();
  std::vector<R> matrix(m * m, zero);
  R const one = uni20::scalar_like(zero, 1), n = uni20::scalar_like(zero, sites);
  for (std::size_t i = 0; i < m; ++i)
  {
    bethe::detail::CompensatedSum<R> diagonal(zero);
    diagonal.add(2 * n / (u[i] * u[i] + one / 4));
    for (std::size_t j = 0; j < m; ++j)
      if (i != j)
      {
        R const d = u[i] - u[j], s = u[i] + u[j];
        R const direct = 2 * one / (one + d * d), reflected = 2 * one / (one + s * s);
        diagonal.add(-direct - reflected);
        matrix[i * m + j] = direct - reflected;
      }
    matrix[i * m + i] = diagonal.value();
  }
  return matrix;
}
} // namespace bethe::experimental::xxx::open::boundary_detail
