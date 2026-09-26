// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz.hpp>

namespace bethe::potts
{
/// Critical ferromagnetic three-state clock chain, periodic, J=1:
/// H=-sum(X+X^dagger+Z_j Z_{j+1}^dagger+Z_j^dagger Z_{j+1}).
/// L=2 counts the bond twice. Only the vacuum and the charged one-hole family
/// below are selected; arbitrary twisted-XXZ states are NOT Potts states.
template <uni20::Real Real> struct State
{
    Real energy = Real{0}, momentum = Real{0}, residual_norm = Real{0};
    std::size_t sites = 0, momentum_index = 0, iterations = 0;
    int charge = 0;
    bool converged = false;
    /// Auxiliary XXZ roots z=tanh(lambda)/tan(pi/12), NOT direct Potts roots.
    std::vector<Real> auxiliary_roots;
    xxz::QuantumNumbers auxiliary_numbers;
    Real auxiliary_twist = Real{0};
};

namespace detail
{
inline void validate_sites(std::size_t sites)
{
  if (sites < 2 || sites > std::size_t(std::numeric_limits<std::int64_t>::max() / 8))
    throw std::invalid_argument("critical Potts chain requires 2<=sites<=INT64_MAX/8");
}

template <uni20::Real Real>
State<Real> solve(std::size_t sites, int charge, std::size_t momentum_index, bethe::SolverOptions<Real> const& options)
{
  validate_sites(sites);
  Real const pi = Real{4} * std::atan(Real{1}), root3 = std::sqrt(Real{3});
  xxz::QuantumNumbers numbers;
  numbers.reserve(sites);
  auto const l = static_cast<std::int64_t>(sites);
  if (charge == 0)
    for (std::int64_t j = 0; j < l; ++j)
      numbers.push_back(uni20::from_twice(2 * j - l + 1));
  else
    // L occupied slots in the L+1-slot interval. The omitted slot h labels
    // the physical Potts momentum 2*pi*h/L. h=L is the duplicate endpoint
    // at momentum zero, not an extra physical state of either charge.
    for (std::size_t j = 0; j <= sites; ++j)
      if (j != momentum_index) numbers.push_back(uni20::from_twice(2 * static_cast<std::int64_t>(j) - l - 1));
  Real const twist = charge == 0 ? pi / Real{3} : pi;
  auto roots = xxz::detail::solve_validated<Real>(2 * sites, root3 / Real{2}, numbers, options, {}, twist);
  State<Real> result;
  result.sites = sites;
  result.charge = charge;
  result.momentum_index = momentum_index;
  result.momentum = Real{2} * pi * (Real(momentum_index) / Real(sites));
  result.energy = Real{2} * root3 * roots.energy + Real(sites) / Real{2};
  if (!uni20::isfinite(result.energy)) throw std::runtime_error("nonfinite Potts energy");
  result.residual_norm = roots.residual_norm;
  result.iterations = roots.iterations;
  result.converged = roots.converged;
  result.auxiliary_roots = std::move(roots.rapidities);
  result.auxiliary_numbers = std::move(roots.quantum_numbers);
  result.auxiliary_twist = twist;
  return result;
}
} // namespace detail

/// The charge-zero, momentum-zero vacuum. Exact finite-size Bethe equations,
/// not a CFT energy formula or a string-center approximation.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> ground_state(std::size_t sites, bethe::SolverOptions<Real> const& options = {})
{
  return detail::solve<Real>(sites, 0, 0, options);
}

/// Selected charged one-hole level. charge=+1 or -1 labels eigenvalue
/// exp(2*pi*i*charge/3) of the global clock rotation; both have equal energy.
/// h=0,...,L-1 gives exactly one member per momentum per charge, NOT all levels.
/// h=0 tends to x=2/15 and h=1,L-1 to its first x=17/15 descendants.
/// No all-spectrum completeness or general sector-minimum claim is made.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> charged_one_hole(std::size_t sites, int charge, std::size_t momentum_index,
                                           bethe::SolverOptions<Real> const& options = {})
{
  detail::validate_sites(sites);
  if (charge != 1 && charge != -1) throw std::invalid_argument("charged Potts branch requires charge=+1 or -1");
  if (momentum_index >= sites) throw std::invalid_argument("Potts momentum index must satisfy 0<=k<L");
  return detail::solve<Real>(sites, charge, momentum_index, options);
}

template <uni20::Real Real = double> [[nodiscard]] Real velocity() { return Real{3} * std::sqrt(Real{3}) / Real{2}; }

template <uni20::Real Real = double> [[nodiscard]] Real bulk_energy_density()
{
  return -Real{4} / Real{3} - std::sqrt(Real{3}) / (Real{2} * std::atan(Real{1}));
}
} // namespace bethe::potts
