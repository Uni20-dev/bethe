// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/hubbard_thermo.hpp>

namespace bethe::hubbard::thermo
{
template <uni20::Real Real> struct TwoSpinonContinuum
{
    Real interaction{}, momentum{};
    std::optional<Real> lower_energy, upper_energy;
    std::array<Real, 2> lower_momenta{}, upper_momenta{};
    // Quadrature energy errors exclude propagation of momentum-inversion error.
    Real lower_energy_quad_error = uni20::numeric_limits<Real>::infinity();
    Real upper_energy_quad_error = uni20::numeric_limits<Real>::infinity();
    Real lower_momentum_error = uni20::numeric_limits<Real>::infinity();
    Real upper_momentum_error = uni20::numeric_limits<Real>::infinity();
    std::size_t evaluations = 0, iterations = 0;
    bool converged = false;
    Status status = Status::precision_limit;
};

/// Half-filled, zero-field U>0 two-spinon scattering-continuum edges.
/// Essler-Korepin (cond-mat/9808018), Eq. (18) and following paragraph:
/// lower edge has one endpoint spinon; upper edge has equal rapidities.
/// Delta N=0, so symmetric/unshifted and Hamiltonian/Fermi energies coincide.
/// The quadrature budget is shared across BOTH edge calculations.
template <uni20::Real Real = double>
[[nodiscard]] TwoSpinonContinuum<Real> two_spinon_continuum(Real interaction, Real momentum,
                                                            Options<Real> const& options = {})
{
  detail::validate(interaction, options, Convention::symmetric);
  Real const pi = detail::pi<Real>();
  if (!uni20::isfinite(momentum) || std::abs(momentum) > pi)
    throw std::invalid_argument("two-spinon total momentum must be in [-pi,pi]");
  TwoSpinonContinuum<Real> out;
  out.interaction = interaction;
  out.momentum = momentum;
  Real const q = std::abs(momentum);
  out.lower_momenta = momentum < Real{0} ? std::array<Real, 2>{pi, pi - q} : std::array<Real, 2>{Real{0}, q};
  Real const upper_p = momentum < Real{0} ? pi - q / Real{2} : q / Real{2};
  out.upper_momenta = {upper_p, upper_p};
  // Reflection lets the native small-momentum inversion retain accuracy.
  auto const lower = dispersion(Branch::spinon, interaction, q, Convention::symmetric, options);
  out.evaluations = lower.evaluations;
  out.iterations = lower.iterations;
  out.status = lower.status;
  if (!lower.converged) return out;
  auto remaining = options;
  remaining.max_evaluations -= lower.evaluations;
  auto const upper = dispersion(Branch::spinon, interaction, q / Real{2}, Convention::symmetric, remaining);
  out.evaluations += upper.evaluations;
  out.iterations += upper.iterations;
  out.status = upper.status;
  if (!upper.converged) return out;
  out.lower_energy = lower.energy;
  out.upper_energy = Real{2} * *upper.energy;
  out.lower_energy_quad_error = lower.energy_error;
  out.upper_energy_quad_error = Real{2} * upper.energy_error;
  out.lower_momentum_error = lower.momentum_error;
  out.upper_momentum_error = Real{2} * upper.momentum_error;
  out.converged = true;
  out.status = Status::converged;
  return out;
}
} // namespace bethe::hubbard::thermo
