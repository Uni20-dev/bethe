// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz_regularity.hpp>
#include <optional>

namespace bethe::xxz::detail
{
enum class PhantomReductionStatus
{
  regular_reduced_equations,
  no_endpoint_factor,
  phase_mismatch,
  unresolved_reduced_state,
  nonfinite
};

/// Necessary mixed-phantom checks, NOT a nonzero-lift/ground-state certificate.
template <uni20::Real Real> struct PhantomReduction
{
    PhantomReductionStatus status = PhantomReductionStatus::nonfinite;
    std::size_t phantom_count = 0;
    int chirality = 0;                     // Scaled endpoint z=chirality*sqrt((1+Delta)/(1-Delta)).
    std::vector<Real> finite_coefficients; // Same affine coordinate as the input.
    std::complex<Real> equation_rotation{1};
    Real reconstruction_error = uni20::numeric_limits<Real>::infinity();
    Real commensurability_error = uni20::numeric_limits<Real>::infinity();
    Real energy_difference = uni20::numeric_limits<Real>::infinity();
    std::optional<PolynomialRegularity<Real>> finite_regularity;
};

template <uni20::Real Real> std::complex<Real> phantom_phase_power(std::complex<Real> phase, std::size_t exponent)
{
  std::complex<Real> result{1};
  for (; exponent; exponent >>= 1, phase *= phase)
    if (exponent & 1) result *= phase;
  return result;
}

/// Remove an explicitly requested, one-sided cluster of p infinite roots.
/// The remaining r=M-p roots must solve a twisted problem with
/// exp(i*phi)=exp(-2*i*chirality*p*gamma), Delta=cos(gamma), while
/// exp(i*(N-2*r)*gamma)=1 is a SEPARATE condition. See Popkov et al. (2021),
/// Eq. (15), and docs/xxz-phantom.md. No nearby coupling is substituted.
template <uni20::Real Real>
PhantomReduction<Real>
reduce_phantom_polynomial(PolynomialBetheSystem<Real> const& system, std::span<Real const> c, Real delta,
                          std::size_t phantom_count, int chirality,
                          Real factor_tolerance = Real{128} * uni20::numeric_limits<Real>::epsilon(),
                          Real residual_tolerance = Real{32} * uni20::numeric_limits<Real>::epsilon())
{
  if (!phantom_count || phantom_count >= system.order)
    throw std::invalid_argument("mixed phantom reduction requires 0 < phantom count < M");
  if (chirality != -1 && chirality != 1) throw std::invalid_argument("phantom chirality must be -1 or +1");
  if (!uni20::isfinite(factor_tolerance) || factor_tolerance <= Real{0} || !uni20::isfinite(residual_tolerance) ||
      residual_tolerance <= Real{0})
    throw std::invalid_argument("phantom reduction tolerances must be finite and positive");
  if (!uni20::isfinite(delta) || delta <= -Real{1} || delta > Real{0})
    throw std::invalid_argument("XXZ phantom reduction requires -1 < Delta <= 0");
  if (c.size() != system.order) throw std::invalid_argument("phantom coefficient count differs from degree");
  for (Real value : c)
    if (!uni20::isfinite(value)) throw std::invalid_argument("nonfinite phantom polynomial coefficient");
  PhantomReduction<Real> out;
  out.phantom_count = phantom_count;
  out.chirality = chirality;
  Real const endpoint = Real(chirality) * std::sqrt((Real{1} + delta) / (Real{1} - delta));
  Real const x = (endpoint - system.center) / system.coordinate_scale;
  if (!uni20::isfinite(x)) return out;
  auto quotient = std::vector<Real>(c.begin(), c.end());
  for (std::size_t count = 0; count < phantom_count; ++count)
  {
    // Monic synthetic division. Keep the quotient in the same coordinates;
    // reconstruction below checks the complete factorization, not just Q(x).
    std::vector<Real> next(quotient.size() - 1);
    Real value = Real{1};
    for (std::size_t j = quotient.size(); j-- > 1;)
    {
      value = quotient[j] + x * value;
      if (!uni20::isfinite(value)) return out;
      next[j - 1] = value;
    }
    quotient = std::move(next);
  }
  out.finite_coefficients = quotient;
  auto reconstructed = quotient;
  reconstructed.push_back(Real{1});
  for (std::size_t count = 0; count < phantom_count; ++count)
  {
    reconstructed.push_back(Real{0});
    for (std::size_t j = reconstructed.size(); j-- > 0;)
    {
      reconstructed[j] = (j ? reconstructed[j - 1] : Real{0}) - x * reconstructed[j];
      if (!uni20::isfinite(reconstructed[j])) return out;
    }
  }
  out.reconstruction_error = Real{0};
  for (std::size_t j = 0; j < c.size(); ++j)
  {
    Real const scale = std::max({Real{1}, std::abs(c[j]), std::abs(reconstructed[j])});
    out.reconstruction_error = std::max(out.reconstruction_error, std::abs(c[j] / scale - reconstructed[j] / scale));
  }
  if (out.reconstruction_error > factor_tolerance)
  {
    out.status = PhantomReductionStatus::no_endpoint_factor;
    return out;
  }
  Real const sine = std::sqrt(Real{1} + delta) * std::sqrt(Real{1} - delta);
  std::complex<Real> const phase{delta, Real(chirality) * sine};
  auto const r = quotient.size();
  out.commensurability_error = std::abs(phantom_phase_power(phase, system.sites - 2 * r) - std::complex<Real>{1});
  if (!uni20::isfinite(out.commensurability_error)) return out;
  if (out.commensurability_error > factor_tolerance)
  {
    out.status = PhantomReductionStatus::phase_mismatch;
    return out;
  }
  out.equation_rotation = phantom_phase_power(std::conj(phase), phantom_count);
  // Remove accumulated radial roundoff, not the physical phase mismatch.
  Real const magnitude = std::abs(out.equation_rotation);
  if (!(magnitude > Real{0}) || !uni20::isfinite(magnitude)) return out;
  out.equation_rotation /= magnitude;
  PolynomialBetheSystem<Real> const reduced(system.sites, r, system.center, system.coordinate_scale);
  try
  {
    out.finite_regularity =
        check_regular_polynomial<Real>(reduced, quotient, delta, residual_tolerance, out.equation_rotation);
    if (out.finite_regularity->status == RegularityStatus::nonfinite) return out;
    out.energy_difference = std::abs(system.energy(c, delta) - reduced.energy(quotient, delta));
    if (!uni20::isfinite(out.energy_difference)) return out;
    out.status = out.finite_regularity->status == RegularityStatus::regular_on_shell
                     ? PhantomReductionStatus::regular_reduced_equations
                     : PhantomReductionStatus::unresolved_reduced_state;
  }
  catch (std::runtime_error const&)
  {} // No physical conclusion from nonfinite arithmetic.
  return out;
}
} // namespace bethe::xxz::detail
