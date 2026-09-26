// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/sine_gordon_particles.hpp>
#include <bethe/sine_gordon_scattering.hpp>
#include <uni20/common/half_int.hpp>

namespace bethe::sine_gordon
{
enum class BetheYangStatus
{
  converged,
  phase_limit,
  iteration_limit,
  precision_limit
};
template <uni20::Real Real> struct BetheYangOptions
{
    Real tolerance = Real{16384} * uni20::numeric_limits<Real>::epsilon(); // absolute counting residual
    std::size_t max_iterations = 256;
    KernelOptions<Real> phase; // tolerance assigned from counting budget
};
template <uni20::Real Real> struct BetheYangResult
{
    std::optional<Real> rapidity, energy; // roots +/-rapidity; excitation energy, NOT E-L*e_bulk
    Real residual = uni20::numeric_limits<Real>::infinity();
    Real phase_error = uni20::numeric_limits<Real>::infinity();
    Real quadrature_error = uni20::numeric_limits<Real>::infinity();
    Real tail_bound = uni20::numeric_limits<Real>::infinity();
    std::size_t iterations = 0, evaluations = 0;
    KernelStatus phase_status = KernelStatus::converged;
    BetheYangStatus status = BetheYangStatus::iteration_limit;
    bool converged = false;
};

/// Same-charge pair ss or anti-s anti-s, topological charge +/-2, P=0.
/// Repulsive/free p>=1, I>0 half-odd-integer. Solve
/// M*L*sinh(theta)+chi(2*theta)=2*pi*I. Only asymptotic Bethe--Yang;
/// numerical convergence does not certify accuracy at small M*L.
template <uni20::Real Real>
BetheYangResult<Real> same_charge_pair(Real mass, Real length, Real p, uni20::half_int number,
                                       BetheYangOptions<Real> options = {})
{
  ParticleSpectrum<Real> const spectrum(mass, p);
  if (!uni20::isfinite(length) || length <= Real{0} || p < Real{1} || number.twice() <= 0 || number.is_integral() ||
      !uni20::isfinite(options.tolerance) || options.tolerance <= Real{0})
    throw std::invalid_argument("same-charge pair requires L>0, p>=1, positive half-odd I and tolerance>0");
  BetheYangResult<Real> out;
  Real const twice_number = Real(number.twice());
  if (twice_number + Real{1} == twice_number)
  {
    out.status = BetheYangStatus::precision_limit;
    return out; // Do not round a half-odd label onto an integral one.
  }
  Real const pi = Real{4} * std::atan(Real{1}), u = mass * length, target = pi * twice_number;
  if (!uni20::isfinite(u) || u == Real{0} || !uni20::isfinite(target))
  {
    out.status = BetheYangStatus::precision_limit;
    return out;
  }
  // chi is nonnegative on theta>=0 in the repulsive regime. The free
  // solution brackets the unique positive root from above.
  Real lo = Real{0}, hi = std::asinh(target / u);
  if (!uni20::isfinite(hi) || hi == Real{0})
  {
    out.status = BetheYangStatus::precision_limit;
    return out;
  }
  options.phase.tolerance = options.tolerance / Real{16};
  for (; out.iterations < options.max_iterations;)
  {
    ++out.iterations;
    Real const theta = p == Real{1} ? hi : lo + (hi - lo) / Real{2};
    auto const phase = soliton_phase(Real{2} * theta, p, options.phase);
    out.evaluations += phase.evaluations;
    out.phase_status = phase.status;
    out.quadrature_error = phase.quadrature_error;
    out.tail_bound = phase.tail_bound;
    out.phase_error = phase.tail_bound + phase.quadrature_error;
    if (!phase.converged)
    {
      out.status = BetheYangStatus::phase_limit;
      return out;
    }
    Real const residual = u * std::sinh(theta) + *phase.phase - target;
    out.residual = std::abs(residual);
    if (!uni20::isfinite(residual))
    {
      out.status = BetheYangStatus::precision_limit;
      return out;
    }
    if (out.residual + out.phase_error <= options.tolerance)
    {
      Real const energy = mass * (Real{2} * std::cosh(theta));
      if (!uni20::isfinite(energy))
      {
        out.status = BetheYangStatus::precision_limit;
        return out;
      }
      out.energy = energy;
      out.rapidity = theta;
      out.status = BetheYangStatus::converged;
      out.converged = true;
      return out;
    }
    if (theta == lo || theta == hi)
    {
      out.status = BetheYangStatus::precision_limit;
      return out;
    }
    if (residual > Real{0})
      hi = theta;
    else
      lo = theta;
  }
  return out;
}
} // namespace bethe::sine_gordon
