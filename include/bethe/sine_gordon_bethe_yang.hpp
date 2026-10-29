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
    uni20::Precision precision() const { return uni20::native_precision<Real>(); }
};
#if UNI20_HAS_MPREAL
template <> struct BetheYangOptions<uni20::mpreal>
{
    explicit BetheYangOptions(uni20::Precision p) : tolerance(default_tolerance(p)), phase(p), precision_(p) {}
    uni20::mpreal tolerance;
    std::size_t max_iterations = 256;
    KernelOptions<uni20::mpreal> phase;
    uni20::Precision precision() const { return precision_; }

  private:
    static uni20::mpreal default_tolerance(uni20::Precision p)
    {
      if (p.is_exact()) throw std::invalid_argument("Bethe-Yang levels require finite working precision");
      return 16384 * uni20::numeric_limits<uni20::mpreal>::epsilon(p);
    }
    uni20::Precision precision_;
};
#endif
template <uni20::Real Real> struct BetheYangResult
{
    std::optional<Real> rapidity, energy; // roots +/-rapidity; excitation energy, NOT E-L*e_bulk
    Real residual{}, phase_error{}, quadrature_error{}, tail_bound{};
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
  Real const zero = uni20::scalar_at_precision<Real>(0, options.precision());
  for (auto const& value : {mass, length, p, options.tolerance})
    bethe::detail::require_same_precision(value, zero);
  if (options.phase.precision() != options.precision())
    throw std::invalid_argument("Bethe-Yang phase and solve controls require matching precision");
  (void)detail::validate_kernel_options(options.phase);
  ParticleSpectrum<Real> const spectrum(mass, p);
  if (!uni20::isfinite(length) || length <= Real{0} || p < Real{1} || number.twice() <= 0 || number.is_integral() ||
      !uni20::isfinite(options.tolerance) || options.tolerance <= Real{0})
    throw std::invalid_argument("same-charge pair requires L>0, p>=1, positive half-odd I and tolerance>0");
  BetheYangResult<Real> out;
  out.residual = out.phase_error = out.quadrature_error = out.tail_bound =
      uni20::scalar_like(zero, std::numeric_limits<double>::infinity());
  Real const twice_number = uni20::scalar_like(zero, number.twice());
  if (twice_number + Real{1} == twice_number)
  {
    out.status = BetheYangStatus::precision_limit;
    return out; // Do not round a half-odd label onto an integral one.
  }
  Real const pi = bethe::detail::pi_like(zero), u = mass * length, target = pi * twice_number;
  if (!uni20::isfinite(u) || u == Real{0} || !uni20::isfinite(target))
  {
    out.status = BetheYangStatus::precision_limit;
    return out;
  }
  // chi is nonnegative on theta>=0 in the repulsive regime. The free
  // solution brackets the unique positive root from above.
  Real lo = zero, hi = uni20::math::asinh(target / u);
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
    Real const residual = u * uni20::math::sinh(theta) + *phase.phase - target;
    out.residual = uni20::math::abs(residual);
    if (!uni20::isfinite(residual))
    {
      out.status = BetheYangStatus::precision_limit;
      return out;
    }
    if (out.residual + out.phase_error <= options.tolerance)
    {
      Real const energy = mass * (Real{2} * uni20::math::cosh(theta));
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
