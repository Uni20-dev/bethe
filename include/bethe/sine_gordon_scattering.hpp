// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/sine_gordon.hpp>

namespace bethe::sine_gordon
{
template <uni20::Real Real> struct PhaseResult
{
    std::optional<Real> phase; // odd chi, S_ss(theta)=-exp(i*chi(theta))
    Real fourier_cutoff{}, tail_bound{}, quadrature_error{};
    std::size_t evaluations = 0, cutoffs = 0;
    bool converged = false;
    KernelStatus status = KernelStatus::cutoff_limit;
};

/// Real-rapidity same-charge scattering, all p>0. chi'=2*pi*G_p;
/// the constant fermionic minus sign is NOT included in the odd phase.
/// Feher--Palmai--Takacs, arXiv:1112.6322, Eq. (2.1). Numerical errors describe
/// quadrature only, not omitted finite-volume wrapping corrections.
template <uni20::Real Real> PhaseResult<Real> soliton_phase(Real theta, Real p, KernelOptions<Real> options = {})
{
  Real const zero = detail::validate_kernel_options(options);
  bethe::detail::require_same_precision(p, zero);
  bethe::detail::require_same_precision(theta, zero);
  if (!uni20::isfinite(p) || p <= Real{0} || !uni20::isfinite(theta) || !uni20::isfinite(options.tolerance) ||
      options.tolerance <= Real{0})
    throw std::invalid_argument("sine-Gordon phase requires finite theta, p>0 and tolerance>0");
  PhaseResult<Real> out;
  out.fourier_cutoff = zero;
  out.tail_bound = out.quadrature_error = uni20::scalar_like(zero, std::numeric_limits<double>::infinity());
  if (p == Real{1} || theta == Real{0})
  {
    out.phase = zero;
    out.tail_bound = out.quadrature_error = zero;
    out.converged = true;
    out.status = KernelStatus::converged;
    return out;
  }
  Real const pi = bethe::detail::pi_like(zero), decay = pi * std::min(p, uni20::scalar_like(zero, 1));
  Real cutoff = Real{1} / decay;
  bool tail_ok = false;
  for (; out.cutoffs < options.max_cutoffs;)
  {
    ++out.cutoffs;
    if (!uni20::isfinite(cutoff))
    {
      out.status = KernelStatus::precision_limit;
      return out;
    }
    out.fourier_cutoff = cutoff;
    // |R(k)|<=2*exp(-decay*k)/(1-exp(-p*pi*k)), and 1/k<=1/K.
    out.tail_bound =
        Real{2} * uni20::math::exp(-decay * cutoff) / (decay * cutoff) / (-uni20::math::expm1(-p * pi * cutoff));
    if (out.tail_bound <= options.tolerance / Real{4})
    {
      tail_ok = true;
      break;
    }
    cutoff *= Real{2};
  }
  if (!tail_ok) return out;
  auto const q = detail::fourier_quadrature<Real, 1>(
      [&](Real k) {
        Real const x = k * theta;
        // sin(k*theta)/k, retaining tiny theta without cancellation.
        Real const sinc =
            uni20::math::abs(x) < uni20::math::sqrt(uni20::numeric_limits<Real>::epsilon(options.precision()))
                ? theta
                : uni20::math::sin(x) / k;
        return std::array<Real, 1>{Real{2} * pi * detail::fourier_components(k, zero, p).real() * sinc};
      },
      cutoff, options);
  out.evaluations = q.evaluations;
  bool const converged = q.status == uni20::quadrature_status::converged;
  if (converged) out.quadrature_error = q.error[0];
  if (!converged || out.quadrature_error + out.tail_bound > options.tolerance)
  {
    out.status = KernelStatus::quadrature_limit;
    return out;
  }
  if (!uni20::isfinite(q.value[0]))
  {
    out.status = KernelStatus::precision_limit;
    return out;
  }
  out.phase = q.value[0];
  out.status = KernelStatus::converged;
  out.converged = true;
  return out;
}
} // namespace bethe::sine_gordon
