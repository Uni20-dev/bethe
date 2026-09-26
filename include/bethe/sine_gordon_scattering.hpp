// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/sine_gordon.hpp>

namespace bethe::sine_gordon
{
template <uni20::Real Real> struct PhaseResult
{
    std::optional<Real> phase; // odd chi, S_ss(theta)=-exp(i*chi(theta))
    Real fourier_cutoff{}, tail_bound = uni20::numeric_limits<Real>::infinity();
    Real quadrature_error = uni20::numeric_limits<Real>::infinity();
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
  if (!uni20::isfinite(p) || p <= Real{0} || !uni20::isfinite(theta) || !uni20::isfinite(options.tolerance) ||
      options.tolerance <= Real{0})
    throw std::invalid_argument("sine-Gordon phase requires finite theta, p>0 and tolerance>0");
  PhaseResult<Real> out;
  if (p == Real{1} || theta == Real{0})
  {
    out.phase = Real{0};
    out.tail_bound = out.quadrature_error = Real{0};
    out.converged = true;
    out.status = KernelStatus::converged;
    return out;
  }
  Real const pi = Real{4} * std::atan(Real{1}), decay = pi * std::min(p, Real{1});
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
    out.tail_bound = Real{2} * std::exp(-decay * cutoff) / (decay * cutoff) / (-std::expm1(-p * pi * cutoff));
    if (out.tail_bound <= options.tolerance / Real{4})
    {
      tail_ok = true;
      break;
    }
    cutoff *= Real{2};
  }
  if (!tail_ok) return out;
  auto const q = bethe::detail::tanh_sinh<Real, 1>(
      [&](Real k) {
        Real const x = k * theta;
        // sin(k*theta)/k, retaining tiny theta without cancellation.
        Real const sinc = std::abs(x) < std::sqrt(uni20::numeric_limits<Real>::epsilon()) ? theta : std::sin(x) / k;
        return std::array<Real, 1>{Real{2} * pi * detail::fourier_components(k, Real{0}, p).real() * sinc};
      },
      Real{0}, cutoff, options.tolerance / Real{4}, out.evaluations, options.max_evaluations, options.max_levels,
      {Real{1}});
  if (q.converged) out.quadrature_error = q.error[0];
  if (!q.converged || out.quadrature_error + out.tail_bound > options.tolerance)
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
