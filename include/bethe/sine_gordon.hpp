// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/quadrature.hpp>
#include <complex>
#include <optional>

namespace bethe::sine_gordon
{
// p=beta_canonical^2/(8*pi-beta_canonical^2). Repulsive p>1,
// attractive 0<p<1; p=1 is the free massive Dirac point.
enum class KernelStatus
{
  converged,
  quadrature_limit,
  cutoff_limit,
  precision_limit
};
template <uni20::Real Real> struct KernelOptions
{
    Real tolerance = Real{1024} * uni20::numeric_limits<Real>::epsilon(); // absolute component error
    std::size_t max_evaluations = 100000, max_levels = 16, max_cutoffs = 64;
};
template <uni20::Real Real> struct KernelResult
{
    std::optional<std::complex<Real>> value;
    Real fourier_cutoff{}, tail_bound = uni20::numeric_limits<Real>::infinity();
    Real quadrature_error = uni20::numeric_limits<Real>::infinity();
    std::size_t evaluations = 0, cutoffs = 0;
    bool converged = false;
    KernelStatus status = KernelStatus::cutoff_limit;
};
namespace detail
{
// cos(k*z) times the Fourier multiplier, with the decaying exponent
// combined before evaluation. Neither sinh(p*pi*k/2) nor cosh(k*Im(z))
// is formed, so a legal complex contour cannot cause spurious overflow.
template <uni20::Real Real> std::complex<Real> fourier_integrand(Real k, std::complex<Real> z, Real p)
{
  Real const pi = Real{4} * std::atan(Real{1});
  if (k == Real{0}) return {(Real{1} - Real{1} / p) / (Real{2} * pi), Real{0}};
  Real const a = pi * k, decay = pi * std::min(p, Real{1});
  Real const ratio = -std::expm1(-std::abs(p - Real{1}) * a) / (-std::expm1(-p * a) * (Real{1} + std::exp(-a)));
  Real const slow = std::exp(-(decay - std::abs(z.imag())) * k);
  Real const fast_ratio = std::exp(-Real{2} * std::abs(z.imag()) * k);
  Real const even = slow * (Real{1} + fast_ratio) / Real{2};
  // exp(-decay*k)*sinh(k*y), without subtracting nearby exponentials.
  Real const odd = std::copysign(slow * (-std::expm1(-Real{2} * std::abs(z.imag()) * k)) / Real{2}, z.imag());
  Real const sign = p > Real{1} ? Real{1} : Real{-1};
  return sign * ratio / pi * std::complex<Real>(std::cos(k * z.real()) * even, -std::sin(k * z.real()) * odd);
}
} // namespace detail

/// G_p(z) = integral_0^infinity dk cos(k*z)/(2*pi)
///          *sinh((p-1)*pi*k/2)/(sinh(p*pi*k/2)*cosh(pi*k/2)).
/// Principal analytic strip: |Im(z)| < pi*min(1,p), even at p=1.
/// The quadrature estimate is not an interval certificate. No kernel value
/// is published unless both that estimate and the analytic tail bound pass.
template <uni20::Real Real>
KernelResult<Real> scattering_kernel(std::complex<Real> z, Real p, KernelOptions<Real> options = {})
{
  Real const pi = Real{4} * std::atan(Real{1}), decay = pi * std::min(p, Real{1});
  if (!uni20::isfinite(p) || !(p > Real{0}) || !uni20::isfinite(z.real()) || !uni20::isfinite(z.imag()) ||
      !(std::abs(z.imag()) < decay) || !uni20::isfinite(options.tolerance) || !(options.tolerance > Real{0}))
    throw std::invalid_argument("sine-Gordon kernel requires finite p>0, |Im(z)|<pi*min(1,p), and positive tolerance");
  KernelResult<Real> out;
  if (p == Real{1})
  {
    out.value = std::complex<Real>{};
    out.tail_bound = out.quadrature_error = Real{0};
    out.converged = true;
    out.status = KernelStatus::converged;
    return out;
  }
  Real const width = decay - std::abs(z.imag());
  Real cutoff = Real{1} / width;
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
    // |cos(k*z)| <= exp(|Im(z)|*k), and for k>=K the denominator
    // 1-exp(-p*pi*k) >= 1-exp(-p*pi*K). This bounds the entire tail.
    out.tail_bound = std::exp(-width * cutoff) / pi / width / (-std::expm1(-p * pi * cutoff));
    if (out.tail_bound <= options.tolerance / Real{4})
    {
      tail_ok = true;
      break;
    }
    cutoff *= Real{2};
  }
  if (!tail_ok) return out;
  auto const q = bethe::detail::tanh_sinh<Real, 2>(
      [&](Real k) {
        auto const value = detail::fourier_integrand(k, z, p);
        return std::array<Real, 2>{value.real(), value.imag()};
      },
      Real{0}, cutoff, options.tolerance / Real{4}, out.evaluations, options.max_evaluations, options.max_levels,
      {Real{1}, Real{1}});
  if (q.converged) out.quadrature_error = std::max(q.error[0], q.error[1]);
  if (!q.converged || out.quadrature_error + out.tail_bound > options.tolerance)
  {
    out.status = KernelStatus::quadrature_limit;
    return out;
  }
  if (!uni20::isfinite(q.value[0]) || !uni20::isfinite(q.value[1]))
  {
    out.status = KernelStatus::precision_limit;
    return out;
  }
  out.value = std::complex<Real>(q.value[0], q.value[1]);
  out.status = KernelStatus::converged;
  out.converged = true;
  return out;
}
} // namespace bethe::sine_gordon
