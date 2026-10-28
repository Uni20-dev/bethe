// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/quadrature.hpp>
#include <optional>
#include <uni20/core/finite_quadrature.hpp>

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
    uni20::Precision precision() const { return uni20::native_precision<Real>(); }
};
#if UNI20_HAS_MPREAL
template <> struct KernelOptions<uni20::mpreal>
{
    explicit KernelOptions(uni20::Precision p) : tolerance(default_tolerance(p)), precision_(p) {}
    uni20::mpreal tolerance;
    std::size_t max_evaluations = 100000, max_levels = 16, max_cutoffs = 64;
    uni20::Precision precision() const { return precision_; }

  private:
    static uni20::mpreal default_tolerance(uni20::Precision p)
    {
      if (p.is_exact()) throw std::invalid_argument("sine-Gordon kernels require finite working precision");
      return 1024 * uni20::numeric_limits<uni20::mpreal>::epsilon(p);
    }
    uni20::Precision precision_;
};
#endif
template <uni20::Real Real> struct KernelResult
{
    std::optional<uni20::complex<Real>> value;
    Real fourier_cutoff{}, tail_bound{}, quadrature_error{};
    std::size_t evaluations = 0, cutoffs = 0;
    bool converged = false;
    KernelStatus status = KernelStatus::cutoff_limit;
};
namespace detail
{
template <uni20::Real Real> Real validate_kernel_options(KernelOptions<Real> const& options)
{
  Real const zero = uni20::scalar_at_precision<Real>(0, options.precision());
  bethe::detail::require_same_precision(options.tolerance, zero);
  if (!uni20::isfinite(options.tolerance) || options.tolerance <= zero)
    throw std::invalid_argument("sine-Gordon kernel tolerance must be finite and positive");
  if (options.max_levels > 24) throw std::invalid_argument("quadrature levels must not exceed 24");
  return zero;
}

// Both Fourier integrals allocate one quarter of their absolute error budget
// to the finite integral; the remaining error includes the analytic tail.
template <uni20::Real Real, std::size_t Components, typename Function>
auto fourier_quadrature(Function const& f, Real cutoff, KernelOptions<Real> const& options)
{
  uni20::tanh_sinh_options<Real, Components> controls(options.precision());
  controls.relative_tolerance = options.tolerance / Real{4};
  controls.absolute_tolerance.fill(controls.relative_tolerance);
  controls.max_evaluations = options.max_evaluations;
  controls.max_levels = options.max_levels;
  Real const zero = uni20::scalar_like(cutoff, 0);
  return uni20::tanh_sinh_integrate<Real, Components>(
      [&](uni20::quadrature_point<Real> const& point) { return f(point.value); }, zero, cutoff, controls);
}
// cos(k*z) times the Fourier multiplier, with the decaying exponent
// combined before evaluation. Neither sinh(p*pi*k/2) nor cosh(k*Im(z))
// is formed, so a legal complex contour cannot cause spurious overflow.
template <uni20::Real Real> uni20::complex<Real> fourier_components(Real k, Real y, Real p)
{
  Real const zero = uni20::scalar_like(p, 0), pi = bethe::detail::pi_like(zero);
  if (k == Real{0}) return {(Real{1} - Real{1} / p) / (Real{2} * pi), zero};
  Real const a = pi * k, decay = pi * std::min(p, Real{1});
  Real const ratio = -uni20::math::expm1(-uni20::math::abs(p - Real{1}) * a) /
                     (-uni20::math::expm1(-p * a) * (Real{1} + uni20::math::exp(-a)));
  Real const slow = uni20::math::exp(-(decay - uni20::math::abs(y)) * k);
  Real const fast_ratio = uni20::math::exp(-Real{2} * uni20::math::abs(y) * k);
  Real const even = slow * (Real{1} + fast_ratio) / Real{2};
  // exp(-decay*k)*sinh(k*y), without subtracting nearby exponentials.
  Real const odd = uni20::math::copysign(slow * (-uni20::math::expm1(-Real{2} * uni20::math::abs(y) * k)) / Real{2}, y);
  Real const sign = uni20::scalar_like(zero, p > Real{1} ? 1 : -1);
  return sign * ratio / pi * uni20::complex<Real>(even, odd);
}
template <uni20::Real Real> uni20::complex<Real> fourier_integrand(Real k, uni20::complex<Real> z, Real p)
{
  auto const factors = fourier_components(k, z.imag(), p);
  return {uni20::math::cos(k * z.real()) * factors.real(), -uni20::math::sin(k * z.real()) * factors.imag()};
}
} // namespace detail

/// G_p(z) = integral_0^infinity dk cos(k*z)/(2*pi)
///          *sinh((p-1)*pi*k/2)/(sinh(p*pi*k/2)*cosh(pi*k/2)).
/// Principal analytic strip: |Im(z)| < pi*min(1,p), even at p=1.
/// The quadrature estimate is not an interval certificate. No kernel value
/// is published unless both that estimate and the analytic tail bound pass.
template <uni20::Real Real>
KernelResult<Real> scattering_kernel(uni20::complex<Real> z, Real p, KernelOptions<Real> options = {})
{
  Real const zero = detail::validate_kernel_options(options);
  for (auto const& value : {p, z.real(), z.imag()})
    bethe::detail::require_same_precision(value, zero);
  Real const pi = bethe::detail::pi_like(zero), decay = pi * std::min(p, uni20::scalar_like(zero, 1));
  if (!uni20::isfinite(p) || !(p > Real{0}) || !uni20::isfinite(z.real()) || !uni20::isfinite(z.imag()) ||
      !(uni20::math::abs(z.imag()) < decay) || !uni20::isfinite(options.tolerance) || !(options.tolerance > Real{0}))
    throw std::invalid_argument("sine-Gordon kernel requires finite p>0, |Im(z)|<pi*min(1,p), and positive tolerance");
  KernelResult<Real> out;
  out.fourier_cutoff = zero;
  out.tail_bound = out.quadrature_error = uni20::scalar_like(zero, std::numeric_limits<double>::infinity());
  if (p == Real{1})
  {
    out.value = uni20::complex<Real>{zero, zero};
    out.tail_bound = out.quadrature_error = zero;
    out.converged = true;
    out.status = KernelStatus::converged;
    return out;
  }
  Real const width = decay - uni20::math::abs(z.imag());
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
    out.tail_bound = uni20::math::exp(-width * cutoff) / pi / width / (-uni20::math::expm1(-p * pi * cutoff));
    if (out.tail_bound <= options.tolerance / Real{4})
    {
      tail_ok = true;
      break;
    }
    cutoff *= Real{2};
  }
  if (!tail_ok) return out;
  auto const q = detail::fourier_quadrature<Real, 2>(
      [&](Real k) {
        auto const value = detail::fourier_integrand(k, z, p);
        return std::array<Real, 2>{value.real(), value.imag()};
      },
      cutoff, options);
  out.evaluations = q.evaluations;
  bool const converged = q.status == uni20::quadrature_status::converged;
  if (converged) out.quadrature_error = std::max(q.error[0], q.error[1]);
  if (!converged || out.quadrature_error + out.tail_bound > options.tolerance)
  {
    out.status = KernelStatus::quadrature_limit;
    return out;
  }
  if (!uni20::isfinite(q.value[0]) || !uni20::isfinite(q.value[1]))
  {
    out.status = KernelStatus::precision_limit;
    return out;
  }
  out.value = uni20::complex<Real>(q.value[0], q.value[1]);
  out.status = KernelStatus::converged;
  out.converged = true;
  return out;
}
} // namespace bethe::sine_gordon
