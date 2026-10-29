// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/spinon.hpp>
#include <optional>
#include <uni20/core/finite_quadrature.hpp>

namespace bethe::xxz
{
enum class BulkStatus
{
  converged,
  evaluation_limit,
  mesh_limit,
  precision_limit
};

template <uni20::Real Real> struct BulkOptions
{
    /// Absolute error target in units J*max(1,|Delta|).
    Real tolerance = Real{256} * uni20::numeric_limits<Real>::epsilon();
    std::size_t max_evaluations = 200000, max_levels = 18;
    uni20::Precision precision() const { return uni20::native_precision<Real>(); }
};
#if UNI20_HAS_MPREAL
template <> struct BulkOptions<uni20::mpreal>
{
    explicit BulkOptions(uni20::Precision p) : tolerance(default_tolerance(p)), precision_(p) {}
    uni20::mpreal tolerance;
    std::size_t max_evaluations = 200000, max_levels = 18;
    uni20::Precision precision() const { return precision_; }

  private:
    static uni20::mpreal default_tolerance(uni20::Precision p)
    {
      if (p.is_exact()) throw std::invalid_argument("XXZ bulk energy requires finite working precision");
      return 256 * uni20::numeric_limits<uni20::mpreal>::epsilon(p);
    }
    uni20::Precision precision_;
};
#endif

template <uni20::Real Real> struct BulkEnergy
{
    /// Absent on budget/precision failure; never publish an unconverged energy.
    std::optional<Real> energy, error;
    std::size_t evaluations = 0;
    bool converged = false;
    BulkStatus status = BulkStatus::precision_limit;
};

/// Zero-field bulk energy/site, including the J*Delta/4 constant. Bethe-root
/// Fourier solution; see the zero-temperature expression in Bortz--Gohmann
/// (2005), with J*Delta/4 restored. The quadrature
/// error is an estimate from successive meshes plus analytic tail bounds,
/// not an interval-arithmetic certificate.
template <uni20::Real Real>
BulkEnergy<Real> bulk_energy_density(Real delta, Real exchange = Real{1}, BulkOptions<Real> const& options = {})
{
  Real const zero = uni20::scalar_at_precision<Real>(0, options.precision()), one = uni20::scalar_like(zero, 1);
  bethe::detail::require_same_precision(delta, zero);
  bethe::detail::require_same_precision(options.tolerance, zero);
#if UNI20_HAS_MPREAL
  if constexpr (uni20::has_runtime_precision_v<Real>)
    if (!exchange.is_exact()) bethe::detail::require_same_precision(exchange, zero);
#endif
  exchange = uni20::scalar_like(zero, exchange);
  bethe::detail::validate_spinon_parameters(zero, exchange);
  if (!uni20::isfinite(delta) || delta <= -Real{1})
    throw std::invalid_argument("XXZ bulk energy requires finite Delta>-1");
  Real const eps = uni20::numeric_limits<Real>::epsilon(options.precision()), pi = bethe::detail::spinon_pi(zero);
  if (!uni20::isfinite(options.tolerance) || options.tolerance < Real{128} * eps ||
      options.tolerance > one / Real{100} || options.max_levels > 30)
    throw std::invalid_argument("bulk tolerance must be in [128 epsilon, 0.01], max-levels<=30");
  BulkEnergy<Real> out;
  Real const scale = std::max(Real{1}, uni20::math::abs(delta)), tol = options.tolerance * scale;
  Real const minimum = [&] {
    if constexpr (uni20::has_runtime_precision_v<Real>)
      return uni20::math::next_up(zero);
    else
      return uni20::numeric_limits<Real>::min();
  }();
  auto integrate = [&](auto const& function, Real const& end) {
    uni20::tanh_sinh_options<Real> controls(options.precision());
    controls.relative_tolerance = tol / Real{8};
    controls.absolute_tolerance[0] = controls.relative_tolerance;
    controls.max_evaluations = options.max_evaluations;
    controls.max_levels = options.max_levels;
    auto result = uni20::tanh_sinh_integrate<Real, 1>(
        [&](uni20::quadrature_point<Real> const& point) { return function(point.lower_distance); }, zero, end,
        controls);
    out.evaluations = result.evaluations;
    if (result.status == uni20::quadrature_status::evaluation_limit)
      out.status = BulkStatus::evaluation_limit;
    else if (result.status == uni20::quadrature_status::refinement_limit)
      out.status = BulkStatus::mesh_limit;
    return result;
  };
  auto finish = [&](Real value, Real error) {
    error += Real{16} * eps * scale;
    if (error <= tol && uni20::isfinite(value * exchange) && uni20::isfinite(error * exchange) &&
        uni20::math::abs(value * exchange) >= minimum)
    {
      out.energy = value * exchange;
      out.error = error * exchange;
      out.converged = true;
      out.status = BulkStatus::converged;
    }
    return out;
  };
  if (delta == Real{1}) return finish(Real{1} / Real{4} - uni20::math::log(one * Real{2}), zero);
  if (delta == Real{0}) return finish(-one / pi, zero);
  if (delta < Real{1})
  {
    Real const gamma = uni20::math::acos(delta), a = pi / gamma, prefactor = uni20::math::sin(gamma) / gamma;
    Real const end = -uni20::math::log(tol / Real{32}) / Real{2};
    // x=gamma*omega in the standard integral. This has a uniform exp(-2x)
    // envelope even as gamma->0; expm1 avoids cancellation at x=0.
    auto integrand = [&](Real x) {
      if (x == zero) return std::array<Real, 1>{(a - one) / a};
      Real const e = uni20::math::exp(-Real{2} * x);
      return std::array<Real, 1>{
          Real{2} * e / (Real{1} + e) *
          (uni20::math::expm1(-Real{2} * (a - Real{1}) * x) / uni20::math::expm1(-Real{2} * a * x))};
    };
    auto const q = integrate(integrand, end);
    if (q.status != uni20::quadrature_status::converged)
    {
      return out;
    }
    return finish(delta / Real{4} - prefactor * q.value[0],
                  prefactor * (q.error[0] + uni20::math::exp(-Real{2} * end)));
  }

  Real const eta = uni20::math::acosh(delta);
  if (eta >= Real{1})
  {
    // e0/J=Delta/4-sinh(eta)*(1/2+2 sum_{n>=1}1/(1+exp(2*n*eta))).
    // The geometric majorant is particularly cheap in the Ising regime.
    Real const r = uni20::math::exp(-Real{2} * eta);
    Real const sh = delta >= Real{2}
                        ? delta * uni20::math::sqrt((Real{1} - Real{1} / delta) * (Real{1} + Real{1} / delta))
                        : uni20::math::sqrt(delta - Real{1}) * uni20::math::sqrt(delta + Real{1});
    if (!uni20::isfinite(sh)) return out;
    bethe::detail::CompensatedSum<Real> sum(zero);
    Real next = r;
    for (;;)
    {
      Real const tail = sh * (Real{2} * next / (Real{1} - r));
      if (tail <= tol / Real{8}) return finish(delta / Real{4} - sh / Real{2} - sh * (Real{2} * sum.value()), tail);
      if (out.evaluations == options.max_evaluations)
      {
        out.status = BulkStatus::evaluation_limit;
        return out;
      }
      ++out.evaluations;
      sum.add(next / (Real{1} + next));
      next *= r;
    }
  }

  // Poisson-resummed density, lambda=eta*x:
  // rho(lambda)=(1/(2*eta))*sum_m sech(pi*(x-m*pi/eta)).
  // Unlike the Fourier sum this remains cheap arbitrarily close to Delta=1.
  Real const period = pi / eta, half = period / Real{2};
  Real const sh = uni20::math::sinh(eta / Real{2}), ch = uni20::math::cosh(eta / Real{2}), bound = Real{2} * ch * ch;
  Real const end = std::min(half, -uni20::math::log(tol / (Real{64} * bound)) / pi);
  Real const geometric = -uni20::math::expm1(-pi * period);
  Real const tail = end == half ? zero : bound * Real{4} / pi * uni20::math::exp(-pi * end) / geometric;
  unsigned images = 0;
  Real image_error;
  for (;; ++images)
  {
    image_error = bound * end * Real{4} * uni20::math::exp(-pi * (Real(images + 1) * period - end)) / geometric;
    if (image_error <= tol / Real{16}) break;
  }
  auto sech = [&](Real x) {
    Real const e = uni20::math::exp(-pi * uni20::math::abs(x));
    return Real{2} * e / (Real{1} + e * e);
  };
  auto integrand = [&](Real x) {
    bethe::detail::CompensatedSum<Real> density(zero);
    density.add(sech(x));
    for (unsigned m = 1; m <= images; ++m)
    {
      density.add(sech(x - Real(m) * period));
      density.add(sech(x + Real(m) * period));
    }
    Real const ratio = uni20::math::sin(eta * x) / sh;
    return std::array<Real, 1>{bound * density.value() / (Real{1} + ratio * ratio)};
  };
  auto const q = integrate(integrand, end);
  if (q.status != uni20::quadrature_status::converged)
  {
    return out;
  }
  return finish(delta / Real{4} - q.value[0], q.error[0] + tail + image_error);
}
} // namespace bethe::xxz
