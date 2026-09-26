// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/quadrature.hpp>
#include <bethe/spinon.hpp>
#include <optional>

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
};

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
  bethe::detail::validate_spinon_parameters(Real{0}, exchange);
  if (!uni20::isfinite(delta) || delta <= -Real{1})
    throw std::invalid_argument("XXZ bulk energy requires finite Delta>-1");
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = bethe::detail::pi<Real>();
  if (!uni20::isfinite(options.tolerance) || options.tolerance < Real{128} * eps || options.tolerance > Real{0.01} ||
      options.max_levels > 30)
    throw std::invalid_argument("bulk tolerance must be in [128 epsilon, 0.01], max-levels<=30");
  BulkEnergy<Real> out;
  Real const scale = std::max(Real{1}, std::abs(delta)), tol = options.tolerance * scale;
  auto finish = [&](Real value, Real error) {
    error += Real{16} * eps * scale;
    if (error <= tol && uni20::isfinite(value * exchange) && uni20::isfinite(error * exchange) &&
        std::abs(value * exchange) >= uni20::numeric_limits<Real>::min())
    {
      out.energy = value * exchange;
      out.error = error * exchange;
      out.converged = true;
      out.status = BulkStatus::converged;
    }
    return out;
  };
  if (delta == Real{1}) return finish(Real{1} / Real{4} - std::log(Real{2}), Real{0});
  if (delta == Real{0}) return finish(-Real{1} / pi, Real{0});
  if (delta < Real{1})
  {
    Real const gamma = std::acos(delta), a = pi / gamma, prefactor = std::sin(gamma) / gamma;
    Real const end = -std::log(tol / Real{32}) / Real{2};
    // x=gamma*omega in the standard integral. This has a uniform exp(-2x)
    // envelope even as gamma->0; expm1 avoids cancellation at x=0.
    auto integrand = [&](Real x) {
      Real const e = std::exp(-Real{2} * x);
      return std::array<Real, 1>{Real{2} * e / (Real{1} + e) *
                                 (std::expm1(-Real{2} * (a - Real{1}) * x) / std::expm1(-Real{2} * a * x))};
    };
    auto const q = bethe::detail::tanh_sinh<Real, 1>(integrand, Real{0}, end, tol / Real{8}, out.evaluations,
                                                     options.max_evaluations, options.max_levels, {Real{1}});
    if (!q.converged)
    {
      out.status = out.evaluations == options.max_evaluations ? BulkStatus::evaluation_limit : BulkStatus::mesh_limit;
      return out;
    }
    return finish(delta / Real{4} - prefactor * q.value[0], prefactor * (q.error[0] + std::exp(-Real{2} * end)));
  }

  Real const eta = std::acosh(delta);
  if (eta >= Real{1})
  {
    // e0/J=Delta/4-sinh(eta)*(1/2+2 sum_{n>=1}1/(1+exp(2*n*eta))).
    // The geometric majorant is particularly cheap in the Ising regime.
    Real const r = std::exp(-Real{2} * eta);
    Real const sh = delta >= Real{2} ? delta * std::sqrt((Real{1} - Real{1} / delta) * (Real{1} + Real{1} / delta))
                                     : std::sqrt(delta - Real{1}) * std::sqrt(delta + Real{1});
    if (!uni20::isfinite(sh)) return out;
    bethe::detail::CompensatedSum<Real> sum;
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
  Real const sh = std::sinh(eta / Real{2}), ch = std::cosh(eta / Real{2}), bound = Real{2} * ch * ch;
  Real const end = std::min(half, -std::log(tol / (Real{64} * bound)) / pi);
  Real const geometric = -std::expm1(-pi * period);
  Real const tail = end == half ? Real{0} : bound * Real{4} / pi * std::exp(-pi * end) / geometric;
  unsigned images = 0;
  Real image_error;
  for (;; ++images)
  {
    image_error = bound * end * Real{4} * std::exp(-pi * (Real(images + 1) * period - end)) / geometric;
    if (image_error <= tol / Real{16}) break;
  }
  auto sech = [&](Real x) {
    Real const e = std::exp(-pi * std::abs(x));
    return Real{2} * e / (Real{1} + e * e);
  };
  auto integrand = [&](Real x) {
    bethe::detail::CompensatedSum<Real> density;
    density.add(sech(x));
    for (unsigned m = 1; m <= images; ++m)
    {
      density.add(sech(x - Real(m) * period));
      density.add(sech(x + Real(m) * period));
    }
    Real const ratio = std::sin(eta * x) / sh;
    return std::array<Real, 1>{bound * density.value() / (Real{1} + ratio * ratio)};
  };
  auto const q = bethe::detail::tanh_sinh<Real, 1>(integrand, Real{0}, end, tol / Real{8}, out.evaluations,
                                                   options.max_evaluations, options.max_levels, {Real{1}});
  if (!q.converged)
  {
    out.status = out.evaluations == options.max_evaluations ? BulkStatus::evaluation_limit : BulkStatus::mesh_limit;
    return out;
  }
  return finish(delta / Real{4} - q.value[0], q.error[0] + tail + image_error);
}
} // namespace bethe::xxz
