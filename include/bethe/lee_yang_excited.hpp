// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/lee_yang.hpp>

namespace bethe::lee_yang
{
template <uni20::Real Real> struct OneParticleOptions : Options<Real>
{
    std::size_t max_root_iterations = 256;
    OneParticleOptions()
    {
      this->tolerance = Real{65536} * uni20::numeric_limits<Real>::epsilon();
      this->max_kernel_products = 1000000000;
    }
};

template <uni20::Real Real> struct OneParticleState : Diagnostics<Real>
{
    std::optional<Real> scaling_function, casimir_energy, beta, pole_displacement;
    Real quantization_residual = uni20::numeric_limits<Real>::infinity();
    Real source_error = uni20::numeric_limits<Real>::infinity();
    std::size_t root_iterations = 0;
};

namespace detail
{
template <uni20::Real Real>
GridResult<Real> one_particle_grid(Real r, Real cutoff, std::size_t n, OneParticleOptions<Real> const& options,
                                   OneParticleState<Real>& work, Real& final_displacement)
{
  GridResult<Real> failed;
  work.quantization_residual = work.source_error = uni20::numeric_limits<Real>::infinity();
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  Real lower = -r - Real{10}, upper = std::log(pi / Real{12});
  Real lower_f{}, upper_f{};
  auto evaluate = [&](Real log_d, Real& f) {
    Real const d = std::exp(log_d), beta = pi / Real{6} + d;
    auto out = grid(r, cutoff, n, static_cast<Options<Real> const&>(options), work, beta);
    if (!out.converged) return out;
    Real const denominator = Real{2} * std::cos(pi / Real{3} + d) * std::sin(d);
    Real const log_s = std::log(denominator + std::sqrt(Real{3})) - std::log(denominator);
    f = r * std::cos(beta) - log_s + out.source_integral;
    out.y += Real{2} * r * std::sin(beta);
    out.error += Real{64} * eps * std::abs(out.y);
    out.source_error += Real{64} * eps * r;
    return out;
  };
  auto a = evaluate(lower, lower_f);
  if (!a.converged) return a;
  auto b = evaluate(upper, upper_f);
  if (!b.converged) return b;
  if (!(lower_f < Real{0} && upper_f > Real{0}))
  {
    failed.status = Status::precision_limit;
    return failed;
  }
  Real previous = lower, previous_f = lower_f, candidate = (lower + upper) / Real{2};
  Real previous_d = std::exp(lower);
  Real slope = Real{1};
  for (std::size_t iteration = 0; iteration < options.max_root_iterations; ++iteration)
  {
    ++work.root_iterations;
    Real const mid = candidate;
    if (mid == lower || mid == upper)
    {
      failed.status = Status::precision_limit;
      return failed;
    }
    Real f{};
    auto out = evaluate(mid, f);
    if (!out.converged) return out;
    Real const d = std::exp(mid);
    // Local secant response, combined with a deliberately conservative
    // beta-to-energy factor. Like the mesh errors this is an estimate, not
    // an interval certificate. Require successive source positions to agree.
    // Do not differentiate residual noise at the inner solve's accuracy.
    if (std::abs(f - previous_f) > Real{4} * out.source_error)
    {
      Real const secant = (f - previous_f) / (mid - previous);
      if (secant > Real{0} && uni20::isfinite(secant)) slope = secant;
    }
    if (!uni20::isfinite(f))
    {
      failed.status = Status::precision_limit;
      return failed;
    }
    Real const root_error = Real{30} * r * (std::abs(d - previous_d) + d * (std::abs(f) + out.source_error) / slope);
    work.quantization_residual = std::abs(f);
    work.source_error = root_error;
    // Reserve 1/24 of the total error target for source quantization. This
    // accommodates its r-scaled rounding floor while leaving room for two
    // adjacent grids inside the controller's 1/8 mesh budget.
    if (root_error <= options.tolerance / Real{24} && std::abs(f) <= options.tolerance)
    {
      out.error += root_error;
      out.source_coordinate = Real{30} * r * d;
      work.nonlinear_error = out.error;
      final_displacement = d;
      return out;
    }
    if (f > Real{0})
    {
      upper = mid;
      upper_f = f;
    }
    else
    {
      lower = mid;
      lower_f = f;
    }
    previous = mid;
    previous_f = f;
    previous_d = d;
    candidate = mid - f / slope;
    if (std::abs(candidate - mid) < Real{4} * eps * std::max(Real{1}, std::abs(mid)))
    {
      Real const step = options.tolerance / (Real{4096} * r * d);
      candidate = mid + (upper - mid > mid - lower ? step : -step);
    }
    if (!(candidate > lower && candidate < upper)) candidate = (lower + upper) / Real{2};
  }
  return failed; // iteration_limit; no physical outputs
}
} // namespace detail

/// Periodic spin-zero one-particle branch for 5<=mL<=30, not UV continuation.
/// E1_C is bulk-subtracted, NOT a gap. Subtract a separately converged E0_C.
template <uni20::Real Real = double>
OneParticleState<Real> one_particle(Real mass, Real length, OneParticleOptions<Real> const& options = {})
{
  Real const r = detail::validate(mass, length, static_cast<Options<Real> const&>(options));
  if (r < Real{5} || r > Real{30})
    throw std::invalid_argument("Lee-Yang regular one-particle solver requires 5<=mL<=30");
  OneParticleState<Real> out;
  Real d{};
  Real const pi = Real{4} * std::atan(Real{1});
  Real const tail = std::exp(-detail::particle_source(Real{0}, pi / Real{4}));
  auto const y = detail::refine(mass, length, static_cast<Options<Real> const&>(options), out, tail,
                                [&](Real scaled, Real cutoff, std::size_t n, auto const&, auto&) {
                                  return detail::one_particle_grid(scaled, cutoff, n, options, out, d);
                                });
  if (y)
  {
    out.scaling_function = *y;
    out.casimir_energy = *y / length;
    out.pole_displacement = d;
    out.beta = pi / Real{6} + d;
  }
  return out;
}
} // namespace bethe::lee_yang
