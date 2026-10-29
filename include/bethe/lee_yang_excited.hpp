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
      requires(!uni20::has_runtime_precision_v<Real>)
    {
      this->tolerance = Real{65536} * uni20::numeric_limits<Real>::epsilon();
      this->max_kernel_products = 1000000000;
    }
    explicit OneParticleOptions(uni20::Precision p)
      requires(uni20::has_runtime_precision_v<Real>)
        : Options<Real>(p)
    {
      this->tolerance = 65536 * uni20::numeric_limits<Real>::epsilon(p);
      this->max_kernel_products = 1000000000;
    }
};

template <uni20::Real Real> struct OneParticleState : Diagnostics<Real>
{
    OneParticleState()
      requires(!uni20::has_runtime_precision_v<Real>)
        : OneParticleState(Real{})
    {}
    explicit OneParticleState(Real const& zero)
        : Diagnostics<Real>(zero), quantization_residual(this->nonlinear_residual),
          source_error(this->nonlinear_residual)
    {}
    std::optional<Real> scaling_function, casimir_energy, beta, pole_displacement;
    Real quantization_residual;
    Real source_error;
    std::size_t root_iterations = 0;
};

namespace detail
{
template <uni20::Real Real> void validate_one_particle(Real mass, Real length, OneParticleOptions<Real> const& options)
{
  Real const r = validate(mass, length, static_cast<Options<Real> const&>(options));
  if (r < Real{5} || r > Real{30})
    throw std::invalid_argument("Lee-Yang regular one-particle solver requires 5<=mL<=30");
}

template <uni20::Real Real>
GridResult<Real> one_particle_grid(Real r, Real cutoff, std::size_t n, OneParticleOptions<Real> const& options,
                                   OneParticleState<Real>& work, Real& final_displacement)
{
  GridResult<Real> failed(uni20::scalar_like(r, 0));
  work.quantization_residual = work.source_error = uni20::scalar_like(r, std::numeric_limits<double>::infinity());
  Real const pi = Real{4} * uni20::math::atan(uni20::scalar_like(r, 1)),
             eps = uni20::numeric_limits<Real>::epsilon(options.precision());
  Real lower = -r - Real{10}, upper = uni20::math::log(pi / Real{12});
  Real lower_f = uni20::scalar_like(r, 0), upper_f = lower_f;
  auto evaluate = [&](Real log_d, Real& f) {
    Real const d = uni20::math::exp(log_d), beta = pi / Real{6} + d;
    auto out = grid(r, cutoff, n, static_cast<Options<Real> const&>(options), work, beta);
    if (!out.converged) return out;
    Real const denominator = Real{2} * uni20::math::cos(pi / Real{3} + d) * uni20::math::sin(d);
    Real const log_s =
        uni20::math::log(denominator + uni20::math::sqrt(uni20::scalar_like(r, 3))) - uni20::math::log(denominator);
    f = r * uni20::math::cos(beta) - log_s + out.source_integral;
    out.y += Real{2} * r * uni20::math::sin(beta);
    out.error += Real{64} * eps * uni20::math::abs(out.y);
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
  Real previous_d = uni20::math::exp(lower);
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
    Real f = uni20::scalar_like(r, 0);
    auto out = evaluate(mid, f);
    if (!out.converged) return out;
    Real const d = uni20::math::exp(mid);
    // Local secant response, combined with a deliberately conservative
    // beta-to-energy factor. Like the mesh errors this is an estimate, not
    // an interval certificate. Require successive source positions to agree.
    // Do not differentiate residual noise at the inner solve's accuracy.
    if (uni20::math::abs(f - previous_f) > Real{4} * out.source_error)
    {
      Real const secant = (f - previous_f) / (mid - previous);
      if (secant > Real{0} && uni20::isfinite(secant)) slope = secant;
    }
    if (!uni20::isfinite(f))
    {
      failed.status = Status::precision_limit;
      return failed;
    }
    Real const root_error =
        Real{30} * r * (uni20::math::abs(d - previous_d) + d * (uni20::math::abs(f) + out.source_error) / slope);
    work.quantization_residual = uni20::math::abs(f);
    work.source_error = root_error;
    // Reserve 1/24 of the total error target for source quantization. This
    // accommodates its r-scaled rounding floor while leaving room for two
    // adjacent grids inside the controller's 1/8 mesh budget.
    if (root_error <= options.tolerance / Real{24} && uni20::math::abs(f) <= options.tolerance)
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
    if (uni20::math::abs(candidate - mid) < Real{4} * eps * std::max(Real{1}, uni20::math::abs(mid)))
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
  detail::validate_one_particle(mass, length, options);
  Real const r = mass * length;
  OneParticleState<Real> out(uni20::scalar_like(r, 0));
  Real d = uni20::scalar_like(r, 0);
  Real const pi = Real{4} * uni20::math::atan(uni20::scalar_like(r, 1));
  Real const tail = uni20::math::exp(-detail::particle_source(uni20::scalar_like(r, 0), pi / Real{4}));
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
