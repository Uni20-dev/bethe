// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/thermal.hpp>
#include <bethe/solver.hpp>
#include <optional>
#include <stdexcept>
#include <vector>

namespace bethe::lee_yang
{
enum class Status
{
  converged,
  iteration_limit,
  mesh_limit,
  cutoff_limit,
  work_limit,
  precision_limit
};
template <uni20::Real Real> struct Options
{
    Real tolerance = Real{8192} * uni20::numeric_limits<Real>::epsilon(); // Absolute error target in Y=L E_C.
    std::size_t initial_intervals = 32, max_intervals = 2048, max_iterations = 1000, max_cutoffs = 3;
    std::size_t max_kernel_products = 200000000;
    std::optional<Real> initial_cutoff{};
};
template <uni20::Real Real> struct Diagnostics
{
    Real mass{}, length{}, scaled_length{}, cutoff{};
    Real nonlinear_residual = uni20::numeric_limits<Real>::infinity();
    Real nonlinear_error = uni20::numeric_limits<Real>::infinity();
    Real mesh_error = uni20::numeric_limits<Real>::infinity();
    Real cutoff_error = uni20::numeric_limits<Real>::infinity();
    Real direct_tail_bound = uni20::numeric_limits<Real>::infinity();
    std::size_t intervals = 0, iterations = 0, cutoffs = 0, kernel_products = 0;
    bool converged = false;
    Status status = Status::mesh_limit;
};
template <uni20::Real Real> struct State : Diagnostics<Real>
{
    std::optional<Real> scaling_function, casimir_energy, effective_central_charge;
};

/// Positive K=-phi/(2*pi); finite real x, including zero. Integral K dx=1.
template <uni20::Real Real> Real kernel(Real x)
{
  if (!uni20::isfinite(x)) throw std::invalid_argument("Lee-Yang kernel requires a finite argument");
  Real const t = std::exp(-std::abs(x)), t2 = t * t, pi = Real{4} * std::atan(Real{1});
  return std::sqrt(Real{3}) / pi * t * (Real{1} + t2) / (Real{1} + t2 + t2 * t2);
}

namespace detail
{
template <uni20::Real Real> struct GridResult
{
    Real y{}, error{}, residual{};
    Real source_integral{}, source_error{};
    Real source_coordinate{}; // Zero for vacuum; excited displacement in error-estimate units.
    bool converged = false;
    Status status = Status::iteration_limit;
};

// Real-axis source and real part of K(i*beta-x), without complex arithmetic
// or overflowing hyperbolic functions. beta=0 denotes the vacuum problem.
template <uni20::Real Real> Real particle_source(Real x, Real beta)
{
  Real const t = std::exp(-x), a = std::sqrt(Real{3}) * t / Real{2};
  Real const u = (Real{1} - t * t) * std::cos(beta) / Real{2};
  Real const v = (Real{1} + t * t) * std::sin(beta) / Real{2};
  return std::log1p(-Real{4} * v * a / (u * u + (v + a) * (v + a)));
}
template <uni20::Real Real> Real continued_kernel(Real x, Real beta)
{
  Real const t = std::exp(-x), s = (Real{1} - t * t) / Real{2}, c = (Real{1} + t * t) / Real{2};
  Real const cb = std::cos(beta), sb = std::sin(beta);
  Real const re = s * s * cb * cb - c * c * sb * sb + Real{3} * t * t / Real{4};
  Real const im = -Real{2} * s * c * cb * sb;
  Real const pi = Real{4} * std::atan(Real{1});
  return std::sqrt(Real{3}) / (Real{2} * pi) * t * (c * cb * re - s * sb * im) / (re * re + im * im);
}

// Even half-line trapezoid. Difference/sum kernel cache is O(N), not a dense
// matrix. Work counts one folded kernel-times-logarithm term per (i,j).
template <uni20::Real Real>
GridResult<Real> grid(Real r, Real cutoff, std::size_t n, Options<Real> const& options, Diagnostics<Real>& work,
                      Real beta = Real{0})
{
  using bethe::detail::CompensatedSum;
  using bethe::detail::thermal_factors;
  GridResult<Real> out;
  work.nonlinear_residual = uni20::numeric_limits<Real>::infinity();
  work.nonlinear_error = uni20::numeric_limits<Real>::infinity();
  if (!uni20::isfinite(Real{2} * cutoff))
  {
    out.status = Status::precision_limit;
    return out;
  }
  auto const points = n + 1, cost = points * points;
  if (cost > options.max_kernel_products - work.kernel_products)
  {
    out.status = Status::work_limit;
    return out;
  }
  Real const pi = Real{4} * std::atan(Real{1}), h = cutoff / Real(n), eps = uni20::numeric_limits<Real>::epsilon();
  Real const logr = std::log(r), log2 = std::log(Real{2});
  std::vector<Real> k(2 * n + 1), drive(points), weights(points), u(points), next(points), logs(points), source(points),
      continued(points);
  CompensatedSum<Real> mass;
  for (std::size_t j = 0; j < k.size(); ++j)
  {
    k[j] = kernel(h * Real(j));
    mass.add((j ? Real{2} : Real{1}) * h * k[j]);
  }
  // This full difference-grid sum bounds every folded row sum from above.
  Real const q = mass.value() * thermal_factors(r + particle_source(Real{0}, beta), Real{1}).filling;
  if (!uni20::isfinite(q) || q >= Real{1})
  {
    out.status = Status::mesh_limit;
    return out;
  }
  CompensatedSum<Real> response;
  CompensatedSum<Real> source_response;
  for (std::size_t j = 0; j < points; ++j)
  {
    Real const x = h * Real(j);
    drive[j] = j == 0 ? r : std::exp(logr + x - log2) + std::exp(logr - x - log2);
    if (!uni20::isfinite(drive[j]))
    {
      out.status = Status::precision_limit;
      return out;
    }
    weights[j] = j == 0 || j == n ? h / Real{2} : h;
    source[j] = particle_source(x, beta);
    continued[j] = beta == Real{0} ? Real{0} : Real{2} * continued_kernel(x, beta);
    Real const filling = thermal_factors(drive[j] + source[j], Real{1}).filling;
    response.add(weights[j] * drive[j] * filling / pi);
    source_response.add(weights[j] * std::abs(continued[j]) * filling);
  }
  for (std::size_t iteration = 0;; ++iteration)
  {
    if (cost > options.max_kernel_products - work.kernel_products)
    {
      out.status = Status::work_limit;
      return out;
    }
    work.kernel_products += cost;
    for (std::size_t j = 0; j < points; ++j)
      logs[j] = weights[j] * thermal_factors(drive[j] + source[j] + u[j], Real{1}).log_weight;
    out.residual = Real{0};
    Real magnitude = Real{1};
    bool changed = false;
    for (std::size_t i = 0; i < points; ++i)
    {
      CompensatedSum<Real> sum;
      for (std::size_t j = 0; j < points; ++j)
        sum.add((k[i > j ? i - j : j - i] + k[i + j]) * logs[j]);
      next[i] = sum.value();
      if (!uni20::isfinite(next[i]))
      {
        out.status = Status::precision_limit;
        return out;
      }
      out.residual = std::max(out.residual, std::abs(next[i] - u[i]));
      magnitude = std::max(magnitude, std::abs(next[i]));
      changed = changed || next[i] != u[i];
    }
    Real const floor = Real{64} * eps * magnitude * response.value() / (Real{1} - q);
    out.error = response.value() * out.residual / (Real{1} - q) + floor;
    out.source_error = source_response.value() * (out.residual + Real{64} * eps * magnitude) / (Real{1} - q);
    work.nonlinear_residual = out.residual;
    work.nonlinear_error = out.error;
    if (out.error <= options.tolerance / Real{64} && out.source_error <= options.tolerance / Real{64})
    {
      CompensatedSum<Real> sum;
      for (std::size_t j = 0; j < points; ++j)
        sum.add(drive[j] * logs[j]);
      out.y = -sum.value() / pi;
      CompensatedSum<Real> integral;
      for (std::size_t j = 0; j < points; ++j)
        integral.add(continued[j] * logs[j]);
      out.source_integral = integral.value();
      out.error += Real{64} * eps * std::abs(out.y);
      work.nonlinear_error = out.error;
      out.converged = true;
      out.status = Status::converged;
      return out;
    }
    if (floor > options.tolerance / Real{64} || !changed)
    {
      out.status = Status::precision_limit;
      return out;
    }
    if (iteration == options.max_iterations) return out;
    ++work.iterations;
    u.swap(next);
  }
}
template <uni20::Real Real> Real validate(Real mass, Real length, Options<Real> const& options)
{
  if (!uni20::isfinite(mass) || !uni20::isfinite(length) || mass <= Real{0} || length <= Real{0})
    throw std::invalid_argument("Lee-Yang mass and length must be finite and positive");
  Real const r = mass * length;
  if (!uni20::isfinite(r) || r <= Real{0})
    throw std::invalid_argument("Lee-Yang mL must be representable and positive");
  if (!uni20::isfinite(options.tolerance) || options.tolerance <= Real{0} || options.tolerance >= Real{1})
    throw std::invalid_argument("Lee-Yang absolute Y tolerance must lie in (0,1)");
  if (options.initial_intervals < 2 || options.max_intervals < options.initial_intervals ||
      options.max_intervals > 8192)
    throw std::invalid_argument("Lee-Yang requires 2<=initial_intervals<=max_intervals<=8192");
  if (options.initial_cutoff && (!uni20::isfinite(*options.initial_cutoff) || *options.initial_cutoff <= Real{0}))
    throw std::invalid_argument("Lee-Yang cutoff must be finite and positive");
  return r;
}

// Shared mesh/cutoff controller; callbacks provide the state-specific grid
// equations. The tail multiplier bounds exp(-source) on the real axis.
template <uni20::Real Real, typename Grid>
std::optional<Real> refine(Real mass, Real length, Options<Real> const& options, Diagnostics<Real>& out,
                           Real tail_multiplier, Grid&& solve_grid)
{
  Real const r = mass * length;
  out.mass = mass;
  out.length = length;
  out.scaled_length = r;
  Real const pi = Real{4} * std::atan(Real{1});
  Real const target = std::max(Real{2}, -std::log(options.tolerance) + std::log(Real{128}));
  out.cutoff = options.initial_cutoff.value_or(std::max(Real{2}, std::log(Real{2} * target) - std::log(r)));
  std::optional<Real> previous_cutoff_y;
  Real previous_cutoff_source{};
  Real previous_cutoff_error{};
  for (std::size_t cut = 0; cut < options.max_cutoffs; ++cut)
  {
    ++out.cutoffs;
    std::optional<Real> previous_y;
    Real previous_error{};
    Real previous_source{};
    detail::GridResult<Real> last;
    unsigned stable = 0;
    bool mesh_ok = false;
    out.mesh_error = uni20::numeric_limits<Real>::infinity();
    for (std::size_t n = options.initial_intervals;;)
    {
      out.intervals = n;
      last = solve_grid(r, out.cutoff, n, options, out);
      if (!last.converged && last.status != Status::mesh_limit)
      {
        out.status = last.status;
        return {};
      }
      if (last.converged)
      {
        if (previous_y)
        {
          out.mesh_error = std::abs(last.y - *previous_y) + last.error + previous_error;
          out.mesh_error += std::abs(last.source_coordinate - previous_source);
          stable = out.mesh_error <= options.tolerance / Real{8} ? stable + 1 : 0;
          if (stable == 2)
          {
            mesh_ok = true;
            break;
          }
        }
        previous_y = last.y;
        previous_error = last.error;
        previous_source = last.source_coordinate;
      }
      if (n > options.max_intervals / 2) break;
      n *= 2;
    }
    if (!mesh_ok)
    {
      out.status = Status::mesh_limit;
      return {};
    }
    Real const drive =
        std::exp(std::log(r) + out.cutoff - std::log(Real{2})) + std::exp(std::log(r) - out.cutoff - std::log(Real{2}));
    out.direct_tail_bound = tail_multiplier * std::exp(-drive) / std::tanh(out.cutoff) / pi;
    if (previous_cutoff_y)
    {
      out.cutoff_error = std::abs(last.y - *previous_cutoff_y) + last.error + out.mesh_error + previous_cutoff_error +
                         out.direct_tail_bound;
      out.cutoff_error += std::abs(last.source_coordinate - previous_cutoff_source);
      if (out.cutoff_error <= options.tolerance / Real{2})
      {
        if (!uni20::isfinite(last.y / length))
        {
          out.status = Status::precision_limit;
          return {};
        }
        out.converged = true;
        out.status = Status::converged;
        return last.y;
      }
    }
    previous_cutoff_y = last.y;
    previous_cutoff_error = last.error + out.mesh_error;
    previous_cutoff_source = last.source_coordinate;
    if (cut + 1 < options.max_cutoffs)
    {
      Real const next = out.cutoff + Real{1};
      if (!uni20::isfinite(next) || next == out.cutoff)
      {
        out.status = Status::precision_limit;
        return {};
      }
      out.cutoff = next;
    }
  }
  out.status = Status::cutoff_limit;
  return {};
}
} // namespace detail

/// Periodic massive scaling Lee-Yang source-free ground-state TBA.
/// Returns BULK-SUBTRACTED E_C, not an absolute extensive vacuum energy.
/// c_eff=-6 L E_C/pi approaches 2/5, whereas the theory has c=-22/5.
/// Mesh/cutoff error estimates are refinement diagnostics, not certificates.
template <uni20::Real Real = double> State<Real> ground_state(Real mass, Real length, Options<Real> const& options = {})
{
  detail::validate(mass, length, options);
  State<Real> out;
  auto const y = detail::refine(mass, length, options, out, Real{1},
                                [](Real r, Real cutoff, std::size_t n, auto const& opt, auto& work) {
                                  return detail::grid(r, cutoff, n, opt, work);
                                });
  if (y)
  {
    Real const pi = Real{4} * std::atan(Real{1}), ceff = -Real{6} * *y / pi;
    if (!uni20::isfinite(ceff))
    {
      out.converged = false;
      out.status = Status::precision_limit;
      return out;
    }
    out.scaling_function = *y;
    out.casimir_energy = *y / length;
    out.effective_central_charge = ceff;
  }
  return out;
}
} // namespace bethe::lee_yang
