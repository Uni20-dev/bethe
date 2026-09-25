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
template <uni20::Real Real> struct State
{
    Real mass{}, length{}, scaled_length{}, cutoff{};
    std::optional<Real> scaling_function, casimir_energy, effective_central_charge;
    Real nonlinear_residual = uni20::numeric_limits<Real>::infinity();
    Real nonlinear_error = uni20::numeric_limits<Real>::infinity();
    Real mesh_error = uni20::numeric_limits<Real>::infinity();
    Real cutoff_error = uni20::numeric_limits<Real>::infinity();
    Real direct_tail_bound = uni20::numeric_limits<Real>::infinity();
    std::size_t intervals = 0, iterations = 0, cutoffs = 0, kernel_products = 0;
    bool converged = false;
    Status status = Status::mesh_limit;
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
    bool converged = false;
    Status status = Status::iteration_limit;
};

// Even half-line trapezoid. Difference/sum kernel cache is O(N), not a dense
// matrix. Work counts one folded kernel-times-logarithm term per (i,j).
template <uni20::Real Real>
GridResult<Real> grid(Real r, Real cutoff, std::size_t n, Options<Real> const& options, State<Real>& work)
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
  std::vector<Real> k(2 * n + 1), drive(points), weights(points), u(points), next(points), logs(points);
  CompensatedSum<Real> mass;
  for (std::size_t j = 0; j < k.size(); ++j)
  {
    k[j] = kernel(h * Real(j));
    mass.add((j ? Real{2} : Real{1}) * h * k[j]);
  }
  // This full difference-grid sum bounds every folded row sum from above.
  Real const q = mass.value() * thermal_factors(r, Real{1}).filling;
  if (!uni20::isfinite(q) || q >= Real{1})
  {
    out.status = Status::mesh_limit;
    return out;
  }
  CompensatedSum<Real> response;
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
    response.add(weights[j] * drive[j] * thermal_factors(drive[j], Real{1}).filling / pi);
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
      logs[j] = weights[j] * thermal_factors(drive[j] + u[j], Real{1}).log_weight;
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
    work.nonlinear_residual = out.residual;
    work.nonlinear_error = out.error;
    if (out.error <= options.tolerance / Real{64})
    {
      CompensatedSum<Real> sum;
      for (std::size_t j = 0; j < points; ++j)
        sum.add(drive[j] * logs[j]);
      out.y = -sum.value() / pi;
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
} // namespace detail

/// Periodic massive scaling Lee-Yang source-free ground-state TBA.
/// Returns BULK-SUBTRACTED E_C, not an absolute extensive vacuum energy.
/// c_eff=-6 L E_C/pi approaches 2/5, whereas the theory has c=-22/5.
/// Mesh/cutoff error estimates are refinement diagnostics, not certificates.
template <uni20::Real Real = double> State<Real> ground_state(Real mass, Real length, Options<Real> const& options = {})
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
  State<Real> out;
  out.mass = mass;
  out.length = length;
  out.scaled_length = r;
  Real const pi = Real{4} * std::atan(Real{1});
  Real const target = std::max(Real{2}, -std::log(options.tolerance) + std::log(Real{128}));
  out.cutoff = options.initial_cutoff.value_or(std::max(Real{2}, std::log(Real{2} * target) - std::log(r)));
  std::optional<Real> previous_cutoff_y;
  Real previous_cutoff_error{};
  for (std::size_t cut = 0; cut < options.max_cutoffs; ++cut)
  {
    ++out.cutoffs;
    std::optional<Real> previous_y;
    Real previous_error{};
    detail::GridResult<Real> last;
    unsigned stable = 0;
    bool mesh_ok = false;
    out.mesh_error = uni20::numeric_limits<Real>::infinity();
    for (std::size_t n = options.initial_intervals;;)
    {
      out.intervals = n;
      last = detail::grid(r, out.cutoff, n, options, out);
      if (!last.converged && last.status != Status::mesh_limit)
      {
        out.status = last.status;
        return out;
      }
      if (last.converged)
      {
        if (previous_y)
        {
          out.mesh_error = std::abs(last.y - *previous_y) + last.error + previous_error;
          stable = out.mesh_error <= options.tolerance / Real{8} ? stable + 1 : 0;
          if (stable == 2)
          {
            mesh_ok = true;
            break;
          }
        }
        previous_y = last.y;
        previous_error = last.error;
      }
      if (n > options.max_intervals / 2) break;
      n *= 2;
    }
    if (!mesh_ok)
    {
      out.status = Status::mesh_limit;
      return out;
    }
    Real const drive =
        std::exp(std::log(r) + out.cutoff - std::log(Real{2})) + std::exp(std::log(r) - out.cutoff - std::log(Real{2}));
    out.direct_tail_bound = std::exp(-drive) / std::tanh(out.cutoff) / pi;
    if (previous_cutoff_y)
    {
      out.cutoff_error = std::abs(last.y - *previous_cutoff_y) + last.error + out.mesh_error + previous_cutoff_error +
                         out.direct_tail_bound;
      if (out.cutoff_error <= options.tolerance / Real{2})
      {
        Real const energy = last.y / length, ceff = -Real{6} * last.y / pi;
        if (!uni20::isfinite(energy) || !uni20::isfinite(ceff))
        {
          out.status = Status::precision_limit;
          return out;
        }
        out.scaling_function = last.y;
        out.casimir_energy = energy;
        out.effective_central_charge = ceff;
        out.converged = true;
        out.status = Status::converged;
        return out;
      }
    }
    previous_cutoff_y = last.y;
    previous_cutoff_error = last.error + out.mesh_error;
    if (cut + 1 < options.max_cutoffs)
    {
      Real const next = out.cutoff + Real{1};
      if (!uni20::isfinite(next) || next == out.cutoff)
      {
        out.status = Status::precision_limit;
        return out;
      }
      out.cutoff = next;
    }
  }
  out.status = Status::cutoff_limit;
  return out;
}
} // namespace bethe::lee_yang
