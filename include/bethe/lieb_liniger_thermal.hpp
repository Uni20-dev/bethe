// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <array>
#include <bethe/detail/continuum_newton.hpp>
#include <bethe/detail/gauss_legendre.hpp>
#include <bethe/detail/newton.hpp>
#include <bethe/detail/thermal.hpp>
#include <optional>

namespace bethe::lieb_liniger::thermal
{
enum class Status
{
  converged,
  iteration_limit,
  mesh_limit,
  cutoff_limit,
  density_limit,
  precision_limit
};
template <uni20::Real Real> struct Options
{
    Real tolerance = Real{16384} * uni20::numeric_limits<Real>::epsilon();
    std::size_t initial_nodes = 16, max_nodes = 256, max_iterations = 512, max_cutoffs = 4;
    std::optional<Real> initial_cutoff{};
};
template <uni20::Real Real> struct State
{
    Real interaction{}, temperature{}, chemical_potential{};
    std::optional<Real> pressure, density, energy_per_length, entropy_per_length;
    Real cutoff{}, nonlinear_residual{}, mesh_error{}, cutoff_error{};
    std::size_t nodes = 0, iterations = 0, cutoffs = 0;
    bool converged = false;
    Status status = Status::cutoff_limit;
};
namespace detail
{
template <uni20::Real Real> struct Mesh
{
    std::vector<Real> k, w, epsilon, rho_total;
    std::array<Real, 4> values{}; // pressure, density, energy/length, entropy/length
    Real residual{};
    bool converged = false;
    Status status = Status::mesh_limit;
};

// Even pseudoenergy: integrate the positive half and use C(k-q)+C(k+q).
// All variables are dimensionless in units sqrt(max(T,abs(mu))).
template <uni20::Real Real>
Mesh<Real> solve_mesh(Real c, Real t, Real mu, Real cutoff, bethe::detail::GaussLegendreRule<Real> const& rule,
                      Options<Real> const& options, std::size_t& iterations)
{
  using bethe::detail::CompensatedSum;
  Real const pi = Real{4} * std::atan(Real{1});
  std::size_t const n = rule.x.size();
  Mesh<Real> out;
  out.k.resize(n);
  out.w.resize(n);
  out.epsilon.resize(n);
  std::vector<Real> bare(n), kernel(n * n);
  for (std::size_t j = 0; j < n; ++j)
  {
    out.k[j] = cutoff * (rule.x[j] + Real{1}) / Real{2};
    out.w[j] = cutoff * rule.w[j] / Real{2};
    bare[j] = out.k[j] * out.k[j] - mu;
    out.epsilon[j] = bare[j];
  }
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      kernel[i * n + j] = out.w[j] / pi *
                          (bethe::detail::rational_scattering_kernel(out.k[i] - out.k[j], c) +
                           bethe::detail::rational_scattering_kernel(out.k[i] + out.k[j], c));
  std::vector<Real> residual(n), jacobian(n * n);
  auto evaluate = [&](std::vector<Real> const& epsilon, bool jac) {
    std::vector<bethe::detail::ThermalFactors<Real>> factors(n);
    for (std::size_t j = 0; j < n; ++j)
    {
      if (!uni20::isfinite(epsilon[j])) return uni20::numeric_limits<Real>::infinity();
      factors[j] = bethe::detail::thermal_factors(epsilon[j], t);
    }
    Real norm{};
    for (std::size_t i = 0; i < n; ++i)
    {
      CompensatedSum<Real> sum;
      sum.add(epsilon[i] - bare[i]);
      for (std::size_t j = 0; j < n; ++j)
      {
        sum.add(kernel[i * n + j] * factors[j].log_weight);
        if (jac) jacobian[i * n + j] = Real(i == j) - kernel[i * n + j] * factors[j].filling;
      }
      residual[i] = sum.value();
      if (!uni20::isfinite(residual[i])) return uni20::numeric_limits<Real>::infinity();
      norm = std::max(norm, std::abs(residual[i]) / std::max({Real{1}, std::abs(epsilon[i]), std::abs(bare[i])}));
    }
    return norm;
  };
  for (;;)
  {
    out.residual = evaluate(out.epsilon, true);
    if (out.residual <= options.tolerance * t / Real{64}) break;
    if (iterations == options.max_iterations)
    {
      out.status = Status::iteration_limit;
      return out;
    }
    std::vector<Real> step = residual;
    for (auto& x : step)
      x = -x;
    if (!bethe::detail::newton_step(jacobian, step)) return out;
    std::vector<Real> trial(n);
    Real damping{1};
    bool accepted = false;
    for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits; ++backtrack)
    {
      for (std::size_t j = 0; j < n; ++j)
        trial[j] = out.epsilon[j] + damping * step[j];
      if (evaluate(trial, false) < out.residual)
      {
        accepted = true;
        out.epsilon = std::move(trial);
        break;
      }
      damping /= Real{2};
    }
    if (!accepted)
    {
      out.status = Status::precision_limit;
      return out;
    }
    ++iterations;
  }
  // -d epsilon/d mu=2*pi*rho_total, using the same dressed linear operator.
  out.rho_total.assign(n, Real{1} / (Real{2} * pi));
  if (!bethe::detail::newton_step(jacobian, out.rho_total)) return out;
  std::array<CompensatedSum<Real>, 4> values;
  for (std::size_t j = 0; j < n; ++j)
  {
    if (!(out.rho_total[j] > Real{0})) return out;
    auto const f = bethe::detail::thermal_factors(out.epsilon[j], t);
    Real const occupied = Real{2} * out.w[j] * f.filling * out.rho_total[j];
    values[0].add(out.w[j] * f.log_weight / pi);
    values[1].add(occupied);
    values[2].add(occupied * out.k[j] * out.k[j]);
    values[3].add(Real{2} * out.w[j] * f.entropy * out.rho_total[j]);
  }
  for (std::size_t j = 0; j < 4; ++j)
  {
    out.values[j] = values[j].value();
    Real const floor = (uni20::numeric_limits<Real>::min() * uni20::numeric_limits<Real>::epsilon()) * Real{8};
    if (!(out.values[j] > Real{0}) || !uni20::isfinite(out.values[j]) || !(out.values[j] * options.tolerance >= floor))
    {
      out.status = Status::precision_limit;
      return out;
    }
  }
  out.converged = true;
  out.status = Status::converged;
  return out;
}
template <uni20::Real Real> Real difference(Mesh<Real> const& a, Mesh<Real> const& b)
{
  Real error{};
  for (std::size_t j = 0; j < 4; ++j)
    error = std::max(error, std::abs(a.values[j] - b.values[j]) / a.values[j]);
  return Real{8} * error;
}
} // namespace detail

namespace detail
{
template <uni20::Real Real>
void validate(Real interaction, Real temperature, Real chemical_potential, Options<Real> const& options)
{
  if (!uni20::isfinite(interaction) || !(interaction > Real{0}) || !uni20::isfinite(temperature) ||
      !(temperature > Real{0}) || !uni20::isfinite(chemical_potential))
    throw std::invalid_argument("Yang-Yang thermodynamics require finite c>0, T>0 and chemical potential");
  if (!uni20::isfinite(options.tolerance) || !(options.tolerance > Real{0}) || options.tolerance >= Real{1} ||
      options.initial_nodes < 4 || options.max_nodes > 512 || options.initial_nodes > options.max_nodes ||
      options.max_cutoffs > 32 ||
      (options.initial_cutoff && (!uni20::isfinite(*options.initial_cutoff) || !(*options.initial_cutoff > Real{0}))))
    throw std::invalid_argument("require 0<tol<1, 4<=initial_nodes<=max_nodes<=512, max_cutoffs<=32, cutoff>0");
}
} // namespace detail

/// Grand-canonical Yang-Yang equilibrium, c,T>0, H=-sum d^2+2c sum delta.
/// No observables are published until both mesh and cutoff refinement pass.
template <uni20::Real Real>
State<Real> equilibrium(Real interaction, Real temperature, Real chemical_potential, Options<Real> options = {})
{
  detail::validate(interaction, temperature, chemical_potential, options);
  State<Real> out;
  out.interaction = interaction;
  out.temperature = temperature;
  out.chemical_potential = chemical_potential;
  Real const scale = std::max(temperature, std::abs(chemical_potential)), momentum_scale = std::sqrt(scale);
  Real const c = interaction / momentum_scale, t = temperature / scale, mu = chemical_potential / scale;
  Real cutoff =
      options.initial_cutoff
          ? *options.initial_cutoff / momentum_scale
          : std::max(Real{2}, std::sqrt(std::max(mu, Real{0}) + t * (Real{16} - std::log(options.tolerance))));
  if (!(c > Real{0}) || !uni20::isfinite(c) || !uni20::isfinite(Real{1} / c) ||
      !(options.tolerance * t / Real{64} > Real{0}) || !(cutoff > Real{0}) || !uni20::isfinite(cutoff))
  {
    out.status = Status::precision_limit;
    return out;
  }
  std::optional<detail::Mesh<Real>> previous_cutoff;
  for (std::size_t domain = 0; domain < options.max_cutoffs; ++domain)
  {
    out.cutoffs = domain + 1;
    out.cutoff = cutoff * momentum_scale;
    if (!uni20::isfinite(out.cutoff))
    {
      out.status = Status::precision_limit;
      return out;
    }
    std::optional<detail::Mesh<Real>> previous_mesh, current;
    for (std::size_t nodes = options.initial_nodes;; nodes *= 2)
    {
      out.nodes = nodes;
      auto candidate =
          detail::solve_mesh(c, t, mu, cutoff, bethe::detail::gauss_legendre<Real>(nodes), options, out.iterations);
      out.nonlinear_residual = candidate.residual;
      if (candidate.status == Status::iteration_limit)
      {
        out.status = candidate.status;
        return out;
      }
      if (candidate.converged && previous_mesh)
      {
        out.mesh_error = detail::difference(candidate, *previous_mesh);
        if (out.mesh_error <= options.tolerance)
        {
          current = std::move(candidate);
          break;
        }
      }
      if (candidate.converged)
        previous_mesh = std::move(candidate);
      else
        previous_mesh.reset();
      if (nodes > options.max_nodes / 2)
      {
        out.status = candidate.status == Status::precision_limit ? Status::precision_limit : Status::mesh_limit;
        return out;
      }
    }
    if (previous_cutoff)
    {
      out.cutoff_error = detail::difference(*current, *previous_cutoff);
      if (out.cutoff_error <= options.tolerance)
      {
        auto v = current->values;
        v[0] = v[0] * momentum_scale * scale;
        v[1] *= momentum_scale;
        v[2] = v[2] * momentum_scale * scale;
        v[3] *= momentum_scale;
        Real const floor = (uni20::numeric_limits<Real>::min() * uni20::numeric_limits<Real>::epsilon()) * Real{8};
        for (auto value : v)
          if (!(value > Real{0}) || !uni20::isfinite(value) || !(value * options.tolerance >= floor))
          {
            out.status = Status::precision_limit;
            return out;
          }
        out.pressure = v[0];
        out.density = v[1];
        out.energy_per_length = v[2];
        out.entropy_per_length = v[3];
        out.converged = true;
        out.status = Status::converged;
        return out;
      }
    }
    previous_cutoff = std::move(current);
    cutoff *= Real{5} / Real{4};
    if (!uni20::isfinite(cutoff))
    {
      out.status = Status::precision_limit;
      return out;
    }
  }
  return out;
}

template <uni20::Real Real> struct DensityOptions
{
    Options<Real> equilibrium{};
    Real tolerance = Real{65536} * uni20::numeric_limits<Real>::epsilon();
    std::size_t max_evaluations = 128;
};
template <uni20::Real Real> struct DensityState
{
    Real requested_density{}, density_error{};
    std::size_t evaluations = 0, iterations = 0;
    std::optional<State<Real>> state;
    bool converged = false;
    Status status = Status::density_limit;
};

/// Canonical equilibrium: invert the monotone n(mu) at fixed c,T.
/// Only successful inner solves participate in the bracket; no failed inner
/// observable is used to infer a sign or silently treated as vacuum.
template <uni20::Real Real>
DensityState<Real> at_density(Real interaction, Real temperature, Real density, DensityOptions<Real> options = {})
{
  detail::validate(interaction, temperature, Real{0}, options.equilibrium);
  if (!(density > Real{0}) || !uni20::isfinite(density) || !(options.tolerance > Real{0}) ||
      !uni20::isfinite(options.tolerance) || options.tolerance >= Real{1})
    throw std::invalid_argument("require finite density>0 and 0<density tolerance<1");
  DensityState<Real> out;
  out.requested_density = density;
  if (options.max_evaluations == 0) return out;
  options.equilibrium.tolerance = std::min(options.equilibrium.tolerance, options.tolerance / Real{8});
  if (!(options.equilibrium.tolerance > Real{0}))
  {
    out.status = Status::precision_limit;
    return out;
  }
  Real const pi = Real{4} * std::atan(Real{1});
  // Classical gas seed, evaluated in log space rather than forming fugacity.
  Real mu = temperature * (std::log(density) + std::log(Real{2} * std::sqrt(pi)) - std::log(temperature) / Real{2});
  Real step = temperature;
  std::optional<Real> lo, hi;
  Real flo{}, fhi{};
  int previous_side = 0;
  for (; out.evaluations < options.max_evaluations;)
  {
    if (!uni20::isfinite(mu))
    {
      out.status = Status::precision_limit;
      return out;
    }
    auto candidate = equilibrium(interaction, temperature, mu, options.equilibrium);
    ++out.evaluations;
    out.iterations += candidate.iterations;
    if (!candidate.converged)
    {
      out.status = candidate.status;
      return out;
    }
    Real const f = (*candidate.density - density) / density;
    out.density_error = std::abs(f);
    if (out.density_error <= options.tolerance / Real{2})
    {
      out.state = std::move(candidate);
      out.converged = true;
      out.status = Status::converged;
      return out;
    }
    if (!uni20::isfinite(f))
    {
      out.status = Status::precision_limit;
      return out;
    }
    // Illinois regula falsi: halve the retained endpoint's interpolation
    // weight on repeated updates to the same side, preventing stagnation.
    int const side = f < Real{0} ? -1 : 1;
    if (side < 0)
    {
      lo = mu;
      flo = f;
      if (previous_side == side) fhi /= Real{2};
    }
    else
    {
      hi = mu;
      fhi = f;
      if (previous_side == side) flo /= Real{2};
    }
    previous_side = side;
    Real next;
    if (lo && hi)
    {
      Real const weight = (-flo) / std::max(-flo, fhi);
      Real const fraction = weight / (weight + fhi / std::max(-flo, fhi));
      next = (Real{1} - fraction) * *lo + fraction * *hi;
      if (!(next > *lo && next < *hi)) next = *lo / Real{2} + *hi / Real{2};
      if (!(next > *lo && next < *hi))
      {
        out.status = Status::precision_limit;
        return out;
      }
    }
    else
    {
      next = mu + (side < 0 ? step : -step);
      step *= Real{2};
    }
    if (next == mu)
    {
      out.status = Status::precision_limit;
      return out;
    }
    mu = next;
  }
  return out;
}
} // namespace bethe::lieb_liniger::thermal
