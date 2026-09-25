// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/continuum_newton.hpp>
#include <bethe/detail/gauss_legendre.hpp>
#include <bethe/detail/newton.hpp>
#include <optional>

namespace bethe::lieb_liniger::thermo
{
enum class Status
{
  converged,
  mesh_limit,
  density_limit,
  precision_limit
};
template <uni20::Real Real> struct Options
{
    Real tolerance = Real{4096} * uni20::numeric_limits<Real>::epsilon();
    std::size_t initial_nodes = 16, max_nodes = 256, max_iterations = 128;
};
template <uni20::Real Real> struct Background
{
    Real interaction{}, density{};
    std::optional<Real> fermi_rapidity, energy_per_length, chemical_potential;
    // Relative mesh/density errors; mesh_error is an estimate, not a bound.
    Real mesh_error{}, density_error{};
    std::size_t nodes = 0, iterations = 0;
    bool converged = false;
    Status status = Status::mesh_limit;
};
namespace detail
{
// Solve at unit density scale: k/n, gamma=c/n, e/n^3, mu/n^2.
// The trial mesh's density is not yet constrained to one.
template <uni20::Real Real> struct Mesh
{
    Real q{}, number{}, energy{}, mu{}, derivative{};
    std::vector<Real> k, w, rho, epsilon;
};
template <uni20::Real Real>
std::optional<Mesh<Real>> mesh(Real gamma, Real q, bethe::detail::GaussLegendreRule<Real> const& rule)
{
  using bethe::detail::CompensatedSum;
  Real const pi = Real{4} * std::atan(Real{1});
  std::size_t const n = rule.x.size();
  Mesh<Real> out;
  out.q = q;
  out.k.resize(n);
  out.w.resize(n);
  out.rho.assign(n, Real{1} / (Real{2} * pi));
  std::vector<Real> a(n * n), h(n);
  for (std::size_t j = 0; j < n; ++j)
  {
    out.k[j] = q * rule.x[j];
    out.w[j] = q * rule.w[j];
    h[j] = out.k[j] * out.k[j];
  }
  auto kernel = [=](Real x) { return bethe::detail::rational_scattering_kernel(x, gamma) / pi; };
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      a[i * n + j] = Real(i == j) - out.w[j] * kernel(out.k[i] - out.k[j]);
  if (!bethe::detail::newton_step(a, out.rho) || !bethe::detail::newton_step(a, h)) return std::nullopt;
  CompensatedSum<Real> number, energy, rq, hq;
  rq.add(Real{1} / (Real{2} * pi));
  hq.add(q * q);
  for (std::size_t j = 0; j < n; ++j)
  {
    if (!(out.rho[j] > Real{0}) || !uni20::isfinite(out.rho[j])) return std::nullopt;
    number.add(out.w[j] * out.rho[j]);
    energy.add(out.w[j] * out.rho[j] * out.k[j] * out.k[j]);
    rq.add(out.w[j] * kernel(q - out.k[j]) * out.rho[j]);
    hq.add(out.w[j] * kernel(q - out.k[j]) * h[j]);
  }
  // Dressed unit charge is 2*pi*rho. epsilon(Q)=0 fixes mu.
  out.number = number.value();
  out.energy = energy.value();
  out.mu = hq.value() / (Real{2} * pi * rq.value());
  out.derivative = Real{4} * pi * rq.value() * rq.value();
  out.epsilon.resize(n);
  for (std::size_t j = 0; j < n; ++j)
    out.epsilon[j] = h[j] - out.mu * Real{2} * pi * out.rho[j];
  if (!uni20::isfinite(out.number) || !uni20::isfinite(out.energy) || !uni20::isfinite(out.mu) ||
      !uni20::isfinite(out.derivative) || !(out.derivative > Real{0}))
    return std::nullopt;
  return out;
}
} // namespace detail

/// Repulsive zero-temperature Lieb equation, H=-sum d^2+2c sum delta.
/// Only mesh-verified observables are populated. Weak coupling can require
/// more nodes than the configured budget; no asymptotic formula is substituted.
template <uni20::Real Real> Background<Real> ground_state(Real interaction, Real density, Options<Real> options = {})
{
  if (!uni20::isfinite(interaction) || !(interaction > Real{0}) || !uni20::isfinite(density) || !(density > Real{0}))
    throw std::invalid_argument("Lieb-Liniger thermodynamics require finite c>0 and density>0");
  if (!uni20::isfinite(options.tolerance) || !(options.tolerance > Real{0}) || options.tolerance >= Real{1} ||
      options.initial_nodes < 4 || options.max_nodes > 512 || options.initial_nodes > options.max_nodes)
    throw std::invalid_argument("require 0<tol<1 and 4<=initial_nodes<=max_nodes<=512");
  Background<Real> out;
  out.interaction = interaction;
  out.density = density;
  Real const gamma = interaction / density, pi = Real{4} * std::atan(Real{1});
  if (!(gamma > Real{0}) || !uni20::isfinite(gamma) || !uni20::isfinite(Real{1} / gamma))
  {
    out.status = Status::precision_limit;
    return out;
  }
  std::optional<detail::Mesh<Real>> previous;
  Real q = std::min(pi / Real{2}, std::sqrt(gamma));
  for (std::size_t nodes = options.initial_nodes;; nodes *= 2)
  {
    out.nodes = nodes;
    auto const rule = bethe::detail::gauss_legendre<Real>(nodes);
    Real lo{0}, hi = pi;
    std::optional<detail::Mesh<Real>> candidate;
    for (;;)
    {
      if (out.iterations == options.max_iterations)
      {
        out.status = Status::density_limit;
        return out;
      }
      ++out.iterations;
      candidate = detail::mesh(gamma, q, rule);
      // A coarse mesh may not resolve the narrow kernel, even losing positivity.
      // Refine instead of treating that invalid discretization as a density root.
      if (!candidate) break;
      Real const residual = candidate->number - Real{1};
      out.density_error = std::abs(residual);
      if (out.density_error <= options.tolerance / Real{16}) break;
      if (residual > Real{0})
        hi = q;
      else
        lo = q;
      Real next = q - residual / candidate->derivative;
      if (!(next > lo && next < hi)) next = lo + (hi - lo) / Real{2};
      if (next == q)
      {
        out.status = Status::precision_limit;
        return out;
      }
      q = next;
    }
    // At least resolve the kernel width before comparing meshes. Otherwise
    // small absolute energies at tiny gamma can give spurious agreement.
    if (candidate && q / gamma / Real(nodes) > Real{1} / (Real{2} * pi)) candidate.reset();
    if (candidate && previous)
    {
      out.mesh_error = Real{8} * std::max({std::abs(candidate->q - previous->q) / candidate->q,
                                           std::abs(candidate->energy - previous->energy) / candidate->energy,
                                           std::abs(candidate->mu - previous->mu) / candidate->mu});
      if (out.mesh_error <= options.tolerance)
      {
        Real const fermi = q * density, energy = candidate->energy * density * density * density,
                   mu = candidate->mu * density * density;
        if (!uni20::isfinite(fermi) || !uni20::isfinite(energy) || !uni20::isfinite(mu) || !(fermi > Real{0}) ||
            !(energy > Real{0}) || !(mu > Real{0}))
        {
          out.status = Status::precision_limit;
          return out;
        }
        out.fermi_rapidity = fermi;
        out.energy_per_length = energy;
        out.chemical_potential = mu;
        out.converged = true;
        out.status = Status::converged;
        return out;
      }
    }
    previous = std::move(candidate);
    if (nodes > options.max_nodes / 2) return out;
  }
}
} // namespace bethe::lieb_liniger::thermo
