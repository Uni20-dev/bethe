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
  momentum_limit,
  precision_limit
};
template <uni20::Real Real> struct Options
{
    Real tolerance = Real{4096} * uni20::numeric_limits<Real>::epsilon();
    std::size_t initial_nodes = 16, max_nodes = 256, max_iterations = 128;
    std::size_t max_momentum_iterations = 160;
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
/// Repulsive zero-temperature Lieb equation, H=-sum d^2+2c sum delta.
/// Only mesh-verified observables are populated. Weak coupling can require
/// more nodes than the configured budget; no asymptotic formula is substituted.
template <uni20::Real Real>
Background<Real> background(Real interaction, Real density, Options<Real> options,
                            std::optional<Mesh<Real>>* fine = nullptr, std::optional<Mesh<Real>>* coarse = nullptr)
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
        if (fine) *fine = std::move(candidate);
        if (coarse) *coarse = std::move(previous);
        return out;
      }
    }
    previous = std::move(candidate);
    if (nodes > options.max_nodes / 2) return out;
  }
}
} // namespace detail

template <uni20::Real Real> Background<Real> ground_state(Real interaction, Real density, Options<Real> options = {})
{
  return detail::background(interaction, density, options);
}

enum class Branch
{
  type_i,
  type_ii
};
template <uni20::Real Real> struct Point
{
    Branch branch = Branch::type_i;
    Real momentum{}, energy_error{}, momentum_error{};
    std::optional<Real> energy, rapidity, edge_distance;
    std::size_t iterations = 0;
    bool converged = false;
    Status status = Status::mesh_limit;
};

/// Fixed-particle-number elementary branches: particle/hole relative to a
/// Fermi-edge particle/hole. Energies are excitation gaps, not absolute energies.
template <uni20::Real Real> class Solver {
  public:
    Solver(Real interaction, Real density, Options<Real> options = {}) : options_(options)
    {
      info_ = detail::background(interaction, density, options, &fine_, &coarse_);
    }
    Background<Real> const& background() const { return info_; }
    std::pair<Real, Real> momentum_range(Branch branch) const
    {
      if (branch == Branch::type_i) return {Real{0}, uni20::numeric_limits<Real>::infinity()};
      if (branch == Branch::type_ii) return {Real{0}, Real{2} * pi() * info_.density};
      throw std::invalid_argument("invalid Lieb-Liniger excitation branch");
    }
    Point<Real> at_momentum(Branch branch, Real momentum) const
    {
      auto const [lo, hi] = momentum_range(branch);
      if (!uni20::isfinite(momentum) || momentum < lo || momentum > hi)
        throw std::invalid_argument("momentum outside Lieb-Liniger branch range");
      Point<Real> out;
      out.branch = branch;
      out.momentum = momentum;
      if (!info_.converged)
      {
        out.status = info_.status;
        return out;
      }
      bool const hole = branch == Branch::type_ii;
      bool const reflected = hole && momentum > hi / Real{2};
      Real const target = (reflected ? hi - momentum : momentum) / info_.density;
      if (momentum == Real{0} || (hole && momentum == hi))
      {
        out.energy = Real{0};
        out.rapidity = reflected ? -*info_.fermi_rapidity : *info_.fermi_rapidity;
        out.edge_distance = Real{0};
        out.converged = true;
        out.status = Status::converged;
        return out;
      }
      if (!(target > Real{0}) || !uni20::isfinite(target) || !(options_.tolerance * target / Real{16} > Real{0}))
      {
        out.status = Status::precision_limit;
        return out;
      }
      auto const f = invert(*fine_, hole, target, out.iterations);
      if (f.status != Status::converged)
      {
        out.status = f.status;
        return out;
      }
      auto const c = invert(*coarse_, hole, target, out.iterations);
      if (c.status != Status::converged)
      {
        out.status = c.status;
        return out;
      }
      Real const rounding = Real{64} * uni20::numeric_limits<Real>::epsilon();
      Real const error = Real{8} * std::abs(f.energy - c.energy) + rounding * f.energy +
                         Real{2} * (f.error * std::abs(f.velocity) + c.error * std::abs(c.velocity));
      Real const subnormal_floor =
          (uni20::numeric_limits<Real>::min() * uni20::numeric_limits<Real>::epsilon()) * Real{8};
      out.energy_error = std::max(error * info_.density * info_.density, subnormal_floor);
      out.momentum_error = std::max(f.error, c.error) * info_.density;
      Real const energy = f.energy * info_.density * info_.density;
      Real const distance = f.distance * info_.density;
      Real const rapidity =
          (reflected ? -Real{1} : Real{1}) * (fine_->q + (hole ? -f.distance : f.distance)) * info_.density;
      if (!(energy > Real{0}) || !(options_.tolerance * energy >= subnormal_floor) || !(distance > Real{0}) ||
          !uni20::isfinite(energy) || !uni20::isfinite(distance) || !uni20::isfinite(rapidity) ||
          !uni20::isfinite(out.energy_error))
      {
        out.status = Status::precision_limit;
        return out;
      }
      if (error > options_.tolerance * f.energy)
      {
        out.status = Status::mesh_limit;
        return out;
      }
      out.energy = energy;
      out.rapidity = rapidity;
      out.edge_distance = distance;
      out.converged = true;
      out.status = Status::converged;
      return out;
    }

  private:
    static Real pi() { return Real{4} * std::atan(Real{1}); }
    struct Evaluation
    {
        Real momentum{}, energy{}, dp{}, de{};
    };
    Evaluation evaluate(detail::Mesh<Real> const& mesh, bool hole, Real distance) const
    {
      using bethe::detail::CompensatedSum;
      Real const gamma = info_.interaction / info_.density, sign = hole ? -Real{1} : Real{1};
      Real const shift = sign * distance;
      CompensatedSum<Real> momentum, energy, dp, de;
      momentum.add(distance);
      energy.add(distance * (Real{2} * mesh.q + shift));
      dp.add(Real{1});
      de.add(Real{2} * (mesh.q + shift));
      for (std::size_t j = 0; j < mesh.k.size(); ++j)
      {
        Real const a = mesh.q - mesh.k[j], b = a + shift;
        Real const scale = std::max({gamma, std::abs(a), std::abs(b)});
        Real const gs = gamma / scale, as = a / scale, bs = b / scale;
        Real const phase = std::atan2((distance / scale) * gs, gs * gs + as * bs);
        momentum.add(Real{2} * mesh.w[j] * mesh.rho[j] * phase);
        Real const cb = bethe::detail::rational_scattering_kernel(b, gamma) / pi();
        Real difference;
        Real const sa = std::max(gamma, std::abs(a));
        if (distance < sa / Real{4})
        {
          // Difference of kernels factored by the displacement, even when
          // Q +/- distance rounds back to Q. No subtraction of dressed energies.
          Real const ga = gamma / sa, aa = a / sa;
          difference = -(shift / sa) * ((Real{2} * a + shift) / sa) / (ga * ga + aa * aa) * cb;
        }
        else
          difference = cb - bethe::detail::rational_scattering_kernel(a, gamma) / pi();
        energy.add(sign * mesh.w[j] * mesh.epsilon[j] * difference);
        dp.add(Real{2} * pi() * mesh.w[j] * mesh.rho[j] * cb);
        Real const sb = std::max(gamma, std::abs(b)), gb = gamma / sb, bb = b / sb;
        Real const derivative = -Real{2} * (bb / sb) / (gb * gb + bb * bb) * cb;
        de.add(mesh.w[j] * mesh.epsilon[j] * derivative);
      }
      return {momentum.value(), energy.value(), dp.value(), de.value()};
    }
    struct Inversion
    {
        Real energy{}, distance{}, error{}, velocity{};
        Status status = Status::precision_limit;
    };
    Inversion invert(detail::Mesh<Real> const& mesh, bool hole, Real target, std::size_t& iterations) const
    {
      Inversion out;
      Real low{0}, high = hole ? mesh.q : target;
      auto const edge = evaluate(mesh, hole, Real{0});
      Real x = std::min(high, target / edge.dp);
      for (;;)
      {
        auto const value = evaluate(mesh, hole, x);
        out.energy = value.energy;
        out.distance = x;
        out.error = std::abs(value.momentum - target);
        out.velocity = value.de / value.dp;
        if (!uni20::isfinite(value.energy) || !uni20::isfinite(out.error) || !uni20::isfinite(out.velocity) ||
            !(value.energy > Real{0}) || !(value.dp > Real{0}))
          return out;
        if (out.error <= options_.tolerance * target / Real{16})
        {
          out.status = Status::converged;
          return out;
        }
        if (iterations == options_.max_momentum_iterations)
        {
          out.status = Status::momentum_limit;
          return out;
        }
        ++iterations;
        if (value.momentum > target)
          high = x;
        else
          low = x;
        Real next = x - (value.momentum - target) / value.dp;
        if (!(next > low && next < high)) next = low + (high - low) / Real{2};
        if (next == x || next == low || next == high) return out;
        x = next;
      }
    }
    Options<Real> options_;
    Background<Real> info_;
    std::optional<detail::Mesh<Real>> fine_, coarse_;
};
} // namespace bethe::lieb_liniger::thermo
