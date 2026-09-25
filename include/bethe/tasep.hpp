// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/newton.hpp>
#include <bethe/solver.hpp>
#include <complex>
#include <optional>
#include <span>

namespace bethe::tasep
{
enum class Status
{
  converged,
  stationary_only,
  seed_limit,
  iteration_limit,
  stalled,
  precision_limit
};
template <uni20::Real Real> struct Options
{
    Real tolerance = Real{128} * uni20::numeric_limits<Real>::epsilon();
    std::size_t max_iterations = 100, max_seed_iterations = 1000, max_sites = 256;
};
template <uni20::Real Real> struct GapState
{
    std::size_t sites = 0, particles = 0, effective_particles = 0, iterations = 0, seed_iterations = 0;
    Real rate{}, residual_norm = uni20::numeric_limits<Real>::infinity();
    std::vector<std::complex<Real>> roots; // Z, in the particle-hole-reduced representation
    std::complex<Real> log_y{};
    std::optional<std::complex<Real>> eigenvalue;
    std::optional<Real> gap;
    bool converged = false;
    Status status = Status::iteration_limit;
};
namespace detail
{
template <uni20::Real Real> bool finite(std::complex<Real> z)
{
  return uni20::isfinite(z.real()) && uni20::isfinite(z.imag());
}
template <uni20::Real Real> Real wrap(Real x)
{
  Real const pi = Real{4} * std::atan(Real{1});
  return std::remainder(x, Real{2} * pi);
}

// Factored Cassini polynomial: coefficients in the monomial basis become
// badly conditioned near half filling. Ehrlich-Aberth here evaluates P/Y
// logarithmically instead, retaining all L roots before selecting a branch.
template <uni20::Real Real>
bool seed(std::size_t l, std::size_t n, std::complex<Real> b, std::vector<std::complex<Real>>& z,
          std::size_t& iterations, Options<Real> const& options)
{
  using C = std::complex<Real>;
  Real const pi = Real{4} * std::atan(Real{1}), tolerance = Real{128} * uni20::numeric_limits<Real>::epsilon();
  z.resize(l);
  for (std::size_t j = 0; j < l; ++j)
    z[j] = Real{2} * std::exp(C(0, Real{2} * pi * (Real(j) + Real{1} / Real{3}) / Real(l)));
  while (iterations < options.max_seed_iterations)
  {
    Real change{};
    for (std::size_t j = 0; j < l; ++j)
    {
      C const p = Real(n) * std::log(Real{1} - z[j]) + Real(l - n) * std::log(Real{1} + z[j]);
      C const difference = b - p;
      C const scale = std::exp(difference.real() > Real{0} ? -difference : difference);
      C const a = difference.real() > Real{0} ? scale - Real{1} : Real{1} - scale;
      C derivative = Real(n) / (z[j] - Real{1}) + Real(l - n) / (z[j] + Real{1});
      // Multiply numerator and denominator by P/Y when |Y/P|>1;
      // exp(log(Y/P)) itself can overflow far from a perfectly finite root.
      if (difference.real() > Real{0}) derivative *= scale;
      bethe::detail::CompensatedSum<C> repulsion;
      for (std::size_t k = 0; k < l; ++k)
        if (k != j)
        {
          if (z[j] == z[k]) return false;
          repulsion.add(Real{1} / (z[j] - z[k]));
        }
      C const step = a / (derivative - a * repulsion.value());
      if (!finite(step)) return false;
      z[j] -= step;
      change = std::max(change, std::abs(step));
    }
    ++iterations;
    if (change <= tolerance)
    {
      // Reject duplicates or incomplete root sets, even if each root has a
      // small residual. At this negative-Y seed all roots are simple.
      for (std::size_t j = 0; j < l; ++j)
      {
        C const f = Real(n) * std::log(Real{1} - z[j]) + Real(l - n) * std::log(Real{1} + z[j]) - b;
        if (!finite(f) || std::abs(C(f.real(), wrap(f.imag()))) > tolerance * Real(l)) return false;
        for (std::size_t k = 0; k < j; ++k)
          if (std::abs(z[j] - z[k]) < std::sqrt(tolerance)) return false;
      }
      return true;
    }
  }
  return false;
}

template <uni20::Real Real> struct System
{
    std::size_t l, n;
    using C = std::complex<Real>;
    struct Evaluation
    {
        std::vector<Real> residual, jacobian;
        Real norm{};
    };
    Evaluation evaluate(std::span<C const> z, C b, bool jacobian) const
    {
      std::size_t const order = 2 * (n + 1);
      Evaluation out{std::vector<Real>(order), jacobian ? std::vector<Real>(order * order) : std::vector<Real>{}};
      auto put = [&](std::size_t row, std::size_t col, C v) {
        if (!jacobian) return;
        out.jacobian[(2 * row) * order + 2 * col] = v.real();
        out.jacobian[(2 * row) * order + 2 * col + 1] = -v.imag();
        out.jacobian[(2 * row + 1) * order + 2 * col] = v.imag();
        out.jacobian[(2 * row + 1) * order + 2 * col + 1] = v.real();
      };
      auto residual = [&](std::size_t row, C f) {
        f.imag(wrap(f.imag()));
        if (!finite(f))
        {
          out.norm = uni20::numeric_limits<Real>::infinity();
          return;
        }
        out.residual[2 * row] = f.real();
        out.residual[2 * row + 1] = f.imag();
        out.norm = std::max(out.norm, std::abs(f) / Real(l));
      };
      bethe::detail::CompensatedSum<C> sum;
      sum.add(C(Real(l) * std::log(Real{2}), Real{4} * std::atan(Real{1})) - b);
      for (std::size_t j = 0; j < n; ++j)
      {
        residual(j, Real(n) * std::log(Real{1} - z[j]) + Real(l - n) * std::log(Real{1} + z[j]) - b);
        put(j, j, Real(n) / (z[j] - Real{1}) + Real(l - n) / (z[j] + Real{1}));
        put(j, n, C(-1));
        sum.add(std::log((z[j] - Real{1}) / (z[j] + Real{1})));
        put(n, j, Real{2} / (z[j] * z[j] - Real{1}));
      }
      residual(n, sum.value());
      put(n, n, C(-1));
      return out;
    }
};
} // namespace detail

/// Unit-direction periodic exclusion: each particle jumps right at rate r
/// when its neighbor is empty. dP/dt=M P; Re(lambda)<=0, gap=-Re(lambda_1).
/// Returns one member (Im>=0) of the leading conjugate pair, not a full spectrum.
template <uni20::Real Real>
GapState<Real> relaxation_gap(std::size_t l, std::size_t particles, Real rate = Real{1}, Options<Real> options = {})
{
  using C = std::complex<Real>;
  if (l < 2 || particles > l || !uni20::isfinite(rate) || !(rate > Real{0}) || !uni20::isfinite(options.tolerance) ||
      !(options.tolerance > Real{0}))
    throw std::invalid_argument("TASEP requires L>=2, 0<=N<=L, finite positive rate and tolerance");
  if (l > options.max_sites) throw std::length_error("TASEP L exceeds max_sites work budget");
  auto const n = std::min(particles, l - particles);
  if (n + 1 > std::numeric_limits<std::size_t>::max() / 2) throw std::length_error("TASEP matrix too large");
  auto const order = 2 * (n + 1);
  if (order > static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real) / order)
    throw std::length_error("TASEP matrix too large");
  GapState<Real> out;
  out.sites = l;
  out.particles = particles;
  out.effective_particles = n;
  out.rate = rate;
  if (n == 0)
  {
    out.converged = true;
    out.status = Status::stationary_only;
    out.residual_norm = Real{0};
    return out;
  }
  Real const pi = Real{4} * std::atan(Real{1});
  if (n == 1)
  {
    Real const angle = Real{2} * pi / Real(l), sine = std::sin(pi / Real(l));
    C const value(-Real{2} * sine * sine, l == 2 ? Real{0} : std::sin(angle));
    out.roots = {Real{1} + Real{2} * value};
    out.log_y = std::log(Real{1} - out.roots[0]) + Real(l - 1) * std::log(Real{1} + out.roots[0]);
    out.eigenvalue = rate * value;
    out.gap = -out.eigenvalue->real();
    out.residual_norm = Real{0};
    if (!detail::finite(*out.eigenvalue) || !(*out.gap > Real{0}))
    {
      out.eigenvalue.reset();
      out.gap.reset();
      out.status = Status::precision_limit;
      return out;
    }
    out.converged = true;
    out.status = Status::converged;
    return out;
  }
  Real const rho = Real(n) / Real(l);
  out.log_y = C(Real(l) * std::log(Real{2}) + Real(n) * std::log(rho) + Real(l - n) * std::log1p(-rho) + Real{3.5}, pi);
  std::vector<C> all;
  if (!detail::seed(l, n, out.log_y, all, out.seed_iterations, options))
  {
    out.status = out.seed_iterations == options.max_seed_iterations ? Status::seed_limit : Status::precision_limit;
    return out;
  }
  std::sort(all.begin(), all.end(), [](C a, C b) { return a.real() > b.real(); });
  out.roots.assign(all.begin(), all.begin() + n);
  // Golinelli-Mallick c1: replace Z_N by Z_(N+1), the upper roots
  // straddling the Cassini pinch. Ignore tiny imaginary noise on real roots.
  Real const threshold = std::sqrt(uni20::numeric_limits<Real>::epsilon());
  std::size_t remove = n, add = l;
  for (std::size_t j = 0; j < n; ++j)
    if (all[j].imag() > threshold) remove = j;
  for (std::size_t j = n; j < l; ++j)
    if (all[j].imag() > threshold)
    {
      add = j;
      break;
    }
  if (remove == n || add == l)
  {
    out.status = Status::precision_limit;
    return out;
  }
  out.roots[remove] = all[add];
  detail::System<Real> system{l, n};
  for (;;)
  {
    auto evaluation = system.evaluate(out.roots, out.log_y, true);
    out.residual_norm = evaluation.norm;
    if (!uni20::isfinite(evaluation.norm))
    {
      out.status = Status::precision_limit;
      return out;
    }
    if (evaluation.norm <= options.tolerance) break;
    if (out.iterations == options.max_iterations)
    {
      out.status = Status::iteration_limit;
      return out;
    }
    for (auto& v : evaluation.residual)
      v = -v;
    if (!bethe::detail::newton_step(std::move(evaluation.jacobian), evaluation.residual))
    {
      out.status = Status::precision_limit;
      return out;
    }
    auto trial = out.roots;
    C b{};
    Real damping{1};
    bool accepted = false;
    for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits; ++backtrack)
    {
      for (std::size_t j = 0; j < n; ++j)
        trial[j] = out.roots[j] + damping * C(evaluation.residual[2 * j], evaluation.residual[2 * j + 1]);
      b = out.log_y + damping * C(evaluation.residual[2 * n], evaluation.residual[2 * n + 1]);
      if (system.evaluate(trial, b, false).norm < evaluation.norm)
      {
        accepted = true;
        break;
      }
      damping /= Real{2};
    }
    if (!accepted)
    {
      out.status = Status::stalled;
      return out;
    }
    out.roots = std::move(trial);
    out.log_y = b;
    ++out.iterations;
  }
  bethe::detail::CompensatedSum<C> energy;
  bethe::detail::CompensatedSum<C> translation;
  for (C z : out.roots)
  {
    energy.add((z - Real{1}) / Real{2});
    translation.add(std::log((Real{1} + z) / Real{2}));
  }
  // Distinct regular roots and a first-harmonic translation phase are
  // necessary checks on the selected relaxation branch, not completeness.
  Real const phase = std::abs(detail::wrap(translation.value().imag()));
  if (!detail::finite(translation.value()) ||
      std::abs(translation.value().real()) > Real{8} * Real(l) * options.tolerance ||
      std::abs(phase - Real{2} * pi / Real(l)) > Real{8} * Real(l) * options.tolerance)
  {
    out.status = Status::precision_limit;
    return out;
  }
  for (std::size_t j = 0; j < n; ++j)
    for (std::size_t k = 0; k < j; ++k)
      if (std::abs(out.roots[j] - out.roots[k]) < threshold)
      {
        out.status = Status::precision_limit;
        return out;
      }
  C value = energy.value();
  if (n == l - n)
  {
    if (std::abs(value.imag()) > Real{8} * Real(l) * options.tolerance)
    {
      out.status = Status::precision_limit;
      return out;
    }
    value.imag(Real{0}); // The half-filled leading eigenvalue is real.
  }
  if (value.imag() < Real{0})
  {
    value = std::conj(value);
    out.log_y = std::conj(out.log_y);
    for (auto& z : out.roots)
      z = std::conj(z);
  }
  value *= rate;
  if (!detail::finite(value) || !(value.real() < Real{0}))
  {
    out.status = Status::precision_limit;
    return out;
  }
  out.eigenvalue = value;
  out.gap = -value.real();
  out.converged = true;
  out.status = Status::converged;
  return out;
}
} // namespace bethe::tasep
