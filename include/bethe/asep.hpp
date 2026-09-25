// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/complex_math.hpp>
#include <bethe/tasep.hpp>

namespace bethe::asep
{
enum class Status
{
  converged,
  stationary_only,
  seed_limit,
  iteration_limit,
  continuation_limit,
  stalled,
  precision_limit
};
template <uni20::Real Real> struct Options
{
    Real tolerance = Real{128} * uni20::numeric_limits<Real>::epsilon();
    std::size_t max_iterations = 100, max_continuation_steps = 256;
    tasep::Options<Real> seed_options{};
};
template <uni20::Real Real> struct GapState
{
    std::size_t sites = 0, particles = 0, effective_particles = 0, wave_index = 0;
    std::size_t iterations = 0, continuation_steps = 0, seed_iterations = 0, seed_newton_iterations = 0;
    Real right_rate{}, left_rate{}, reached_ratio{}, residual_norm = uni20::numeric_limits<Real>::infinity();
    std::complex<Real> wave_base{}; // z0-1, exact first harmonic
    // z_j=1+(1-x)*v_j, except z_wave=1+wave_base+(1-x)*v_wave.
    // Analytic symmetric/one-particle results need no regular Bethe roots.
    std::vector<std::complex<Real>> scaled_roots;
    std::optional<std::complex<Real>> eigenvalue;
    std::optional<Real> gap;
    bool converged = false, analytic = false;
    Status status = Status::iteration_limit;
};
namespace detail
{
template <uni20::Real Real> struct System
{
    using C = std::complex<Real>;
    std::size_t l, n, wave;
    Real x, delta;
    C base;
    struct Evaluation
    {
        std::vector<Real> residual, jacobian;
        Real norm{};
    };
    Evaluation evaluate(std::span<C const> v, bool jacobian = false) const
    {
      auto const order = 2 * n;
      Evaluation out{std::vector<Real>(order), jacobian ? std::vector<Real>(order * order) : std::vector<Real>{}};
      auto put = [&](std::size_t i, std::size_t j, C value) {
        if (!jacobian) return;
        out.jacobian[(2 * i) * order + 2 * j] = value.real();
        out.jacobian[(2 * i) * order + 2 * j + 1] = -value.imag();
        out.jacobian[(2 * i + 1) * order + 2 * j] = value.imag();
        out.jacobian[(2 * i + 1) * order + 2 * j + 1] = value.real();
      };
      auto u = [&](std::size_t j) { return delta * v[j] + (j == wave ? base : C{}); };
      for (std::size_t i = 0; i < n; ++i)
      {
        C const ui = u(i), reference = i == wave ? Real{1} + base : C(1);
        bethe::detail::CompensatedSum<C> f, diagonal;
        f.add(Real(l) * bethe::detail::complex_log1p(delta * v[i] / reference) / delta);
        diagonal.add(Real(l) / (Real{1} + ui));
        for (std::size_t j = 0; j < n; ++j)
          if (j != i)
          {
            C const uj = u(j);
            C a, b, ai, aj, bi, bj;
            if (i == wave)
            {
              a = ui + delta * v[j] - Real{2} * x * ui * v[j];
              b = x * ui * (Real{1} + delta * v[j]) - delta * v[j];
              ai = delta * (Real{1} - Real{2} * x * v[j]);
              aj = delta - Real{2} * x * ui;
              bi = x * delta * (Real{1} + delta * v[j]);
              bj = delta * (x * ui - Real{1});
            }
            else if (j == wave)
            {
              a = uj + delta * v[i] - Real{2} * x * v[i] * uj;
              b = delta * x * v[i] * (uj + Real{1}) - uj;
              ai = delta - Real{2} * x * uj;
              aj = delta * (Real{1} - Real{2} * x * v[i]);
              bi = delta * x * (uj + Real{1});
              bj = delta * (delta * x * v[i] - Real{1});
            }
            else
            {
              a = v[i] + v[j] - Real{2} * x * v[i] * v[j];
              b = x * delta * v[i] * v[j] + x * v[i] - v[j];
              ai = Real{1} - Real{2} * x * v[j];
              aj = Real{1} - Real{2} * x * v[i];
              bi = x * delta * v[j] + x;
              bj = x * delta * v[i] - Real{1};
            }
            C const s = a / b;
            f.add(-bethe::detail::complex_log1p(delta * s) / delta);
            diagonal.add(-(ai - s * bi) / (b * (Real{1} + delta * s)));
            put(i, j, -(aj - s * bj) / (b * (Real{1} + delta * s)));
          }
        C value = f.value();
        value.imag(tasep::detail::wrap(value.imag() * delta) / delta);
        if (!tasep::detail::finite(value))
        {
          out.norm = uni20::numeric_limits<Real>::infinity();
          return out;
        }
        put(i, i, diagonal.value());
        out.residual[2 * i] = value.real();
        out.residual[2 * i + 1] = value.imag();
        out.norm = std::max(out.norm, std::abs(value) / Real(l));
      }
      return out;
    }
};

template <uni20::Real Real>
Status correct(System<Real> const& system, std::vector<std::complex<Real>>& v, Options<Real> const& options,
               std::size_t& total_iterations)
{
  using C = std::complex<Real>;
  for (std::size_t iteration = 0;; ++iteration)
  {
    auto e = system.evaluate(v, true);
    if (!uni20::isfinite(e.norm)) return Status::precision_limit;
    if (e.norm <= options.tolerance) return Status::converged;
    if (iteration == options.max_iterations) return Status::iteration_limit;
    for (auto& r : e.residual)
      r = -r;
    if (!bethe::detail::newton_step(std::move(e.jacobian), e.residual)) return Status::precision_limit;
    auto trial = v;
    Real damping{1};
    bool accepted = false;
    for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits; ++backtrack)
    {
      for (std::size_t j = 0; j < v.size(); ++j)
        trial[j] = v[j] + damping * C(e.residual[2 * j], e.residual[2 * j + 1]);
      if (system.evaluate(trial).norm < e.norm)
      {
        accepted = true;
        break;
      }
      damping /= Real{2};
    }
    if (!accepted) return Status::stalled;
    v = std::move(trial);
    ++total_iterations;
  }
}
} // namespace detail

/// Leading periodic ASEP relaxation mode. Rates need not sum to one.
/// Reflection/particle-hole symmetry selects an Im(lambda)>=0 representative.
template <uni20::Real Real>
GapState<Real> relaxation_gap(std::size_t l, std::size_t particles, Real right_rate, Real left_rate,
                              Options<Real> options = {})
{
  using C = std::complex<Real>;
  if (l < 2 || particles > l || !uni20::isfinite(right_rate) || !uni20::isfinite(left_rate) || right_rate < Real{0} ||
      left_rate < Real{0} || !(std::max(right_rate, left_rate) > Real{0}) || !uni20::isfinite(options.tolerance) ||
      !(options.tolerance > Real{0}))
    throw std::invalid_argument(
        "ASEP requires L>=2, 0<=N<=L, finite nonnegative rates with one positive, and positive tolerance");
  if (l > options.seed_options.max_sites) throw std::length_error("ASEP L exceeds max_sites work budget");
  GapState<Real> out;
  out.sites = l;
  out.particles = particles;
  out.effective_particles = std::min(particles, l - particles);
  out.right_rate = right_rate;
  out.left_rate = left_rate;
  auto const n = out.effective_particles;
  Real const scale = std::max(right_rate, left_rate), small = std::min(right_rate, left_rate), x = small / scale;
  Real const delta = (scale - small) / scale, pi = Real{4} * std::atan(Real{1}), sine = std::sin(pi / Real(l));
  auto publish = [&](C value) {
    value *= scale;
    if (value.imag() < Real{0})
    {
      value = std::conj(value);
      out.wave_base = std::conj(out.wave_base);
      for (auto& root : out.scaled_roots)
        root = std::conj(root);
    }
    if (!tasep::detail::finite(value) || !(value.real() < Real{0}))
    {
      out.status = Status::precision_limit;
      return;
    }
    out.eigenvalue = value;
    out.gap = -value.real();
    out.converged = true;
    out.status = Status::converged;
  };
  if (n == 0)
  {
    out.converged = true;
    out.analytic = true;
    out.status = Status::stationary_only;
    out.residual_norm = Real{0};
    out.reached_ratio = x;
    return out;
  }
  if (n == 1 || right_rate == left_rate)
  {
    out.analytic = true;
    out.residual_norm = Real{0};
    out.reached_ratio = x;
    publish(C(-Real{2} * (Real{1} + x) * sine * sine,
              l == 2 || right_rate == left_rate ? Real{0} : delta * std::sin(Real{2} * pi / Real(l))));
    return out;
  }
  if ((small > Real{0} && x == Real{0}) || !(delta > Real{0}))
  {
    out.status = Status::precision_limit;
    return out;
  }
  auto const seed = tasep::relaxation_gap(l, particles, Real{1}, options.seed_options);
  out.seed_iterations = seed.seed_iterations;
  out.seed_newton_iterations = seed.iterations;
  out.residual_norm = seed.residual_norm;
  if (!seed.converged)
  {
    out.status = seed.status == tasep::Status::seed_limit || seed.status == tasep::Status::iteration_limit
                     ? Status::seed_limit
                     : Status::precision_limit;
    return out;
  }
  std::vector<C> z;
  bethe::detail::CompensatedSum<C> translation;
  for (C root : seed.roots)
  {
    z.push_back(Real{2} / (Real{1} + root));
    translation.add(std::log(z.back()));
  }
  Real const theta = std::copysign(Real{2} * pi / Real(l), tasep::detail::wrap(translation.value().imag()));
  out.wave_base = C(-Real{2} * sine * sine, std::sin(theta));
  for (std::size_t j = 0; j < n; ++j)
    if (std::abs(z[j] - Real{1}) > std::abs(z[out.wave_index] - Real{1})) out.wave_index = j;
  for (std::size_t j = 0; j < n; ++j)
    out.scaled_roots.push_back(z[j] - Real{1} - (j == out.wave_index ? out.wave_base : C{}));
  if (x == Real{0})
  {
    out.residual_norm = seed.residual_norm;
    publish(*seed.eigenvalue);
    return out;
  }
  Real step = Real{1} / Real{20};
  while (out.reached_ratio < x)
  {
    if (out.continuation_steps == options.max_continuation_steps)
    {
      out.status = Status::continuation_limit;
      return out;
    }
    Real const next = std::min(x, out.reached_ratio + step);
    if (next == out.reached_ratio)
    {
      out.status = Status::precision_limit;
      return out;
    }
    detail::System<Real> system{l, n, out.wave_index, next, next == x ? delta : Real{1} - next, out.wave_base};
    auto trial = out.scaled_roots;
    ++out.continuation_steps;
    auto const status = detail::correct(system, trial, options, out.iterations);
    if (status != Status::converged)
    {
      if (options.max_iterations == 0)
      {
        out.status = status;
        return out;
      }
      step /= Real{2};
      continue;
    }
    out.scaled_roots = std::move(trial);
    out.reached_ratio = next;
    out.residual_norm = system.evaluate(out.scaled_roots).norm;
    step = std::min(Real{1} / Real{20}, step * Real{1.5});
  }
  // Retain tiny drift frequencies: evaluate the finite wave analytically,
  // then its correction, rather than subtracting nearly conjugate O(1) terms.
  C const z0 = Real{1} + out.wave_base;
  // Test non-wave separation in scaled coordinates: physical z separations
  // legitimately vanish with the bias and must not be compared to sqrt(eps).
  Real const separation = Real{32} * uni20::numeric_limits<Real>::epsilon();
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < i; ++j)
    {
      C difference = out.scaled_roots[i] - out.scaled_roots[j];
      if (i == out.wave_index || j == out.wave_index)
        difference = delta * difference + (i == out.wave_index ? out.wave_base : -out.wave_base);
      if (!tasep::detail::finite(difference) || std::abs(difference) <= separation)
      {
        out.status = Status::precision_limit;
        return out;
      }
    }
  bethe::detail::CompensatedSum<C> eigen, phase;
  eigen.add(C(-Real{2} * (Real{1} + x) * sine * sine, -delta * std::sin(theta)));
  for (std::size_t j = 0; j < n; ++j)
  {
    C const u = delta * out.scaled_roots[j];
    if (j == out.wave_index)
    {
      eigen.add(u * (x - Real{1} / (z0 * (z0 + u))));
      phase.add(bethe::detail::complex_log1p(u / z0));
    }
    else
    {
      eigen.add(delta * delta * out.scaled_roots[j] * (x * out.scaled_roots[j] - Real{1}) / (Real{1} + u));
      phase.add(bethe::detail::complex_log1p(u));
    }
  }
  if (std::abs(phase.value()) > Real{8} * Real(l) * delta * options.tolerance)
  {
    out.status = Status::precision_limit;
    return out;
  }
  C value = eigen.value();
  if (n == l - n)
  {
    if (std::abs(value.imag()) > Real{8} * Real(l) * delta * options.tolerance)
    {
      out.status = Status::precision_limit;
      return out;
    }
    value.imag(Real{0});
  }
  publish(value);
  return out;
}
} // namespace bethe::asep
