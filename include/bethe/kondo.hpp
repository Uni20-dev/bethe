// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/quadrature.hpp>
#include <optional>

namespace bethe::kondo
{
enum class Status
{
  converged,
  series_limit,
  quadrature_limit,
  tail_limit,
  precision_limit
};
template <uni20::Real Real> struct Options
{
    // Absolute estimates in M and Delta E/|b|, not absolute physical energy.
    Real tolerance = Real{1048576} * uni20::numeric_limits<Real>::epsilon();
    std::size_t max_series_terms = 256, max_lobes = 256;
    std::size_t max_evaluations = 200000, max_quadrature_levels = 12;
};
template <uni20::Real Real> struct State
{
    Real field{}, scale{};
    std::optional<Real> magnetization, energy_change, zero_field_susceptibility;
    Real magnetization_error = uni20::numeric_limits<Real>::infinity();
    Real scaled_energy_error = uni20::numeric_limits<Real>::infinity();
    std::size_t series_terms = 0, lobes = 0, evaluations = 0;
    bool converged = false;
    Status status = Status::series_limit;
};
namespace detail
{
template <uni20::Real Real> Real pi() { return Real{4} * std::atan(Real{1}); }
template <uni20::Real Real> struct Estimate
{
    Real value{}, error{};
};

// Gamma(w+1/2)*exp(w)*w^(-w), w>=0. This restricted scaled amplitude
// avoids a general Gamma dependency and stays native for binary128.
// Stirling expansion/remainder: DLMF 5.11.1 and 5.11(ii), positive real z.
template <uni20::Real Real> Estimate<Real> amplitude(Real w)
{
  Real const eps = uni20::numeric_limits<Real>::epsilon(), half_log = std::log(Real{2} * pi<Real>()) / Real{2};
  if (w == Real{0}) return {std::sqrt(pi<Real>()), Real{16} * eps};
  Real z = w + Real{0.5}, base{}, magnitude{};
  if (w >= Real{32})
  {
    base = w * std::log1p(Real{0.5} / w) - Real{0.5} + half_log;
    magnitude = Real{4};
  }
  else
  {
    bethe::detail::CompensatedSum<Real> logs;
    while (z < Real{32})
    {
      logs.add(std::log(z));
      z += Real{1};
    }
    Real const main = (z - Real{0.5}) * std::log(z), wl = w * std::log(w);
    bethe::detail::CompensatedSum<Real> sum;
    for (Real term : {main, -z, half_log, -logs.value(), -wl, w})
      sum.add(term);
    base = sum.value();
    magnitude = std::abs(main) + z + std::abs(logs.value()) + std::abs(wl) + w + Real{4};
  }
  std::array<Real, 17> const c{Real{1} / Real{12},
                               -Real{1} / Real{360},
                               Real{1} / Real{1260},
                               -Real{1} / Real{1680},
                               Real{1} / Real{1188},
                               -Real{691} / Real{360360},
                               Real{1} / Real{156},
                               -Real{3617} / Real{122400},
                               Real{43867} / Real{244188},
                               -Real{174611} / Real{125400},
                               Real{77683} / Real{5796},
                               -Real{236364091} / Real{1506960},
                               Real{657931} / Real{300},
                               -Real{3392780147LL} / Real{93960},
                               Real{1723168255201LL} / Real{2492028},
                               -Real{7709321041217LL} / Real{505920},
                               Real{151628697551LL} / Real{396}};
  Real power = Real{1} / z, step = power * power;
  bethe::detail::CompensatedSum<Real> correction;
  for (std::size_t n = 0; n + 1 < c.size(); ++n)
  {
    correction.add(c[n] * power);
    power *= step;
  }
  Real const value = std::exp(base + correction.value());
  Real const log_error = std::abs(c.back() * power) + Real{32} * eps * magnitude;
  return {value, value * std::expm1(log_error)};
}

template <uni20::Real Real> struct Pair
{
    std::array<Real, 2> value{}, error{};
};
// Euler transformation of an alternating sequence. Halve at each difference
// to avoid exponentially large intermediate values; propagate input errors.
template <uni20::Real Real> Pair<Real> euler(std::vector<Pair<Real>> values)
{
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  std::array<bethe::detail::CompensatedSum<Real>, 2> sum, errors;
  for (std::size_t count = values.size(); count; --count)
    for (std::size_t j = 0; j < 2; ++j)
    {
      sum[j].add(values[0].value[j] / Real{2});
      errors[j].add(values[0].error[j] / Real{2} + Real{4} * eps * std::abs(values[0].value[j]));
      for (std::size_t i = 0; i + 1 < count; ++i)
      {
        values[i].error[j] = (values[i].error[j] + values[i + 1].error[j]) / Real{2} +
                             eps * (std::abs(values[i].value[j]) + std::abs(values[i + 1].value[j]));
        values[i].value[j] = (values[i].value[j] - values[i + 1].value[j]) / Real{2};
      }
    }
  return {{sum[0].value(), sum[1].value()}, {errors[0].value(), errors[1].value()}};
}
template <uni20::Real Real> bool resolved(Pair<Real>& now, Pair<Real> const& before, Real tolerance)
{
  bool good = true;
  for (std::size_t j = 0; j < 2; ++j)
  {
    now.error[j] += std::abs(now.value[j] - before.value[j]);
    good = good && uni20::isfinite(now.error[j]) && now.error[j] <= tolerance;
  }
  return good;
}
// Values: M/x and -(Delta E/|b|)/x; normalized to retain small-field accuracy.
template <uni20::Real Real> std::optional<Pair<Real>> low(Real x, Options<Real> const& options, State<Real>& work)
{
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real a = Real{1} / std::sqrt(Real{2} * pi<Real>() * std::exp(Real{1})), power = Real{1};
  std::vector<Pair<Real>> terms;
  Pair<Real> previous;
  unsigned stable = 0;
  for (std::size_t n = 0; n < options.max_series_terms; ++n)
  {
    if (n)
    {
      Real const t = Real(n) - Real{0.5};
      a *= t / Real(n) * std::exp(t * std::log1p(Real{1} / t) - Real{1});
      power *= x * x;
    }
    Real const v = a * power, e = Real{64} * eps * Real(n + 1) * std::abs(v);
    terms.push_back({{v, v / Real(2 * n + 2)}, {e, e / Real(2 * n + 2)}});
    ++work.series_terms;
    if (terms.size() % 16) continue;
    auto current = euler(terms);
    stable = terms.size() > 16 && resolved(current, previous, options.tolerance / Real{8}) ? stable + 1 : 0;
    previous = current;
    if (stable == 2) return current;
  }
  work.status = Status::series_limit;
  return {};
}

// Values are M and Delta E/|b|; log_x>=0 also supports unrepresentable x.
template <uni20::Real Real>
std::optional<Pair<Real>> high(Real log_x, Pair<Real> const& anchor, Options<Real> const& options, State<Real>& work)
{
  Real const p = pi<Real>(), c = Real{1} / (Real{2} * p * std::sqrt(p)), inverse = std::exp(-log_x);
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  std::vector<Pair<Real>> lobes;
  Pair<Real> previous;
  unsigned stable = 0;
  for (std::size_t n = 0; n < options.max_lobes; ++n)
  {
    auto f = [&](Real t) {
      Real const w = Real(n) + t, exponent = (Real{1} - Real{2} * w) * log_x;
      Real const damping = std::exp(-Real{2} * w * log_x);
      Real integral;
      if (w == Real{0.5})
        integral = inverse * log_x;
      else if (std::abs(exponent) < Real{0.5})
        integral = inverse * std::expm1(exponent) / (Real{1} - Real{2} * w);
      else
        integral = (damping - inverse) / (Real{1} - Real{2} * w);
      auto const a = amplitude(w);
      Real const factor = c * std::sin(p * t) / w;
      Real const m = factor * a.value * damping, e = factor * a.value * integral;
      Real const uncertainty = a.error / a.value + Real{64} * eps * (Real{1} + w * log_x);
      return std::array<Real, 4>{m, e, std::abs(m) * uncertainty, std::abs(e) * uncertainty};
    };
    auto const q = bethe::detail::tanh_sinh<Real, 4>(
        f, Real{0}, Real{1}, options.tolerance / Real{64}, work.evaluations, options.max_evaluations,
        options.max_quadrature_levels, {Real{1}, Real{1}, Real{1}, Real{1}});
    ++work.lobes;
    if (!q.converged)
    {
      work.status = Status::quadrature_limit;
      return {};
    }
    lobes.push_back({{q.value[0], q.value[1]}, {q.error[0] + q.value[2], q.error[1] + q.value[3]}});
    if (lobes.size() % 16) continue;
    auto current = euler(lobes);
    stable = lobes.size() > 16 && resolved(current, previous, options.tolerance / Real{4}) ? stable + 1 : 0;
    previous = current;
    if (stable == 2)
    {
      current.value[0] = Real{0.5} - current.value[0];
      current.value[1] += -anchor.value[1] * inverse + std::expm1(-log_x) / Real{2};
      current.error[0] += Real{8} * eps;
      current.error[1] += anchor.error[1] * inverse + Real{16} * eps;
      return current;
    }
  }
  work.status = Status::tail_limit;
  return {};
}
} // namespace detail

/// Universal T=0 impurity response, single-channel antiferromagnetic spin-1/2
/// Kondo scaling limit. Full Zeeman splitting b, scale T_B=2*T1 (Barcza 2020).
template <uni20::Real Real = double>
[[nodiscard]] State<Real> ground_response(Real field, Real scale, Options<Real> options = {})
{
  if (!uni20::isfinite(field) || !uni20::isfinite(scale) || !(scale > Real{0}) || !uni20::isfinite(options.tolerance) ||
      !(options.tolerance > Real{0}))
    throw std::invalid_argument("Kondo field must be finite, scale and tolerance finite and positive");
  if (options.max_series_terms > 4096 || options.max_lobes > 4096)
    throw std::length_error("Kondo series/lobe work ceiling is 4096");
  State<Real> out;
  out.field = field;
  out.scale = scale;
  Real const a0 = Real{1} / std::sqrt(Real{2} * detail::pi<Real>() * std::exp(Real{1})), chi = a0 / scale;
  Real const magnitude = std::abs(field);
  detail::Pair<Real> answer;
  if (magnitude != Real{0})
  {
    if (magnitude <= scale)
    {
      Real const x = magnitude / scale;
      if (x == Real{0})
      {
        out.status = Status::precision_limit;
        return out;
      }
      auto const low = detail::low(x, options, out);
      if (!low) return out;
      answer = *low;
      answer.value[0] *= x;
      answer.value[1] *= -x;
      for (auto& e : answer.error)
        e *= x;
    }
    else
    {
      auto const anchor = detail::low(Real{1}, options, out);
      if (!anchor) return out;
      Real const ratio = magnitude / scale;
      Real const log_x = uni20::isfinite(ratio) ? std::log(ratio) : std::log(magnitude) - std::log(scale);
      auto const high = detail::high(log_x, *anchor, options, out);
      if (!high) return out;
      answer = *high;
    }
  }
  Real const energy = magnitude * answer.value[1];
  out.magnetization_error = answer.error[0];
  out.scaled_energy_error = answer.error[1];
  if (!uni20::isfinite(chi) || !(chi > Real{0}) || !uni20::isfinite(energy) || !uni20::isfinite(answer.value[0]) ||
      (magnitude > Real{0} && (answer.value[0] == Real{0} || energy == Real{0})) ||
      answer.error[0] > options.tolerance || answer.error[1] > options.tolerance)
  {
    out.status = Status::precision_limit;
    return out;
  }
  out.magnetization = std::copysign(answer.value[0], field);
  out.energy_change = energy;
  out.zero_field_susceptibility = chi;
  out.converged = true;
  out.status = Status::converged;
  return out;
}
} // namespace bethe::kondo
