// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/quadrature.hpp>
#include <optional>
#include <uni20/core/finite_quadrature.hpp>

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
    uni20::Precision precision() const { return uni20::native_precision<Real>(); }
};

#if UNI20_HAS_MPREAL
template <> struct Options<uni20::mpreal>
{
    explicit Options(uni20::Precision p) : tolerance(default_tolerance(p)), precision_(p) {}
    uni20::mpreal tolerance;
    std::size_t max_series_terms = 1024, max_lobes = 1024;
    std::size_t max_evaluations = 200000, max_quadrature_levels = 12;
    uni20::Precision precision() const { return precision_; }

  private:
    static uni20::mpreal default_tolerance(uni20::Precision p)
    {
      if (p.is_exact()) throw std::invalid_argument("Kondo response requires finite working precision");
      return 1048576 * uni20::numeric_limits<uni20::mpreal>::epsilon(p);
    }
    uni20::Precision precision_;
};
#endif
template <uni20::Real Real> struct State
{
    State()
      requires(!uni20::has_runtime_precision_v<Real>)
        : State(uni20::native_precision<Real>())
    {}
    explicit State(uni20::Precision p)
        : field(uni20::scalar_at_precision<Real>(0, p)), scale(field),
          magnetization_error(uni20::scalar_like(field, uni20::numeric_limits<double>::infinity())),
          scaled_energy_error(magnetization_error)
    {}
    Real field, scale;
    std::optional<Real> magnetization, energy_change, zero_field_susceptibility;
    Real magnetization_error, scaled_energy_error;
    std::size_t series_terms = 0, lobes = 0, evaluations = 0;
    bool converged = false;
    Status status = Status::series_limit;
};
namespace detail
{
namespace math = uni20::math;
template <uni20::Real Real> Real pi(Real const& one) { return Real{4} * math::atan(one); }
template <uni20::Real Real> void validate(Real const& field, Real const& scale, Options<Real> const& options)
{
  Real const zero = uni20::scalar_at_precision<Real>(0, options.precision());
  for (auto const* x : {&field, &scale, &options.tolerance})
    bethe::detail::require_same_precision(*x, zero);
  if (!uni20::isfinite(field) || !uni20::isfinite(scale) || !(scale > Real{0}) || !uni20::isfinite(options.tolerance) ||
      !(options.tolerance > Real{0}))
    throw std::invalid_argument("Kondo field must be finite, scale and tolerance finite and positive");
  if (options.max_series_terms > 4096 || options.max_lobes > 4096)
    throw std::length_error("Kondo series/lobe work ceiling is 4096");
}
template <uni20::Real Real> Real pi() { return Real{4} * std::atan(Real{1}); }
template <uni20::Real Real> struct Estimate
{
    Real value{}, error{};
};

// Gamma(w+1/2)*exp(w)*w^(-w), w>=0. This restricted scaled amplitude
// avoids a general Gamma dependency and stays native for binary128.
// Stirling expansion/remainder: DLMF 5.11.1 and 5.11(ii), positive real z.
template <uni20::Real Real> Estimate<Real> native_amplitude(Real w)
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

// Runtime precision uses Uni20's MPFR log-Gamma facade. The logarithmic
// form avoids overflow from Gamma(w+1/2) and exp(w) separately. This amplitude
// is sampled only for 0<=w<=max_lobes<=4096; the roundoff estimate includes the
// cancellation in log-Gamma - w*log(w) + w.
template <uni20::Real Real> Estimate<Real> amplitude(Real w)
{
  if constexpr (uni20::has_runtime_precision_v<Real>)
  {
    Real const zero = bethe::detail::precision_zero(w), one = uni20::scalar_like(w, 1);
    Real const eps = uni20::numeric_limits<Real>::epsilon(uni20::precision_of(w));
    if (!uni20::isfinite(w) || w < zero) throw std::invalid_argument("Kondo amplitude requires finite w>=0");
    if (w == zero) return {math::sqrt(pi(one)), 16 * eps};
    Real const g = math::lgamma(w + one / 2), wl = w * math::log(w);
    bethe::detail::CompensatedSum<Real> sum(zero);
    sum.add(g);
    sum.add(-wl);
    sum.add(w);
    Real const value = math::exp(sum.value());
    Real const log_error = 32 * eps * (math::abs(g) + math::abs(wl) + w + 4);
    return {value, value * math::expm1(log_error)};
  }
  else
    return native_amplitude(w);
}

template <uni20::Real Real> struct Pair
{
    std::array<Real, 2> value{}, error{};
};
// Euler transformation of an alternating sequence. Halve at each difference
// to avoid exponentially large intermediate values; propagate input errors.
template <uni20::Real Real> Pair<Real> euler(std::vector<Pair<Real>> values)
{
  if (values.empty()) throw std::invalid_argument("Euler transformation requires a nonempty sequence");
  Real const zero = bethe::detail::precision_zero(values.front().value[0]);
  Real const eps = uni20::numeric_limits<Real>::epsilon(uni20::precision_of(zero));
  std::array<bethe::detail::CompensatedSum<Real>, 2> sum, errors;
  for (std::size_t j = 0; j < 2; ++j)
    sum[j] = errors[j] = bethe::detail::CompensatedSum<Real>(zero);
  for (std::size_t count = values.size(); count; --count)
    for (std::size_t j = 0; j < 2; ++j)
    {
      sum[j].add(values[0].value[j] / Real{2});
      errors[j].add(values[0].error[j] / Real{2} + Real{4} * eps * uni20::math::abs(values[0].value[j]));
      for (std::size_t i = 0; i + 1 < count; ++i)
      {
        values[i].error[j] = (values[i].error[j] + values[i + 1].error[j]) / Real{2} +
                             eps * (uni20::math::abs(values[i].value[j]) + uni20::math::abs(values[i + 1].value[j]));
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
    now.error[j] += uni20::math::abs(now.value[j] - before.value[j]);
    good = good && uni20::isfinite(now.error[j]) && now.error[j] <= tolerance;
  }
  return good;
}
// Values: M/x and -(Delta E/|b|)/x; normalized to retain small-field accuracy.
template <uni20::Real Real> std::optional<Pair<Real>> low(Real x, Options<Real> const& options, State<Real>& work)
{
  Real const one = uni20::scalar_like(x, 1);
  Real const eps = uni20::numeric_limits<Real>::epsilon(options.precision());
  Real a = one / uni20::math::sqrt(2 * pi(one) * uni20::math::exp(one)), power = one;
  std::vector<Pair<Real>> terms;
  Pair<Real> previous;
  unsigned stable = 0;
  for (std::size_t n = 0; n < options.max_series_terms; ++n)
  {
    if (n)
    {
      Real const t = uni20::scalar_like(one, n) - one / 2;
      a *= t / Real(n) * uni20::math::exp(t * uni20::math::log1p(Real{1} / t) - Real{1});
      power *= x * x;
    }
    Real const v = a * power, e = Real{64} * eps * Real(n + 1) * uni20::math::abs(v);
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
  Real const zero = uni20::scalar_like(log_x, 0), one = uni20::scalar_like(log_x, 1);
  Real const p = pi(one), c = one / (2 * p * uni20::math::sqrt(p)), inverse = uni20::math::exp(-log_x);
  Real const eps = uni20::numeric_limits<Real>::epsilon(options.precision());
  std::vector<Pair<Real>> lobes;
  Pair<Real> previous;
  unsigned stable = 0;
  for (std::size_t n = 0; n < options.max_lobes; ++n)
  {
    auto f = [&](Real t, Real sine) {
      Real const w = Real(n) + t, exponent = (Real{1} - Real{2} * w) * log_x;
      Real const damping = uni20::math::exp(-Real{2} * w * log_x);
      Real integral;
      if (w == (one / 2))
        integral = inverse * log_x;
      else if (uni20::math::abs(exponent) < (one / 2))
        integral = inverse * uni20::math::expm1(exponent) / (Real{1} - Real{2} * w);
      else
        integral = (damping - inverse) / (Real{1} - Real{2} * w);
      auto const a = amplitude(w);
      Real const factor = c * sine / w;
      Real const m = factor * a.value * damping, e = factor * a.value * integral;
      Real const uncertainty = a.error / a.value + Real{64} * eps * (Real{1} + w * log_x);
      return std::array<Real, 4>{m, e, uni20::math::abs(m) * uncertainty, uni20::math::abs(e) * uncertainty};
    };
    auto const q = [&] {
      if constexpr (uni20::has_runtime_precision_v<Real>)
      {
        uni20::tanh_sinh_options<Real, 4> controls(options.precision());
        controls.relative_tolerance = options.tolerance / 64;
        controls.absolute_tolerance.fill(controls.relative_tolerance);
        controls.max_evaluations = options.max_evaluations - work.evaluations;
        controls.max_levels = options.max_quadrature_levels;
        auto const integral = uni20::tanh_sinh_integrate<Real, 4>(
            [&](uni20::quadrature_point<Real> const& point) {
              // Preserve sin(pi*t) near either endpoint even when t rounds
              // to one. The remaining factors are regular there.
              Real const distance = std::min(point.lower_distance, point.upper_distance);
              return f(point.lower_distance, uni20::math::sin(p * distance));
            },
            zero, one, controls);
        work.evaluations += integral.evaluations;
        return bethe::detail::QuadratureResult<Real, 4>{integral.value, integral.error,
                                                        integral.status == uni20::quadrature_status::converged};
      }
      else
        return bethe::detail::tanh_sinh<Real, 4>(
            [&](Real t) { return f(t, uni20::math::sin(p * t)); }, zero, one, options.tolerance / Real{64},
            work.evaluations, options.max_evaluations, options.max_quadrature_levels, {one, one, one, one});
    }();
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
      current.value[0] = (one / 2) - current.value[0];
      current.value[1] += -anchor.value[1] * inverse + uni20::math::expm1(-log_x) / Real{2};
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
  detail::validate(field, scale, options);
  Real const zero = uni20::scalar_at_precision<Real>(0, options.precision()), one = uni20::scalar_like(zero, 1);
  State<Real> out(options.precision());
  out.field = field;
  out.scale = scale;
  Real const a0 = one / uni20::math::sqrt(2 * detail::pi(one) * uni20::math::exp(one)), chi = a0 / scale;
  Real const magnitude = uni20::math::abs(field);
  detail::Pair<Real> answer{{zero, zero}, {zero, zero}};
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
      auto const anchor = detail::low(one, options, out);
      if (!anchor) return out;
      Real const ratio = magnitude / scale;
      Real const log_x =
          uni20::isfinite(ratio) ? uni20::math::log(ratio) : uni20::math::log(magnitude) - uni20::math::log(scale);
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
  out.magnetization = uni20::math::copysign(answer.value[0], field);
  out.energy_change = energy;
  out.zero_field_susceptibility = chi;
  out.converged = true;
  out.status = Status::converged;
  return out;
}
} // namespace bethe::kondo
