// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <array>
#include <bethe/solver.hpp>
#include <complex>
#include <cstdint>
#include <stdexcept>

namespace bethe::detail
{
/// theta_kind(pi*u | i*t), t>0. Derivatives are with respect to u, NOT pi*u.
/// All three entries have the SAME real exponential scale:
/// actual derivative[r] = exp(log_scale) * derivative[r].
/// Keep the scaled representation for ratios and logarithmic derivatives.
template <uni20::Real Real> struct ThetaJet
{
    std::array<std::complex<Real>, 3> derivative{};
    Real log_scale{};
};
namespace theta_detail
{
template <uni20::Real Real> bool finite(std::complex<Real> z)
{
  return uni20::isfinite(z.real()) && uni20::isfinite(z.imag());
}

template <uni20::Real Real> struct Sum
{
    CompensatedSum<Real> real, imag;
    void add(std::complex<Real> z)
    {
      real.add(z.real());
      imag.add(z.imag());
    }
    std::complex<Real> value() const { return {real.value(), imag.value()}; }
};

// Fundamental rectangle |Re(u)|<=1/2, |Im(u)|<=t/2. Both representations
// converge geometrically faster than exp(-pi*n*n) in their selected regimes.
template <uni20::Real Real> ThetaJet<Real> theta1_reduced(std::complex<Real> u, Real t, bool gaussian)
{
  using Complex = std::complex<Real>;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  auto const terms = static_cast<unsigned>(std::ceil(std::sqrt((-std::log(eps) + Real{32}) / pi))) + 3;
  Real const x = u.real(), y = u.imag();
  ThetaJet<Real> result;
  // Isolate the largest exponential before evaluating any series term.
  if (gaussian)
    result.log_scale = pi * (y * y - (Real{0.5} - std::abs(x)) * (Real{0.5} - std::abs(x))) / t - std::log(t) / Real{2};
  else
    result.log_scale = -pi * (t / Real{4} - std::abs(y));
  std::array<Sum<Real>, 3> sums;
  for (unsigned n = 0; n < terms; ++n)
  {
    Real const a = Real(n) + Real{0.5}, sign = n % 2 ? -Real{1} : Real{1};
    Complex value, first, second;
    if (!gaussian)
    {
      Real const omega = Real{2} * pi * a, fall = pi * t * (Real(n) * Real(n + 1));
      if (std::abs(omega * y) <= Real{1})
      {
        Real const weight = std::exp(-fall - pi * std::abs(y));
        value = Real{2} * sign * weight * std::sin(omega * u);
        first = Real{2} * sign * omega * weight * std::cos(omega * u);
      }
      else
      {
        Complex const positive = std::exp(Complex(-fall - omega * y - pi * std::abs(y), omega * x));
        Complex const negative = std::exp(Complex(-fall + omega * y - pi * std::abs(y), -omega * x));
        value = Complex(0, -sign) * (positive - negative);
        first = sign * omega * (positive + negative);
      }
      second = -omega * omega * value;
    }
    else
    {
      Real const omega = Real{2} * pi * a / t, fall = Real(n) * Real(n + 1);
      Complex difference, sum;
      if (std::abs(omega * x) <= Real{1})
      {
        Complex const weight = std::exp(Complex(-pi * (fall + std::abs(x)) / t, -Real{2} * pi * x * (y / t)));
        difference = Real{2} * weight * std::sinh(omega * u);
        sum = Real{2} * weight * std::cosh(omega * u);
      }
      else
      {
        Complex const positive =
            std::exp(Complex(-pi * (fall + std::abs(x) - Real{2} * a * x) / t, -Real{2} * pi * (x - a) * (y / t)));
        Complex const negative =
            std::exp(Complex(-pi * (fall + std::abs(x) + Real{2} * a * x) / t, -Real{2} * pi * (x + a) * (y / t)));
        difference = positive - negative;
        sum = positive + negative;
      }
      Complex const center = Real{2} * pi * (u / t);
      value = sign * difference;
      first = sign * (-center * difference + omega * sum);
      second =
          sign * ((center * center + omega * omega - Real{2} * pi / t) * difference - Real{2} * center * omega * sum);
    }
    sums[0].add(value);
    sums[1].add(first);
    sums[2].add(second);
  }
  for (unsigned r = 0; r < 3; ++r)
    result.derivative[r] = sums[r].value();
  if (u == Complex{})
  {
    result.derivative[0] = Complex{};
    result.derivative[2] = Complex{};
  }
  return result;
}

template <uni20::Real Real> ThetaJet<Real> even_reduced(unsigned kind, std::complex<Real> u, Real t)
{
  using Complex = std::complex<Real>;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  // Translate theta2 near its real zero without subtracting rounded cosines.
  if (kind == 2 && std::abs(u.real()) > Real{0.25})
  {
    bool const positive = u.real() > Real{0};
    auto result = theta1_reduced(u + (positive ? -Real{0.5} : Real{0.5}), t, false);
    if (positive)
      for (auto& d : result.derivative)
        d = -d;
    return result;
  }
  ThetaJet<Real> result;
  bool const half = kind == 2;
  result.log_scale = half ? -pi * (t / Real{4} - std::abs(u.imag())) : Real{0};
  std::array<Sum<Real>, 3> sums;
  if (!half) sums[0].add(Complex(1, 0));
  auto const terms = static_cast<unsigned>(std::ceil(std::sqrt((-std::log(eps) + Real{32}) / pi))) + 3;
  for (unsigned n = 0; n < terms; ++n)
  {
    Real const a = Real(n) + (half ? Real{0.5} : Real{1}), omega = Real{2} * pi * a;
    Real const sign = kind == 4 && n % 2 == 0 ? -Real{1} : Real{1};
    Real const fall = pi * t * (half ? Real(n) * Real(n + 1) : a * a), shift = half ? pi * std::abs(u.imag()) : Real{0};
    Complex value, first;
    if (std::abs(omega * u.imag()) <= Real{1})
    {
      Real const weight = Real{2} * sign * std::exp(-fall - shift);
      value = weight * std::cos(omega * u);
      first = -weight * omega * std::sin(omega * u);
    }
    else
    {
      Complex const aterm = std::exp(Complex(-fall - shift - omega * u.imag(), omega * u.real()));
      Complex const bterm = std::exp(Complex(-fall - shift + omega * u.imag(), -omega * u.real()));
      value = sign * (aterm + bterm);
      first = Complex(0, sign * omega) * (aterm - bterm);
    }
    sums[0].add(value);
    sums[1].add(first);
    sums[2].add(-omega * omega * value);
  }
  for (unsigned d = 0; d < 3; ++d)
    result.derivative[d] = sums[d].value();
  if (std::abs(u.imag()) == t / Real{2} &&
      ((kind == 3 && std::abs(u.real()) == Real{0.5}) || (kind == 4 && u.real() == Real{0})))
    result.derivative[0] = Complex{};
  return result;
}

template <uni20::Real Real> ThetaJet<Real> theta1(std::complex<Real> u, Real t, unsigned kind = 1)
{
  using Complex = std::complex<Real>;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  Real const limit = std::min(Real{1} / (Real{16} * eps), Real(std::numeric_limits<std::int64_t>::max() / 4));
  Real const vertical = u.imag() / t;
  if (!uni20::isfinite(vertical) || std::abs(u.real()) > limit || std::abs(vertical) > limit)
    throw std::overflow_error("theta argument reduction exceeds scalar precision");
  Real const nr = std::floor(u.real() + Real{0.5}), mr = std::floor(vertical + Real{0.5});
  auto const n = static_cast<std::int64_t>(nr), m = static_cast<std::int64_t>(mr);
  Complex const reduced(u.real() - nr, u.imag() - mr * t);
  auto result = kind == 1 ? theta1_reduced(reduced, t, t < Real{1}) : even_reduced(kind, reduced, t);
  // Quasi-periodicity, differentiated before reconstructing the common scale.
  Complex const slope(0, -Real{2} * pi * mr);
  result.derivative[2] += Real{2} * slope * result.derivative[1] + slope * slope * result.derivative[0];
  result.derivative[1] += slope * result.derivative[0];
  auto const parity = kind == 1 ? n + m : kind == 2 ? n : kind == 4 ? m : 0;
  Complex const phase = Real(parity % 2 ? -1 : 1) * std::exp(Complex(0, -Real{2} * pi * mr * reduced.real()));
  for (auto& d : result.derivative)
    d *= phase;
  result.log_scale += pi * mr * (mr * t + Real{2} * reduced.imag());
  return result;
}
} // namespace theta_detail

/// Rectangular nome only: tau=i*t, t finite and positive. No double conversion.
/// Throws for invalid inputs or unrepresentable scaled derivatives/reduction.
template <uni20::Real Real> ThetaJet<Real> elliptic_theta(unsigned kind, std::complex<Real> u, Real t)
{
  using Complex = std::complex<Real>;
  if (kind < 1 || kind > 4 || !theta_detail::finite(u) || !uni20::isfinite(t) || !(t > Real{0}))
    throw std::invalid_argument("theta requires kind 1..4, finite u, and finite t>0");
  Real const pi = Real{4} * std::atan(Real{1});
  if (kind > 1 && t < Real{1} && u != Complex{} &&
      std::abs(u) < std::sqrt(uni20::numeric_limits<Real>::epsilon()) * t / Real{8})
  {
    // Avoid losing tiny arguments when adding a half-period to an even function.
    auto result = elliptic_theta(kind, Complex{}, t);
    result.derivative[0] += result.derivative[2] * u * u / Real{2};
    result.derivative[1] = result.derivative[2] * u;
    return result;
  }
  // In the Gaussian regime, translate the other kinds through theta1.
  Complex shifted = u;
  if (t < Real{1} && (kind == 2 || kind == 3)) shifted += Real{0.5};
  if (t < Real{1} && (kind == 3 || kind == 4)) shifted += Complex(0, t / Real{2});
  auto result = theta_detail::theta1(shifted, t, t < Real{1} ? 1u : kind);
  if (t < Real{1} && (kind == 3 || kind == 4))
  {
    Complex const slope(0, pi);
    result.derivative[2] += Real{2} * slope * result.derivative[1] + slope * slope * result.derivative[0];
    result.derivative[1] += slope * result.derivative[0];
    Complex const phase = (kind == 4 ? Complex(0, -1) : Complex(1, 0)) * std::exp(Complex(0, pi * u.real()));
    for (auto& d : result.derivative)
      d *= phase;
    result.log_scale -= pi * (t / Real{4} + u.imag());
  }
  if (kind > 1 && u == Complex{}) result.derivative[1] = Complex{};
  if (!uni20::isfinite(result.log_scale)) throw std::overflow_error("theta logarithmic scale exceeds precision");
  for (auto d : result.derivative)
    if (!theta_detail::finite(d)) throw std::overflow_error("theta scaled derivatives exceed precision");
  return result;
}
/// Ratio of two scaled theta derivatives, reconstructing only their quotient.
template <uni20::Real Real>
std::complex<Real> theta_ratio(ThetaJet<Real> const& a, unsigned da, ThetaJet<Real> const& b, unsigned db)
{
  if (da > 2 || db > 2) throw std::invalid_argument("theta derivative order must be 0..2");
  Real const numerator = std::abs(a.derivative[da]), denominator = std::abs(b.derivative[db]);
  if (!(denominator > Real{0})) throw std::domain_error("theta ratio has a zero denominator");
  if (numerator == Real{0}) return {};
  Real const shift = a.log_scale - b.log_scale;
  if (std::abs(shift) <= Real{1})
  {
    auto const direct = (a.derivative[da] / b.derivative[db]) * std::exp(shift);
    if (theta_detail::finite(direct) && std::abs(direct) > Real{0}) return direct;
  }
  Real const magnitude = std::exp(a.log_scale - b.log_scale + std::log(numerator) - std::log(denominator));
  if (!uni20::isfinite(magnitude) || !(magnitude > Real{0}))
    throw std::overflow_error("theta ratio exceeds scalar range");
  return magnitude * (a.derivative[da] / numerator) * std::conj(b.derivative[db] / denominator);
}
} // namespace bethe::detail
