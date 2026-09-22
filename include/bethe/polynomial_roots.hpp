// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <algorithm>
#include <bethe/solver.hpp>
#include <cmath>
#include <complex>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace bethe::detail
{
// Error-free real transforms under round-to-nearest, without fast-math
// reassociation and away from under/overflow. Keep every operation in Real;
// no double/long-double fallback or extra scalar backend is introduced.
template <uni20::Real Real> std::pair<Real, Real> root_two_sum(Real a, Real b)
{
  Real const sum = a + b, part = sum - a;
  return {sum, (a - (sum - part)) + (b - part)};
}

template <uni20::Real Real> std::pair<Real, Real> root_two_product(Real a, Real b)
{
  Real const splitter = std::ldexp(Real{1}, (uni20::numeric_limits<Real>::digits + 1) / 2) + Real{1};
  Real const ca = splitter * a, cb = splitter * b;
  Real const ah = ca - (ca - a), al = a - ah, bh = cb - (cb - b), bl = b - bh;
  Real const product = a * b;
  Real const error1 = product - ah * bh;
  Real const error2 = error1 - al * bh;
  Real const error3 = error2 - ah * bl;
  return {product, al * bl - error3};
}

// Rounded complex a*b+c and its first-order rounding correction. This is
// the elementary building block of compensated Horner evaluation.
template <uni20::Real Real>
std::pair<std::complex<Real>, std::complex<Real>> root_multiply_add(std::complex<Real> a, std::complex<Real> b,
                                                                    std::complex<Real> c)
{
  auto const [rr, err] = root_two_product(a.real(), b.real());
  auto const [ii, eii] = root_two_product(-a.imag(), b.imag());
  auto const [ri, eri] = root_two_product(a.real(), b.imag());
  auto const [ir, eir] = root_two_product(a.imag(), b.real());
  auto const [r, er] = root_two_sum(rr, ii);
  auto const [i, ei] = root_two_sum(ri, ir);
  auto const [real, erc] = root_two_sum(r, c.real());
  auto const [imag, eic] = root_two_sum(i, c.imag());
  return {{real, imag}, {((err + eii) + er) + erc, ((eri + eir) + ei) + eic}};
}

enum class PolynomialRootStatus
{
  resolved,
  iteration_limit,
  unresolved_cluster,
  stalled,
  nonfinite_or_unrepresentable
};

template <uni20::Real Real> struct PolynomialRootOptions
{
    Real tolerance = Real{128} * uni20::numeric_limits<Real>::epsilon();
    std::size_t max_iterations = 200;
    std::size_t max_degree = 64;
};

/// Numerical diagnostics, not interval enclosures or a proof of simplicity.
/// Roots use the original polynomial coordinate; their order is unspecified.
template <uni20::Real Real> struct PolynomialRoots
{
    std::vector<std::complex<Real>> roots;
    PolynomialRootStatus status = PolynomialRootStatus::nonfinite_or_unrepresentable;
    std::size_t iterations = 0; // Attempted full-root sweeps; no hidden retries.
    Real residual_norm = uni20::numeric_limits<Real>::infinity();
    Real reconstruction_error = uni20::numeric_limits<Real>::infinity();
    Real max_root_uncertainty = uni20::numeric_limits<Real>::infinity();
    Real separation_ratio = Real{0}; // min |x_i-x_j|/(estimated radii_i+radii_j)
};

/// Recover the roots of monic Q(x)=x^n+sum_(k<n) c[k]*x^k, with REAL c.
/// Native-precision Ehrlich-Aberth iteration (Gauss-Seidel sweeps), no explicit
/// deflation. A resolved result has small backward errors and numerically
/// separated roots; it does NOT promise tolerance-sized forward errors.
/// See docs/polynomial-roots.md for scaling, diagnostics, and limitations.
template <uni20::Real Real>
PolynomialRoots<Real> recover_polynomial_roots(std::span<Real const> coefficients,
                                               PolynomialRootOptions<Real> options = {},
                                               std::span<std::complex<Real> const> initial = {})
{
  using C = std::complex<Real>;
  Real const infinity = uni20::numeric_limits<Real>::infinity();
  auto finite = [](C z) { return uni20::isfinite(z.real()) && uni20::isfinite(z.imag()); };
  auto const n = coefficients.size();
  if (!uni20::isfinite(options.tolerance) || options.tolerance <= Real{0})
    throw std::invalid_argument("polynomial root tolerance must be finite and positive");
  if (n > options.max_degree) throw std::length_error("polynomial root degree exceeds the work budget");
  if (!initial.empty() && initial.size() != n)
    throw std::invalid_argument("polynomial root seed count differs from degree");
  for (Real c : coefficients)
    if (!uni20::isfinite(c)) throw std::invalid_argument("nonfinite polynomial coefficient");
  for (C x : initial)
    if (!finite(x)) throw std::invalid_argument("nonfinite polynomial root seed");
  PolynomialRoots<Real> out;
  // Use a power-of-two coordinate scale, so scaling does not perturb the
  // represented coefficients. Reject subnormal scaled coefficients rather
  // than silently lose input bits. The usual bound is <=1; if rounding the
  // scale upward overflows, the last finite power of two is retained.
  Real scale{1};
  for (std::size_t k = 0; k < n; ++k)
    scale = std::max(scale, std::pow(std::abs(coefficients[k]), Real{1} / Real(n - k)));
  if (!uni20::isfinite(scale)) return out;
  int exponent = 0;
  Real const mantissa = std::frexp(scale, &exponent);
  scale = std::ldexp(Real{1}, exponent - 1);
  if (mantissa > Real{1} / Real{2} && uni20::isfinite(Real{2} * scale)) scale *= Real{2};
  std::vector<Real> c(coefficients.begin(), coefficients.end());
  for (std::size_t k = 0; k < n; ++k)
  {
    for (std::size_t j = k; j < n; ++j)
      c[k] /= scale;
    if (!uni20::isfinite(c[k]) || (coefficients[k] != Real{0} && std::abs(c[k]) < uni20::numeric_limits<Real>::min()))
      return out;
  }
  std::vector<C> roots(n);
  Real const pi = Real{4} * std::atan(Real{1});
  for (std::size_t j = 0; j < n; ++j)
  {
    if (!initial.empty())
      roots[j] = initial[j] / scale;
    else if (n == 1)
      roots[j] = C{-c[0]};
    else
    {
      Real const angle = Real{2} * pi * (Real(j) + Real{1} / Real{3}) / Real(n);
      roots[j] = Real{4} * C{std::cos(angle), std::sin(angle)};
    }
    if (!finite(roots[j]) || (!initial.empty() && initial[j] != C{} && roots[j] == C{})) return out;
  }
  struct Evaluation
  {
      C value{1}, derivative{};
      Real magnitude{1}, derivative_magnitude{};
  };
  auto evaluate = [&](C x) {
    Evaluation e;
    C correction{};
    Real const radius = std::abs(x);
    for (std::size_t k = n; k-- > 0;)
    {
      e.derivative = e.derivative * x + e.value;
      auto const [value, error] = root_multiply_add(e.value, x, C{c[k]});
      correction = correction * x + error;
      e.value = value;
      e.derivative_magnitude = e.derivative_magnitude * radius + e.magnitude;
      e.magnitude = e.magnitude * radius + std::abs(c[k]);
    }
    e.value += correction;
    return e;
  };
  bool diagnostics_current = false;
  auto finish = [&](PolynomialRootStatus status) {
    out.status = status;
    out.roots = roots;
    for (C& x : out.roots)
    {
      x *= scale;
      if (!finite(x))
      {
        out.status = PolynomialRootStatus::nonfinite_or_unrepresentable;
        diagnostics_current = false;
      }
    }
    if (!diagnostics_current)
    {
      out.residual_norm = out.reconstruction_error = out.max_root_uncertainty = infinity;
      out.separation_ratio = Real{0};
    }
    return out;
  };
  Real const roundoff = Real{16} * Real(n) * uni20::numeric_limits<Real>::epsilon();
  for (;;)
  {
    diagnostics_current = false;
    out.residual_norm = Real{0};
    out.max_root_uncertainty = Real{0};
    out.separation_ratio = infinity;
    std::vector<Real> radii(n);
    for (std::size_t j = 0; j < n; ++j)
    {
      auto const e = evaluate(roots[j]);
      if (!finite(e.value) || !finite(e.derivative) || !uni20::isfinite(e.magnitude) ||
          !uni20::isfinite(e.derivative_magnitude))
        return finish(PolynomialRootStatus::nonfinite_or_unrepresentable);
      Real const residual =
          e.magnitude > Real{0} ? std::abs(e.value) / e.magnitude : (e.value == C{} ? Real{0} : infinity);
      out.residual_norm = std::max(out.residual_norm, residual);
      Real const derivative = std::abs(e.derivative) - roundoff * e.derivative_magnitude;
      radii[j] = derivative > Real{0} ? (std::abs(e.value) + roundoff * e.magnitude) / derivative : infinity;
      out.max_root_uncertainty = std::max(out.max_root_uncertainty, radii[j] * scale);
    }
    bool collision = false;
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = i + 1; j < n; ++j)
      {
        Real const distance = std::abs(roots[i] - roots[j]), radius = radii[i] + radii[j];
        collision = collision || distance == Real{0};
        Real const ratio = distance == Real{0} ? Real{0} : (radius > Real{0} ? distance / radius : infinity);
        out.separation_ratio = std::min(out.separation_ratio, ratio);
      }
    // Leja-style ordering distributes roots before multiplying factors;
    // neighboring roots on a circle are a poorly conditioned product order.
    std::vector<C> ordered = roots;
    std::vector<Real> distances(n, Real{1});
    for (std::size_t j = 0; j < n; ++j)
    {
      auto best = j;
      for (std::size_t k = j + 1; k < n; ++k)
        if (distances[k] > distances[best]) best = k;
      std::swap(ordered[j], ordered[best]);
      std::swap(distances[j], distances[best]);
      for (std::size_t k = j + 1; k < n; ++k)
        distances[k] *= std::abs(ordered[k] - ordered[j]);
    }
    std::vector<C> reconstructed{C{1}}, corrections{C{0}};
    for (C x : ordered)
    {
      reconstructed.push_back(C{});
      corrections.push_back(C{});
      for (std::size_t k = reconstructed.size(); k-- > 0;)
      {
        auto const [value, error] = root_multiply_add(reconstructed[k], -x, k ? reconstructed[k - 1] : C{});
        corrections[k] = (k ? corrections[k - 1] : C{}) - x * corrections[k] + error;
        reconstructed[k] = value;
        if (!finite(value) || !finite(corrections[k]))
          return finish(PolynomialRootStatus::nonfinite_or_unrepresentable);
      }
    }
    out.reconstruction_error = Real{0};
    for (std::size_t k = 0; k < n; ++k)
      out.reconstruction_error =
          std::max(out.reconstruction_error,
                   std::abs((reconstructed[k] - c[k]) + corrections[k]) / std::max(Real{1}, std::abs(c[k])));
    diagnostics_current = true;
    bool const backward = out.residual_norm <= options.tolerance && out.reconstruction_error <= options.tolerance;
    if (backward)
      return finish(out.separation_ratio > Real{8} && uni20::isfinite(out.max_root_uncertainty)
                        ? PolynomialRootStatus::resolved
                        : PolynomialRootStatus::unresolved_cluster);
    if (out.iterations == options.max_iterations) return finish(PolynomialRootStatus::iteration_limit);
    if (collision) return finish(PolynomialRootStatus::unresolved_cluster);
    ++out.iterations;
    bool changed = false;
    for (std::size_t j = 0; j < n; ++j)
    {
      auto const e = evaluate(roots[j]);
      CompensatedSum<Real> real, imag;
      for (std::size_t k = 0; k < n; ++k)
        if (j != k)
        {
          if (roots[j] == roots[k]) return finish(PolynomialRootStatus::unresolved_cluster);
          C const reciprocal = Real{1} / (roots[j] - roots[k]);
          real.add(reciprocal.real());
          imag.add(reciprocal.imag());
        }
      // Algebraically equivalent to Newton/(1-Newton*sum), also defined
      // at Q'=0 when the Ehrlich-Aberth denominator is nonzero.
      C const denominator = e.derivative - e.value * C{real.value(), imag.value()};
      if (!finite(denominator)) return finish(PolynomialRootStatus::nonfinite_or_unrepresentable);
      if (denominator == C{}) return finish(PolynomialRootStatus::stalled);
      C const next = roots[j] - e.value / denominator;
      if (!finite(next)) return finish(PolynomialRootStatus::nonfinite_or_unrepresentable);
      if (next != roots[j])
      {
        changed = true;
        diagnostics_current = false;
      }
      roots[j] = next;
    }
    if (!changed) return finish(PolynomialRootStatus::stalled);
  }
}
} // namespace bethe::detail
