// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz.hpp>
#include <complex>
#include <uni20/tensor/tensor.hpp>

namespace bethe::xxz::detail
{
// A directional derivative of complex polynomial arithmetic with respect to
// one real coefficient. No finite differences or complex transcendental
// functions are needed, including for fp128.
template <uni20::Real Real> struct PolynomialJet
{
    using Complex = std::complex<Real>;
    Complex value{}, tangent{};
    friend PolynomialJet operator+(PolynomialJet a, PolynomialJet b)
    {
      return {a.value + b.value, a.tangent + b.tangent};
    }
    friend PolynomialJet operator-(PolynomialJet a, PolynomialJet b)
    {
      return {a.value - b.value, a.tangent - b.tangent};
    }
    friend PolynomialJet operator*(PolynomialJet a, PolynomialJet b)
    {
      return {a.value * b.value, a.tangent * b.value + a.value * b.tangent};
    }
};

/// Coefficient equations for periodic XXZ, -1 < Delta <= 0.
/// Q(z)=z^M+sum_(k=0)^(M-1) c[k]*z^k has REAL coefficients, so it can
/// describe both real roots and complex conjugate pairs without changing
/// variables at collisions. This is an equation/observable building block,
/// NOT a ground-state solver or a certification of physical Bethe states.
template <uni20::Real Real> class PolynomialBetheSystem {
  public:
    using Complex = std::complex<Real>;
    using Jet = PolynomialJet<Real>;
    using Polynomial = std::vector<Jet>;

    /// Coefficients may instead describe Q(x), with physical z=center+scale*x.
    /// The default preserves the original physical-z coefficient convention.
    PolynomialBetheSystem(std::size_t sites, std::size_t roots, Real center = Real{0}, Real scale = Real{1})
        : sites(sites), order(roots), center(center), coordinate_scale(scale)
    {
      checked_sites(sites);
      if (roots > sites / 2) throw std::invalid_argument("XXZ polynomial degree requires M <= N/2");
      auto const elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Jet);
      if (roots > elements / 2) throw std::length_error("XXZ polynomial workspace is too large");
      if (!uni20::isfinite(center) || !uni20::isfinite(scale) || scale <= Real{0})
        throw std::invalid_argument("XXZ polynomial coordinates require finite center and positive finite scale");
    }

    struct Evaluation
    {
        /// Raw real or imaginary coefficients of (1+i*z)^N K_-(z) modulo Q.
        std::vector<Real> residual;
        /// max|residual| / max|complex coefficient|; not a rootwise phase error.
        Real norm = Real{0};
        Real scale = Real{1};
    };

    /// Jacobian entries differentiate the RAW residual, not the normalization.
    Evaluation evaluate(std::span<Real const> coefficients, Real delta,
                        uni20::DenseMatrix<Real>* jacobian = nullptr) const
    {
      validate(coefficients);
      validate_delta(delta);
      if (jacobian && (std::size_t(jacobian->extent(0)) != order || std::size_t(jacobian->extent(1)) != order))
        throw std::invalid_argument("XXZ polynomial Jacobian has wrong shape");
      Evaluation out{.residual = std::vector<Real>(order)};
      if (!order) return out;
      bool const imaginary = (sites - order - 1) % 2 == 0;
      auto component = [&](Complex z) { return imaginary ? z.imag() : z.real(); };
      Real scale = Real{0}, maximum = Real{0};
      // One evaluation for values only, or M directional evaluations for the
      // full analytic Jacobian. The derivative includes reduction modulo Q:
      // the quotient algebra itself depends on the unknown coefficients.
      for (std::size_t j = 0; j < (jacobian ? order : 1); ++j)
      {
        auto const f = direction(coefficients, delta, jacobian ? j : order);
        for (std::size_t k = 0; k < order; ++k)
        {
          if (!finite(f[k].value) || (jacobian && !finite(f[k].tangent)))
            throw std::overflow_error("nonfinite XXZ polynomial equation or derivative");
          if (j == 0)
          {
            out.residual[k] = component(f[k].value);
            maximum = std::max(maximum, std::abs(out.residual[k]));
            scale = std::max(scale, std::abs(f[k].value));
          }
          if (jacobian) (*jacobian)[k, j] = component(f[k].tangent);
        }
      }
      // Zero can result from a singular/exact-string polynomial, not evidence
      // that a normalized physical-state residual is zero.
      if (!(scale > Real{0}) || !uni20::isfinite(scale))
        throw std::runtime_error("XXZ polynomial normalization vanished or overflowed");
      out.scale = scale;
      out.norm = maximum / scale;
      return out;
    }

    /// Energy and momentum are rational functions of Q at i, so roots need
    /// not be extracted (an ill-conditioned operation at multiple roots).
    Real energy(std::span<Real const> coefficients, Real delta) const
    {
      validate(coefficients);
      validate_delta(delta);
      auto const [q, derivative] = at_i(coefficients);
      Real const scale = std::max(std::abs(q.real()), std::abs(q.imag()));
      Real const result = Real(sites) * delta / Real{4} + Real(order) * (Real{1} - delta) +
                          Real{2} * ((derivative / scale) / (q / scale)).imag();
      if (!uni20::isfinite(result)) throw std::overflow_error("nonfinite XXZ polynomial energy");
      return result;
    }

    Complex momentum_phase(std::span<Real const> coefficients) const
    {
      validate(coefficients);
      auto const [q, derivative] = at_i(coefficients);
      Real const scale = std::max(std::abs(q.real()), std::abs(q.imag()));
      Complex const normalized = q / scale;
      return std::conj(normalized) / normalized;
    }

    /// A NECESSARY physical-state check, not a sufficient admissibility test.
    /// A polynomial remainder can be tiny even for nonquantized momentum.
    Real momentum_defect(std::span<Real const> coefficients) const
    {
      Complex power{1}, phase = momentum_phase(coefficients);
      for (std::size_t n = sites; n; n >>= 1)
      {
        if (n & 1) power *= phase;
        phase *= phase;
      }
      return std::abs(power - Complex{1});
    }

    std::size_t sites, order;
    Real center, coordinate_scale;

  private:
    static bool finite(Complex z) { return uni20::isfinite(z.real()) && uni20::isfinite(z.imag()); }
    static void validate_delta(Real delta)
    {
      if (!uni20::isfinite(delta) || delta <= -Real{1} || delta > Real{0})
        throw std::invalid_argument("XXZ polynomial equations require -1 < Delta <= 0");
    }
    void validate(std::span<Real const> c) const
    {
      if (c.size() != order) throw std::invalid_argument("XXZ polynomial coefficient count differs from degree");
      for (Real x : c)
        if (!uni20::isfinite(x)) throw std::invalid_argument("nonfinite XXZ polynomial coefficient");
    }

    std::pair<Complex, Complex> at_i(std::span<Real const> c) const
    {
      if (center != Real{0} || coordinate_scale != Real{1}) return shifted_at_i(c);
      // Powers of i cycle exactly. Compensated real and imaginary sums avoid
      // repeated Horner cancellation when evaluating the observables.
      bethe::detail::CompensatedSum<Real> qr, qi, dr, di;
      for (std::size_t k = 0; k <= order; ++k)
      {
        Real const value = k == order ? Real{1} : c[k];
        Real const signed_value = k % 4 < 2 ? value : -value;
        (k % 2 ? qi : qr).add(signed_value);
        if (k)
        {
          Real const d = Real(k) * value;
          ((k - 1) % 2 ? di : dr).add((k - 1) % 4 < 2 ? d : -d);
        }
      }
      Complex const q{qr.value(), qi.value()}, derivative{dr.value(), di.value()};
      if (!finite(q) || !finite(derivative) || q == Complex{})
        throw std::runtime_error("XXZ polynomial observable has a pole or nonfinite coefficients");
      return {q, derivative};
    }

    std::pair<Complex, Complex> shifted_at_i(std::span<Real const> c) const
    {
      // R(z)=scale^M Q((z-center)/scale). Evaluate its reciprocal polynomial
      // W(t)=1+c[M-1]*t+...+c[0]*t^M, t=scale/(i-center), without large powers
      // of (i-center)/scale. The returned pair is R(i),R'(i) divided by the
      // SAME positive scale |i-center|^M, preserving both ratio and phase.
      Complex const denominator{-center, Real{1}}, t = coordinate_scale / denominator;
      Complex unit = denominator / std::abs(denominator), phase{1}, power{1};
      for (std::size_t n = order; n; n >>= 1)
      {
        if (n & 1) phase *= unit;
        unit *= unit;
      }
      bethe::detail::CompensatedSum<Real> wr, wi, dr, di;
      wr.add(Real{1});
      for (std::size_t j = 1; j <= order; ++j)
      {
        power *= t;
        Complex const value = c[order - j] * power;
        wr.add(value.real());
        wi.add(value.imag());
        dr.add(Real(j) * value.real());
        di.add(Real(j) * value.imag());
      }
      Complex const w{wr.value(), wi.value()}, dw{dr.value(), di.value()};
      Complex const q = phase * w, derivative = phase * ((Real(order) * w - dw) / denominator);
      if (!finite(q) || !finite(derivative) || q == Complex{})
        throw std::runtime_error("XXZ polynomial observable has a pole or nonfinite coefficients");
      return {q, derivative};
    }

    // Multiply by alpha+beta*z and reduce immediately. Keeping all
    // polynomials below degree M avoids forming a degree N+2M expression.
    Polynomial linear(Polynomial const& p, Jet alpha, Jet beta, Polynomial const& q) const
    {
      Polynomial out(order);
      Jet const highest = beta * p.back();
      for (std::size_t k = 0; k < order; ++k)
        out[k] = alpha * p[k] + (k ? beta * p[k - 1] : Jet{}) - highest * q[k];
      return out;
    }

    Polynomial multiply(Polynomial const& a, Polynomial const& b, Polynomial const& q) const
    {
      Polynomial out(2 * order - 1);
      for (std::size_t i = 0; i < order; ++i)
        for (std::size_t j = 0; j < order; ++j)
          out[i + j] = out[i + j] + a[i] * b[j];
      for (std::size_t j = out.size(); j-- > order;)
        for (std::size_t k = 0; k < order; ++k)
          out[j - order + k] = out[j - order + k] - out[j] * q[k];
      out.resize(order);
      return out;
    }

    Polynomial direction(std::span<Real const> c, Real delta, std::size_t column) const
    {
      Polynomial q(order), k(order), z(order), power(order);
      for (std::size_t j = 0; j < order; ++j)
        q[j] = {Complex{c[j]}, Complex{j == column ? Real{1} : Real{0}}};
      z[0] = power[0] = Jet{Complex{1}};
      Jet const zero{}, one{Complex{1}};
      // A=1+Delta-i*Delta*z, B=(1-Delta)*z-i*Delta.
      // K=(B^M Q(A/B)-B^M Q(z))/(A-B*z) is a POLYNOMIAL.
      // Divided Horner recurrence removes the self-scattering factor
      // A-B*z=1+Delta-(1-Delta)*z^2 before imposing divisibility by Q.
      // Leaving that factor in would introduce spurious roots at +/-s.
      Real const p = Real{1} + delta, v = Real{1} - delta;
      Jet const a0{Complex{p - v * center * center}};
      Jet const a1{Complex{-v * center * coordinate_scale, -delta * coordinate_scale}};
      Jet const b0{Complex{v * center * coordinate_scale, -delta * coordinate_scale}};
      Jet const b1{Complex{v * coordinate_scale * coordinate_scale}};
      for (std::size_t r = 0; r < order; ++r)
      {
        auto const product = multiply(power, z, q);
        k = linear(k, a0, a1, q);
        for (std::size_t j = 0; j < order; ++j)
          k[j] = k[j] + product[j];
        z = linear(z, zero, one, q);
        z[0] = z[0] + q[order - r - 1];
        power = linear(power, b0, b1, q);
      }
      for (std::size_t j = 0; j < sites; ++j)
        k = linear(k, Jet{Complex{Real{1}, center}}, Jet{Complex{Real{0}, coordinate_scale}}, q);
      return k;
    }
};
} // namespace bethe::xxz::detail
