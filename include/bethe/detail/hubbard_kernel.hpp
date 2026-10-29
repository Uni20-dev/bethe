// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/solver.hpp>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <uni20/core/complex_log_gamma.hpp>
#include <utility>

namespace bethe::detail
{
// Immutable coefficients and explicit arithmetic precision for the zero-field
// Hubbard kernel. Here u=U/4; meshes and concurrent queries share one context.
template <uni20::Real Real> class HubbardKernel {
  public:
    explicit HubbardKernel(Real u)
      requires(!uni20::has_runtime_precision_v<Real>)
        : HubbardKernel(std::move(u), uni20::native_precision<Real>())
    {}

    HubbardKernel(Real u, uni20::Precision precision, std::optional<std::size_t> max_terms = {})
        : u_(std::move(u)), zero_(checked_zero(precision)), one_(uni20::scalar_at_precision<Real>(1, precision)),
          pi_(Real{4} * uni20::math::atan(one_)), epsilon_(uni20::numeric_limits<Real>::epsilon(precision)),
          max_terms_(max_terms.value_or(2 * uni20::numeric_limits<Real>::digits(precision) + 32)), gamma_(precision)
    {
      require_same_precision(u_, zero_);
      if (!uni20::isfinite(u_) || !(u_ > zero_)) throw std::invalid_argument("Hubbard kernel requires finite u=U/4>0");
    }

    Real const& u() const { return u_; }

    // R(x) = integral_0^infinity cos(w*x)/(1+exp(2*u*w)) dw/pi.
    // Euler's beta-integral expansion avoids subtracting nearly equal digammas.
    Real value(Real const& x) const
    {
      this->validate(x);
      // In a distant tail, include 1/u in the terms before accumulation.
      // Otherwise their real parts can underflow even when R(x) is representable.
      bool const scaled = uni20::math::abs(x) * epsilon_ > u_;
      Real const scale = scaled ? u_ : one_;
      uni20::complex<Real> const a{scale, scaled ? x / Real{2} : this->scaled_argument(x, Real{2})};
      uni20::complex<Real> term = (one_ / Real{2}) / a;
      CompensatedSum<Real> sum(zero_);
      for (std::size_t j = 0; j < max_terms_; ++j)
      {
        sum.add(term.real());
        Real const numerator = scale * Real(j + 1);
        term *= (numerator / (a + numerator)) / Real{2};
        // Successive moduli shrink by a factor strictly below 1/2.
        // The exact-arithmetic remainder is therefore at most twice the next
        // modulus. This controls truncation, not accumulated rounding error.
        if (sum.value() > zero_ &&
            uni20::math::abs(term.real()) + uni20::math::abs(term.imag()) <= epsilon_ * sum.value() / Real{16})
        {
          Real const result = (scaled ? sum.value() : sum.value() / u_) / (Real{2} * pi_);
          if (!uni20::isfinite(result)) throw std::overflow_error("Hubbard kernel exceeds arithmetic range");
          return result;
        }
      }
      throw std::runtime_error("Hubbard kernel exhausted its Euler-series term budget");
    }

    // Integral_0^x R(t)dt = Im log[Gamma(1+i*x/(4u))/Gamma(1/2+i*x/(4u))]/pi.
    // Uni20 evaluates the paired ratio directly, with precision-dependent
    // coefficients and recurrence. Subtracting two large log-Gamma phases loses
    // the small difference in the rapidity tail.
    Real primitive(Real const& x) const
    {
      this->validate(x);
      if (x == zero_) return zero_;
      Real const b = this->scaled_argument(uni20::math::abs(x), Real{4});
      // The omitted correction is O(1/b), below working roundoff here. This
      // also handles an overflowing ratio of otherwise finite inputs.
      if (!uni20::isfinite(b) || b > one_ / epsilon_) return x < zero_ ? -one_ / Real{4} : one_ / Real{4};
      Real const half = one_ / Real{2};
      auto const ratio = gamma_.log_gamma_ratio_right(uni20::complex<Real>{half, b}, half);
      Real const phase = ratio.value.imag() / pi_;
      return x < zero_ ? -phase : phase;
    }

  private:
    static Real checked_zero(uni20::Precision precision)
    {
      if (precision.is_exact()) throw std::invalid_argument("Hubbard kernel requires finite working precision");
      return uni20::scalar_at_precision<Real>(0, precision);
    }

    void validate(Real const& x) const
    {
      require_same_precision(x, zero_);
      if (!uni20::isfinite(x)) throw std::invalid_argument("Hubbard kernel requires finite rapidity");
    }

    Real scaled_argument(Real const& x, Real factor) const
    {
      Real const denominator = factor * u_;
      return uni20::isfinite(denominator) ? x / denominator : (x / u_) / factor;
    }

    Real u_, zero_, one_, pi_, epsilon_;
    std::size_t max_terms_;
    uni20::log_gamma_context<Real> gamma_;
};

// Native convenience functions. Repeated evaluations should retain a context;
// runtime-precision callers always construct it with an explicit precision.
template <uni20::Real Real>
  requires(!uni20::has_runtime_precision_v<Real>)
Real hubbard_r(Real x, Real u)
{
  return HubbardKernel<Real>(u).value(x);
}
template <uni20::Real Real>
  requires(!uni20::has_runtime_precision_v<Real>)
Real hubbard_r_primitive(Real x, Real u)
{
  return HubbardKernel<Real>(u).primitive(x);
}
} // namespace bethe::detail
