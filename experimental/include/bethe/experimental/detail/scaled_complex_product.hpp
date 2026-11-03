// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <algorithm>
#include <complex>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <uni20/core/math.hpp>

namespace bethe::experimental::detail
{
/// Product arithmetic with a separate binary exponent. Never reconstruct an
/// enormous q-product before its prefactor cancels it. No logarithmic branch
/// choices and no conversion to double; exact powers of two do the scaling.
template <uni20::Real R> class ScaledComplexProduct {
  public:
    using C = uni20::complex<R>;
    explicit ScaledComplexProduct(C value = C{R{1}}) : mantissa_(value) { normalize(); }
    ScaledComplexProduct& operator*=(ScaledComplexProduct const& b)
    {
      mantissa_ *= b.mantissa_;
      add_exponent(b.exponent_);
      normalize();
      return *this;
    }
    ScaledComplexProduct& operator/=(ScaledComplexProduct const& b)
    {
      if (b.mantissa_ == C{R{0}}) throw std::domain_error("zero scaled-product denominator");
      mantissa_ /= b.mantissa_;
      add_exponent(-b.exponent_);
      normalize();
      return *this;
    }
    ScaledComplexProduct& operator*=(C b) { return *this *= ScaledComplexProduct(b); }
    ScaledComplexProduct& operator/=(C b) { return *this /= ScaledComplexProduct(b); }
    C value() const
    {
      C const result{uni20::math::ldexp(mantissa_.real(), exponent_), uni20::math::ldexp(mantissa_.imag(), exponent_)};
      if (!uni20::isfinite(result.real()) || !uni20::isfinite(result.imag()) ||
          (result == C{R{0}} && mantissa_ != C{R{0}}))
        throw std::overflow_error("scaled amplitude outside scalar range; increase precision");
      return result;
    }

  private:
    void add_exponent(std::int64_t increment)
    {
      constexpr auto limit = std::numeric_limits<std::int64_t>::max() / 4;
      if (increment > limit || increment < -limit || (increment > 0 && exponent_ > limit - increment) ||
          (increment < 0 && exponent_ < -limit - increment))
        throw std::overflow_error("scaled-product exponent budget exceeded");
      exponent_ += increment;
    }
    void normalize()
    {
      R const size = std::max(uni20::math::abs(mantissa_.real()), uni20::math::abs(mantissa_.imag()));
      if (!uni20::isfinite(size)) throw std::overflow_error("invalid scaled-product factor");
      if (size == R{0})
      {
        exponent_ = 0;
        return;
      }
      auto const shift = uni20::math::frexp(size).exponent;
      add_exponent(shift);
      mantissa_ = {uni20::math::ldexp(mantissa_.real(), -shift), uni20::math::ldexp(mantissa_.imag(), -shift)};
    }
    C mantissa_;
    std::int64_t exponent_ = 0;
};
} // namespace bethe::experimental::detail
