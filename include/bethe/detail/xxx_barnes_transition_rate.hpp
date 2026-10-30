// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/spinon_barnes.hpp>

namespace bethe::heisenberg::detail
{
/// Independent Barnes-G representation of log(exp(-I(rho))). Does not replace
/// the existing quadrature path. Errors refer to the logarithm, not the rate.
template <uni20::Real Real> class XXXBarnesTransitionRate {
  public:
    using Complex = uni20::complex<Real>;
    template <class... Args> explicit XXXBarnesTransitionRate(Args... args) : functions_(args...)
    {
      constant_ = functions_.log_abs_barnes_g_right(Complex(Real{1} / Real{2}, Real{0}));
      one_ = (constant_.value - constant_.value) + Real{1};
      using std::atan;
      pi_ = Real{4} * atan(one_);
    }

    bethe::detail::SpinonFunctionEstimate<Real> operator()(Real rho) const
    {
      if (!uni20::isfinite(rho)) throw std::invalid_argument("transition rate requires finite rho");
      // Exact constants acquire the context precision; mismatched approximate
      // inputs must be explicitly converted by the caller.
      rho = uni20::math::abs(rho) * one_;
      Real zero = one_ - one_;
      if (rho == zero) return {uni20::math::log(zero), zero, zero, 0, 0};
      auto a = functions_.log_abs_barnes_g_right(Complex(one_, rho));
      auto b = functions_.log_abs_barnes_g_right(Complex(one_ / Real{2}, rho));
      Real x = pi_ * rho, hyperbolic;
      Real eps = uni20::numeric_limits<Real>::epsilon(one_);
      std::size_t terms = 0;
      if (x <= one_)
      {
        // sinh(x)/x and cosh(x), without subtracting almost equal exponentials.
        Real s = one_, c = one_, st = one_, ct = one_, square = x * x;
        do
        {
          ++terms;
          st *= square / (Real(2 * terms) * Real(2 * terms + 1));
          ct *= square / (Real(2 * terms - 1) * Real(2 * terms));
          s += st;
          c += ct;
          if (terms == 4096)
            throw bethe::detail::SpinonFunctionError(bethe::detail::SpinonFunctionFailure::budget,
                                                     "hyperbolic series budget exhausted");
        }
        while (st + ct > eps);
        hyperbolic = uni20::math::log(x) + uni20::math::log(s) + uni20::math::log(c);
      }
      else
        hyperbolic =
            Real{2} * x - uni20::math::log(Real{4} * one_) + uni20::math::log(one_ - uni20::math::exp(-Real{4} * x));
      Real elementary = uni20::math::log(rho) + hyperbolic - Real{2} * uni20::math::log(pi_);
      Real value = -Real{4} * constant_.value + elementary + Real{4} * (a.value - b.value);
      Real truncation = Real{4} * (constant_.truncation_bound + a.truncation_bound + b.truncation_bound);
      Real roundoff =
          Real{4} * (constant_.roundoff_indicator + a.roundoff_indicator + b.roundoff_indicator) +
          Real{64} * eps * (one_ + uni20::math::abs(elementary) + uni20::math::abs(value) + uni20::math::abs(x));
      if (!uni20::isfinite(value) || !uni20::isfinite(roundoff))
        throw bethe::detail::SpinonFunctionError(bethe::detail::SpinonFunctionFailure::arithmetic,
                                                 "nonfinite Barnes transition rate");
      return {value, truncation, roundoff, a.terms + b.terms + terms, a.shifts + b.shifts};
    }

  private:
    bethe::detail::SpinonBarnesContext<Real> functions_;
    bethe::detail::SpinonFunctionEstimate<Real> constant_{};
    Real one_, pi_;
};
} // namespace bethe::heisenberg::detail
