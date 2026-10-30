// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/spinon_numerics.hpp>
#include <uni20/core/complex_log_gamma.hpp>
#include <uni20/core/scalar_precision.hpp>
#include <vector>

namespace bethe::detail
{
/// Internal right-half-plane log|Barnes G|. Recurrence/asymptotic expansion
/// with a complex remainder majorant; numerical bounds are not intervals.
/// Gamma is delegated to Uni20, and no ambient MPFR precision is changed.
template <uni20::Real Real> class SpinonBarnesContext {
  public:
    using Complex = uni20::complex<Real>;
    explicit SpinonBarnesContext(uni20::log_gamma_options options = {})
      requires(!uni20::has_runtime_precision_v<Real>)
        : SpinonBarnesContext(uni20::native_precision<Real>(), options)
    {}

    explicit SpinonBarnesContext(uni20::Precision precision, uni20::log_gamma_options options = {})
        : output_(precision), working_(working_precision(precision)), options_(options), gamma_(working_, options)
    {
      std::size_t count = 16;
#if UNI20_HAS_MPREAL
      if constexpr (uni20::has_runtime_precision_v<Real>)
      {
        count = static_cast<std::size_t>(working_.bit_count() / 6) + 8;
        if (count > options.max_terms) throw std::invalid_argument("Barnes precision exceeds coefficient budget");
        // Shared exact Bernoulli recurrence; internal dependency is pinned with
        // Uni20. No finite table caps the MPFR precision.
        uni20::detail::exact_bernoulli_sequence sequence(2 * count);
        for (std::size_t n = 0; n <= 2 * count; ++n)
        {
          auto b = sequence.next();
          if (n >= 2 && n % 2 == 0)
            coefficients_.push_back(uni20::scalar_at_precision<Real>(
                b / (uni20::exact_rational{n} * uni20::exact_rational{n - 1}), working_));
        }
      }
      else
#endif
      {
        if (options.max_terms < count) throw std::invalid_argument("Barnes requires sixteen native coefficients");
        auto table = uni20::detail::native_stirling_coefficients<Real>();
        coefficients_.assign(table.begin(), table.end());
      }
      log_a_ = barnes_core(Complex(r(1), r(0)));
    }

    SpinonFunctionEstimate<Real> log_abs_barnes_g_right(Complex const& input) const
    {
      if (!uni20::isfinite(input) || input.real() <= 0) throw std::invalid_argument("Barnes G requires finite Re(z)>0");
      if constexpr (uni20::has_runtime_precision_v<Real>)
        if ((!input.real().is_exact() && input.real().precision() != output_) ||
            (!input.imag().is_exact() && input.imag().precision() != output_))
          throw std::invalid_argument("Barnes G input must match context precision");
      auto z = uni20::scalar_at_precision<Complex>(input, working_);
      auto result = barnes_core(z);
      result.value -= log_a_.value;
      result.truncation_bound += log_a_.truncation_bound;
      result.roundoff_indicator += log_a_.roundoff_indicator + eps() * uni20::math::abs(log_a_.value);
      auto output_epsilon = uni20::scalar_at_precision<Real>(uni20::numeric_limits<Real>::epsilon(output_), working_);
      result.roundoff_indicator += output_epsilon * (r(1) + uni20::math::abs(result.value));
      for (auto const& x : {result.value, result.truncation_bound, result.roundoff_indicator})
        if (!uni20::isfinite(x))
          throw SpinonFunctionError(SpinonFunctionFailure::arithmetic, "nonfinite Barnes evaluation");
      result.value = uni20::scalar_at_precision<Real>(result.value, output_);
      result.truncation_bound = uni20::scalar_at_precision<Real>(result.truncation_bound, output_);
      result.roundoff_indicator = uni20::scalar_at_precision<Real>(result.roundoff_indicator, output_);
      return result;
    }

  private:
    static uni20::Precision working_precision(uni20::Precision p)
    {
      if (p.is_exact()) throw std::invalid_argument("Barnes context requires finite precision");
      if constexpr (uni20::has_runtime_precision_v<Real>)
      {
#if UNI20_HAS_MPREAL
        if (p.bit_count() < 24 || p.bit_count() > MPFR_PREC_MAX - 128)
          throw std::invalid_argument("Barnes precision needs at least 24 bits and room for guard bits");
        return uni20::Precision::bits(p.bit_count() + 64);
#endif
      }
      else if (p != uni20::native_precision<Real>())
        throw std::invalid_argument("Barnes native context requires native precision");
      return p;
    }
    Real r(std::size_t n) const { return uni20::scalar_at_precision<Real>(n, working_); }
    Real eps() const { return uni20::numeric_limits<Real>::epsilon(working_); }
    uni20::Precision output_, working_;
    uni20::log_gamma_options options_;
    uni20::log_gamma_context<Real> gamma_;
    std::vector<Real> coefficients_;
    SpinonFunctionEstimate<Real> log_a_;

    Real power(Real x, std::size_t n) const
    {
      Real y = this->r(1);
      while (n)
      {
        if (n & 1) y *= x;
        n >>= 1;
        if (n) x *= x;
      }
      return y;
    }

    // DLMF 5.11(ii) and Nemes (2014), Eq. (1.10). Secant factors
    // bound complex remainders; a small first omitted term alone is insufficient.
    Real remainder(Complex const& z) const
    {
      Real magnitude = uni20::math::abs(z);
      Real sec2 = this->r(2) / (this->r(1) + z.real() / magnitude);
      Real inverse = this->r(1) / magnitude;
      Real ratio = sec2 * inverse * inverse;
      auto k = coefficients_.size();
      auto n = k - 1;
      return uni20::math::abs(coefficients_.back()) / this->r(2 * n) * this->power(ratio, n) * uni20::math::sqrt(sec2);
    }

    // DLMF 5.17.5, without -log A. Shift w=z+s and use
    // sum_j log Gamma(z+j) = s log Gamma(w) - sum_j (j+1) log(z+j).
    // This needs ONE Gamma evaluation, not s separate evaluations.
    SpinonFunctionEstimate<Real> barnes_core(Complex z) const
    {
      Complex w = z;
      SpinonSum<Real> recurrence(this->r(0));
      Real mass = this->r(0), tail;
      std::size_t shifts = 0;
      while (w.real() <= this->r(1) || (tail = this->remainder(w - this->r(1))) > this->eps() / this->r(16))
      {
        if (shifts == options_.max_shifts)
          throw SpinonFunctionError(SpinonFunctionFailure::budget, "Barnes G shift budget exhausted");
        Real term = this->r(shifts + 1) * uni20::math::log(w).real();
        recurrence.add(term);
        mass += uni20::math::abs(term);
        w += this->r(1);
        ++shifts;
      }
      Complex v = w - this->r(1);
      auto lg = gamma_.log_gamma_right(w);
      Complex first = v * v / this->r(4), second = (z - this->r(1)) * lg.value;
      Complex third = -(v * w / this->r(2) + this->r(1) / this->r(12)) * uni20::math::log(v);
      SpinonSum<Real> sum(recurrence.value());
      for (auto const& term : {first, second, third})
      {
        sum.add(term.real());
        mass += uni20::math::abs(term);
      }
      Complex inverse = this->r(1) / v, square = inverse * inverse, p = square;
      for (std::size_t k = 1; k + 1 < coefficients_.size(); ++k)
      {
        Complex term = coefficients_[k] / this->r(2 * k) * p;
        sum.add(term.real());
        mass += uni20::math::abs(term);
        p *= square;
      }
      Real factor = uni20::math::abs(z - this->r(1));
      return {sum.value(), tail + factor * lg.truncation_bound,
              factor * lg.roundoff_indicator + this->r(64) * this->eps() * (this->r(1) + mass),
              coefficients_.size() - 2, shifts + lg.shifts};
    }
};
} // namespace bethe::detail
