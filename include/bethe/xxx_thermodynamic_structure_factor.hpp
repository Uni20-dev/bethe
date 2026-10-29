// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/xxx_transition_rate.hpp>
#include <bethe/spinon.hpp>

namespace bethe::heisenberg
{
enum class SpectralDensityStatus
{
  converged,
  outside_continuum,
  lower_threshold,
  numerical_failure
};

template <uni20::Real Real> struct SpectralDensity
{
    // Szz itself (not Szz/(2*pi)), in inverse energy units. The divergent
    // lower threshold has no finite value; exact upper-edge limit is zero.
    std::optional<Real> value, error;
    std::size_t evaluations = 0;
    SpectralDensityStatus status = SpectralDensityStatus::numerical_failure;
    bool converged() const { return status != SpectralDensityStatus::numerical_failure; }
};

/// Infinite XXX chain, J>0, zero field and temperature; exact two-spinon
/// contribution only. No broadening, finite-N roots, or sum-rule rescaling.
/// Caux--Hagemans (2006), Eqs. (9)--(12). Immutable/thread-safe evaluator.
template <uni20::Real Real> class ThermodynamicTwoSpinonStructureFactor {
  public:
    explicit ThermodynamicTwoSpinonStructureFactor(Real exchange = Real{1}, StructureFactorOptions<Real> options = {})
        : exchange_(exchange), kernel_(options)
    {
      bethe::detail::require_same_precision(exchange, options.tolerance);
      if (!uni20::isfinite(exchange) || exchange <= Real{0})
        throw std::invalid_argument("XXX structure factor requires finite J>0");
    }
    std::pair<Real, Real> boundaries(Real q) const
    {
      bethe::detail::require_same_precision(q, exchange_);
      Real const pi = Real{4} * uni20::math::atan(uni20::scalar_like(exchange_, 1));
      if (!uni20::isfinite(q) || q < Real{0} || q > Real{2} * pi)
        throw std::invalid_argument("momentum must be in [0,2*pi]");
      q = std::min(q, Real{2} * pi - q);
      Real const upper = exchange_ * (pi * uni20::math::sin(q / Real{2}));
      Real const lower = q == pi ? uni20::scalar_like(exchange_, 0) : upper * uni20::math::cos(q / Real{2});
      if (!uni20::isfinite(upper)) throw std::overflow_error("continuum energy exceeds scalar range");
      if (q > Real{0} && upper < positive_minimum())
        throw std::underflow_error("continuum energy is below the normal scalar range");
      return {lower, upper};
    }
    SpectralDensity<Real> operator()(Real q, Real omega) const
    {
      bethe::detail::require_same_precision(omega, exchange_);
      if (!uni20::isfinite(omega)) throw std::invalid_argument("frequency must be finite");
      auto const [lower, upper] = boundaries(q);
      SpectralDensity<Real> out;
      Real const zero = uni20::scalar_like(exchange_, 0);
      if (upper == Real{0} || omega < lower || omega > upper)
        return {zero, zero, 0, SpectralDensityStatus::outside_continuum};
      if (upper == lower) return out; // unresolved narrow continuum
      if (omega == lower) return {{}, {}, 0, SpectralDensityStatus::lower_threshold};
      if (omega == upper) return {zero, zero, 0, SpectralDensityStatus::converged};
      // Factored differences and asinh avoid subtracting 1 inside acosh
      // close to the upper edge. Normalize by upper to avoid energy squares.
      Real const above = (omega - lower) / upper * (omega / upper + lower / upper);
      Real const below = (upper - omega) / upper * (Real{1} + omega / upper);
      Real const rho = uni20::math::asinh(uni20::math::sqrt(below) / uni20::math::sqrt(above)) /
                       (Real{4} * uni20::math::atan(uni20::scalar_like(exchange_, 1)));
      if (!(rho > Real{0}) || !uni20::isfinite(rho)) return out;
      auto const rate = kernel_(rho);
      out.evaluations = rate.evaluations;
      if (!rate.log_rate) return out;
      Real const log_value = *rate.log_rate - uni20::math::log(uni20::scalar_like(exchange_, 2)) -
                             uni20::math::log(upper) - uni20::math::log(below) / Real{2};
      Real const value = uni20::math::exp(log_value);
      if (!uni20::isfinite(value) || value < positive_minimum()) return out;
      out.value = value;
      out.error =
          value * (*rate.error + Real{32} * uni20::numeric_limits<Real>::epsilon(uni20::precision_of(exchange_)) *
                                     (Real{1} + uni20::math::abs(log_value)));
      out.status = SpectralDensityStatus::converged;
      return out;
    }

  private:
    Real positive_minimum() const
    {
      if constexpr (uni20::has_runtime_precision_v<Real>)
        return uni20::math::next_up(uni20::scalar_like(exchange_, 0));
      else
        return uni20::numeric_limits<Real>::min();
    }
    Real exchange_;
    detail::XXXTransitionRate<Real> kernel_;
};
} // namespace bethe::heisenberg
