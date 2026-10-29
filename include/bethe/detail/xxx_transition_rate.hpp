// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <algorithm>
#include <bethe/detail/gauss_legendre.hpp>
#include <bethe/solver.hpp>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace bethe::heisenberg
{
template <uni20::Real Real> struct StructureFactorOptions
{
    // Absolute target for I(rho), hence approximately relative for exp(-I).
    Real tolerance = Real{4096} * uni20::numeric_limits<Real>::epsilon();
    std::size_t max_evaluations = 200000, max_order = 64;
    uni20::Precision precision() const { return uni20::native_precision<Real>(); }
};

#if UNI20_HAS_MPREAL
template <> struct StructureFactorOptions<uni20::mpreal>
{
    explicit StructureFactorOptions(uni20::Precision p) : tolerance(default_tolerance(p)), precision_(p) {}
    uni20::mpreal tolerance;
    std::size_t max_evaluations = 200000, max_order = 256;
    uni20::Precision precision() const { return precision_; }

  private:
    static uni20::mpreal default_tolerance(uni20::Precision p)
    {
      if (p.is_exact()) throw std::invalid_argument("a structure factor requires finite working precision");
      return 4096 * uni20::numeric_limits<uni20::mpreal>::epsilon(p);
    }
    uni20::Precision precision_;
};
#endif

namespace detail
{
template <uni20::Real Real> struct TransitionRate
{
    std::optional<Real> log_rate, error;
    std::size_t evaluations = 0;
};

// Caux--Hagemans (2006), Eq. (12). Subtract
//   2 (1-exp(-2x)) cos(a*x)/x, a=4*rho,
// whose Abel integral is log(1+4/a^2). The remainder is exponentially
// decaying, not conditionally convergent. For x>=1 its absolute tail
// beyond L is bounded by 8 exp(-2L)/L. No tabulated fp64 constants.
// Even rapidity-difference kernel, independent of continuum kinematics.
template <uni20::Real Real> class XXXTransitionRate {
  public:
    explicit XXXTransitionRate(StructureFactorOptions<Real> options = {}) : options_(options)
    {
      Real const zero = uni20::scalar_at_precision<Real>(0, options.precision());
      bethe::detail::require_same_precision(options.tolerance, zero);
      Real const eps = uni20::numeric_limits<Real>::epsilon(options.precision());
      if (!uni20::isfinite(options.tolerance) || options.tolerance < Real{256} * eps ||
          options.tolerance > uni20::scalar_like(zero, 1) / Real{100})
        throw std::invalid_argument("transition-rate tolerance must be in [256 epsilon, 0.01]");
      end_ = -uni20::math::log(options.tolerance / Real{64}) / Real{2};
      if (options.max_order < 32 || options.max_order > 2048 || (options.max_order & (options.max_order - 1)))
        throw std::invalid_argument("maximum quadrature order must be a power of two in [32,2048]");
      for (std::size_t n = 16; n <= options.max_order; n *= 2)
        rules_.push_back(bethe::detail::gauss_legendre<Real>(n, options.precision()));
    }

    TransitionRate<Real> operator()(Real rho) const
    {
      Real const zero = uni20::scalar_at_precision<Real>(0, options_.precision()), one = uni20::scalar_like(zero, 1);
      bethe::detail::require_same_precision(rho, zero);
      if (!uni20::isfinite(rho) || rho == Real{0})
        throw std::invalid_argument("transition-rate kernel requires finite nonzero rho");
      TransitionRate<Real> out;
      Real const a = Real{4} * uni20::math::abs(rho), pi = Real{4} * uni20::math::atan(one);
      Real const panels_real = uni20::math::ceil(end_ * std::max(Real{1}, a / pi));
      // Bound before converting to an integer (also rejects overflow).
      if (!uni20::isfinite(panels_real) || panels_real > Real(options_.max_evaluations / 48)) return out;
      auto const panels = static_cast<std::size_t>(panels_real);
      Real const width = end_ / Real(panels);
      Real const analytic = a < Real{2} ? uni20::math::log1p(a * a / Real{4}) +
                                              Real{2} * (uni20::math::log(Real{2} * one) - uni20::math::log(a))
                                        : uni20::math::log1p(Real{4} / (a * a));
      Real previous = zero;
      for (auto const& rule : rules_)
      {
        if (panels > (options_.max_evaluations - out.evaluations) / rule.x.size()) return out;
        bethe::detail::CompensatedSum<Real> sum(zero), absolute(zero);
        for (std::size_t p = 0; p < panels; ++p)
          for (std::size_t j = 0; j < rule.x.size(); ++j)
          {
            Real const x = width * (Real(p) + (Real{1} + rule.x[j]) / Real{2});
            Real const u = uni20::math::exp(-Real{2} * x), cosine = uni20::math::cos(a * x);
            Real f;
            if (x < Real{0.25})
            {
              Real const t = uni20::math::tanh(x), sine = uni20::math::sin(a * x / Real{2});
              f = (-Real{2} * (Real{1} + t) / uni20::math::tanh(Real{2} * x) * sine * sine + t + t * t +
                   Real{2} * uni20::math::expm1(-Real{2} * x) * cosine) /
                  x;
            }
            else
              f = (Real{2} * u * u * (Real{3} - u * u) * cosine - Real{4} * u) /
                  ((Real{1} - u) * (Real{1} + u) * (Real{1} + u) * x);
            Real const term = width / Real{2} * rule.w[j] * f;
            sum.add(term);
            absolute.add(uni20::math::abs(term));
          }
        out.evaluations += panels * rule.x.size();
        Real const value = sum.value();
        Real const error = uni20::math::abs(value - previous) + Real{8} * uni20::math::exp(-Real{2} * end_) / end_ +
                           Real{32} * uni20::numeric_limits<Real>::epsilon(options_.precision()) *
                               (Real{1} + absolute.value() + uni20::math::abs(analytic));
        if (rule.x.size() > 16 && error <= options_.tolerance && uni20::isfinite(value))
        {
          out.log_rate = -analytic - value; // log(exp(-I)), not the extra 1/2 in Szz
          out.error = error;
          return out;
        }
        previous = value;
      }
      return out;
    }

  private:
    StructureFactorOptions<Real> options_;
    Real end_{};
    std::vector<bethe::detail::GaussLegendreRule<Real>> rules_;
};
} // namespace detail
} // namespace bethe::heisenberg
