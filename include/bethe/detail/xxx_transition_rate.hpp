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
    std::size_t max_evaluations = 200000;
};

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
// This even rapidity-difference kernel is also needed for four spinons.
template <uni20::Real Real> class XXXTransitionRate {
  public:
    explicit XXXTransitionRate(StructureFactorOptions<Real> options = {}) : options_(options)
    {
      Real const eps = uni20::numeric_limits<Real>::epsilon();
      if (!uni20::isfinite(options.tolerance) || options.tolerance < Real{256} * eps || options.tolerance > Real{0.01})
        throw std::invalid_argument("transition-rate tolerance must be in [256 epsilon, 0.01]");
      end_ = -std::log(options.tolerance / Real{64}) / Real{2};
      for (std::size_t n : {16, 32, 64})
        rules_.push_back(bethe::detail::gauss_legendre<Real>(n));
    }

    TransitionRate<Real> operator()(Real rho) const
    {
      if (!uni20::isfinite(rho) || rho == Real{0})
        throw std::invalid_argument("transition-rate kernel requires finite nonzero rho");
      TransitionRate<Real> out;
      Real const a = Real{4} * std::abs(rho), pi = Real{4} * std::atan(Real{1});
      Real const panels_real = std::ceil(end_ * std::max(Real{1}, a / pi));
      // Bound before converting to an integer (also rejects overflow).
      if (!uni20::isfinite(panels_real) || panels_real > Real(options_.max_evaluations / 48)) return out;
      auto const panels = static_cast<std::size_t>(panels_real);
      Real const width = end_ / Real(panels);
      Real const analytic = a < Real{2} ? std::log1p(a * a / Real{4}) + Real{2} * (std::log(Real{2}) - std::log(a))
                                        : std::log1p(Real{4} / (a * a));
      Real previous{};
      for (auto const& rule : rules_)
      {
        if (panels > (options_.max_evaluations - out.evaluations) / rule.x.size()) return out;
        bethe::detail::CompensatedSum<Real> sum, absolute;
        for (std::size_t p = 0; p < panels; ++p)
          for (std::size_t j = 0; j < rule.x.size(); ++j)
          {
            Real const x = width * (Real(p) + (Real{1} + rule.x[j]) / Real{2});
            Real const u = std::exp(-Real{2} * x), cosine = std::cos(a * x);
            Real f;
            if (x < Real{0.25})
            {
              Real const t = std::tanh(x), sine = std::sin(a * x / Real{2});
              f = (-Real{2} * (Real{1} + t) / std::tanh(Real{2} * x) * sine * sine + t + t * t +
                   Real{2} * std::expm1(-Real{2} * x) * cosine) /
                  x;
            }
            else
              f = (Real{2} * u * u * (Real{3} - u * u) * cosine - Real{4} * u) /
                  ((Real{1} - u) * (Real{1} + u) * (Real{1} + u) * x);
            Real const term = width / Real{2} * rule.w[j] * f;
            sum.add(term);
            absolute.add(std::abs(term));
          }
        out.evaluations += panels * rule.x.size();
        Real const value = sum.value();
        Real const error =
            std::abs(value - previous) + Real{8} * std::exp(-Real{2} * end_) / end_ +
            Real{32} * uni20::numeric_limits<Real>::epsilon() * (Real{1} + absolute.value() + std::abs(analytic));
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
